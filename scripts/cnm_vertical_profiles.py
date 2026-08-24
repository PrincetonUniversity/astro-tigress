#!/usr/bin/env python
"""Compute vertical profiles of the Cold Neutral Medium (CNM) for TIGRESS R8_2pc.

This standalone script reads the released ``MHD`` snapshots of a TIGRESS model
(default ``R8_2pc``) through :class:`astro_tigress.Model`, selects the CNM by a
temperature threshold (``T < Tcnm``, default 500 K), and reduces the 3D data to
one-dimensional *vertical* profiles by averaging horizontally (over the ``x`` and
``y`` axes) at each height ``z``.

For every snapshot it writes one self-describing NetCDF file holding the profiles
of

    * area filling fraction of the CNM                      ``f_A``
    * mean hydrogen number density                          ``nH``       [cm^-3]
    * mean thermal pressure  P/k_B                           ``Pth``      [K cm^-3]
    * mean temperature (mass- and volume-weighted)          ``T``, ``T_vw`` [K]
    * mean bulk velocities                                  ``vx,vy,vz`` [km/s]
    * velocity dispersions (mass-weighted)                  ``sigma_x/y/z`` [km/s]
    * thermal sound speed                                   ``cs``       [km/s]
    * effective (thermal+kinetic) vertical dispersion       ``sigma_eff_z`` [km/s]
    * 1D turbulent / effective dispersions                  ``sigma_turb_1D``,
                                                            ``sigma_eff_1D`` [km/s]

See ``docstring`` of :func:`compute_profiles` for the exact mathematical
definitions; the companion notebook ``notebooks/cnm_vertical_profiles.ipynb``
reproduces them in LaTeX.

Usage
-----
    python scripts/cnm_vertical_profiles.py \
        --data-root /projects/EOSTRIKE/TIGRESS_data_release \
        --model R8_2pc --outdir cnm_profiles

Run ``python scripts/cnm_vertical_profiles.py --help`` for all options.
"""

import os
import os.path as osp
import sys
import gc
import time
import argparse

import numpy as np

# Make the in-repo astro_tigress package importable when run from anywhere.
_REPO_ROOT = osp.dirname(osp.dirname(osp.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import astro_tigress as at  # noqa: E402
from astro_tigress import const  # noqa: E402
from astro_tigress.ytathena import muH  # mean molecular weight per H  # noqa: E402


# ----------------------------------------------------------------------------
# Physical constants (CGS) -- taken from astro_tigress.const for consistency
# ----------------------------------------------------------------------------
KB = const.kb  # Boltzmann constant [erg/K]
MH = const.mh  # hydrogen mass [g]
KMS = 1.0e5  # cm/s per km/s


def compute_profiles(model, ivtk, Tcnm=500.0, Omega=28.0e-3, qshear=1.0):
    r"""Compute CNM vertical profiles for a single snapshot.

    The snapshot is loaded onto a level-0 covering grid of shape
    ``(Nx, Ny, Nz)`` (axes = ``x, y, z``).  At each height ``z`` we average over
    the horizontal plane, but only over cells belonging to the CNM,

    .. math:: \mathrm{CNM}(z) \equiv \{(x,y) : T(x,y,z) < T_{\rm cnm}\}.

    Let :math:`w_i = n_{{\rm H},i}` be the (volume-element-normalised) mass weight
    of cell ``i`` in the plane (the cell volume is constant on a uniform grid, so
    weighting by :math:`n_{\rm H}` is mass weighting).  The profiles are

    * Area filling fraction  :math:`f_A = N_{\rm CNM}/(N_x N_y)`.
    * Volume-weighted density :math:`\langle n_{\rm H}\rangle
      = \frac{1}{N_{\rm CNM}}\sum_i n_{{\rm H},i}`.
    * Volume-weighted thermal pressure
      :math:`\langle P_{\rm th}/k_B\rangle = \frac{1}{N_{\rm CNM}}\sum_i P_i/k_B`.
    * Mass-weighted temperature
      :math:`\langle T\rangle_M = \sum_i w_i T_i / \sum_i w_i`.
    * Mass-weighted mean velocity
      :math:`\bar v_\alpha = \sum_i w_i v_{\alpha,i}/\sum_i w_i` and dispersion
      :math:`\sigma_\alpha^2 = \sum_i w_i v_{\alpha,i}^2/\sum_i w_i - \bar v_\alpha^2`
      for :math:`\alpha \in \{x,y,z\}`.  For the azimuthal component the
      background shear is removed, :math:`\delta v_y = v_y + q\Omega x`.
    * Thermal sound speed (mass-weighted)
      :math:`c_s^2 = \langle P/\rho\rangle_M = \sum_i P_i / \sum_i \rho_i`.
    * Effective vertical dispersion
      :math:`\sigma_{\rm eff,z} = \sqrt{\sigma_z^2 + c_s^2}`.

    Parameters
    ----------
    model : astro_tigress.Model
        Model handle (already knows the data directory layout).
    ivtk : int
        VTK output number of the snapshot.
    Tcnm : float, optional
        CNM temperature threshold in K (default 500).
    Omega : float, optional
        Galactic angular velocity in km/s/pc (default 0.028 = 28 km/s/kpc).
    qshear : float, optional
        Shear parameter :math:`q=-d\ln\Omega/d\ln R` (default 1, flat rotation).

    Returns
    -------
    dict
        Mapping of profile name -> 1D ``numpy`` array (length ``Nz``), plus the
        ``z`` coordinate and scalar metadata.
    """
    model.load(ivtk, "MHD")
    grid = model.MHD.grid

    # --- coordinates (1D) -------------------------------------------------
    # Covering-grid arrays are indexed [i, j, k] = (x, y, z).
    x1d = grid["index", "x"][:, 0, 0].to("pc").d.astype(np.float64)  # (Nx,)
    z1d = grid["index", "z"][0, 0, :].to("pc").d.astype(np.float64)  # (Nz,)

    # Sanity check on axis ordering (cheap; only std over the reduced axes).
    assert np.allclose(grid["index", "z"][:, :, 0].std(), 0.0), "z not along axis 2"
    assert np.allclose(grid["index", "x"][0, :, :].std(), 0.0), "x not along axis 0"

    Nx, Ny, Nz = grid.ActiveDimensions
    Nxy = float(Nx * Ny)

    # --- fields as float32 numpy (memory friendly, sums done in float64) ---
    T = grid["gas", "temperature"].to("K").d.astype(np.float32)
    nH = grid["gas", "nH"].to("cm**-3").d.astype(np.float32)
    pok = grid["gas", "pok"].to("K/cm**3").d.astype(np.float32)
    vx = grid["gas", "velocity_x"].to("km/s").d.astype(np.float32)
    vy = grid["gas", "velocity_y"].to("km/s").d.astype(np.float32)
    vz = grid["gas", "velocity_z"].to("km/s").d.astype(np.float32)

    # Remove the background shear from the azimuthal velocity:
    #   v_y,bg = -q*Omega*x   ->   delta v_y = v_y + q*Omega*x
    dvy = vy + np.float32(qshear * Omega) * x1d[:, None, None]
    del vy

    # --- CNM mask and masked weighted sums over the horizontal plane ------
    ax = (0, 1)  # reduce over x and y, keep z
    m = T < np.float32(Tcnm)  # boolean CNM mask
    wm = np.where(m, nH, np.float32(0.0))  # mass weight, zero outside CNM

    def psum(a):
        """Plane sum over (x, y) accumulated in float64."""
        return a.sum(axis=ax, dtype=np.float64)

    N = m.sum(axis=ax, dtype=np.float64)  # cell count per plane
    S_w = psum(wm)  # sum n_H over CNM   (=Sum rho / (muH mh))
    S_nH_all = psum(nH)  # sum n_H over ALL cells (for mass fraction)
    S_pok = psum(np.where(m, pok, np.float32(0.0)))
    S_T = psum(np.where(m, T, np.float32(0.0)))
    S_wT = psum(wm * T)

    S_wvx = psum(wm * vx)
    S_wvx2 = psum(wm * vx * vx)
    S_wvy = psum(wm * dvy)
    S_wvy2 = psum(wm * dvy * dvy)
    S_wvz = psum(wm * vz)
    S_wvz2 = psum(wm * vz * vz)

    # free the big arrays early
    del T, nH, pok, vx, vz, dvy, wm, m
    gc.collect()

    # --- assemble profiles (guard against empty planes) -------------------
    empty = N == 0
    Nsafe = np.where(empty, np.nan, N)
    Wsafe = np.where(S_w <= 0, np.nan, S_w)

    def variance(S_wv2, S_wv):
        mean = S_wv / Wsafe
        var = S_wv2 / Wsafe - mean * mean
        return mean, np.clip(var, 0.0, None)  # guard tiny negative from roundoff

    vxm, sx2 = variance(S_wvx2, S_wvx)
    vym, sy2 = variance(S_wvy2, S_wvy)
    vzm, sz2 = variance(S_wvz2, S_wvz)

    # thermal sound speed: cs^2 = Sum P / Sum rho = kB*Sum(pok) / (muH*mh*Sum(nH))
    cs2_cgs = KB * S_pok / (muH * MH * Wsafe)  # [cm^2/s^2]
    cs2 = cs2_cgs / (KMS * KMS)  # [km^2/s^2]

    sig_turb_1d2 = (sx2 + sy2 + sz2) / 3.0

    prof = {
        "z": z1d,
        "f_A": N / Nxy,  # CNM area filling fraction
        "f_M": S_w / np.where(S_nH_all > 0, S_nH_all, np.nan),  # CNM mass fraction
        "nH": S_w / Nsafe,  # volume-weighted mean density
        "Pth": S_pok / Nsafe,  # volume-weighted mean P/kB
        "T": S_wT / Wsafe,  # mass-weighted temperature
        "T_vw": S_T / Nsafe,  # volume-weighted temperature
        "vx": vxm,
        "vy": vym,
        "vz": vzm,
        "sigma_x": np.sqrt(sx2),
        "sigma_y": np.sqrt(sy2),
        "sigma_z": np.sqrt(sz2),
        "cs": np.sqrt(cs2),
        "sigma_eff_z": np.sqrt(sz2 + cs2),
        "sigma_turb_1D": np.sqrt(sig_turb_1d2),
        "sigma_eff_1D": np.sqrt(sig_turb_1d2 + cs2),
    }
    # blank out fully empty planes for the count-normalised profiles
    prof["nH"][empty] = np.nan
    prof["Pth"][empty] = np.nan
    prof["T_vw"][empty] = np.nan

    meta = dict(
        ivtk=int(ivtk),
        t_Myr=float(model.t_Myr[np.where(model.ivtks == ivtk)[0][0]]),
        Tcnm=float(Tcnm),
        Omega_km_s_kpc=float(Omega * 1e3),
        qshear=float(qshear),
        muH=float(muH),
        Nx=int(Nx),
        Ny=int(Ny),
        Nz=int(Nz),
        model=model.model_id,
    )

    # release the yt dataset / covering grid before the next snapshot
    del grid
    del model.MHD
    gc.collect()

    return prof, meta


def to_dataset(prof, meta):
    """Pack the profile dict into an :class:`xarray.Dataset`."""
    import xarray as xr

    z = prof["z"]
    data_vars = {}
    units = {
        "f_A": "",
        "f_M": "",
        "nH": "cm**-3",
        "Pth": "K cm**-3",
        "T": "K",
        "T_vw": "K",
        "vx": "km/s",
        "vy": "km/s",
        "vz": "km/s",
        "sigma_x": "km/s",
        "sigma_y": "km/s",
        "sigma_z": "km/s",
        "cs": "km/s",
        "sigma_eff_z": "km/s",
        "sigma_turb_1D": "km/s",
        "sigma_eff_1D": "km/s",
    }
    for k, v in prof.items():
        if k == "z":
            continue
        data_vars[k] = (
            "z",
            np.asarray(v, dtype=np.float64),
            {"units": units.get(k, "")},
        )
    ds = xr.Dataset(
        data_vars,
        coords={"z": ("z", z, {"units": "pc", "long_name": "height"})},
        attrs=meta,
    )
    # add a time coordinate for convenient concatenation later
    ds = ds.assign_coords(t_Myr=meta["t_Myr"], ivtk=meta["ivtk"])
    return ds


def save_netcdf(ds, fname):
    """Write ``ds`` to NetCDF, trying available engines in order."""
    last = None
    for engine in ("netcdf4", "h5netcdf", "scipy"):
        try:
            ds.to_netcdf(fname, engine=engine)
            return engine
        except Exception as e:  # engine missing / unsupported -> try next
            last = e
    raise RuntimeError(f"Could not write {fname}: {last}")


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--data-root",
        default="/projects/EOSTRIKE/TIGRESS_data_release",
        help="master directory containing <model>/ subdirectories",
    )
    p.add_argument("--model", default="R8_2pc", help="model id (default R8_2pc)")
    p.add_argument(
        "--ivtks", default="all", help="comma-separated vtk numbers, or 'all' (default)"
    )
    p.add_argument(
        "--outdir",
        default=None,
        help="output directory (default: <data-root>/<model>/cnm_profiles, "
        "i.e. beside the model data; pass an explicit path if that is read-only)",
    )
    p.add_argument(
        "--Tcnm",
        type=float,
        default=500.0,
        help="CNM temperature threshold in K (default 500)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="recompute even if the output file exists",
    )
    args = p.parse_args(argv)

    # default: write beside the model data (per-model, collision-free)
    outdir = args.outdir or osp.join(args.data_root, args.model, "cnm_profiles")
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as e:
        raise SystemExit(
            f"cannot create output directory {outdir!r}: {e}\n"
            "The model directory is likely read-only (e.g. the public data "
            "release). Re-run with --outdir <writable path>."
        )
    if not os.access(outdir, os.W_OK):
        raise SystemExit(
            f"output directory {outdir!r} is not writable; "
            "re-run with --outdir <writable path>."
        )
    args.outdir = outdir
    model = at.Model(args.model, args.data_root)

    if args.ivtks == "all":
        ivtks = list(model.ivtks)
    else:
        ivtks = [int(s) for s in args.ivtks.split(",")]

    print(f"model={args.model}  Tcnm={args.Tcnm:g} K  outdir={args.outdir}")
    print(f"processing {len(ivtks)} snapshot(s): {ivtks}")

    for ivtk in ivtks:
        fname = osp.join(args.outdir, f"{args.model}.cnm_zprof.{ivtk:04d}.nc")
        if osp.isfile(fname) and not args.overwrite:
            print(f"[{ivtk:04d}] exists, skip -> {fname}")
            continue
        t0 = time.time()
        prof, meta = compute_profiles(model, ivtk, Tcnm=args.Tcnm)
        ds = to_dataset(prof, meta)
        engine = save_netcdf(ds, fname)
        dt = time.time() - t0
        print(
            f"[{ivtk:04d}] t={meta['t_Myr']:.1f} Myr  "
            f"<f_A>_mid={float(ds['f_A'].sel(z=0, method='nearest')):.3f}  "
            f"wrote {osp.basename(fname)} ({engine}) in {dt:.1f}s"
        )

    print("done.")


if __name__ == "__main__":
    main()
