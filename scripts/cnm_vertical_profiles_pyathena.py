#!/usr/bin/env python
"""CNM vertical profiles using the ``pyathena`` backend (any TIGRESS model).

This is the general-backend port of ``scripts/cnm_vertical_profiles.py``.  It
produces the **same** per-snapshot NetCDF artifacts (identical variable names,
units and schema, so the same combine notebook works) but reads the simulation
output through `pyathena <https://github.com/jeonggyukim/pyathena>`_ instead of
``astro_tigress``.  See ``docs/cnm_vertical_profiles_spec.md`` for the full
specification.

Given a model **base directory**, for every VTK snapshot it writes one
self-describing NetCDF holding the horizontally-averaged vertical profiles of the
cold neutral medium (``T < Tcnm``, default 500 K): area/mass fractions, number
density, thermal pressure, temperature, per-component and effective
(thermal+kinetic) velocity dispersions, and sound speed, with the background
shear removed from ``vy``.

It works on both

* the **public TIGRESS data release** (classic Athena 4.2 VTK, no stored
  temperature -> derived from the cooling function; per-snapshot subdirectory
  layout that ``load_vtk`` can't address -> read the discovered file directly), and
* standard pyathena model trees (``load_vtk`` + native ``T`` when present).

Output layout
-------------
By default the per-snapshot NetCDF files are written **beside the model data**,
into ``<basedir>/cnm_profiles/`` (the companion figure script then writes to
``<basedir>/cnm_profiles/figures/``).  Because each model has its own
``basedir`` this is automatically per-model and collision-free.  When the model
directory is read-only (e.g. the public data release), pass ``--outdir`` to a
writable location instead.

Usage
-----
    # writable model tree: outputs land in <basedir>/cnm_profiles/
    python scripts/cnm_vertical_profiles_pyathena.py --basedir /path/to/MODEL

    # read-only public release: choose a writable --outdir
    python scripts/cnm_vertical_profiles_pyathena.py \
        --basedir /projects/EOSTRIKE/TIGRESS_data_release/R8_2pc \
        --outdir cnm_profiles/R8_2pc

Run with ``--help`` for all options.
"""

import os
import os.path as osp
import sys
import time
import argparse

import numpy as np
import xarray as xr

# Physical constants (CGS) and unit helpers -- match the astro_tigress reference.
KB = 1.380658e-16   # Boltzmann constant [erg/K]
MH = 1.673534e-24   # hydrogen mass [g]
KMS = 1.0e5         # cm/s per km/s

_REPO_ROOT = osp.dirname(osp.dirname(osp.abspath(__file__)))
# default location of the general backend (../pyathena_master next to this repo)
_DEFAULT_PYATHENA = osp.join(osp.dirname(_REPO_ROOT), "pyathena_master")

UNITS = {
    "f_A": "", "f_M": "", "nH": "cm**-3", "Pth": "K cm**-3", "T": "K",
    "T_vw": "K", "vx": "km/s", "vy": "km/s", "vz": "km/s",
    "sigma_x": "km/s", "sigma_y": "km/s", "sigma_z": "km/s", "cs": "km/s",
    "sigma_eff_z": "km/s", "sigma_turb_1D": "km/s", "sigma_eff_1D": "km/s",
}


def import_pyathena(pyathena_path):
    """Import ``pyathena``, preferring an explicit checkout on ``sys.path``."""
    if pyathena_path and osp.isdir(pyathena_path):
        if pyathena_path not in sys.path:
            sys.path.insert(0, pyathena_path)
    import pyathena as pa  # noqa: E402

    return pa


def load_cube(pa, s, num, fmap, zmax=None):
    """Return ``(dd, t_Myr)`` for one snapshot.

    ``dd`` is an xarray Dataset with fields ``nH`` [cm^-3], ``pok`` [K cm^-3],
    ``vx, vy, vz`` [km/s] and ``T`` [K], on coordinates ``x, y, z`` [pc].  The
    temperature is the native ``T`` field when the backend can build it,
    otherwise it is derived from the classic two-phase cooling function
    (``T1 = pok / (nH * muH)`` -> ``coolftn.get_temp``).  ``zmax`` (pc) restricts
    the read to ``|z| < zmax`` -- useful for tall boxes.
    """
    # Prefer the framework loader; fall back to a direct read of the discovered
    # file (the data-release layout puts each snapshot in its own subdirectory,
    # which load_vtk cannot address).
    try:
        ds = s.load_vtk(num)
    except Exception:
        ds = pa.AthenaDataSet(fmap[num], units=s.u, dfi=s.dfi)

    le = re = None
    if zmax is not None:  # read only the |z| < zmax slab
        le, re = list(ds.domain["le"]), list(ds.domain["re"])
        le[2], re[2] = -abs(zmax), abs(zmax)

    base_fields = ["nH", "pok", "velocity"]
    try:
        dd = ds.get_field(base_fields + ["T"], le=le, re=re)
        temperature = dd["T"]
    except Exception:  # native T unavailable -> classic cooling function
        dd = ds.get_field(base_fields, le=le, re=re)
        from pyathena.classic.cooling import coolftn

        T1 = dd["pok"] / (dd["nH"] * s.u.muH)
        temperature = xr.apply_ufunc(coolftn().get_temp, T1)

    dd = dd.rename({"velocity1": "vx", "velocity2": "vy", "velocity3": "vz"})
    dd["T"] = temperature
    t_Myr = float(ds.domain["time"]) * s.u.Myr
    return dd, t_Myr


def compute_profiles(dd, t_Myr, num, Omega, qshear, muH, Tcnm=500.0):
    """Horizontally-averaged CNM vertical profiles for one snapshot.

    All reductions are labeled xarray sums over the horizontal dims ``x, y``,
    restricted to cold cells (``T < Tcnm``).  Densities/pressures are
    volume-weighted; temperature and velocity moments are ``nH`` (mass) weighted.
    ``vy`` is measured in the shearing frame, ``dvy = vy + qshear*Omega*x``.
    """
    hor = ["x", "y"]
    cnm = dd["T"] < Tcnm
    w = dd["nH"].where(cnm, 0.0)                       # mass weight, 0 outside CNM
    dvy = dd["vy"] + qshear * Omega * dd["x"]          # shear-subtracted vy

    N = cnm.sum(hor)                                   # cold cell count per plane
    S_w = w.sum(hor)                                   # sum nH over CNM
    S_all = dd["nH"].sum(hor)                          # sum nH over all cells
    S_pok = dd["pok"].where(cnm, 0.0).sum(hor)
    S_T = dd["T"].where(cnm, 0.0).sum(hor)
    S_wT = (w * dd["T"]).sum(hor)
    S_w_safe = S_w.where(S_w > 0)

    def moments(v):
        mean = (w * v).sum(hor) / S_w_safe
        var = (w * v * v).sum(hor) / S_w_safe - mean * mean
        return mean, np.maximum(var, 0.0)             # clip roundoff negatives

    vxm, sx2 = moments(dd["vx"])
    vym, sy2 = moments(dvy)
    vzm, sz2 = moments(dd["vz"])

    cs2 = (KB * S_pok / (muH * MH * S_w_safe)) / (KMS * KMS)     # km^2/s^2
    sig_turb_1d2 = (sx2 + sy2 + sz2) / 3.0

    Nx, Ny = dd.sizes["x"], dd.sizes["y"]
    empty = N == 0
    N_safe = N.where(~empty)

    prof = xr.Dataset(
        dict(
            f_A=N / (Nx * Ny),
            f_M=S_w / S_all.where(S_all > 0),
            nH=(S_w / N_safe),
            Pth=(S_pok / N_safe),
            T=(S_wT / S_w_safe),
            T_vw=(S_T / N_safe),
            vx=vxm, vy=vym, vz=vzm,
            sigma_x=np.sqrt(sx2), sigma_y=np.sqrt(sy2), sigma_z=np.sqrt(sz2),
            cs=np.sqrt(cs2),
            sigma_eff_z=np.sqrt(sz2 + cs2),
            sigma_turb_1D=np.sqrt(sig_turb_1d2),
            sigma_eff_1D=np.sqrt(sig_turb_1d2 + cs2),
        ),
        coords=dict(z=dd["z"]),
    )
    for k, u in UNITS.items():
        if k in prof:
            prof[k].attrs["units"] = u
    prof["z"].attrs.update(units="pc", long_name="height")
    prof.attrs.update(
        ivtk=int(num), t_Myr=float(t_Myr), Tcnm=float(Tcnm),
        Omega_km_s_kpc=float(Omega * 1e3), qshear=float(qshear),
        muH=float(muH), Nx=int(Nx), Ny=int(Ny),
    )
    return prof.assign_coords(t_Myr=float(t_Myr), ivtk=int(num))


def resolve_outdir(outdir):
    """Create ``outdir`` (must be writable) or exit with actionable guidance."""
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
    return outdir


def save_netcdf(ds, fname):
    """Write ``ds`` to NetCDF, trying available engines in order."""
    last = None
    for engine in ("netcdf4", "h5netcdf", "scipy"):
        try:
            ds.to_netcdf(fname, engine=engine)
            return engine
        except Exception as e:
            last = e
    raise RuntimeError(f"Could not write {fname}: {last}")


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--basedir", required=True, help="model base directory")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: <basedir>/cnm_profiles)")
    p.add_argument("--model", default=None,
                   help="model name for filenames (default: basename of basedir)")
    p.add_argument("--nums", default="all",
                   help="comma-separated snapshot numbers, or 'all' (default)")
    p.add_argument("--every", type=int, default=1,
                   help="process every Nth snapshot (stride on the sorted list)")
    p.add_argument("--zmax", type=float, default=None,
                   help="restrict reads to |z| < zmax pc (for tall boxes)")
    p.add_argument("--Tcnm", type=float, default=500.0,
                   help="CNM temperature threshold in K (default 500)")
    p.add_argument("--pyathena", default=_DEFAULT_PYATHENA,
                   help="path to the pyathena checkout to use")
    p.add_argument("--overwrite", action="store_true",
                   help="recompute even if the output file exists")
    args = p.parse_args(argv)

    pa = import_pyathena(args.pyathena)
    print(f"pyathena: {pa.__file__}")

    s = pa.LoadSim(args.basedir, verbose=False)
    model = args.model or osp.basename(osp.normpath(args.basedir))
    # default: write beside the model data (per-model, collision-free)
    outdir = args.outdir or osp.join(args.basedir, "cnm_profiles")
    outdir = resolve_outdir(outdir)

    Omega = float(s.par["problem"]["Omega"])     # km/s/pc (code units)
    qshear = float(s.par["problem"]["qshear"])
    muH = float(s.u.muH)
    fmap = dict(zip(s.nums, s.files.get("vtk", [])))  # only used for the direct-read fallback

    nums = list(s.nums) if args.nums == "all" else [int(x) for x in args.nums.split(",")]
    nums = nums[:: args.every]
    print(f"model={model}  Tcnm={args.Tcnm:g} K  Omega={Omega:g} km/s/pc  "
          f"qshear={qshear:g}  muH={muH:g}")
    print(f"outdir={outdir}  processing {len(nums)} snapshot(s): {nums}")

    for num in nums:
        fname = osp.join(outdir, f"{model}.cnm_zprof.{num:04d}.nc")
        if osp.isfile(fname) and not args.overwrite:
            print(f"[{num:04d}] exists, skip")
            continue
        t0 = time.time()
        dd, t_Myr = load_cube(pa, s, num, fmap, zmax=args.zmax)
        prof = compute_profiles(dd, t_Myr, num, Omega, qshear, muH, Tcnm=args.Tcnm)
        engine = save_netcdf(prof, fname)
        fmid = float(prof["f_A"].sel(z=0, method="nearest"))
        print(f"[{num:04d}] t={t_Myr:.1f} Myr  <f_A>_mid={fmid:.3f}  "
              f"wrote {osp.basename(fname)} ({engine}) in {time.time() - t0:.1f}s")

    print("done.")


if __name__ == "__main__":
    main()
