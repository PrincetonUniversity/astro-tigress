#!/usr/bin/env python
"""Mock HI-emission/absorption + 3D-dust observer experiment (TIGRESS-NCR).

Place an observer at the midplane of a simulation box and look out along many
``(l, b)`` sightlines.  Along each ray build, from the 3D cube,

  * the 21 cm **emission** ``T_B(v)`` and **absorption** ``tau(v)`` profiles
    (Kim, Ostriker & Kim 2014 radiative transfer, spin temperature ``T_s = T_kin``;
    the forward model matches ``astro_tigress.synthetic_HI``), and
  * the **3D dust** differential extinction ``dA_V/dd(d)`` assuming a constant
    dust-to-gas ratio ``N_H / A_V = 1.87e21 cm^-2 mag^-1``.

Select sightlines with a *single* CNM absorption component in velocity **and** a
*single* dust structure in distance, Gaussian-fit each, and record a per-cloud
catalog:

  HI  : v0, dv (FWHM), tau0, N_CNM, T_s, T_kmax(=21.86 dv^2)
  dust: d0, dd (FWHM), A_V
  geom: l, b, z = z_obs + d0 * sin(b)

The companion notebook compares these observer-recovered quantities vs z against
the true ``cnm_zprof`` vertical profiles (ground truth).

Emission/absorption pair: for an isolated CNM component the spin temperature
``T_s = T_B(v0) / (1 - exp(-tau0))`` recovers the true kinetic temperature, while
``T_kmax`` from the linewidth is the (larger) turbulent upper bound.
``N_CNM = 1.823e18 * T_s * integral(tau dv)``.
"""

import os
import os.path as osp
import sys
import argparse

import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

_REPO_ROOT = osp.dirname(osp.dirname(osp.abspath(__file__)))
_DEFAULT_PYATHENA = osp.join(osp.dirname(_REPO_ROOT), "pyathena_master")

PC = 3.085677581467192e18          # cm
NH_AV = 1.87e21                    # N_H / A_V  [cm^-2 mag^-1]
NHI_C = 1.823e18                   # N_HI = NHI_C * T_s * int(tau dv)  [cm^-2]
TKMAX_C = 21.86                    # T_kmax = TKMAX_C * (FWHM[km/s])^2  = m_H/(8 k ln2)
SIG2FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


def import_pyathena(path):
    if path and osp.isdir(path) and path not in sys.path:
        sys.path.insert(0, path)
    import pyathena as pa

    return pa


def load_cube(pa, basedir, num):
    """Load one snapshot: nH, nHI, T, and shear-subtracted velocities (z,y,x)."""
    s = pa.LoadSim(basedir, verbose=False)
    ds = s.load_vtk(num)
    dd = ds.get_field(["nH", "T", "velocity", "xHI"])
    Om = float(s.par["problem"]["Omega"])
    q = float(s.par["problem"]["qshear"])
    cube = dict(
        nH=dd["nH"].values.astype(np.float32),
        nHI=(dd["nH"] * dd["xHI"]).values.astype(np.float32),
        T=dd["T"].values.astype(np.float32),
        vx=dd["velocity1"].values.astype(np.float32),
        vy=(dd["velocity2"] + q * Om * dd["x"]).values.astype(np.float32),
        vz=dd["velocity3"].values.astype(np.float32),
    )
    x, y, z = dd["x"].values, dd["y"].values, dd["z"].values
    dx = float(x[1] - x[0])
    geom = dict(
        dx=dx, Nz=cube["T"].shape[0], Ny=cube["T"].shape[1], Nx=cube["T"].shape[2],
        xmin=float(x[0] - dx / 2), ymin=float(y[0] - dx / 2), zmin=float(z[0] - dx / 2),
        zmax=float(z[-1] + dx / 2),
        t_Myr=float(ds.domain["time"]) * s.u.Myr,
    )
    return cube, geom


def sample_ray(cube, g, obs, l, b, dmax, ds=None):
    """Sample the cube along a ray to distance dmax; periodic wrap in x,y."""
    ds = ds or g["dx"]
    nhat = np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])
    t = np.arange(0.0, dmax, ds)
    p = obs[None, :] + t[:, None] * nhat[None, :]
    ix = (((p[:, 0] - g["xmin"]) / g["dx"]).astype(int)) % g["Nx"]
    iy = (((p[:, 1] - g["ymin"]) / g["dx"]).astype(int)) % g["Ny"]
    iz = ((p[:, 2] - g["zmin"]) / g["dx"]).astype(int)
    ok = (iz >= 0) & (iz < g["Nz"])
    ix, iy, iz, t = ix[ok], iy[ok], iz[ok], t[ok]
    vlos = (cube["vx"][iz, iy, ix] * nhat[0]
            + cube["vy"][iz, iy, ix] * nhat[1]
            + cube["vz"][iz, iy, ix] * nhat[2])
    return dict(t=t, ds=ds, nHI=cube["nHI"][iz, iy, ix], nH=cube["nH"][iz, iy, ix],
                T=cube["T"][iz, iy, ix], vlos=vlos)


def hi_profiles(ray, vchan):
    """Vectorized 21 cm RT along the ray -> (tau(v), T_B(v)); dust dA_V/dd."""
    n, Tk, vlos = ray["nHI"], ray["T"], ray["vlos"]
    vL = 0.21394414 * np.sqrt(Tk)                       # thermal FWHM [km/s]
    dss = ray["ds"] * PC
    phi = 0.00019827867 / vL * np.exp(
        -(1.6651092223153954 * (vchan[:, None] - vlos[None, :]) / vL[None, :]) ** 2)
    tl = (2.6137475e-15 * n / Tk * phi) * dss           # (Nv, Ncell) per-cell tau
    tau = np.nansum(tl, axis=1)
    tc = np.cumsum(tl, axis=1) - tl                     # foreground tau (exclusive)
    TB = np.nansum(Tk[None, :] * (1 - np.exp(-tl)) * np.exp(-tc), axis=1)
    dAv = ray["nH"] * dss / NH_AV                        # per-cell extinction
    return tau, TB, dAv


def _gauss(x, A, x0, sig):
    return A * np.exp(-0.5 * ((x - x0) / sig) ** 2)


def single_component(y, xgrid, ymin, sep_frac=0.35):
    """Return (is_single, peak_index) for a 1D profile with one dominant peak."""
    pk, props = find_peaks(y, height=ymin)
    if len(pk) == 0:
        return False, None
    order = np.argsort(props["peak_heights"])[::-1]
    top = pk[order[0]]
    if len(pk) == 1:
        return True, top
    second = props["peak_heights"][order[1]]
    return (second < sep_frac * y[top]), top


def fit_gauss(x, y, i0, w0):
    try:
        p, _ = curve_fit(_gauss, x, y, p0=[y[i0], x[i0], w0],
                         maxfev=4000)
        return abs(p[0]), p[1], abs(p[2])
    except Exception:
        return None


def process_sightline(cube, g, obs, l, b, vchan, dgrid, dmax,
                      tau_min=0.1, av_min=0.02):
    """Return (catalog dict, profiles) for one cloud, or (None, None)."""
    ray = sample_ray(cube, g, obs, l, b, dmax)
    if ray["nHI"].size < 5:
        return None, None
    tau, TB, dAv = hi_profiles(ray, vchan)
    d = ray["t"]                                          # distance along ray [pc]

    ok_hi, ipk = single_component(tau, vchan, tau_min)
    if not ok_hi:
        return None, None
    hi = fit_gauss(vchan, tau, ipk, 2.0)
    if hi is None:
        return None, None
    tau0, v0, sig_v = hi
    dv = SIG2FWHM * sig_v
    int_tau = np.trapz(tau, vchan)
    # spin temperature from the emission/absorption pair at line center
    TB0 = np.interp(v0, vchan, TB)
    Ts = TB0 / (1.0 - np.exp(-max(tau0, 1e-3)))
    N_CNM = NHI_C * Ts * int_tau
    Tkmax = TKMAX_C * dv * dv

    # dust structure in distance (resample onto a regular grid for fitting)
    dA = np.interp(dgrid, d, dAv, left=0, right=0)
    ok_du, jpk = single_component(dA, dgrid, av_min)
    if not ok_du:
        return None, None
    du = fit_gauss(dgrid, dA, jpk, 40.0)
    if du is None:
        return None, None
    Av0, d0, sig_d = du
    dd_fwhm = SIG2FWHM * sig_d
    A_V = np.trapz(dA, dgrid)
    if d0 <= 0:
        return None, None

    z = obs[2] + d0 * np.sin(b)
    cat = dict(l=np.rad2deg(l), b=np.rad2deg(b), v0=v0, dv=dv, tau0=tau0,
               N_CNM=N_CNM, Ts=Ts, Tkmax=Tkmax, d0=d0, dd=dd_fwhm, A_V=A_V, z=z)
    prof = dict(tau=tau, TB=TB, dA=dA, v0=v0, sig_v=sig_v, tau0=tau0,
                d0=d0, sig_d=sig_d)
    return cat, prof


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--basedir", required=True)
    p.add_argument("--num", type=int, required=True)
    p.add_argument("--pyathena", default=_DEFAULT_PYATHENA)
    p.add_argument("--nobs", type=int, default=8, help="number of observer positions")
    p.add_argument("--nl", type=int, default=48, help="samples in longitude l")
    p.add_argument("--nb", type=int, default=16, help="samples in latitude b")
    p.add_argument("--bmin", type=float, default=20.0)
    p.add_argument("--bmax", type=float, default=80.0)
    p.add_argument("--dmax", type=float, default=500.0,
                   help="max LOS distance from the observer in pc (default 500)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--outfile", default=None)
    args = p.parse_args(argv)

    pa = import_pyathena(args.pyathena)
    print(f"pyathena: {pa.__file__}")
    cube, g = load_cube(pa, args.basedir, args.num)
    print(f"cube {cube['T'].shape}  t={g['t_Myr']:.1f} Myr  zmax={g['zmax']:.0f} pc")

    rng = np.random.default_rng(args.seed)
    Lx = g["Nx"] * g["dx"]
    obs_list = np.column_stack([
        g["xmin"] + rng.uniform(0.2, 0.8, args.nobs) * Lx,
        g["ymin"] + rng.uniform(0.2, 0.8, args.nobs) * (g["Ny"] * g["dx"]),
        np.zeros(args.nobs),
    ])
    lgrid = np.linspace(0, 2 * np.pi, args.nl, endpoint=False)
    bgrid = np.deg2rad(np.concatenate([
        np.linspace(args.bmin, args.bmax, args.nb // 2),
        -np.linspace(args.bmin, args.bmax, args.nb // 2)]))
    vchan = np.linspace(-50, 50, 401)
    dgrid = np.arange(0, args.dmax, g["dx"])          # distance grid, capped at dmax

    cat, examples = [], []
    ntot = args.nobs * len(lgrid) * len(bgrid)
    k = 0
    for obs in obs_list:
        for l in lgrid:
            for b in bgrid:
                k += 1
                c, prof = process_sightline(cube, g, obs, l, b, vchan, dgrid, args.dmax)
                if c is not None:
                    cat.append(c)
                    if len(examples) < 3:
                        examples.append(prof)
        print(f"  observer {int(obs[0]):+5d},{int(obs[1]):+5d}: "
              f"{len(cat)} clouds so far ({k}/{ntot} rays)")

    keys = ["l", "b", "z", "v0", "dv", "tau0", "N_CNM", "Ts", "Tkmax", "d0", "dd", "A_V"]
    arrs = {kk: np.array([c[kk] for c in cat]) for kk in keys}
    ex = {}
    for j, e in enumerate(examples):                 # a few full profiles for the notebook
        for kk in ("tau", "TB", "dA", "v0", "sig_v", "tau0", "d0", "sig_d"):
            ex[f"ex{j}_{kk}"] = e[kk]
    out = args.outfile or osp.join(
        _REPO_ROOT, "cnm_profiles", "mock_obs",
        f"{osp.basename(osp.normpath(args.basedir))}.mockobs.{args.num:04d}.npz")
    os.makedirs(osp.dirname(out), exist_ok=True)
    np.savez(out, t_Myr=g["t_Myr"], dmax=args.dmax, vchan=vchan, dgrid=dgrid,
             n_examples=len(examples), **arrs, **ex)
    print(f"selected {len(cat)} single-component clouds -> {out}")
    if cat:
        z = arrs["z"] / 1e3
        print(f"  z range {z.min():+.2f}..{z.max():+.2f} kpc ; "
              f"median T_s={np.median(arrs['Ts']):.0f} K, "
              f"T_kmax={np.median(arrs['Tkmax']):.0f} K, dv={np.median(arrs['dv']):.1f} km/s")


if __name__ == "__main__":
    main()
