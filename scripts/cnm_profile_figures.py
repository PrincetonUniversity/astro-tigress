#!/usr/bin/env python
"""Combine per-snapshot CNM profiles in a directory and write the figure set.

This is the headless / batch version of
``notebooks/cnm_vertical_profiles.ipynb``.  Given a ``--datadir`` of per-snapshot
NetCDF files (written by ``cnm_vertical_profiles.py`` or
``cnm_vertical_profiles_pyathena.py``), it performs the *same* calculations as
the notebook — combine over snapshots, time median/mean and 16-84 percentiles,
locate the CNM layer (``f_M > 1%``) — and writes the same six figures to PNG.

The notebook remains the interactive example; this script reproduces its figures
for any model's ``datadir`` in one command.  Point ``--datadir`` at the folder
of per-snapshot files (by default ``<basedir>/cnm_profiles`` next to the model
data); figures land in ``<datadir>/figures`` alongside them:

    python scripts/cnm_profile_figures.py --datadir /path/to/MODEL/cnm_profiles

Figures written (to ``--figdir``, default ``<datadir>/figures``):
    cnm_thermodynamic_profiles.png   cnm_layer_extent_vs_time.png
    cnm_velocity_dispersions.png     cnm_layer_thickness.png
    cnm_filling_mass_fraction.png    cnm_layer_profiles.png
"""

import os
import os.path as osp
import glob
import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

plt.switch_backend("Agg")
plt.rcParams.update({"figure.dpi": 100, "font.size": 11})

FMTHR = 0.01  # mass-fraction threshold defining the CNM layer


# ---------------------------------------------------------------------------
# data handling
# ---------------------------------------------------------------------------
def load_combined(datadir, model=None):
    """Concatenate all per-snapshot profile files along the snapshot axis."""
    pat = f"{model}.cnm_zprof.*.nc" if model else "*.cnm_zprof.*.nc"
    files = sorted(glob.glob(osp.join(datadir, pat)))
    if not files:
        raise SystemExit(f"no files matching '{pat}' in {datadir}")
    # if the model was not given, make sure there is only one in the directory
    prefixes = {osp.basename(f).split(".cnm_zprof.")[0] for f in files}
    if model is None and len(prefixes) > 1:
        raise SystemExit(f"multiple models in {datadir}: {sorted(prefixes)}; "
                         "pass --model to choose one")
    comb = xr.concat([xr.open_dataset(f) for f in files], dim="ivtk")
    return comb, files, prefixes.pop() if len(prefixes) == 1 else model


def add_effective_temperature(comb):
    """Add the kinetic temperature and effective pressure from existing vars.

    ``T_kmax = T (sigma_eff_1D / cs)^2`` is the temperature a purely thermal gas
    would need to reproduce the *full* 1D effective (thermal+turbulent)
    dispersion -- the maximum kinetic temperature.  ``nT_kmax = Pth
    (sigma_eff_1D / cs)^2`` is the matching effective (thermal+turbulent)
    pressure ``n T_kmax``.  Both reduce to ``T`` and ``Pth`` when turbulence
    vanishes (``sigma_eff_1D -> cs``); no reprocessing of the raw data is needed.
    """
    ratio = (comb["sigma_eff_1D"] / comb["cs"]) ** 2
    comb["T_kmax"] = comb["T"] * ratio
    comb["nT_kmax"] = comb["Pth"] * ratio
    return comb


def cnm_layer(comb, fMthr=FMTHR):
    """Per-snapshot outermost heights (kpc) where f_M exceeds the threshold."""
    fM = comb["f_M"].values                 # (n_snap, n_z)
    zk = comb["z"].values / 1e3             # kpc
    tt = np.asarray(comb["t_Myr"].values)   # Myr
    above = fM > fMthr
    zbot = np.full(fM.shape[0], np.nan)
    ztop = np.full(fM.shape[0], np.nan)
    for i in range(fM.shape[0]):
        zi = zk[above[i]]
        if zi.size:
            zbot[i], ztop[i] = zi.min(), zi.max()
    return zbot, ztop, ztop - zbot, tt


# ---------------------------------------------------------------------------
# plotting helpers (identical to the notebook)
# ---------------------------------------------------------------------------
def plot_band(ax, stat, name, zkpc, color="C0", logy=False, label=None):
    # median line + shaded 16-84 percentile band vs height (z on the x-axis)
    med = stat[name].sel(quantile=0.50)
    lo = stat[name].sel(quantile=0.16)
    hi = stat[name].sel(quantile=0.84)
    ax.fill_between(zkpc, lo, hi, color=color, alpha=0.25, lw=0)
    ax.plot(zkpc, med, color=color, lw=2, label=label)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(r"$z\ [{\rm kpc}]$")
    ax.grid(alpha=0.3)
    return ax


def plot_layer(ax, stat, meands, zkpc, zsel, zL, zH, series,
               logy=False, title="", ylabel="", ymin=None):
    # one or more series in one panel; each is (name, color, label) drawn as
    # median (solid) + mean (dashed) + 16-84 band, restricted to the mean layer
    chunks = []
    for name, color, label in series:
        med = stat[name].sel(quantile=0.50)
        lo = stat[name].sel(quantile=0.16)
        hi = stat[name].sel(quantile=0.84)
        mean = meands[name]
        ax.fill_between(zkpc, lo, hi, color=color, alpha=0.18, lw=0)
        ax.plot(zkpc, med, color=color, lw=2, label=label)
        ax.plot(zkpc, mean, color=color, lw=1.5, ls="--")
        chunks += [lo.values[zsel], hi.values[zsel], mean.values[zsel]]
    ax.set_xlim(zL, zH)
    if logy:
        ax.set_yscale("log")
    vals = np.concatenate(chunks)
    vals = vals[np.isfinite(vals)]
    if vals.size:
        vmax = float(vals.max())
        vmin = ymin if ymin is not None else float(vals.min())
        if logy:
            ax.set_ylim(ymin if ymin is not None else vmin / 1.4, vmax * 1.4)
        else:
            pad = 0.10 * (vmax - vmin) + 1e-12
            ax.set_ylim(ymin if ymin is not None else vmin - pad, vmax + pad)
    ax.set_xlabel(r"$z\ [{\rm kpc}]$")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    if len(series) > 1:
        ax.legend(fontsize=8, loc="best")


# ---------------------------------------------------------------------------
# figures (each returns a Figure)
# ---------------------------------------------------------------------------
def fig_thermo(stat, zkpc):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True)
    plot_band(axes[0], stat, "nH", zkpc, color="C0", logy=True)
    axes[0].set_ylabel(r"$n_{\rm H}\ [{\rm cm^{-3}}]$")
    axes[0].set_title("CNM number density")
    plot_band(axes[1], stat, "Pth", zkpc, color="C3", logy=True)
    axes[1].set_ylabel(r"$P_{\rm th}/k_B\ [{\rm K\,cm^{-3}}]$")
    axes[1].set_title("CNM thermal pressure")
    plot_band(axes[2], stat, "T", zkpc, color="C1")
    axes[2].set_ylabel(r"$T\ [{\rm K}]$")
    axes[2].set_title("CNM temperature (mass-weighted)")
    fig.suptitle(r"CNM vertical profiles ($T<500$ K): median and 16-84%", y=1.02)
    fig.tight_layout()
    return fig


def fig_velocity(stat, zkpc):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True, sharey=True)
    for name, color, lab in [("sigma_x", "C0", r"$\sigma_x$"),
                             ("sigma_y", "C2", r"$\sigma_y$"),
                             ("sigma_z", "C3", r"$\sigma_z$")]:
        plot_band(axes[0], stat, name, zkpc, color=color, label=lab)
    axes[0].set_ylabel(r"$\sigma\ [{\rm km\,s^{-1}}]$")
    axes[0].set_title("Turbulent velocity dispersions")
    axes[0].legend()
    plot_band(axes[1], stat, "sigma_z", zkpc, color="C3", label=r"$\sigma_z$ (kinetic)")
    plot_band(axes[1], stat, "cs", zkpc, color="C1", label=r"$c_s$ (thermal)")
    plot_band(axes[1], stat, "sigma_eff_z", zkpc, color="k",
              label=r"$\sigma_{\rm eff,z}=\sqrt{\sigma_z^2+c_s^2}$")
    axes[1].set_ylabel(r"$\sigma\ [{\rm km\,s^{-1}}]$")
    axes[1].set_title("Vertical support of the CNM")
    axes[1].legend()
    plot_band(axes[2], stat, "sigma_turb_1D", zkpc, color="C0",
              label=r"$\sigma_{\rm turb,1D}$")
    plot_band(axes[2], stat, "sigma_eff_1D", zkpc, color="k",
              label=r"$\sigma_{\rm eff,1D}$")
    axes[2].set_ylabel(r"$\sigma\ [{\rm km\,s^{-1}}]$")
    axes[2].set_title("1D effective dispersion")
    axes[2].legend()
    for ax in axes:
        ax.set_ylim(0, 10)
    fig.tight_layout()
    return fig


def fig_fractions(stat, zkpc):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    plot_band(axes[0], stat, "f_A", zkpc, color="C4")
    axes[0].set_ylabel(r"$f_A$")
    axes[0].set_title("CNM area filling fraction")
    axes[0].set_ylim(0, None)
    plot_band(axes[1], stat, "f_M", zkpc, color="C5")
    axes[1].set_ylabel(r"$f_M$")
    axes[1].set_title("CNM mass fraction")
    axes[1].set_ylim(0, None)
    fig.suptitle("Where is the CNM?", y=1.02)
    fig.tight_layout()
    return fig


def fig_layer_extent(zbot, ztop, thick, tt):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    axes[0].vlines(tt, zbot, ztop, color="C5", lw=7, alpha=0.4)
    axes[0].plot(tt, ztop, "o-", color="C3", label=r"$z_{\rm top}$")
    axes[0].plot(tt, zbot, "o-", color="C0", label=r"$z_{\rm bot}$")
    axes[0].axhline(0, color="k", lw=0.7, ls=":")
    axes[0].set_xlabel(r"$t\ [{\rm Myr}]$")
    axes[0].set_ylabel(r"$z\ [{\rm kpc}]$")
    axes[0].set_title(r"CNM layer extent ($f_M>1\%$)")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].plot(tt, thick, "s-", color="C1")
    axes[1].axhline(np.nanmedian(thick), color="k", ls="--",
                    label=f"median = {np.nanmedian(thick):.2f} kpc")
    axes[1].set_xlabel(r"$t\ [{\rm Myr}]$")
    axes[1].set_ylabel(r"$\Delta z_{\rm CNM}\ [{\rm kpc}]$")
    axes[1].set_title("CNM layer thickness")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    return fig


def fig_layer_summary(stat, zkpc, thick):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    med = stat["f_M"].sel(quantile=0.50)
    axes[0].fill_between(zkpc, stat["f_M"].sel(quantile=0.16),
                         stat["f_M"].sel(quantile=0.84), color="C5", alpha=0.25, lw=0)
    axes[0].plot(zkpc, med, color="C5", lw=2)
    axes[0].axhline(FMTHR, color="k", ls="--", label=r"$f_M=1\%$")
    zk_arr = np.asarray(zkpc)
    regm = np.asarray(med) > FMTHR
    if regm.any():
        axes[0].axvspan(zk_arr[regm].min(), zk_arr[regm].max(), color="C5", alpha=0.12)
    axes[0].set_xlabel(r"$z\ [{\rm kpc}]$")
    axes[0].set_ylabel(r"$f_M$")
    axes[0].set_title(r"Median $f_M$ and the $>1\%$ layer")
    axes[0].set_ylim(0, None)
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].hist(thick[np.isfinite(thick)],
                 bins=np.linspace(0, np.nanmax(thick) * 1.1 + 1e-6, 9),
                 color="C1", edgecolor="k", alpha=0.8)
    axes[1].axvline(np.nanmedian(thick), color="k", ls="--",
                    label=f"median = {np.nanmedian(thick):.2f} kpc")
    axes[1].set_xlabel(r"$\Delta z_{\rm CNM}\ [{\rm kpc}]$")
    axes[1].set_ylabel("number of snapshots")
    axes[1].set_title("CNM layer thickness distribution")
    axes[1].legend()
    fig.tight_layout()
    return fig


def fig_inlayer(stat, meands, zkpc, zbot, ztop):
    zL, zH = float(np.nanmean(zbot)), float(np.nanmean(ztop))
    zsel = (np.asarray(zkpc) >= zL) & (np.asarray(zkpc) <= zH)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    a = (stat, meands, zkpc, zsel, zL, zH)
    # top row: density; thermal T with kinetic T_kmax; thermal & effective pressure
    plot_layer(axes[0, 0], *a, [("nH", "C0", "")], logy=True, ymin=1.0,
               title="number density", ylabel=r"$n_{\rm H}\ [{\rm cm^{-3}}]$")
    plot_layer(axes[0, 1], *a,
               [("T", "C1", r"$T$ (thermal)"),
                ("T_kmax", "C4", r"$T_{k,\max}$")],
               logy=True, title="temperature", ylabel=r"$T\ [{\rm K}]$")
    plot_layer(axes[0, 2], *a,
               [("Pth", "C3", "thermal"),
                ("nT_kmax", "C4", r"$n\,T_{k,\max}$")],
               logy=True, ymin=5.0e2,
               title="pressure", ylabel=r"$P/k_B\ [{\rm K\,cm^{-3}}]$")
    # bottom row: sound speed + 1D turbulence, 1D effective, volume + mass fractions
    plot_layer(axes[1, 0], *a,
               [("cs", "C4", r"$c_s$"),
                ("sigma_turb_1D", "C2", r"$\sigma_{\rm turb,1D}$")],
               title="sound speed & turbulence", ylabel=r"$[{\rm km\,s^{-1}}]$")
    plot_layer(axes[1, 1], *a, [("sigma_eff_1D", "k", "")],
               title="1D effective dispersion",
               ylabel=r"$\sigma_{\rm eff,1D}\ [{\rm km\,s^{-1}}]$")
    plot_layer(axes[1, 2], *a,
               [("f_A", "C4", r"$f_A$ (volume)"),
                ("f_M", "C5", r"$f_M$ (mass)")],
               title="CNM volume & mass fractions", ylabel="fraction", ymin=0.0)
    fig.suptitle(r"CNM profiles within the mean CNM layer ($f_M>1\%$) — "
                 r"solid: median, dashed: mean, shaded: 16-84%", y=1.01)
    fig.tight_layout()
    return fig, zL, zH


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--datadir", required=True,
                   help="directory with <model>.cnm_zprof.*.nc files")
    p.add_argument("--model", default=None,
                   help="model name (inferred from filenames if omitted)")
    p.add_argument("--figdir", default=None,
                   help="output directory for figures (default: <datadir>/figures)")
    args = p.parse_args(argv)

    comb, files, model = load_combined(args.datadir, args.model)
    comb = add_effective_temperature(comb)
    figdir = args.figdir or osp.join(args.datadir, "figures")
    os.makedirs(figdir, exist_ok=True)

    stat = comb.quantile([0.16, 0.50, 0.84], dim="ivtk", skipna=True)
    meands = comb.mean(dim="ivtk", skipna=True)
    zkpc = comb["z"] / 1e3
    zbot, ztop, thick, tt = cnm_layer(comb)

    print(f"model={model}  {len(files)} snapshots  "
          f"t={float(comb.t_Myr.min()):.1f}-{float(comb.t_Myr.max()):.1f} Myr")

    saved = []

    def save(fig, name):
        path = osp.join(figdir, name)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)

    save(fig_thermo(stat, zkpc), "cnm_thermodynamic_profiles.png")
    save(fig_velocity(stat, zkpc), "cnm_velocity_dispersions.png")
    save(fig_fractions(stat, zkpc), "cnm_filling_mass_fraction.png")
    save(fig_layer_extent(zbot, ztop, thick, tt), "cnm_layer_extent_vs_time.png")
    save(fig_layer_summary(stat, zkpc, thick), "cnm_layer_thickness.png")
    fig, zL, zH = fig_inlayer(stat, meands, zkpc, zbot, ztop)
    save(fig, "cnm_layer_profiles.png")

    print(f"mean CNM layer: z = {zL:+.3f} .. {zH:+.3f} kpc "
          f"(thickness {zH - zL:.3f} kpc; median {np.nanmedian(thick):.3f} kpc)")
    print(f"wrote {len(saved)} figures to {figdir}")


if __name__ == "__main__":
    main()
