# CNM vertical-profile analysis — porting specification

This document specifies the **Cold Neutral Medium (CNM) vertical-profile analysis**
so it can be reproduced for *any* TIGRESS-like model, using a general reading
backend (here [`pyathena`](../../pyathena_master)) in place of the
`astro_tigress` package.

The goal is a script that, **given the base directory of a model**, produces the
same *CNM-profile artifacts* that were built for `R8_2pc`:

1. one **NetCDF file per snapshot** holding the horizontally-averaged vertical
   profiles of the CNM, and
2. a **notebook / figure set** that combines all snapshots and shows time
   statistics (mean, median, 16–84 percentiles), including the profiles zoomed
   into the mean CNM layer.

The reference implementation for `R8_2pc` lives in
[`scripts/cnm_vertical_profiles.py`](../scripts/cnm_vertical_profiles.py) and
[`notebooks/cnm_vertical_profiles.ipynb`](../notebooks/cnm_vertical_profiles.ipynb)
(built by [`scripts/_build_cnm_notebook.py`](../scripts/_build_cnm_notebook.py)).
This spec generalizes them.

---

## 1. Physics and definitions

All quantities are computed **per height `z`** by averaging over the horizontal
(`x`, `y`) plane, restricted to the cold gas.

### 1.1 CNM selection

At each height the CNM is the set of cold cells

$$
\mathcal{C}(z) = \{(x,y) : T(x,y,z) < T_{\rm cnm}\},\qquad T_{\rm cnm}=500\ \mathrm{K},
$$

with $N_{\rm CNM}(z)=|\mathcal{C}(z)|$ cold cells out of $N_x N_y$. $T_{\rm cnm}$
is a parameter (500 K = the classic two-phase cold threshold).

### 1.2 Weighting

The grid is uniform, so cell volume $\Delta V$ is constant and weighting a plane
sum by $n_{\rm H}$ is **mass weighting** ($\mathrm{d}m=\mu_{\rm H}m_{\rm H}n_{\rm H}\Delta V$).
Two horizontal averages are used over $i\in\mathcal{C}(z)$:

$$
\langle q\rangle_V = \frac{1}{N_{\rm CNM}}\sum_i q_i \quad\text{(volume-weighted)},\qquad
\langle q\rangle_M = \frac{\sum_i n_{{\rm H},i}\,q_i}{\sum_i n_{{\rm H},i}} \quad\text{(mass-weighted)}.
$$

### 1.3 Profiles (functions of `z`)

| Quantity | Symbol | Definition | Weighting |
|---|---|---|---|
| Area filling fraction | `f_A` | $N_{\rm CNM}/(N_xN_y)$ | — |
| Mass fraction | `f_M` | $\sum_{\mathcal C} n_{\rm H} / \sum_{\rm all} n_{\rm H}$ | — |
| Number density | `nH` | $\langle n_{\rm H}\rangle_V$ | volume |
| Thermal pressure | `Pth` | $\langle P/k_B\rangle_V$ | volume |
| Temperature | `T` (`T_vw`) | $\langle T\rangle_M$ ($\langle T\rangle_V$) | mass (volume) |
| Bulk velocities | `vx,vy,vz` | $\langle v_\alpha\rangle_M$ | mass |
| Velocity dispersions | `sigma_x/y/z` | $\sigma_\alpha^2=\langle v_\alpha^2\rangle_M-\langle v_\alpha\rangle_M^2$ | mass |
| Sound speed | `cs` | $c_s^2=\langle P/\rho\rangle_M=\dfrac{k_B\sum_{\mathcal C}(P/k_B)}{\mu_{\rm H}m_{\rm H}\sum_{\mathcal C}n_{\rm H}}$ | mass |
| Effective vertical | `sigma_eff_z` | $\sqrt{\sigma_z^2+c_s^2}$ | — |
| 1D turbulent / effective | `sigma_turb_1D`, `sigma_eff_1D` | $\sqrt{(\sigma_x^2+\sigma_y^2+\sigma_z^2)/3}$, $\sqrt{\sigma_{\rm turb,1D}^2+c_s^2}$ | — |

### 1.4 Shear frame (critical)

TIGRESS is a **local shearing box**: the azimuthal velocity `vy` contains the
background shear $v_{y,\rm bg}=-q\,\Omega\,x$. It **must** be removed before
computing `vy`/`sigma_y`, otherwise the box-wide shear (~$\Omega L_x$) inflates
the dispersion:

$$
\delta v_y = v_y + q\,\Omega\,x .
$$

Use the model's own $\Omega$ and $q=q_{\rm shear}$ (see §2.3). `vx` and `vz` need
no correction.

### 1.5 Empty planes

Where a plane has no CNM, the count-normalized profiles (`nH`, `Pth`, `T_vw`) and
the mass-weighted quantities are `NaN`; `f_A` and `f_M` are `0`. Downstream time
statistics skip NaNs (`skipna=True`).

---

## 2. Backend interface

The analysis needs only a small set of capabilities from the reading backend.
Any framework that can provide the **adapter** below will work.

### 2.1 Required capabilities (abstract)

| Need | Description |
|---|---|
| `list_snapshots()` | available VTK output numbers for the model |
| `time_Myr(num)` | physical time of a snapshot in Myr |
| `get_cube(num)` | an xarray Dataset on the uniform grid with fields `nH` [cm⁻³], `T` [K], `pok` [K cm⁻³], `vx,vy,vz` [km s⁻¹], and coordinates `x,y,z` [pc] |
| `Nx, Ny` | horizontal cell counts (for `f_A`) |
| `Omega, qshear` | shear parameters ($\Omega$ in km s⁻¹ pc⁻¹, i.e. code units) |
| `muH` | mean molecular weight per H |

### 2.2 pyathena mapping

```python
import pyathena as pa

s = pa.LoadSim(basedir)                 # or pa.LoadSimTIGRESSNCR(basedir, verbose=False)

nums   = s.nums                         # list_snapshots()
Omega  = s.par['problem']['Omega']      # km/s/pc (code units) — see §2.3
qshear = s.par['problem']['qshear']
muH    = s.u.muH                        # mean molecular weight per H
Nx, Ny = s.domain['Nx'][0], s.domain['Nx'][1]

ds  = s.load_vtk(num)                   # one snapshot
t   = ds.domain['time'] * s.u.Myr       # time_Myr(num)
dd  = ds.get_field(['nH', 'T', 'pok', 'velocity'])   # xarray Dataset
# the vector 'velocity' expands to velocity1/2/3; rename to abstract names
dd = dd.rename({'velocity1': 'vx', 'velocity2': 'vy', 'velocity3': 'vz'})
```

`get_field` returns an **xarray Dataset** with dims `(z, y, x)` and coordinates
`x, y, z` in pc, so every reduction below is a *labeled* reduction over `x` and
`y` — no axis-index bookkeeping, and the shear term `qshear*Omega*dd['x']`
broadcasts automatically.

### 2.3 `Omega` units — verify per model

In Athena-TIGRESS code units (length = pc, velocity = km s⁻¹), $\Omega$ has units
km s⁻¹ pc⁻¹, so `q*Omega*x_pc` is in km s⁻¹ directly (for R8, `Omega ≈ 0.028`,
`qshear = 1`). **Check** `s.par['problem']['Omega']` for each model; if a
framework stores $\Omega$ differently (e.g. km s⁻¹ kpc⁻¹ or 1/Myr), convert so
that `Omega * x[pc]` yields km s⁻¹. A quick check: without the correction
`sigma_y` should be $\sim\Omega L_x/\sqrt{12}$ larger than with it.

### 2.4 astro_tigress ↔ pyathena equivalents

| abstract | astro_tigress (reference) | pyathena |
|---|---|---|
| load model | `Model(id, dir_master)` | `pa.LoadSim(basedir)` |
| snapshots | `model.ivtks` | `s.nums` |
| time | `model.t_Myr` | `ds.domain['time']*s.u.Myr` |
| load cube | `model.load(ivtk,'MHD'); grid=model.MHD.grid` | `ds=s.load_vtk(num); dd=ds.get_field([...])` |
| fields | `grid['gas','nH'\|'pok'\|'temperature'\|'velocity_x'...]` (numpy `[x,y,z]`) | `dd['nH'\|'pok'\|'T'\|'vx'...]` (xarray `(z,y,x)`) |
| Omega,q,muH | hard-coded 0.028, 1, 1.4271 | `s.par['problem'][...]`, `s.u.muH` |

---

## 3. Per-snapshot algorithm (xarray reference)

```python
import numpy as np
import xarray as xr

KB, MH, KMS = 1.380658e-16, 1.673534e-24, 1.0e5   # cgs, cm/s per km/s

def compute_profiles(dd, Nx, Ny, t_Myr, num, Omega, qshear, muH, Tcnm=500.0):
    hor = ['x', 'y']                                   # reduce over the plane
    cnm = dd['T'] < Tcnm                               # boolean CNM mask
    w   = dd['nH'].where(cnm, 0.0)                     # mass weight (0 outside)
    dvy = dd['vy'] + qshear * Omega * dd['x']          # shear-subtracted vy

    N      = cnm.sum(hor)                              # cells per plane
    S_w    = w.sum(hor)                                # sum nH over CNM
    S_all  = dd['nH'].sum(hor)                         # sum nH over all cells
    S_pok  = dd['pok'].where(cnm, 0.0).sum(hor)
    S_T    = dd['T'].where(cnm, 0.0).sum(hor)
    S_wT   = (w * dd['T']).sum(hor)

    def moments(v):                                    # mass-weighted mean, var
        m  = (w * v).sum(hor) / S_w
        v2 = (w * v * v).sum(hor) / S_w
        return m, np.maximum(v2 - m * m, 0.0)

    vxm, sx2 = moments(dd['vx'])
    vym, sy2 = moments(dvy)
    vzm, sz2 = moments(dd['vz'])

    cs2 = (KB * S_pok / (muH * MH * S_w)) / (KMS * KMS)     # km^2/s^2
    empty = N == 0

    out = xr.Dataset(
        dict(
            f_A=N / (Nx * Ny),
            f_M=S_w / S_all.where(S_all > 0),
            nH=(S_w / N).where(~empty),
            Pth=(S_pok / N).where(~empty),
            T=S_wT / S_w.where(S_w > 0),
            T_vw=(S_T / N).where(~empty),
            vx=vxm, vy=vym, vz=vzm,
            sigma_x=np.sqrt(sx2), sigma_y=np.sqrt(sy2), sigma_z=np.sqrt(sz2),
            cs=np.sqrt(cs2),
            sigma_eff_z=np.sqrt(sz2 + cs2),
            sigma_turb_1D=np.sqrt((sx2 + sy2 + sz2) / 3),
            sigma_eff_1D=np.sqrt((sx2 + sy2 + sz2) / 3 + cs2),
        ),
        coords=dict(z=dd['z']),
        attrs=dict(num=int(num), t_Myr=float(t_Myr), Tcnm=float(Tcnm),
                   Omega=float(Omega), qshear=float(qshear), muH=float(muH),
                   Nx=int(Nx), Ny=int(Ny)),
    )
    return out.assign_coords(t_Myr=t_Myr)
```

This is the exact numeric content of the reference `astro_tigress` script,
re-expressed with labeled xarray reductions. Attach `units` attributes per
variable (see §4) before writing.

---

## 4. Output artifact — NetCDF per snapshot

Write one file per snapshot to `<outdir>/<model>.cnm_zprof.<num:04d>.nc`.

- **coordinate:** `z` [pc] (and scalar `t_Myr`, `num`).
- **variables** (all along `z`), with `units` attrs:

  | var | units | var | units |
  |---|---|---|---|
  | `f_A` | (—) | `sigma_x/y/z` | km/s |
  | `f_M` | (—) | `cs` | km/s |
  | `nH` | cm⁻³ | `sigma_eff_z` | km/s |
  | `Pth` | K cm⁻³ | `sigma_turb_1D` | km/s |
  | `T`, `T_vw` | K | `sigma_eff_1D` | km/s |
  | `vx/vy/vz` | km/s | | |

- **attrs:** `num`, `t_Myr`, `Tcnm`, `Omega`, `qshear`, `muH`, `Nx`, `Ny`, `model`.
- **engine:** try `netcdf4`, then `h5netcdf`, then `scipy` (NetCDF3 fallback).
- Skip snapshots whose output already exists unless `--overwrite`.

> The reference used field name `ivtk`; keep whatever the backend calls the
> snapshot number, but be consistent between the writer and the notebook.

---

## 5. CNM layer and combined figures (notebook)

Combine all per-snapshot files (`xr.concat` along the snapshot dim, carrying
`t_Myr`), compute `stat = comb.quantile([0.16,0.5,0.84], dim=<snap>, skipna=True)`
and `meands = comb.mean(dim=<snap>, skipna=True)`, then produce:

1. **Thermodynamic profiles** — `nH`, `Pth` (log y), `T`; z on the x-axis, median
   line + 16–84 band.
2. **Velocity dispersions** — `sigma_x/y/z`; `sigma_z`,`cs`,`sigma_eff_z`;
   `sigma_turb_1D`,`sigma_eff_1D`. Clip y to `0–10 km/s`.
3. **Fractions** — `f_A` and `f_M` vs z.
4. **CNM layer extent** (§5.1) — per-snapshot `[z_bot, z_top]` and thickness vs time.
5. **Layer summary** — median `f_M(z)` with the 1% line + shaded layer, and the
   thickness histogram.
6. **In-layer profiles** — every quantity zoomed to the **mean CNM layer**
   $[\langle z_{\rm bot}\rangle,\langle z_{\rm top}\rangle]$, showing **mean
   (dashed) + median (solid) + 16–84 band**, with y-limits auto-scaled to the
   data inside the layer (fixed floors `nH ≥ 1`, `Pth ≥ 5e2`).

Every figure is also saved to `figures/<name>.png`.

### 5.1 CNM layer definition

Per snapshot, with threshold $f_{M,\rm thr}=1\%$:

$$
z_{\rm bot}=\min\{z:f_M(z)>f_{M,\rm thr}\},\quad
z_{\rm top}=\max\{z:f_M(z)>f_{M,\rm thr}\},\quad
\Delta z_{\rm CNM}=z_{\rm top}-z_{\rm bot}.
$$

The **mean CNM layer** used for the zoomed profiles is
$[\langle z_{\rm bot}\rangle,\langle z_{\rm top}\rangle]$ averaged over snapshots.

---

## 6. Running for a new model

```
python cnm_vertical_profiles.py --basedir /path/to/MODEL --outdir cnm_profiles/MODEL
python cnm_vertical_profiles.py --basedir /path/to/MODEL --nums all --Tcnm 500
```

Then point the combine notebook at `cnm_profiles/MODEL` and run it.

Suggested CLI: `--basedir`, `--nums` (`all` or comma list), `--outdir`,
`--Tcnm` (default 500), `--zrange` (optional, restrict `get_field` to
`|z| < zmax` for tall boxes to save memory), `--overwrite`.

---

## 7. Validation checklist

Run one snapshot first and confirm, near the midplane, CNM-like values:

- `nH ~ 10–50 cm⁻³`, `T ~ 50–300 K`, `Pth/kB ~ 10³–10⁴ K cm⁻³`.
- `cs` matches the analytic isothermal sound speed for the mass-weighted `T`
  (e.g. `T≈200 K → cs≈1.1 km/s`); this validates units end-to-end.
- **Shear removed:** `sigma_y` is *not* inflated to $\sim\Omega L_x$ — recompute
  without the `+qΩx` term and confirm it drops back to the turbulent value.
- Grid orientation: confirm the reduction is over the two *horizontal* axes and
  `z` is the vertical coordinate (in pyathena, reduce over dims `x,y`).
- High-|z| planes with almost no CNM are noisy/`NaN` by construction — expected;
  they are skipped by the percentile reduction and `f_A,f_M → 0`.

---

## 8. Notes / choices to keep consistent

- Weighting: densities & pressures volume-weighted; temperature and all velocity
  moments mass-weighted (both `T` weightings are stored).
- `muH` comes from the model (`s.u.muH`), **not** hard-coded.
- `Tcnm=500 K` is the default cold threshold; expose it as a parameter.
- Keep the NetCDF schema (names/units/attrs) identical across models so a single
  notebook works for all of them.
