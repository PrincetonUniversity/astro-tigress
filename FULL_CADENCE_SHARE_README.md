# TIGRESS full-domain, high-cadence MHD data

This is the README draft for the collaborator-facing Globus collection at:

```text
/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share
```

Before publication, replace this paragraph with the collection UUID, release
date, content version, maintainer contact, and final manifest link.

## Included data

| Model | Shared outputs | Resolution | Full domain |
|---|---:|---:|---|
| `R8_2pc` | `0285-0448` | 2 pc | `[-512,512] x [-512,512] x [-3584,3584] pc` |
| `R8_4pc` | `0200-0674` currently; policy limit `0700` | 4 pc | `[-512,512] x [-512,512] x [-3584,3584] pc` |

Each output contains the complete native MHD VTK processor set, matching
`starpar.vtk`, and canonical `.par`, `.hst`, and `.sn` files. The `R8_4pc`
source currently ends at `0674`; outputs through `0700` will be added only if
they become available and pass completeness validation.

Restart files, plots, FITS and pickle products, projections, slices, maps,
scripts, logs, executables, and other analysis or administrative files are
intentionally excluded.

## Difference from the public release

The DOI-backed [TIGRESS public data release](https://doi.org/10.34770/ackh-7y71)
is the stable dataset supported by the
[astro-tigress documentation](https://princetonuniversity.github.io/astro-tigress/).
It has selected snapshots at roughly 10 Myr intervals and curated products such
as chemistry, CO emission, and ionizing-radiation results.

Its public `MHD` product is one joined central-cube file:

```text
MODEL/NNNN/MHD/MODEL.NNNN.vtk
```

That cube covers `x,y,z = [-512,512] pc`, with `512^3` cells for `R8_2pc` and
`256^3` cells for `R8_4pc`.

This collection provides the denser native MHD sequence over the full vertical
domain. It preserves 56 processor pieces for `R8_2pc` and 28 for `R8_4pc`:

```text
R8_2pc/snapshots/0300/
├── R8_2pc_rst.par
├── hst/
│   ├── R8_2pc_rst.hst
│   └── R8_2pc_rst.sn
├── starpar/R8_2pc_rst.0300.starpar.vtk
├── id0/R8_2pc_rst.0300.vtk
├── id1/R8_2pc_rst-id1.0300.vtk
├── ...
└── id55/R8_2pc_rst-id55.0300.vtk
```

`R8_4pc` has the same layout through `id27`. One directory per snapshot lets
Globus users download only the outputs they need.

## Download selected snapshots

Select one or more four-digit directories under:

```text
R8_2pc/snapshots/
R8_4pc/snapshots/
```

Transfer recursively and preserve the internal structure. Approximate VTK sizes
are 42 GiB per `R8_2pc` snapshot and 5.25 GiB per `R8_4pc` snapshot.

## Read with pyathena

[pyathena](https://github.com/jeonggyukim/pyathena) is recommended for this
native layout:

```python
import pyathena as pa

snapshot_dir = "/path/to/R8_2pc/snapshots/0300"
sim = pa.LoadSim(snapshot_dir)
ds = sim.load_vtk(num=300, id0=True, load_method="xarray")
starpar = sim.load_starpar_vtk(num=300)

print(ds.domain["le"], ds.domain["re"])
print(ds.domain["Nx"], ds.domain["dx"])
```

Expected `R8_2pc` values are:

```text
le = [-512, -512, -3584] pc
re = [ 512,  512,  3584] pc
Nx = [512, 512, 3584]
dx = [2, 2, 2] pc
```

For `R8_4pc`, the edges are identical, with `Nx = [256,256,1792]` and
`dx = [4,4,4] pc`.

The public-release readers in `astro-tigress` assume the public
snapshot/product layout. Use pyathena here unless you explicitly rearrange or
join the native files.

pyathena also reads its per-snapshot uncompressed tar format, but tar archives
are not stored here because they would duplicate the VTK data. Users who need a
single archive may create one after download.

## Citation

For MHD data, cite Kim & Ostriker (2017),
[ApJ 846, 133](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract),
and the public release DOI,
[10.34770/ackh-7y71](https://doi.org/10.34770/ackh-7y71).

See the public documentation for citations associated with its separate
chemistry, CO-line, and ionizing-radiation products.
