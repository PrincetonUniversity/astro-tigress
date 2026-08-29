# TIGRESS full-domain, high-cadence MHD data

This Globus collection provides collaborator access to raw TIGRESS outputs.

## Included data

| Model | Shared outputs | Resolution | Full domain |
|---|---:|---:|---|
| `R8_2pc` | `0285-0448` | 2 pc | `[-512,512] x [-512,512] x [-3584,3584] pc` |
| `R8_4pc` | `0200-0650` | 4 pc | `[-512,512] x [-512,512] x [-3584,3584] pc` |

Each model contains the complete native MHD VTK processor files for every
listed output, matching `starpar.vtk` files, and one canonical `.par`, `.hst`,
and `.sn` file. Restart files and personal analysis products are excluded.
The entries in this collection are allowlisted symbolic links to immutable
source files; the data bytes are not duplicated.

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

This collection instead provides the high-cadence native MHD sequence over the
full vertical domain. It preserves 56 processor pieces for `R8_2pc` and 28 for
`R8_4pc`:

```text
R8_2pc/
├── R8_2pc_rst.par
├── hst/
│   ├── R8_2pc_rst.hst
│   └── R8_2pc_rst.sn
├── starpar/
│   ├── R8_2pc_rst.0285.starpar.vtk
│   └── ...
├── id0/
│   ├── R8_2pc_rst.0285.vtk
│   └── ...
├── id1/
│   ├── R8_2pc_rst-id1.0285.vtk
│   └── ...
└── ... id55/
```

`R8_4pc` has the same layout through `id27`.

## Download selected snapshots

To download only output `NNNN`, select:

- `<run>.NNNN.vtk` from `id0`;
- `<run>-idN.NNNN.vtk` from every nonzero `idN` directory; and
- `<run>.NNNN.starpar.vtk` from `starpar`.

Preserve the `idN` and `starpar` directory structure. The `.par`, `.hst`, and
`.sn` files are small and may be downloaded once. Approximate VTK sizes are
42 GiB per `R8_2pc` snapshot and 5.25 GiB per `R8_4pc` snapshot.

### Globus batch transfer

For repeated or multi-snapshot transfers, use the
[Globus CLI](https://docs.globus.org/cli/) to create a batch list. Each line has
a source-relative path followed by its destination-relative path. For example,
the following Bash commands create a list for `R8_2pc` outputs `0300` and
`0310`:

```bash
model=R8_2pc
run=R8_2pc_rst
last_rank=55
snapshots=(0300 0310)
batch="${model}-selected.batch"

: > "${batch}"
add_path() {
    printf '%s\t%s\n' "$1" "$1" >> "${batch}"
}

add_path "${run}.par"
add_path "hst/${run}.hst"
add_path "hst/${run}.sn"

for snapshot in "${snapshots[@]}"; do
    add_path "starpar/${run}.${snapshot}.starpar.vtk"
    add_path "id0/${run}.${snapshot}.vtk"
    for rank in $(seq 1 "${last_rank}"); do
        add_path "id${rank}/${run}-id${rank}.${snapshot}.vtk"
    done
done
```

For `R8_4pc`, use `model=R8_4pc`, `run=R8_4pc_newacc`, and
`last_rank=27`. Change the `snapshots` array to the desired output numbers.
Inspect the generated list, then submit it using collection UUIDs and paths
appropriate to your accounts:

```bash
globus login
globus transfer \
    --batch "${batch}" \
    --verify-checksum \
    --preserve-mtime \
    "SOURCE_COLLECTION_UUID:/${model}/" \
    "DESTINATION_COLLECTION_UUID:/desired/path/${model}/"
```

The source and destination prefixes are prepended to each path in the batch
file. The repeated relative path in each row preserves the native `idN`,
`starpar`, and `hst` layout at the destination. Use `globus task show TASK_ID`
to monitor the returned task ID.

The collection uses file symlinks. Globus follows individually transferred
file symlinks and creates ordinary files at the destination, but recursive
transfers generally skip symlinks. Select the required files explicitly or use
a batch list; do not recursively transfer an `idN` directory and assume every
link will be followed. Report any `Path not allowed` error to the collection
maintainer because endpoint path restrictions may need administrator changes.

## Read with pyathena

[pyathena](https://github.com/jeonggyukim/pyathena) is recommended for this
native layout. Point `LoadSim` at the model directory; pyathena automatically
detects the available snapshot numbers from `id0`:

```python
import pyathena as pa

model_dir = "/path/to/downloaded/R8_2pc"
sim = pa.LoadSim(model_dir)

print(sim.nums_vtk)
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

## Citation

For MHD data, cite Kim & Ostriker (2017),
[ApJ 846, 133](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract),
and the public release DOI,
[10.34770/ackh-7y71](https://doi.org/10.34770/ackh-7y71).

See the public documentation for citations associated with its separate
chemistry, CO-line, and ionizing-radiation products.
