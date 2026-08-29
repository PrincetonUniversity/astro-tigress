# Full-cadence TIGRESS data sharing plan

Status: approved for initial build

Last reviewed: 2026-08-29

## Decisions

Build the collaborator-facing data at:

```text
/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share
```

Include:

- all raw VTK outputs currently available for `R8_2pc` (`0285-0448`);
- raw VTK outputs `0200-0650` inclusive for `R8_4pc`;
- the matching `*.starpar.vtk` for every shared snapshot; and
- one canonical `.par`, `.hst`, and `.sn` file per model.

Preserve the native processor-oriented layout. Do not rearrange VTK files into
per-snapshot directories and do not create tar archives. Populate the clean
share tree with file symlinks selected by a strict allowlist, then expose only
the share tree through a read-only Globus guest collection.

This choice:

- lets pyathena automatically discover all available snapshot numbers;
- stores `.par`, `.hst`, and `.sn` only once per model;
- does not duplicate any simulation or star-particle data blocks;
- keeps personal analysis files out of the shared namespace; and
- is reversible by unlinking the generated symlinks, without changing the
  source pathnames or data.

The DOI-backed release at `/projects/EOSTRIKE/TIGRESS_data_release` remains the
canonical public release. This collection is a raw, full-domain, high-cadence
supplement for collaborators.

## Source inventory

| Shared model | Native run | Source | Shared VTK range | Processor directories |
|---|---|---|---:|---:|
| `R8_2pc` | `R8_2pc_rst` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-R8/R8_2pc_rst` | `0285-0448` | `id0-id55` |
| `R8_4pc` | `R8_4pc_newacc` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-local/R8_4pc_newacc` | `0200-0650` | `id0-id27` |

The source inventories are contiguous over both selected ranges, and every
selected `id0` output has a matching star-particle file. The builder must also
verify every nonzero processor file before creating the public tree.

Although the sources and target parent report the same client filesystem device
ID, the NFS server returns `EIO` for hard links, even within a source directory.
The no-copy user-level implementation must therefore use symlinks. If the
Globus endpoint rejects those links, an administrator-provided namespace view
is required; copying is not an approved fallback.

## Allowlisted data

For every approved snapshot `NNNN`, include exactly:

- `id0/<native_run>.NNNN.vtk`;
- every expected `idN/<native_run>-idN.NNNN.vtk`; and
- `starpar/<native_run>.NNNN.starpar.vtk`.

Include once per model:

- `<native_run>.par`;
- `hst/<native_run>.hst`; and
- `hst/<native_run>.sn`.

Exclude restart/checkpoint data; FITS, pickle, PNG, PDF, map, projection, slice,
and phase products; scheduler output, scripts, executables, CSV products, and
logs; merged or auxiliary VTK products; and public-release chemistry, CO-line,
or radiation post-processing.

This must be an allowlist. The native roots contain analysis directories, and
raw-looking directories are mixed too: `id0` contains a derived FITS file,
`starpar` contains PNGs, and `hst` contains pickles. Sharing a source root would
expose personal files.

## Shared tree

```text
/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share/
├── README.md
├── RELEASE.json
├── R8_2pc/
│   ├── MANIFEST.tsv
│   ├── R8_2pc_rst.par
│   ├── hst/
│   │   ├── R8_2pc_rst.hst
│   │   └── R8_2pc_rst.sn
│   ├── starpar/
│   │   ├── R8_2pc_rst.0285.starpar.vtk
│   │   └── ...
│   ├── id0/
│   │   ├── R8_2pc_rst.0285.vtk
│   │   └── ...
│   ├── id1/
│   │   ├── R8_2pc_rst-id1.0285.vtk
│   │   └── ...
│   └── ... id55/
└── R8_4pc/
    ├── MANIFEST.tsv
    ├── R8_4pc_newacc.par
    ├── hst/
    │   ├── R8_4pc_newacc.hst
    │   └── R8_4pc_newacc.sn
    ├── starpar/
    │   ├── R8_4pc_newacc.0200.starpar.vtk
    │   └── ... through 0650
    ├── id0/
    │   ├── R8_4pc_newacc.0200.vtk
    │   └── ... through 0650
    └── ... id27/
```

`FULL_CADENCE_SHARE_README.md` is copied to the root as `README.md` during the
build. Each `MANIFEST.tsv` records the shared relative path, size, UTC
modification time, device, and inode of the symlink target.

## VTK layout and pyathena

[pyathena](https://github.com/jeonggyukim/pyathena) scans `id0` for the primary
VTK outputs, combines those numbers with other supported VTK representations,
and loads matching files from sibling `idN` directories. Pointing `LoadSim` at
the model root therefore exposes the complete `nums_vtk` list automatically:

```python
import pyathena as pa

sim = pa.LoadSim("/path/to/TIGRESS-full-data-share/R8_2pc")
print(sim.nums_vtk)
ds = sim.load_vtk(num=300, id0=True, load_method="xarray")
starpar = sim.load_starpar_vtk(num=300)
```

The behavior was reviewed against pyathena revision
`08b509f1375cf41a4da2e7d735877bb470c6dd08` (2026-07-22). A runtime smoke test
with the installed collaborator environment remains a publication gate.

Users who want only selected snapshots must select the matching filename in
each `idN` directory plus the matching `starpar` file while preserving the
directory structure. With 56 or 28 processor directories, this is manageable
manually or with a Globus batch list. The collection does not store an easier
second representation because per-snapshot tar archives would duplicate about
9 TiB for the previously considered range.

Globus follows a symlink when an individual file is transferred and writes an
ordinary file at the destination. Recursive transfers generally skip symlinks,
and collection path restrictions can reject links whose targets are outside the
allowed namespace. Test an individual VTK link and a complete batch transfer on
the actual guest collection before inviting collaborators. If the test fails,
ask the endpoint administrator for a permitted namespace/bind view or an
appropriate server policy; do not broaden access to the native source roots.

## Domain and public-release comparison

Both shared full-domain models cover:

```text
x = [-512, 512] pc
y = [-512, 512] pc
z = [-3584, 3584] pc
```

| Model | Cell size | Full-domain cells | Processor pieces |
|---|---:|---:|---:|
| `R8_2pc` | 2 pc | `512 x 512 x 3584` | 56 |
| `R8_4pc` | 4 pc | `256 x 256 x 1792` | 28 |

The public release's joined `MHD` files are central cubes with
`x,y,z = [-512,512] pc`: `512^3` cells for `R8_2pc` and `256^3` for `R8_4pc`.
The shared data keep the same horizontal extent and resolution but include the
full vertical domain. The public release has selected snapshots at roughly
10 Myr intervals and separately curated post-processing products. This
collection has the denser native MHD sequence, star particles, and run metadata
only.

The README distinguishes:

- public: `MODEL/NNNN/MHD/<model>.NNNN.vtk`, one joined central-cube file;
- shared: `MODEL/idN/<native-run>[-idN].NNNN.vtk`, full-domain processor files;
- the public readers here assume the public snapshot/product tree; and
- pyathena is recommended for the shared native tree.

## Build and validation procedure

Use `scripts/build_full_data_share.py`. Its default mode performs source
inventory checks and prints the intended build without writing. `--execute`
creates a uniquely named sibling staging directory, symlinks only exact
allowlisted files, validates it, and atomically renames it to the target. It
refuses to overwrite an existing target. `--validate-only` audits the published
tree against the sources.

The builder must reject:

- a snapshot without exactly 56 (`R8_2pc`) or 28 (`R8_4pc`) VTK pieces;
- a missing or duplicate matching star-particle file;
- output numbers outside the approved ranges;
- symbolic links, non-regular files, paths outside the approved sources, and
  duplicate destinations; and
- any source symlink or non-regular source file.

Validation requires every export symlink to have the exact expected absolute
target, every target to be a regular file, manifests to contain the expected
counts, and no unexpected pathname to exist.

After the filesystem audit, load an early, middle, and late snapshot for both
models with the README example. Confirm `nums_vtk`, domain edges, grid count,
resolution, fields, simulation time, and star-particle time. Then transfer one
snapshot through Globus with checksum verification and test access with a
non-owner identity.

## Globus and reversibility

- Root the guest collection at the clean share root, never at a native run or a
  broader `TIGRESS-classic` directory.
- Grant `r` to a dedicated Globus group; do not grant collaborators `rw`,
  anonymous access, or management roles.
- Record collection UUID, group, maintainer, content version, build time,
  ranges, and validation result in `RELEASE.json`.

Globus guest permissions are directory-based and additive. Transfer filters
help users select files but are not access controls. The clean export root makes
a root-level read permission safe.

The share is reversible because all large data entries are symlinks. Removing
the generated share directory unlinks only those references; it does not remove
the target files. Deletion must still be performed only after revoking
Globus access, comparing the target to the committed manifest, and confirming
the exact target path. The build script intentionally provides no removal
command.

## Checklist

- [x] Isolate changes on the `full-data-share` branch.
- [x] Choose target, model ranges, native layout, and allowlist.
- [x] Include matching star particles and one `.par`, `.hst`, `.sn` per model.
- [x] Record full and public-release extents.
- [x] Run the builder dry run.
- [x] Create and validate the symlink share tree.
- [ ] Run pyathena smoke tests in a complete environment.
- [ ] Create the read-only Globus collection and group.
- [ ] Test with a non-owner identity and record the collection UUID.

## References

- [pyathena repository](https://github.com/jeonggyukim/pyathena)
- [pyathena documentation](https://jeonggyukim.github.io/pyathena/intro.html)
- [Sharing with a Globus guest collection](https://docs.globus.org/how-to/share-files/)
- [Globus guest permission behavior](https://docs.globus.org/api/transfer/permissions/)
- [Globus transfer filters and checksums](https://docs.globus.org/cli/reference/transfer/)
