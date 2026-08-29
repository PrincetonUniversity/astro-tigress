# Full-cadence TIGRESS data sharing plan

Status: proposed

Last reviewed: 2026-08-28

## Decisions

Build the collaborator-facing data at:

```text
/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share
```

Include:

- all currently available raw VTK outputs for `R8_2pc` (`0285-0448`);
- raw VTK outputs `0200-0700` for `R8_4pc`, adding only snapshots that exist
  and pass validation (`0200-0674` are currently available);
- the matching `*.starpar.vtk` for every shared snapshot;
- the canonical `.par`, `.hst`, and `.sn` files; and
- documentation and manifests, but no personal analysis products.

The canonical VTK representation will be a **self-contained directory per
snapshot**, populated with hard links to the native processor files. Version 1
will not store tar archives.

This lets a collaborator select one snapshot directory in Globus without
duplicating VTK data blocks. The collection will be read-only to a dedicated
Globus group. The DOI-backed release at
`/projects/EOSTRIKE/TIGRESS_data_release` remains the canonical public release;
this collection is a raw, full-domain, high-cadence supplement for
collaborators.

## VTK layout decision

[pyathena](https://github.com/jeonggyukim/pyathena) supports uncompressed
per-snapshot archives named `vtk/<problem_id>.<NNNN>.tar`. This is convenient,
but building those archives writes a second complete copy of the VTK bytes.

| Selection | Files | Approximate VTK size |
|---|---:|---:|
| One `R8_2pc` snapshot (`0300`) | 56 | 42.00 GiB |
| All current `R8_2pc` snapshots (`0285-0448`) | 9,184 | 6.73 TiB |
| One `R8_4pc` snapshot (`0300`) | 28 | 5.25 GiB |
| Current `R8_4pc` selection (`0200-0674`) | 13,300 | 2.44 TiB |

A canonical tar tier would require roughly another 9.17 TiB for the current
selection. Compression is not a good default for dense binary VTK and would
depart from pyathena's expected simple `.tar` naming.

The ordinary pyathena reader already loads an `id0` VTK and discovers matching
files in sibling `idN` directories. Each exported snapshot will therefore look
like a small, self-contained native run. This provides direct pyathena support,
one-folder Globus selection, and no duplicate data blocks.

Tar remains an optional downstream format. A collaborator may create a tar
after transfer. Add a server-side tar tier only if later usage justifies its
explicit storage and maintenance cost.

This plan checked pyathena revision
`08b509f1375cf41a4da2e7d735877bb470c6dd08` (2026-07-22). Revalidate the loading
example against the installed version before publication.

## Source inventory

| Shared model | Native run | Source | Processor directories |
|---|---|---|---:|
| `R8_2pc` | `R8_2pc_rst` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-R8/R8_2pc_rst` | `id0-id55` |
| `R8_4pc` | `R8_4pc_newacc` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-local/R8_4pc_newacc` | `id0-id27` |

The strict `id0` inventory is contiguous:

| Model | Raw VTK range | VTK count | Matching star-particle count |
|---|---:|---:|---:|
| `R8_2pc` | `0285-0448` | 164 | 164 |
| `R8_4pc` | `0000-0674` | 675 | 675 |

The `R8_4pc` sharing policy is `0200-0700`, but `0675-0700` do not currently
exist. Publish `0200-0674` first. Add later outputs only in a new validated
content version; never create empty snapshot directories or undocumented gaps.

Both sources and the target parent currently report filesystem device ID 47,
which permits hard links. Verify this at build time. The target directory does
not yet exist.

## Allowlisted data

For every approved snapshot `NNNN`, include exactly:

- `id0/<native_run>.NNNN.vtk`;
- every expected `idN/<native_run>-idN.NNNN.vtk`;
- `starpar/<native_run>.NNNN.starpar.vtk`;
- `<native_run>.par`;
- `hst/<native_run>.hst`; and
- `hst/<native_run>.sn`.

Hard-link `.par`, `.hst`, and `.sn` into every snapshot directory. These
repeated directory entries make each download self-contained for pyathena but
refer to the same source inodes.

Exclude restart/checkpoint data; FITS, pickle, PNG, PDF, map, projection, slice,
and phase products; scheduler output, scripts, executables, CSV products, and
logs; merged or auxiliary VTK products; and public-release chemistry, CO-line,
or radiation post-processing.

This must be an allowlist. The source roots contain analysis directories, and
raw-looking directories are mixed too: `id0` contains a derived FITS file,
`starpar` contains PNGs, and `hst` contains pickles. Sharing a native root would
expose personal files.

## Shared tree

```text
/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share/
├── README.md
├── RELEASE.json
├── R8_2pc/
│   ├── MANIFEST.tsv
│   └── snapshots/
│       ├── 0285/
│       │   ├── R8_2pc_rst.par
│       │   ├── hst/
│       │   │   ├── R8_2pc_rst.hst
│       │   │   └── R8_2pc_rst.sn
│       │   ├── starpar/R8_2pc_rst.0285.starpar.vtk
│       │   ├── id0/R8_2pc_rst.0285.vtk
│       │   ├── id1/R8_2pc_rst-id1.0285.vtk
│       │   ├── ...
│       │   └── id55/R8_2pc_rst-id55.0285.vtk
│       ├── 0286/
│       └── ...
└── R8_4pc/
    ├── MANIFEST.tsv
    └── snapshots/
        ├── 0200/
        │   ├── R8_4pc_newacc.par
        │   ├── hst/
        │   │   ├── R8_4pc_newacc.hst
        │   │   └── R8_4pc_newacc.sn
        │   ├── starpar/R8_4pc_newacc.0200.starpar.vtk
        │   ├── id0/R8_4pc_newacc.0200.vtk
        │   ├── id1/R8_4pc_newacc-id1.0200.vtk
        │   ├── ...
        │   └── id27/R8_4pc_newacc-id27.0200.vtk
        ├── 0201/
        └── ...
```

`FULL_CADENCE_SHARE_README.md` is the share-root README draft. Install it as
`README.md` after filling in the collection UUID, contact, release date, and
final inventory.

## Domain and public-release comparison

Both full-domain models cover:

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

The collection README must distinguish:

- public: `MODEL/NNNN/MHD/<model>.NNNN.vtk`, one joined central-cube file;
- shared: `MODEL/snapshots/NNNN/idN/...`, multiple full-domain processor files;
- the public readers here assume the public snapshot/product tree; and
- pyathena is recommended for the shared native tree.

## Build procedure

### 1. Generate and review manifests

Generate strict manifests containing at least:

```text
export_path<TAB>source_relative_path<TAB>size_bytes<TAB>mtime_utc
```

Reject a snapshot without exactly 56 (`R8_2pc`) or 28 (`R8_4pc`) processor
files; a missing or duplicate star-particle file; output numbers outside the
policy; symlinks and non-regular files; paths outside the two approved roots;
duplicate destinations; and filenames that match only a broad suffix such as
`*.vtk` rather than the complete native-run pattern.

Review and commit the rules and manifests before publication. Published
manifests should not contain unrelated internal absolute paths.

### 2. Build in staging

- Create a uniquely named staging tree beside the target.
- Create ordinary snapshot, `idN`, `starpar`, and `hst` directories.
- Hard-link only files in the reviewed manifest.
- Add README and release metadata as ordinary small files.
- Do not `chmod` or `chown` a hard-linked file; that changes the source inode.
- Rename staging into place only after validation.

Treat a published version as immutable. In-place source edits are visible
through hard links, while atomic source replacements are not reflected in old
links. Rebuild a new version rather than silently changing published content.

### 3. Validate

For every snapshot:

- match file count and total bytes to the manifest;
- confirm export/source device and inode numbers are identical;
- confirm the continuous `idN` range and matching output number;
- confirm the star-particle output number matches;
- reject unexpected file types and suffixes; and
- load gas and star particles with pyathena.

Test early, middle, and late snapshots for both models with the exact README
example. Confirm domain edges, grid count, resolution, fields, simulation time,
and star-particle time. Transfer one snapshot through Globus with checksum
verification. A non-owner test identity must not browse native or analysis
paths.

## Globus configuration and lifecycle

- Root the guest collection at the clean share root, never at a native run or a
  broader `TIGRESS-classic` directory.
- Grant `r` to a dedicated Globus group; do not grant collaborators `rw`,
  anonymous access, or management roles.
- Record collection UUID, group, maintainer, content version, build time,
  source ranges, and validation result in `RELEASE.json`.

Globus guest permissions are directory-based and additive. Transfer filters
help users select files but are not access controls. The clean export root makes
a root-level read permission safe.

When `R8_4pc` outputs `0675-0700` appear, run the full manifest and validation
workflow and publish a new version. Do not expose partially written processor
sets.

Removing an export hard link does not remove the source. Deleting the source
pathname does not free data while an export link remains. To retire a version,
revoke its Globus permission, retain its manifest, remove only the reviewed
export tree, and inspect link counts and storage afterward.

## Checklist

- [x] Choose target path and model ranges.
- [x] Include matching star particles, `.par`, `.hst`, and `.sn`.
- [x] Choose hard-linked snapshot directories over canonical tar.
- [x] Record full and public-release extents.
- [ ] Write inventory/build/validation scripts with dry-run support.
- [ ] Generate and review manifests.
- [ ] Create target and staging trees.
- [ ] Build and validate the first content version.
- [ ] Finalize and install the collection README.
- [ ] Create the read-only Globus collection and group.
- [ ] Test with a non-owner identity and record the collection UUID.

## References

- [pyathena repository](https://github.com/jeonggyukim/pyathena)
- [pyathena documentation](https://jeonggyukim.github.io/pyathena/intro.html)
- [Sharing with a Globus guest collection](https://docs.globus.org/how-to/share-files/)
- [Globus guest permission behavior](https://docs.globus.org/api/transfer/permissions/)
- [Globus handling of symlinks](https://docs.globus.org/faq/transfer-sharing/)
- [Globus transfer filters and checksums](https://docs.globus.org/cli/reference/transfer/)
