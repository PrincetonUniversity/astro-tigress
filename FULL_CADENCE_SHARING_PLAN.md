# Full-cadence TIGRESS data sharing plan

Status: proposed

Last reviewed: 2026-08-28

## Decision summary

Create a versioned, raw-only export tree on the same `/tigerdata` filesystem as
the simulation outputs. Populate it from an explicitly reviewed manifest using
**hard links**, then expose only that tree through a read-only Globus guest
collection.

This is the best near-term option because it:

- does not duplicate the data blocks;
- does not require moving or renaming the working simulation trees;
- exposes a simple root that can safely receive a read permission;
- selects individual files, which is necessary because raw and personal-analysis
  files are mixed within some source directories; and
- gives collaborators a stable, documented data contract independent of the
  owner's working layout.

The DOI-backed public release at
`/projects/EOSTRIKE/TIGRESS_data_release` remains the canonical curated release.
The new collection should be described as a collaborator-facing, native-format,
full-cadence supplement, not as a replacement or silent extension of the public
release.

## Goals and boundaries

The shared dataset should:

- provide every approved raw snapshot for `R8_2pc` and `R8_4pc`;
- include the minimum run metadata needed to interpret those snapshots;
- exclude plots, pickles, FITS analysis products, scripts, logs, job files,
  executables, restart variants, and other personal products by default;
- be read-only to collaborators;
- preserve native filenames and processor-directory structure; and
- be reproducible from a reviewed, version-controlled selection manifest.

This plan does not reorganize the DOI data release, promise that the native
full-cadence tree has the same layout as that release, or grant access to all
files simply because they live in a directory with raw output.

## Current layout and risk

The source trees are:

| Shared model name | Native run name | Source |
|---|---|---|
| `R8_2pc` | `R8_2pc_rst` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-R8/R8_2pc_rst` |
| `R8_4pc` | `R8_4pc_newacc` | `/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-local/R8_4pc_newacc` |

An inventory on 2026-08-28 found many clearly derived top-level directories. It
also found derived files inside directories that otherwise contain raw data:

- `R8_2pc_rst/id0` contains raw VTK snapshots and at least one derived FITS file;
- `starpar` contains `*.starpar.vtk` files and derived PNG files; and
- `hst` contains plain history/event files and derived pickle files.

The `R8_2pc` native output has processor directories `id0` through `id55`; the
`R8_4pc` output has `id0` through `id27`. Both source roots and
`/tigerdata/EOSTRIKE` reported filesystem device ID 47, so a hard-link export can
live elsewhere under that same filesystem without duplicating file contents.
The final export location must be checked again before building it.

A direct read permission on either source root is therefore unsafe. Globus guest
collection permissions apply to directories and are additive: a broad read rule
cannot be reduced by adding a narrower rule below it. Globus transfer
`--include`/`--exclude` filters control what a particular transfer selects; they
are not access controls and must not be used to hide personal files.

## Proposed shared layout

Ask the storage or Globus administrator to create an owner-writable location on
the same filesystem, for example:

```text
/tigerdata/EOSTRIKE/TIGRESS-share/full-cadence/v1/
├── README.md
├── SCOPE.md
├── R8_2pc/
│   ├── MANIFEST.tsv
│   ├── R8_2pc_rst.par
│   ├── hst/
│   ├── starpar/
│   ├── id0/
│   ├── ...
│   └── id55/
└── R8_4pc/
    ├── MANIFEST.tsv
    ├── R8_4pc_newacc.par
    ├── athinput.R8
    ├── hst/
    ├── starpar/
    ├── id0/
    ├── ...
    └── id27/
```

The model directory names match this repository's public model IDs. Native run
filenames remain unchanged so their provenance is unambiguous. Do not add a
`current` symlink; create a separate collection or explicit version path when a
new selection is released.

`README.md` should contain the public-release DOI and documentation link, model
and run-name mapping, native-layout loading example, cadence/snapshot range,
creation date, contact, citation requirements, and a clear statement that the
collection is read-only and collaborator-facing. `SCOPE.md` should describe the
approved file classes and exclusions. Each manifest should record at least:

```text
relative_path<TAB>size_bytes<TAB>mtime_utc
```

Checksums may be added incrementally if useful. Globus normally verifies
transfer integrity with checksums, so computing a second complete checksum
catalog need not delay the first share.

## Raw-data contract

Use an allowlist, not an exclusion list. The initial candidate scope is:

1. Main MHD VTK output for every approved output number in each `idN`
   directory.
2. Matching raw `*.starpar.vtk` files.
3. Plain-text run history and event files needed for time conversion, such as
   the canonical `.hst` and `.sn` files.
4. The exact input/parameter files needed to interpret the run, including the
   canonical `.par` file and, after review, `athinput.R8` for `R8_4pc`.

The following are out of scope unless a later manifest explicitly adds them:

- restart/checkpoint data (`rst` and its refine, split, or degrade variants);
- pickles, FITS files, PNGs, figures, maps, projections, slices, PDFs, and
  post-processed diagnostics;
- job scripts, scheduler output/error logs, executables, and build products;
- merged or transformed data products; and
- chemistry, line-emission, and radiation post-processing already governed by
  the curated public release or another separately reviewed release.

The filename rules must be anchored to the complete native run name and an
approved output-number range. For example, a rule should express “the canonical
four-digit MHD VTK for run `R8_2pc_rst` in the matching `idN` directory,” not
the broad pattern `*.vtk`. Before implementation, review an automatically
generated candidate list and explicitly decide whether restart files and any
auxiliary VTK components are scientifically necessary.

## Build and publication procedure

### 1. Freeze the selection

- Agree on the included output-number range for each model.
- Confirm that an output is complete in every expected `idN` directory and has
  the corresponding star-particle output when applicable.
- Generate candidate manifests from strict filename rules.
- Review unexpected filenames and gaps manually; never broaden a rule merely to
  make the counts agree.
- Commit the selection rules and reviewed manifests to this repository before
  publishing the collection. The manifests should use relative paths and should
  not contain secrets or unrelated internal paths.

### 2. Build outside the exposed collection

- Create a uniquely named staging directory on device ID 47.
- Create normal directories, then hard-link each approved regular source file
  into its manifest path.
- Reject symbolic links, non-regular files, paths outside the two approved source
  roots, duplicate destination paths, and files not present in the manifest.
- Do not run `chmod` or `chown` on a hard-linked file: hard links share the same
  inode, so this would also change the source file. Directory modes in the export
  tree may be set independently.
- Write the README, scope statement, manifests, and a machine-readable release
  metadata file.

Build a new version in staging and rename it into place only after validation.
Do not mutate a published manifest in place.

### 3. Validate the export

For each model, require all of the following:

- every exported data file appears exactly once in its reviewed manifest;
- source and export have the same filesystem device and inode numbers;
- file count and total byte count match the manifest;
- no symlink, socket, device, FIFO, or unexpected suffix is present;
- every approved output number is represented in every expected `idN`;
- a representative early, middle, and late snapshot can be loaded from `id0`
  with all expected grids and fields;
- a collaborator can transfer a small test snapshot and Globus checksum
  verification succeeds; and
- a separate test identity cannot browse any source or personal-analysis path.

The native parallel VTK loading command should be validated and then copied
verbatim into the collection README. The public-release `astro_tigress.Model`
layout is snapshot-oriented and should not be claimed to work against the native
processor-oriented export without this test.

### 4. Configure Globus

- Create a guest collection rooted at the validated version directory, not at a
  source simulation directory or at `/tigerdata/EOSTRIKE`.
- Grant `r` permission to a dedicated Globus group for collaborators.
- Do not grant `rw`, anonymous access, or an access-manager/administrator role to
  collaborators. Globus management roles can carry implicit read-write access.
- Keep the owner permission limited to the collection maintainer and use the
  group for membership changes.
- Test listing and transferring with a real non-owner identity before sending
  invitations.
- Record the collection UUID, version, group, maintainer, creation date, and
  review date in the release metadata and this plan.

### 5. Operate it as an immutable release

Raw source files that have been linked into a published version should be
treated as immutable. In-place edits are immediately visible through every hard
link; atomic replacement of a source pathname is not visible through an old
hard link. If the selection or a file changes, build and validate `v2` rather
than silently altering `v1`.

Removing a hard link from the export does not remove the source file. Conversely,
deleting the source pathname does not free its data while an export hard link
remains. Before retiring a version, revoke its Globus permission, retain its
manifest according to project policy, remove only the reviewed export tree, and
check link counts and storage usage afterward.

## Longer-term organization

New simulation output and new analysis should never be written into the same
namespace. For future runs, use sibling trees such as:

```text
RUN_ROOT/
├── raw/                 # simulation-owned, immutable after completion
├── analysis/
│   ├── <person>/
│   └── <project>/
└── admin/               # logs, job files, executable, restart operations
```

If reorganizing these legacy runs later becomes acceptable, move analysis files
within the same filesystem into an unshared sibling tree. That is a metadata
operation rather than a data copy, but it can break existing scripts and should
be driven by a reviewed move manifest. Once each raw directory is genuinely
clean, direct directory sharing or read-only bind mounts become simpler than a
hard-link export.

## Alternatives considered

| Approach | Assessment |
|---|---|
| Direct guest collection on a source root | Reject: broad read access exposes intermingled personal files. |
| Globus include/exclude transfer filters | Reject as a security control: users choose their own transfers. |
| Symlink farm | Avoid: Globus path restrictions commonly block links leaving the allowed namespace; recursive transfers do not reliably preserve/follow all symlinks. |
| Hard-link export from a reviewed manifest | Recommended now: exact file-level selection, stable layout, no duplicated data blocks. |
| Move all analysis out of the source trees | Good long-term design, but disruptive to existing personal analysis paths. |
| Administrator-managed mapped collection with GCS path restrictions | Good institutional alternative if Princeton Research Computing will own and test a strict filename/path allowlist. It removes the hard-link namespace but requires endpoint-level administration and careful policy validation. |
| Bind-mounted export | Useful only after raw files are cleanly separated by directory; the present file-level mixing would otherwise require an impractical number of mounts. |

## Implementation checklist

- [ ] Decide snapshot ranges and whether star particles are required.
- [ ] Decide whether canonical history/event and input files may be shared.
- [ ] Confirm restart/checkpoint files are excluded.
- [ ] Ask the storage/Globus administrator to approve and create the export root.
- [ ] Write strict inventory and hard-link build scripts with dry-run output.
- [ ] Review and commit manifests before publication.
- [ ] Build and validate `v1` in staging.
- [ ] Add and test the native-layout loading example.
- [ ] Create the read-only Globus guest collection and collaborator group.
- [ ] Test with a non-owner identity.
- [ ] Record the collection UUID and send the collection README to collaborators.
- [ ] Schedule an annual access and content review.

## Globus references

- [How to share data using a guest collection](https://docs.globus.org/how-to/share-files/)
- [Guest collection permission behavior](https://docs.globus.org/api/transfer/permissions/)
- [Globus Connect Server data access and path restrictions](https://docs.globus.org/globus-connect-server/v5/data-access-guide/)
- [Globus handling of symlinks](https://docs.globus.org/faq/transfer-sharing/)
- [Globus CLI transfer filters and checksum options](https://docs.globus.org/cli/reference/transfer/)
