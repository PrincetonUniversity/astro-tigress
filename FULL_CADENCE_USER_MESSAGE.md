# Full-cadence TIGRESS data now available

Hello everyone,

We have made the full-domain, high-cadence TIGRESS MHD outputs available through a collaborator Globus collection.

Please note that this collection differs from the public data release. The public release contains selected snapshots, joined central-cube MHD files, and curated post-processing products. The collaborator collection instead contains raw MHD outputs over the full vertical domain in the native parallel `idN` layout, together with the corresponding star-particle files and run metadata. The available snapshot ranges are `0285`–`0448` for `R8_2pc` and `0200`–`0650` for `R8_4pc`.

Before downloading, please read the `README.md` at the root of the Globus collection. It describes the directory structure and domain extent, provides examples for analysis with [pyathena](https://github.com/jeonggyukim/pyathena), and explains how to use a Globus batch list to download selected snapshots. When selecting a snapshot, include its files from every `idN` directory as well as the matching `starpar.vtk` file. Because the collection is assembled using file links, please select files explicitly or use the documented batch-list method rather than recursively transferring entire directories.

Please let us know if you encounter a `Path not allowed` error or have questions about the data layout.
