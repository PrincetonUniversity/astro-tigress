# TIGRESS public data release [<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/solid/book.svg" width=30 height=30>](https://princetonuniversity.github.io/astro-tigress/)


<!--
[:fontawesome-brands-github:](https://github.com/PrincetonUniversity/astro-tigress/tree/packaging)
[:fontawesome-solid-book:](https://princetonuniversity.github.io/astro-tigress/)
-->

**Welcome to the TIGRESS data release!**

This repo contains a series of python scripts and example jupyter notebooks that will help you read and analyse the TIGRESS simulation data.


## Data Location

Currently, the data is located at [this
webpage](https://tigress-web.princeton.edu/~munan/astro-tigress/). You can
browse and download individual files. However, you will need to preserve the
file structure (folders and sub-folders) for the scripts provided here to work
properly.

Using the internal downloader, each snapshot can be downloaded
``` py
import sys
sys.insert.path(0,"PATH-TO-astro-tigress/module/")

import tigress_read
model=tigress_read.Model("R8_2pc")
model.download(300)
```

Currently, this data release only contains the solar neighborhood models, `R8`,
for selected snapshots with an interval of about 10 Myr after a quasi-steady
state is reached. Full MHD snapshots for every 1 Myr are available upon request.

## Data types

* `MHD`: The original TIGRESS output for MHD variables including density, pressure, velocity (vector), and magnetic field (vector). Athena VTK format.
* `chem`: The chemistry post-processing output for molecular hydrogen and CO abundnaces. Athena++ HDF5 format.
* `CO_lines`: The CO lines, J(1-0) and J(2-1). RADMC-3D binary format.
* (to be added) `rad`: The radiation post-processing output. Athena VTK format.

Note that the chemistry post-processing data is only available for `R8-2pc`.

## Quick Start

Start with playing the jupyter notebooks **analysis/**, which contains examples for reading and analysing the data.

* [read_data_1-MHD.ipynb](https://github.com/PrincetonUniversity/astro-tigress/blob/master/analysis/read_data_1-MHD.ipynb) shows you how to work with the original TIGRESS MHD output data ([Kim & Ostriker 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract))
* [read_data_2-chem-CO_lines.ipynb](https://github.com/PrincetonUniversity/astro-tigress/blob/master/analysis/read_data_2-chem-CO_lines.ipynb) gives examples working with the chemistry and CO line post-processing data ([Gong et al 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract), [Gong et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract)).

## License and Attribution

If you make use of

* only `MHD` data, please cite
  `Kim & Ostriker (2017)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/exportcitation))
* `chem` and/or `CO_lines` data, please also cite
  `Gong et al. (2018)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/exportcitation))
  and `Gong et al. (2020)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/exportcitation))
