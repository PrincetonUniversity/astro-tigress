# TIGRESS Data Release

## Data Location

Currently, the data is located at [this webpage](https://tigress-web.princeton.edu/~munan/astro-tigress/). You can browse and download individual files. However, you will need to preserve the file structure (folders and sub-folders) for the scripts provided here to work properly (you can also modify the file names in the scripts if you don't want to use the same file structure).

## Installation

Clone the repository.
```
git clone git@github.com:PrincetonUniversity/astro-tigress.git
```

We use `poetry` for packaging and dependency handling. `poetry` can be installed with conda

```
conda install poetry
```

or with `pip`

```
pip install -U poetrry
```

Then, this will install all requirements.

```
poetry install
```

Otherwise, you can manually install following packages.

## Requirements
### Python 3
The scripts for reading the outputs are written in Python 3. I manage my python installation, versions, and packages with [Anaconda](https://docs.anaconda.com/anaconda/install/).

### yt
`yt` can be used to read the MHD simulation outputs from [Athena++](https://princetonuniversity.github.io/athena/). It allows part of the simulation domain to be accessed without loading the whole simulation output in the memory, which improves the time and memory efficiency of reading the output data. Since the molecular gas is largely located at the galactic mid-plane, we only need to access the MHD simulation output near the mid-plane region.

You can install yt here: https://yt-project.org/#getyt

I recommend installing yt using [Anaconda](https://docs.anaconda.com/anaconda/install/). This will guarantee that all the packages (such as `yt`, `numpy`, `scipy` etc) work together. You may still need to install `h5py`.

### Astropy
The [Astropy](https://www.astropy.org/) package is needed for converting the sythetic CO emission line maps made by [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) to FITS files. Astropy can also be installed using [Anaconda](https://docs.anaconda.com/anaconda/install/).

## Quick start
Start with playing the jupyter notebooks **analysis/**, which contains examples for reading and analysing the data.

* [read_data_1-MHD.ipynb](https://github.com/PrincetonUniversity/astro-tigress/blob/master/analysis/read_data_1-MHD.ipynb) shows you how to work with the original TIGRESS MHD output data ([Kim & Ostriker (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract))
* [read_data_2-chem-CO_lines.ipynb](https://github.com/PrincetonUniversity/astro-tigress/blob/master/analysis/read_data_2-chem-CO_lines.ipynb) gives examples working with the chemistry and CO line post-processing data ([Gong et al 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract), [Gong et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract)).

## File structure
There are three different folders:

* [analysis](https://github.com/PrincetonUniversity/astro-tigress/tree/master/analysis): Jupyter notebooks for data reading and analysis

* [module](https://github.com/PrincetonUniversity/astro-tigress/tree/master/module): Python modules that are used in the Jupyter notebooks

* [miscellaneous](https://github.com/PrincetonUniversity/astro-tigress/tree/master/miscellaneous): Other stuff. [copy_data.ipynb](https://github.com/PrincetonUniversity/astro-tigress/blob/master/miscellaneous/copy_data.ipynb) for example, is used by the developer to move simulation data into the data release storage.
