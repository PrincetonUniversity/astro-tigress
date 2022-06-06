# Installation

## Anaconda + poetry

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
pip install -U poetry
```

Then, this will install all requirements.

```
poetry install
```

Otherwise, you can manually install following packages.

```
conda install -c conda-forge python=3.8 yt h5py astropy
```
