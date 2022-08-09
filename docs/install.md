# Installation

Installation of basic reader and analysis scripts:
```
git clone https://github.com/PrincetonUniversity/astro-tigress.git
cd astro-tigress
conda create --name astro-tigress python=3.8 poetry jupyterlab xarray
conda activate astro-tigress
poetry install
```

Quick data download for tutorials:
```
python prepare_tutorial.py
```
