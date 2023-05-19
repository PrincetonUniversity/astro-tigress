# TIGRESS public data release

[<img src="https://img.shields.io/badge/DOI-10.34770%2Fackh--7y71-blue">](https://doi.org/10.34770/ackh-7y71)
[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://princetonuniversity.github.io/astro-tigress/)

## Data Download

Full data can be downloaded via [Globus](https://app.globus.org/file-manager?origin_id=dc43f461-0ca7-4203-848c-33a9fc00a464&origin_path=%2Fackh-7y71%2F).

A few data files used in the tutorials is also available at [this
webpage](https://tigress-web.princeton.edu/~munan/astro-tigress/).

You can
browse and download whole datasets as well as individual files. However, you will need to preserve the
file structure (folders and sub-folders) for the scripts provided here to work
properly.

Currently, this data release only contains the solar neighborhood models, `R8`,
for selected snapshots with an interval of about 10 Myr after a quasi-steady
state is reached. Full MHD snapshots for every 1 Myr are available upon request.

## Installation of scripts

The python scripts and example notebooks are available at the [GitHub](https://github.com/PrincetonUniversity/astro-tigress) repository.

```sh
git clone git@github.com:PrincetonUniversity/astro-tigress.git
```

For document building, we use `jupyter-book`

```sh
PYTHONPATH=./ jupyter-book build docs
```

## Quickstart
``` py
import astro_tigress
dir_master = ... # master path to your data
model_id = "R8_2pc"
model=astro_tigress.Model(model_id,dir_master)
model.load(300,"MHD")

# make plots using yt
import yt
slc=yt.SlicePlot(model.MHD.ytds,'z',fields=('gas','nH'))
slc.annotate_magnetic_field()
```

## Data types

* `MHD`: Small box (z=+-512 pc) original TIGRESS output for MHD variables including density, pressure, velocity (vector), and magnetic field (vector). Athena VTK format.
* `MHD_PI`: Full box (z=+-3072 pc) MHD data plus electron abundance as a result of UV radiation post-processing. Athena VTK format.
* `chem`: The chemistry post-processing small box output for molecular hydrogen and CO abundances. Athena++ HDF5 format.
* `CO_lines`: The CO lines, J(1-0) and J(2-1). RADMC-3D binary format.

Note that the chemistry post-processing data is only available for `R8-2pc`.

## License and Attribution

If you make use of

* only `MHD` data, please cite
  {cite:t}`2017ApJ...846..133K`
  <!-- `Kim & Ostriker (2017)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/exportcitation)) -->
* `chem` and/or `CO_lines` data, please also cite
  {cite:t}`2018ApJ...858...16G,2020ApJ...903..142G`
  <!-- `Gong et al. (2018)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/exportcitation))
  and `Gong et al. (2020)`
  ([ADS](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract),
   [BibTex](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/exportcitation)) -->
* `MHD_PI`, please also cite
  {cite:t}`2020ApJ...897..143K`
  <!-- ([ADS](https://ui.adsabs.harvard.edu/abs/2020ApJ...897..143K/abstract), -->
   <!-- [BibTex](https://ui.adsabs.harvard.edu/abs/2020ApJ...897..143K/exportcitation)) -->


## Examples

```{tableofcontents}
```

## References

```{bibliography}
```