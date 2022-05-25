# TIGRESS data release

Welcome to the TIGRESS data release! This is a seriew of python scripts that will help to
to read and analyse the data.

See the **wiki pages** for a detailed guide.

Currently, this data release contains the solar neighborhood TIGRESS simulations. The
outputs are selected during the time frame when the simulation has reached steady state.
The released snapshots has a time interval of about 10 Myr. 

We the original TIGRESS MHD output data ([Kim & Ostriker (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract)),
as well as the chemistry and CO line post-processing data ([Gong et al
2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract), [Gong et al.
2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract)).

There are two sets of simulation models:

**R8_2pc**: solar nieghborhood environment with 2pc resolution. This set contains MHD,
chemistry and CO line data.

**R8_4pc**: solar neighborhood environment with 4pc resolution. This set only contains
the MHD data. We do not include the chemistry and CO line data, as the molecular clouds
are not well-resolved in this case ([Gong et al
2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract)).

