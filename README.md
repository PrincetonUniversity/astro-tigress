# TIGRESS data release

Welcome to the TIGRESS data release! This repo contains a series of python scripts and example jupyter notebooks that will help you read and analyse the TIGRESS simulation data.

See the **wiki pages** for a detailed guide.

Currently, this data release contains the solar neighborhood TIGRESS simulations. The outputs are for a time period when the simulation has reached a steady state. The snapshots are separated by a time interval of about 10 Myr.

We the original TIGRESS MHD output data ([Kim & Ostriker (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...846..133K/abstract)), as well as the chemistry and CO line post-processing data ([Gong et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract), [Gong et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...903..142G/abstract)).

There are two sets of simulation models:

**R8_2pc**: Solar neighborhood environment with 2pc resolution. This set contains MHD, chemistry and CO line data.

**R8_4pc**: Solar neighborhood environment with 4pc resolution. This set only contains the MHD data. We do not include the chemistry and CO line data, as the molecular clouds are not well-resolved in this case ([Gong et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...858...16G/abstract)).
