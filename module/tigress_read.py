import os
import numpy as np
import athena_read
import const
import ytathena as ya
import radmc

class Model:
    """Class containing the simulation model information.

    Parameters
    ----------
    model_id: str
        name of the model, e.g. `R8_2pc`.
    """

    def __init__(self, dir_master, model_id):
        self.dir_master = dir_master
        self.model_id = model_id
        self.dir_model = dir_master + model_id + "/"
        subdirs = os.listdir(self.dir_model)
        ivtks = [int(x) for x in subdirs if x.isnumeric()]
        ivtks.sort()
        self.ivtks = np.array(ivtks) #list of vtk numbers of outputs
        #contained data sets, e.g. ["MHD", "chem", "CO_lines"]
        self.data_sets = {}
        for ivtk in self.ivtks:
            dir_ivtk = "{:s}{:04d}/".format(self.dir_model, ivtk)
            self.data_sets[ivtk] = os.listdir(dir_ivtk)
        #history file
        fn_hst = "{:s}history/MHD/{:s}.hst".format(self.dir_model,self.model_id)
        if os.path.isfile(fn_hst):
            self.hst = athena_read.hst(fn_hst)
            self.fn_hst = fn_hst
        #read time output interval in input file
        fn_input_MHD = "{:s}input/MHD/{:s}.par".format(
                          self.dir_model,self.model_id)
        if os.path.isfile(fn_input_MHD):
            self.MHD_input = athena_read.athinput(fn_input_MHD)
            self.fn_input_MHD = fn_input_MHD
        for k in self.MHD_input.keys():
            if k.startswith("output"):
                if self.MHD_input[k]["out_fmt"] == "vtk":
                    self.dt_code = self.MHD_input[k]["dt"]
                    self.dt_myr = self.dt_code * const.pc/const.km/const.myr
                    self.t_myr = self.ivtks * self.dt_myr
        self.CO_lines = {}
        return

    def load(self, ivtk, dataset="MHD", Z=1., iline=1, Tdect=0.):
        """Load simulation data from a certain data set.

        Parameters
        ----------
        ivtk: int
            vtk output number
        dataset: srt, ["MHD", "chem", "CO_lines"]
            name of the dataset.
        Z: float, [0.5, 1, 2]
            metallicity for chemistry post-processing.
        iline: int, [1, 2]
            uppper level of CO rotational line. 1 for J=1-0, 2 for J=2-1
        Tdect: float
            detection limit for atenna temperature in K
            for each velocity channel.
        """
        dir_ivtk = "{:s}{:04d}/".format(self.dir_model, ivtk)
        dict_Z = {0.5:"Z05", 1.:"Z1", 2.:"Z2"}
        Z_id = dict_Z[Z]
        if dataset == "MHD":
            fn_MHD = "{:s}MHD/{:s}.{:04d}.vtk".format(
                      dir_ivtk, self.model_id, ivtk)
            self.MHD = DataMHD(fn_MHD)
        elif dataset == "chem":
            fn_chem = "{:s}chem/{:s}/{:s}-{:s}.{:04d}.athdf".format(
                      dir_ivtk, Z_id, self.model_id, Z_id, ivtk)
            self.chem = DataChem(fn_chem)
        elif dataset == "CO_lines":
            fn_CO = "{:s}CO_lines/{:s}/il{:d}/{:s}-{:s}.il{:d}.{:04d}.bout".format(
                      dir_ivtk, Z_id, iline, self.model_id, Z_id, iline, ivtk)
            self.CO_lines[iline] = DataCO(fn_CO, iline, Tdect)
        else:
            msg = "ERROR: Model.load(): dataset name not recogonozed.\n"
            msg += 'dataset: {"MHD", "chem", "CO_lines"}'
            raise RuntimeError(msg)
        return

class DataMHD:
    """MHD data in one time snapshot.

    Parameters
    ----------
    fn: str
        filename.
    """
    def __init__(self, fn):
        self.fn = fn
        self.ytds = ya.load_athena4p2(fn) #yt data
        self.grid = ya.get_covering_grid(self.ytds) #yt covering grid
        return

class DataChem:
    """Chemistry data in one time snapshot.

    Parameters
    ----------
    fn: str
        filename.
    """
    def __init__(self, fn):
        self.fn = fn
        self.ytds = ya.load_athenapp(fn) #yt data
        self.grid = ya.get_covering_grid(self.ytds) #yt covering grid
        return
    def get_col(self, spec, axis=2):
        """Calculate the column density of species.

        Parameters
        ----------
        spec: str
            name of species, e.g. "H2".
        axis: int, [0, 1, 2]
            axis where the column is calculated along.
            0:x, 1:y, 2:z.

        Returns
        -------
        col: YTarray
            column density, e.g. N_H2.
        """
        dx = self.grid["dx"][0, 0, 0]
        nH = self.grid["nH"]
        abd = self.grid[spec]
        ns = abd*nH
        col = ns.sum(axis=2)*dx
        return col

class DataCO:
    """CO emission lines PPV in one time snapshot.

    Parameters
    ----------
    fn: str
        filename
    iline: int, [1, 2]
        index for CO rotaional line upper level.
    Tdect: float
        detection limit for atenna temperature in K
        for each velocity channel.
    """
    def __init__(self, fn, iline, Tdect=0.):
        if iline == 1:
            self.lam0 = 2600.75763346
        elif iline == 2:
            self.lam0 = 1300.4036558
        else:
            msg = "ERROR: DataCO.__init__(): "
            msg += "Line index {:d} not recognized (iline={1,2}).".format(
                    iline)
            raise RuntimeError(msg)
        self.fn = fn
        self.Tdect = Tdect
        self.img = radmc.Images(fn, lam0=self.lam0)  #PPV cube
        self.WCO = self.img.getimageWCO(Tdect=Tdect) #CO line intensity in K*km/s
        self.Tpeak = self.img.getTpeak(Tdect=Tdect)  #peak atenna temperature in K
        #brightness weighted mean velocity and velocity dispersion in km/s
        self.vz, self.sigmav= self.img.getimage_velocities(Tdect=Tdect)
        return

