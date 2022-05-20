import numpy as np
import os
import ytathena as ya
import radmc

class Data:
    """Class containing the simulation data for a time snapshot."""
    def __init__(self, dir_ath, dir_radmc, model_id, time, 
                 lines, zlim_pc):
        """Initialization.
        Inputs:
            dir_ath: string, directory for athena++ MHD simulation outputs.
            dir_radmc, string, directory for radmc post-processing outputs.
            model_id: string, name of the model. e.g. "R4-Z1".
            time: int, certain time in the simulation in Myr.
            lines: list of integers, CO line(s) to read. e.g. [1,2]
            zlim_pc: float, simulation is loaded at z=[-zlim_pc,zlim_pc], if
                     zlim_pc is None, read the whole simulation.
        """
        self.fn_ath = dir_ath + "{:s}.{:04d}.athdf".format(model_id, time)
        if not os.path.isfile(self.fn_ath):
            msg = "ERROR: Data.__init__(): "
            msg += "Athena++ output file {:s} do not exist.".format(
                    self.fn_ath)
            raise RuntimeError(msg)
        self.fn_radmc = []
        for iline in lines:
            fn_radmc = dir_radmc + "il{:d}/{:s}.{:04d}.il{:d}.bout".format(
                                    iline, model_id, time, iline)
            if not os.path.isfile(fn_radmc):
                msg = "ERROR: Data.__init__(): "
                msg += "Radmc output file {:s} do not exist.".format(fn_radmc)
                raise RuntimeError(msg)
            self.fn_radmc.append(fn_radmc)
        self.model_id = model_id
        self.lines = lines
        self.time = time
        #read MHD data
        self.data_athena = ya.load_athenapp(self.fn_ath)
        self.dx_pc = ( ( self.data_athena.domain_right_edge[0] 
                - self.data_athena.domain_left_edge[0] )
                     /self.data_athena.domain_dimensions[0] ).value
        left_edge = self.data_athena.domain_left_edge.value
        dims = self.data_athena.domain_dimensions
        xleft = left_edge[0]
        yleft = left_edge[1]
        if zlim_pc is None:
            self.zlim_pc = np.abs(left_edge[2])
        else:
            self.zlim_pc = zlim_pc
        left_edge = [xleft, yleft, -self.zlim_pc]
        xdim = dims[0]
        ydim = dims[1]
        dims = [xdim, ydim, int(self.zlim_pc*2/self.dx_pc)]
        self.left_edge = left_edge
        self.dims = dims
        self.grid = self.data_athena.covering_grid(level=0,
                        left_edge=left_edge, dims=dims)
        #column densities
        self.NH = self.get_col("Htot").in_units("cm**-2")
        self.NH2 = self.get_col("H2").in_units("cm**-2")
        self.NCO = self.get_col("CO").in_units("cm**-2")
        self.Sigma_tot = self.get_col("density").in_units("Msun/pc**2")
        #read RADMC synthetic observation data
        self.radmc = {}
        for i in np.arange(len(self.lines)):
            fn = self.fn_radmc[i]
            iline = self.lines[i]
            self.radmc[iline] = self.read_radmc(fn, iline)
            #XCO
            NH20 = self.NH2.in_units("cm**-2").value
            WCO = self.radmc[iline].WCO
            XCO_20 = NH20/WCO/1e20
            self.radmc[iline].XCO_20 = XCO_20
        return
    def read_radmc(self, fn, iline):
        """Read the radmc output."""
        if iline == 1:
            lam0 = 2600.75763346
        elif iline == 2:
            lam0 = 1300.4036558
        else:
            msg = "ERROR: Data.read_radmc(): "
            msg += "Line index {:d} not recognized (iline={1,2}).".format(
                    iline)
            raise RuntimeError(msg)
        Tdect = 0. #detection limit for each channel
        img = radmc.Images(fn, lam0=lam0)
        WCO = img.getimageWCO(Tdect=Tdect)
        Tpeak = img.getTpeak(Tdect=Tdect)
        _, sigmav= img.getimage_velocities(Tdect=Tdect)
        radmc_output = RadmcOutput(img, WCO, Tpeak, sigmav)
        return radmc_output
    def get_col(self, spec):
        dx = self.grid["dx"][0, 0, 0]
        nH = self.grid["nH"]
        if spec == "Htot":
            ns = nH
        elif spec == "density":
            ns = self.grid["density"]
        else:
            abd = self.grid[spec]
            ns = abd*nH
        col = ns.sum(axis=2)*dx
        return col

class RadmcOutput:
    """Class for radmc output data."""
    def __init__(self, img, WCO, Tpeak, sigmav):
        """Initialization.
        Inputs:
            img: class, full information of the radmc image
            WCO: 2D array, WCO in K km/s.
            Tpeak: peak brightness tempreature in K
            sigmav: velocity dispersion in km/s
        """
        self.img = img
        self.WCO = WCO
        self.Tpeak = Tpeak
        self.sigmav = sigmav
        return

class Model:
    """Model class for a set of simulation outputs."""
    def __init__(self, dir_ath, dir_radmc, model_mhd, model_chem, 
                 lines, zlim_pc, time_snapshots):
        """Initialization.
        Inputs:
            dir_ath: string, directory for athena++ MHD simulation outputs.
            dir_radmc, string, directory for radmc post-processing outputs.
            model_mhd: string, name of the MHD model. e.g. "R4".
            model_chem: string, name of the chemistry model. e.g. "Z1".
            lines: list of integers, CO line(s) to read. e.g. [1,2]
            zlim_pc: float, simulation is loaded at z=[-zlim_pc,zlim_pc]
            time_snapshots: arrray of integers, all time snapshots in Myr.
        """
        self.dir_ath = dir_ath
        self.dir_radmc = dir_radmc
        self.model_mhd = model_mhd
        self.model_chem = model_chem
        self.model_id = model_mhd + "-" + model_chem
        self.lines = lines
        self.zlim_pc = zlim_pc
        self.time_snapshots = time_snapshots
        return
    def get_data(self, time):
        """Load output data for a given time snapshot.
        input:
            time: int, certain time in the simulation in Myr.
                  Must be in time_snapshots.
        output:
            ds: instance of the class Data, containing simulation data.
        """
        time = int(time)
        ds = Data(self.dir_ath, self.dir_radmc, self.model_id, time,
                  self.lines, self.zlim_pc)
        return ds

def load(model_id, lines, dir_master, zlim_pc="auto"):
    """Load the simulation and synthetic observation data in
    Gong, Ostriker, Kim & Kim (2020).
    input: 
        model_id: string, Model ID in Table 2 of the paper, e.g. "R4-Z1".
        lines: list of integers, CO line(s) to read.
               e.g. lines=[1, 2] will read the CO (1-0) and (2-1) line. 
               The indexing 1 corresponds to the CO (1-0) line 
               and 2 corresponds to the CO (2-1) line, etc. 
               lines=[2] will read only the (2-1) line. 
        dir_master: string, master directory of the data storage.
        zlim_pc: float, simulation is loaded at z=[-zlim_pc,zlim_pc]
                 One can also set "zlim_pc=None", to read the whole simulation
                 box. The default is "zlim_pc='auto'", which adjust zlim_pc to
                 automatically match the mid-plane region where the CO line PPV
                 cubes are produced from RADMC-3D.
    output:
        model: instance of the Model class."""
    model_id_split = model_id.split("-")
    model_mhd = model_id_split[0]
    model_chem = model_id_split[1]
    dir_ath = "{:s}athena/{:s}/{:s}/".format(
                dir_master, model_mhd, model_chem)
    dir_radmc = "{:s}radmc/{:s}/{:s}/".format(
                dir_master, model_mhd, model_chem)
    if not os.path.isdir(dir_ath):
        msg = "ERROR: load(): In model_id="
        msg += "{:s}, athena++ directory {:s} do not exist.".format(
                model_id, dir_ath)
        raise RuntimeError(msg)
    if not os.path.isdir(dir_radmc):
        msg = "ERROR: load(): "
        msg += "In model_id={:s}, radmc directory {:s} do not exist.".format(
                model_id, dir_radmc)
        raise RuntimeError(msg)
    if model_mhd == "R2":
        time_snapshots = np.arange(40, 84, 4)
    elif model_mhd == "R2B2":
        time_snapshots = np.arange(40, 60, 4)
    elif model_mhd == "R4":
        time_snapshots = np.arange(50, 160, 10)
    elif model_mhd == "R8":
        time_snapshots = np.arange(300, 400, 10)
    else:
        msg = "ERROR: load(): "
        msg += "In model_id={:s}, MHD model {:s} do not exist".format(
                model_id, model_mhd)
        raise RuntimeError(msg)
    if zlim_pc == "auto":
        if model_mhd == "R2":
            zlim_pc = 300.
        elif model_mhd == "R2B2":
            zlim_pc = 400.
        elif model_mhd == "R4":
            zlim_pc = 512.
        elif model_mhd == "R8":
            zlim_pc = 300
    model = Model(dir_ath, dir_radmc, model_mhd, model_chem,
                  lines, zlim_pc, time_snapshots)
    return model
