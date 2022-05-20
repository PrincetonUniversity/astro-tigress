import os
import numpy as np
import athena_read
import const
import ytathena as ya
import radmc

class Model:
    """Class containing the simulation model information."""
    def __init__(self, dir_master, model_id):
        """Initialization.
        Inputs:
            model_id: string, name of the model, e.g. "R8_2pc"."""
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
        return
    def load(self, ivtk, dataset="MHD"):
        """Load simulation data from a certain data set.
        input:
            ivtk: int, vtk output number
            dataset: string, name of the data set. {"MHD", "chem", "CO+lines"}
        """
        dir_ivtk = "{:s}{:04d}/".format(self.dir_model, ivtk)
        if dataset == "MHD":
            fn_MHD = "{:s}MHD/{:s}.{:04d}.vtk".format(
                      dir_ivtk, self.model_id, ivtk)
            self.MHD = DataMHD(fn_MHD)
        else:
            msg = "ERROR: Model.load(): dataset name not recogonozed.\n"
            msg += 'dataset: {"MHD"}'
            raise RuntimeError(msg)
        return

class DataMHD:
    """MHD data in one time snapshot."""
    def __init__(self, fn):
        self.ytds = ya.load_athena4p2(fn)
        self.dx = ( ( self.ytds.domain_right_edge[0] 
                   - self.ytds.domain_left_edge[0] )
                     /self.ytds.domain_dimensions[0] ).value
        self.grid = ya.get_covering_grid(self.ytds)
        return
