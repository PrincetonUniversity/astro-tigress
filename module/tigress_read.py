import os
import os.path as osp
import numpy as np
import athena_read
import const
import ytathena as ya
import radmc
import urllib.request

class Model:
    """Class containing the simulation model information.

    Parameters
    ----------
    model_id: str
        name of the model, e.g. `R8_2pc`.
    """

    def __init__(self, model_id, dir_master="../data/"):
        self.dir_master = dir_master
        self.model_id = model_id
        self.dir_model = osp.join(dir_master, model_id, "")
        if osp.isdir(self.dir_model):
            subdirs = os.listdir(self.dir_model)
            ivtks = sorted([int(x) for x in subdirs if x.isnumeric()])
            self.ivtks = np.array(ivtks) #list of vtk numbers of outputs

            #contained data sets, e.g. ["MHD", "chem", "CO_lines"]
            self.data_sets = {}
            for ivtk in self.ivtks:
                dir_ivtk = "{:s}{:04d}/".format(self.dir_model, ivtk)
                self.data_sets[ivtk] = os.listdir(dir_ivtk)
        # check history file
        fn_hst = "{:s}history/MHD/{:s}.hst".format(self.dir_model,self.model_id)
        if osp.isfile(fn_hst):
            self.hst = athena_read.hst(fn_hst)
            self.fn_hst = fn_hst

        # read time output interval in input file
        fn_input_MHD = "{:s}input/MHD/{:s}.par".format(
                          self.dir_model,self.model_id)
        if osp.isfile(fn_input_MHD):
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

    def show(self):
        """Print available data list"""
        print("Avalable data:")
        if hasattr(self,'fn_hst'):
            print("  history")
        if hasattr(self,'fn_input_MHD'):
            print("  input")
        for i in self.ivtks:
            print("  ivtk = {:d}".format(i), end=" ")
            for d in self.data_sets[i]:
                print(d, end=" ")
            print(" ")

    def download(self, ivtk, dataset="MHD", Z=1., iline=1):
        """Download simulation data from a certain data set.

        Parameters
        ----------
        ivtk: int
            vtk output number
        dataset: str, ["MHD", "chem", "CO_lines", "history", "input", "all"]
            name of the dataset or download all.
        Z: float, [0.5, 1, 2]
            metallicity for chemistry post-processing.
        iline: int, [1, 2]
            uppper level of CO rotational line. 1 for J=1-0, 2 for J=2-1
        """

        fn = self._set_filename(ivtk, Z=Z, iline=iline)
        if dataset == "all":
            for d,f in fn.items():
                self._download_file(f)
        else:
            self._download_file(fn[dataset])


    def load(self, ivtk, dataset="MHD", Z=1., iline=1, Tdect=0.):
        """Load simulation data from a certain data set.

        Parameters
        ----------
        ivtk: int
            vtk output number
        dataset: str, ["MHD", "chem", "CO_lines", "all"]
            name of the dataset or load all.
        Z: float, [0.5, 1, 2]
            metallicity for chemistry post-processing.
        iline: int, [1, 2]
            uppper level of CO rotational line. 1 for J=1-0, 2 for J=2-1
        Tdect: float
            detection limit for atenna temperature in K
            for each velocity channel.
        """
        fn = self._set_filename(ivtk, Z=Z, iline=iline, add_master=True)
        if dataset == "all":
            for d in self.data_sets[ivtk]:
                self._loadone(fn[d], d, iline=iline, Tdect=Tdect)
        else:
            self._loadone(fn[dataset], dataset, iline=iline, Tdect=Tdect)

    def _loadone(self,f,d,iline=1,Tdect=0.):
        if d == "MHD":
            self.MHD = DataMHD(f)
        elif d == "chem":
            self.chem = DataChem(f)
        elif d == "CO_lines":
            self.CO_lines[iline] = DataCO(f, iline, Tdect)
        else:
            msg = "ERROR: Model.load(): dataset name not recogonozed.\n"
            msg += 'dataset: {"MHD", "chem", "CO_lines"}'
            raise RuntimeError(msg)

    def _set_filename(self, ivtk, Z=1., iline=1, add_master=False):
        source_dir_ivtk = "{:s}/{:04d}/".format(self.model_id, ivtk)
        if add_master: source_dir_ivtk = osp.join(self.dir_master, source_dir_ivtk)
        dict_Z = {0.5:"Z05", 1.:"Z1", 2.:"Z2"}
        Z_id = dict_Z[Z]
        fn = dict()
        fn['MHD'] = "{:s}MHD/{:s}.{:04d}.vtk".format(
                      source_dir_ivtk, self.model_id, ivtk)
        fn['chem'] = "{:s}chem/{:s}/{:s}-{:s}.{:04d}.athdf".format(
                      source_dir_ivtk, Z_id, self.model_id, Z_id, ivtk)
        fn['CO_lines'] = "{:s}CO_lines/{:s}/il{:d}/{:s}-{:s}.il{:d}.{:04d}.bout".format(
                      source_dir_ivtk, Z_id, iline, self.model_id, Z_id, iline, ivtk)
        fn['history'] = "{0:s}/history/MHD/{0:s}.hst".format(self.model_id)
        fn['input'] = "{0:s}/input/MHD/{0:s}.par".format(self.model_id)
        return fn

    def _download_file(self, f):
        import urllib.request
        from urllib.error import URLError
        import shutil

        target = osp.join(self.dir_master,f)
        if osp.isfile(target):
            print("{} with size of {}GB already exists".format(f,
                osp.getsize(target)/2**30))
            while True:
                answer = input("overwrite? [y/n]:")
                if answer.lower() in ['y','n']:
                    break
            if answer.lower() == 'n': return
        os.makedirs(osp.dirname(target),exist_ok=True)

        url = 'https://tigress-web.princeton.edu/~munan/astro-tigress/'
        source = url+f
        req = urllib.request.Request(source)

        try:
            response = urllib.request.urlopen(req)
        except URLError as e:
            if hasattr(e, 'reason'):
                print('We failed to reach a server.')
                print('Reason: ', e.reason)
            elif hasattr(e, 'code'):
                print('The server couldn\'t fulfill the request.')
                print('Error code: ', e.code)
        else:
            print("We reached ", url)
            print("  downloading...",osp.basename(source),end=" ")
            with urllib.request.urlopen(source) as response, \
                    open(target, 'wb') as target:
                shutil.copyfileobj(response, target)
            #urllib.request.urlretrieve(source,target)
            print(" complete!")

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

