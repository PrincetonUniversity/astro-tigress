import numpy as np
import struct
import os
from subprocess import Popen, PIPE, STDOUT
import copy
import math
import yt
import yt.units as u
from astropy.io import fits

from . import const
from . import ytathena as ya
from . import map_circular_beam as mcb


class YTData:
    def __init__(self, filename, left_edge=None, dims=None):
        self.fn = filename
        fn_ext = filename.split(".")[-1]
        if fn_ext == "athdf":
            self.ytds = ya.load_athenapp(filename)
        else:
            raise RuntimeError("file extension not recogonized: " + fn_ext)
        if left_edge is None:
            self.left_edge = self.ytds.domain_left_edge
        else:
            self.left_edge = left_edge
        if dims is None:
            self.dims = self.ytds.domain_dimensions
        else:
            self.dims = dims
        self.ytdata = self.ytds.covering_grid(
            level=0, left_edge=self.left_edge, dims=self.dims
        )
        self.shape_xyz = self.dims
        self.nx = self.shape_xyz[0]
        self.ny = self.shape_xyz[1]
        self.nz = self.shape_xyz[2]
        self.dx = (self.ytdata["dx"][0, 0, 0].in_cgs()).value
        nH0 = (
            (self.ytdata["density"] / 2.38858753789e-24 / u.g * u.cm**3)
            .in_cgs()
            .to_ndarray()
        )
        self.nH = np.swapaxes(nH0, 0, 2)
        # velocity in km/s
        velx0 = (self.ytdata["velocity_x"] / (u.km / u.s)).in_cgs().to_ndarray()
        vely0 = (self.ytdata["velocity_y"] / (u.km / u.s)).in_cgs().to_ndarray()
        velz0 = (self.ytdata["velocity_z"] / (u.km / u.s)).in_cgs().to_ndarray()
        vel = np.zeros([self.nx, self.ny, self.nz, 3])
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    vel[ix, iy, iz, 0] = velx0[ix, iy, iz]
                    vel[ix, iy, iz, 1] = vely0[ix, iy, iz]
                    vel[ix, iy, iz, 2] = velz0[ix, iy, iz]
        self.vel_kms = np.swapaxes(vel, 0, 2)
        # abd: CO, H2, H, T
        self.abd = {}
        for s in ["CO", "HCO+", "H2", "H", "e"]:
            self.abd[s] = np.swapaxes(self.ytdata[s].to_ndarray(), 0, 2)
        self.abd["T"] = np.swapaxes(
            (self.ytdata["temperature_chem"] / u.K).in_cgs().to_ndarray(), 0, 2
        )
        self.abd["Tinit"] = np.swapaxes(
            (self.ytdata["temperature_init"] / u.K).in_cgs().to_ndarray(), 0, 2
        )
        return


class CoreData:
    def __init__(self, fn, left_edge=None, dims=None, n0=1.0e3, T0=12.0):
        self.n0 = n0
        self.T0 = T0
        self.mu = 1.4  # mean molecular weight including He
        self.cs0_kms = 0.19 * (T0 / 10.0) ** 0.5
        L0 = 0.87 * const.pc * (n0 / 1.0e3) ** (-0.5) * (T0 / 10.0) ** 0.5
        L0_pc = L0 / const.pc
        M0 = L0**3 * n0 * const.mh * self.mu
        unit_base = {
            "length_unit": (L0_pc, "pc"),
            "time_unit": (L0_pc / self.cs0_kms, "s*pc/km"),
            "mass_unit": (M0, "g"),
            "velocity_unit": (self.cs0_kms, "km/s"),
        }
        # load data
        self.ytds = yt.load(fn, units_override=unit_base)
        if left_edge is None:
            self.left_edge = self.ytds.domain_left_edge
        else:
            self.left_edge = left_edge
        if dims is None:
            self.dims = self.ytds.domain_dimensions
        else:
            self.dims = dims
        self.ytdata = self.ytds.covering_grid(
            level=0, left_edge=self.left_edge, dims=self.dims
        )
        self.shape_xyz = self.dims
        self.nx = self.shape_xyz[0]
        self.ny = self.shape_xyz[1]
        self.nz = self.shape_xyz[2]
        self.dx = (self.ytdata["dx"][0, 0, 0].in_cgs()).value
        rho = self.ytdata["density"].in_cgs().to_ndarray()
        nH0 = rho / const.mh / self.mu
        self.nH = np.swapaxes(nH0, 0, 2)
        # velocity in km/s
        velx0 = (self.ytdata["velocity_x"] / (u.km / u.s)).in_cgs().to_ndarray()
        vely0 = (self.ytdata["velocity_y"] / (u.km / u.s)).in_cgs().to_ndarray()
        velz0 = (self.ytdata["velocity_z"] / (u.km / u.s)).in_cgs().to_ndarray()
        vel = np.zeros([self.nx, self.ny, self.nz, 3])
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    vel[ix, iy, iz, 0] = velx0[ix, iy, iz]
                    vel[ix, iy, iz, 1] = vely0[ix, iy, iz]
                    vel[ix, iy, iz, 2] = velz0[ix, iy, iz]
        self.vel_kms = np.swapaxes(vel, 0, 2)
        # abd: temprature, H2, NH3
        self.abd = {}
        self.abd["T"] = np.zeros(self.nH.shape) + self.T0
        self.abd["H2"] = np.zeros(self.nH.shape) + 0.5
        xNH3 = np.zeros(self.nH.shape)
        n_th = 6.0e3
        xNH3[self.nH > n_th] = 10 ** (-8.5)
        self.abd["pNH3"] = xNH3
        return


class Data:
    """Read from simulation output, and functions to write RADMC-3D input
    files, run radmc3d to produce outputs."""

    def __init__(self, fn, fmrt="TIGRESS", left_edge=None, dims=None, **keys):
        """fmrt: input format, ["TIGRESS", "core"]"""
        self.fn = fn
        self.fmrt = fmrt
        if self.fmrt == "TIGRESS":
            self.data = YTData(self.fn, left_edge=left_edge, dims=dims, **keys)
            for s in ["CO", "HCO+", "H2", "H", "e"]:
                indx = self.data.abd[s] < 0.0
                self.data.abd[s][indx] = 0.0
        elif self.fmrt == "core":
            self.data = CoreData(self.fn, left_edge=left_edge, dims=dims, **keys)
        else:
            raise RuntimeError("file format not recogonized: " + self.fmrt)
        self.dx = self.data.dx
        self.nx = self.data.nx
        self.ny = self.data.ny
        self.nz = self.data.nz
        self.shape_xyz = self.data.shape_xyz
        self.linemodes = {}
        self.linemodes["LTE"] = "1"
        self.linemodes["LVG"] = "3"
        self.linemodes["Thin"] = "4"
        return

    def get_NH2(self, char_axis):
        """return in [nx, ny, nz] shape"""
        if char_axis == "x":
            nax = 0
        elif char_axis == "y":
            nax = 1
        elif char_axis == "z":
            nax = 2
        else:
            raise RuntimeError("char_axis has to be {x, y, z}")
        nH = np.swapaxes(self.data.nH, 0, 2)
        xH2 = np.swapaxes(self.data.abd["H2"], 0, 2)
        return (self.dx * nH * xH2).sum(axis=nax)

    def prepare(
        self,
        work_dir,
        image_dir,
        clean_outfile=True,
        clean_infile=True,
        species_coll=None,
        moldata_dir="/tigress/munan/chemistry/scripts/radmc3d/prob/moldata/",
        line_mode="LVG",
        output_binary=True,
        escprob=None,
        escprob_max_pc=100.0,
        gas_velocity=None,
        micro_turbulence=None,
        gas_temperature="data",
        gas_temperature_max=200.0,
        temperature_bg=2.73,
        gas_temperature_cold=None,
        slow_lvg_conv=None,
        species_main="co",
        minimum_abundance=None,
        ortho2paraH2=3.0,
        iline=1,
        widthkms=1,
        vkms=0,
        linenlam=40,
        isloadlambda=False,
        external_tbg=None,
    ):
        """Prepare for input paramters, create simple input files.
        This can be done multiple times to the Data class (work_dir).
        moldata_dir: store molecular line data for co.
        line_mode: LVG, LTE, Thin
        escprob, gas_veolcity, micro_turbulence, gas_temperature:
        {"data" (in Data class), None (not specified), float_number (constant), np.array([nx,ny,nz]) }
        gas_temperature_max: cap the gas temperature at some maximum value. Default 200. Kelvin.
        temperature_bg: back ground temperature. Default 2.73K for cmb.
        gas_temperature_cold: upper liit of temperature of cold gas. If not
        None, set species abundance in any hotter gas to be zero.
        escprob_max_pc: cap the escape probability maximu, value in pc. Default 100 pc.
        gas_temperature_max, gas_temperature_cold and escprob_max_pc
        are only used when the choice 'data' is specified.
        slow_lvg_conv: if sepecified, add slow lvg and the convergence
        criteria in radmc3d.inp"""
        # create working directory if it does not exist
        self.work_dir = work_dir
        self.image_dir = image_dir
        self.species_coll = species_coll
        self.species_main = species_main
        self.minimum_abundance = minimum_abundance
        if species_coll is not None:
            self.ncoll = len(species_coll)  # number of collisional species
        else:
            self.ncoll = 0
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        else:
            print("{} already exists. Using the directory.".format(work_dir))
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)
        else:
            print("{} already exists. Using the directory.".format(image_dir))
        # clean up directory
        if clean_outfile:
            print("Cleaning radmc3d output files...")
            cmd = "rm " + work_dir + "*\.*out " + work_dir + "*\.*dat"
            p = Popen(
                cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
            )
            print(p.stdout.read())
        if clean_infile:
            print("Cleaning radmc3d input files...")
            cmd = "rm " + work_dir + "*\.*inp "
            p = Popen(
                cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
            )
            print(p.stdout.read())
        # copy molecular line data into directory
        print("Copy molecular data file into directory..")
        fn_mol = "molecule_" + self.species_main + ".inp"
        cmd = "cp " + moldata_dir + fn_mol + " " + work_dir
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
        )
        print(p.stdout.read())
        # parameters for the line module
        self.iline = iline
        self.widthkms = widthkms
        self.vkms = vkms
        self.linenlam = linenlam
        self.isloadlambda = isloadlambda
        self.external_tbg = external_tbg
        # calculate wavelength points
        if self.isloadlambda:
            line_nu0 = read_line_freq(moldata_dir + fn_mol)
            self.wavelength_micron = get_wavelength_micron(
                line_nu0, iline, widthkms, vkms, linenlam
            )
        # external radiation requires loadlambda to be turned on
        if external_tbg is not None:
            if not isloadlambda:
                raise RuntimeError(
                    "ERORR:"
                    + " external radiation option requires loadlambda to be turned on."
                )
            external_intensity = const.bb_inu(
                temp_k=external_tbg, lam_cm=self.wavelength_micron / 1.0e4, units=False
            )
            print("Creating external_source.inp...")
            fw = open(work_dir + "external_source.inp", "w")
            fw.write("2\n{:d}".format(self.linenlam))
            for li in self.wavelength_micron:
                fw.write("\n{:f}".format(li))
            for ii in external_intensity:
                fw.write("\n{:.6e}".format(ii))
            fw.close()
        # create camera_wavelength_micron.inp if loadlambda is used
        if self.isloadlambda:
            print(
                "Creating camera_wavelength_micron.inp required by loadlambda option..."
            )
            fw = open(work_dir + "camera_wavelength_micron.inp", "w")
            fw.write("{:d}".format(self.linenlam))
            for li in self.wavelength_micron:
                fw.write("\n{:f}".format(li))
            fw.close()
            # copy to wavelength_micron.inp, used to check background radiation
            cmd = (
                "cp "
                + work_dir
                + "camera_wavelength_micron.inp "
                + (work_dir + "wavelength_micron.inp")
            )
            p = Popen(
                cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
            )
            print(p.stdout.read())
        else:
            # create wavelength_micron.inp necessary for file
            print("Creating wavelength_micron.inp required by radmc3d..")
            fw = open(work_dir + "wavelength_micron.inp", "w")
            fw.write("2\n0.1\n100\n\n")
            fw.close()
        # create lines.inp
        print("Writing lines.inp...")
        fl = open(work_dir + "lines.inp", "w")
        fl.write(
            "2\n1\n"
            + self.species_main
            + "    leiden    0    0    "
            + str(self.ncoll)
            + "\n"
        )
        if self.species_coll is not None:
            for s in species_coll:
                fl.write(s)
                fl.write("\n")
        fl.write("\n")
        fl.close()
        # create radmc3d.inp
        print("Writing radmc3d.inp...")
        fr = open(work_dir + "radmc3d.inp", "w")
        if output_binary:
            fr.write("writeimage_unformatted = 1\n")
            fr.write("lines_mode = " + self.linemodes[line_mode])
            fr.write("\n")
            fr.write("lines_tbg = " + str(temperature_bg))
            fr.write("\n")
        else:
            fr.write("lines_mode = " + self.linemodes[line_mode])
            fr.write("\n")
        if slow_lvg_conv is not None:
            fr.write(
                "lines_nonlte_convcrit = "
                + str(slow_lvg_conv)
                + "\n"
                + "lines_slowlvg_as_alternative=1\n"
            )
        fr.write("\n\n")
        fr.close()
        # store input parameters
        self.line_mode = line_mode
        self.output_binary = output_binary
        self.escprob = escprob
        self.escprob_max_pc = escprob_max_pc
        self.gas_velocity = gas_velocity
        self.micro_turbulence = micro_turbulence
        self.gas_temperature = gas_temperature
        self.gas_temperature_max = gas_temperature_max
        self.temperature_bg = temperature_bg
        self.gas_temperature_cold = gas_temperature_cold
        self.ortho2paraH2 = ortho2paraH2
        # give summary of parameters.
        print("\nSummary, radmc input parameters:")
        print(
            "species:{}, minimum abundance:{}".format(
                self.species_main, self.minimum_abundance
            )
        )
        print("Working directory: {}".format(self.work_dir))
        print(
            "{} collisional species: {}, oH2/pH2={}".format(
                self.ncoll, self.species_coll, self.ortho2paraH2
            )
        )
        # print("Note that the collisional speceis in molecule_co.inp is pH2 and oH2!")
        print(
            "line_mode={}, output_binary={}".format(self.line_mode, self.output_binary)
        )
        print(
            "escprob={}, escprob_max_pc={}, gas_velocity={}, micro_turbulence={}".format(
                self.escprob,
                self.escprob_max_pc,
                self.gas_velocity,
                self.micro_turbulence,
            )
        )
        print(
            "gas_temperature={}, gas_temperature_max={}, gas_temperature_cold={}".format(
                self.gas_temperature,
                self.gas_temperature_max,
                self.gas_temperature_cold,
            )
        )
        print(
            "iline={}, widthkms={}, vkms={}, linenlam={}, isloadlambda={}".format(
                iline, widthkms, vkms, linenlam, isloadlambda
            )
        )
        return

    def write_inputs(self):
        """Write main radmc3d input files *.binp"""
        print("\nWriting main input files for radmc3d...")
        # amr_grid.binp
        xmin = -self.dx * self.nx / 2
        ymin = -self.dx * self.ny / 2
        zmin = -self.dx * self.nz / 2
        self._write_grid_binp(
            self.work_dir, self.nx, self.ny, self.nz, self.dx, xmin, ymin, zmin
        )
        # numberdens_species.binp
        fn_ns = self.work_dir + "numberdens_" + self.species_main + ".binp"
        if self.species_main == "co":
            ns = np.swapaxes(self.data.abd["CO"] * self.data.nH, 0, 2)
        elif self.species_main == "hcop":
            ns = np.swapaxes(self.data.abd["HCO+"] * self.data.nH, 0, 2)
        elif self.species_main == "pnh3":
            ns = np.swapaxes(self.data.abd["pNH3"] * self.data.nH, 0, 2)
        elif self.species_main == "HI":
            ns = np.swapaxes(self.data.abd["H"] * self.data.nH, 0, 2)
        else:
            raise RuntimeError("species {} not implemented".format(self.species_main))
        if self.gas_temperature_cold is not None:
            Tg1 = np.swapaxes(self.data.abd["Tinit"], 0, 2)
            indx_Thot = np.where(Tg1 > self.gas_temperature_cold)
            ns[indx_Thot] = 0.0
        if self.minimum_abundance is not None:
            indx_empty = np.where(ns < self.minimum_abundance)
            ns[indx_empty] = 0.0
        self._write_binp(ns, self.nx, self.ny, self.nz, fn_ns)
        # numberdens_***.binp for collisional species
        # o-H2 and p-H2
        fo = float(self.ortho2paraH2) / float(1.0 + self.ortho2paraH2)
        fp = 1.0 - fo
        if self.species_coll is not None:
            for si in range(self.ncoll):
                s = self.species_coll[si]
                if s == "oH2":
                    fn_oh2 = self.work_dir + "numberdens_oH2.binp"
                    h20 = np.swapaxes(self.data.abd["H2"] * self.data.nH, 0, 2)
                    oh20 = h20 * fo
                    self._write_binp(oh20, self.nx, self.ny, self.nz, fn_oh2)
                elif s == "pH2":
                    fn_ph2 = self.work_dir + "numberdens_pH2.binp"
                    h20 = np.swapaxes(self.data.abd["H2"] * self.data.nH, 0, 2)
                    ph20 = h20 * fp
                    self._write_binp(ph20, self.nx, self.ny, self.nz, fn_ph2)
                else:
                    fn_s = self.work_dir + "numberdens_" + s + ".binp"
                    s0 = np.swapaxes(self.data.abd[s] * self.data.nH, 0, 2)
                    self._write_binp(s0, self.nx, self.ny, self.nz, fn_s)
                    # co_coll = ["pH2", "oH2"]
                    # print("Warning: species {} used as {}".format(s, co_coll[si]))
        # gas_temperature.binp
        fn_Tg = self.work_dir + "gas_temperature.binp"
        if self.gas_temperature == "data":
            Tg0 = np.swapaxes(self.data.abd["T"], 0, 2)
            indx_Tghigh = np.where(Tg0 > self.gas_temperature_max)
            Tg0[indx_Tghigh] = self.gas_temperature_max
        elif isinstance(self.gas_temperature, float):
            Tg0 = np.zeros([self.nx, self.ny, self.nz]) + self.gas_temperature
        elif isinstance(self.gas_temperature, np.ndarray):
            if self.gas_temperature.shape == self.shape_xyz:
                Tg0 = self.gas_temperature
            else:
                raise RuntimeError(
                    "gas_temperature array has shape {} != grid shape {}".format(
                        self.gas_temperature.shape, self.shape_xyz
                    )
                )
        else:
            raise RuntimeError(
                "gas_temperature: {}, not recognized.".format(self.gas_temperature)
            )
        self._write_binp(Tg0, self.nx, self.ny, self.nz, fn_Tg)

        # escprob_lengthscale.binp
        fn_Lesc = self.work_dir + "escprob_lengthscale.binp"
        if self.escprob == "data":
            grad = np.gradient(np.swapaxes(self.data.nH, 0, 2), self.dx)
            Lesc0 = 1.0 / np.sqrt(grad[0] ** 2 + grad[1] ** 2 + grad[2] ** 2)
            indx = Lesc0 > self.escprob_max_pc * const.pc
            Lesc0[indx] = self.escprob_max_pc * const.pc
        elif self.escprob is None:
            print("No escprob_lengthscale.binp, Lesc = inf.")
        elif self.escprob == "jeans":
            # use jeans length with temperature cap at 40K
            nH0 = self.data.nH
            mu = const.mh * 1.4 / (1 - self.data.abd["H2"] + self.data.abd["e"] + 0.1)
            Tcs = np.copy(self.data.abd["T"])
            indxTcs = Tcs > 40.0
            Tcs[indxTcs] = 40.0
            cs = np.sqrt(const.kb * Tcs / mu)
            LJ = cs * np.sqrt(math.pi / (const.g * nH0 * mu))
            Lesc0 = np.swapaxes(LJ, 0, 2)
        elif isinstance(self.escprob, float):
            Lesc0 = np.zeros(self.shape_xyz) + self.escprob
        elif isinstance(self.escprob, np.ndarray):
            if self.escprob.shape == self.shape_xyz:
                Lesc0 = self.escprob
            else:
                raise RuntimeError(
                    "escprob array has shape {} != grid shape {}".format(
                        self.escprob.shape, self.shape_xyz
                    )
                )
        else:
            raise RuntimeError("escprob: {}, not recognized.".format(self.escprob))
        if self.escprob is not None:
            self._write_binp(Lesc0, self.nx, self.ny, self.nz, fn_Lesc)

        # micro_turbulence.binp
        fn_mt = self.work_dir + "microturbulence.binp"
        if self.micro_turbulence == "data":
            mt0 = np.zeros(self.shape_xyz) + const.km * np.sqrt(self.dx / const.pc)
        elif self.micro_turbulence is None:
            print("No micro_turbulence.binp, vturb = 0.")
        elif isinstance(self.micro_turbulence, float):
            mt0 = np.zeros(self.shape_xyz) + self.micro_turbulence
        elif isinstance(self.micro_turbulence, np.ndarray):
            if self.micro_turbulence.shape == self.shape_xyz:
                mt0 = self.micro_turbulence
            else:
                raise RuntimeError(
                    "micro_turbulence array has shape {} != grid shape {}".format(
                        self.micro_turbulence.shape, self.shape_xyz
                    )
                )
        else:
            raise RuntimeError(
                "micro_turbulence: {}, not recognized.".format(self.micro_turbulence)
            )
        if self.micro_turbulence is not None:
            self._write_binp(mt0, self.nx, self.ny, self.nz, fn_mt)

        # gas_velocity.binp
        vel_shape = (self.nx, self.ny, self.nz, 3)
        if self.gas_velocity is None:
            print("No gas_velocity.binp, vgas = 0.")
        elif self.gas_velocity == "data":
            fn_vel = self.work_dir + "gas_velocity.binp"
            vel = self.data.vel_kms
            vel0 = np.swapaxes(vel, 0, 2) * 1.0e5
        elif isinstance(self.gas_velocity, float):
            vel0 = np.zeros(vel_shape) + self.gas_velocity
        elif isinstance(self.gas_velocity, np.ndarray):
            if self.gas_velocity.shape == vel_shape:
                vel0 = self.gas_velocity
            else:
                raise RuntimeError(
                    "gas_velocity array has shape {} != required shape {}".format(
                        self.gas_velocity.shape, vel_shape
                    )
                )
        else:
            raise RuntimeError(
                "gas_velocity: {}, not recognized.".format(self.gas_velocity)
            )
        if self.gas_velocity is not None:
            self._write_binp(vel0, self.nx, self.ny, self.nz, fn_vel, ndim=3)
            return

    def run(
        self,
        image_new_name=None,
        levelpop_new_name=None,
        tausurf=False,
        cmd_tail="",
        incl=0,
        phi=0,
    ):
        """Run radmc3d in working directory, move image.dat to a new file if image_new_name is specified."""
        if levelpop_new_name is not None:
            writepop = "writepop"
        else:
            writepop = ""
        if tausurf:
            cmd_image = "tausurf"
        else:
            cmd_image = "image"
        if self.isloadlambda:
            cmd_lambda = " loadlambda"
        else:
            cmd_lambda = (
                " iline {:d} widthkms ".format(self.iline)
                + str(self.widthkms)
                + " vkms "
                + str(self.vkms)
                + " linenlam "
                + str(self.linenlam)
            )
        cmd = (
            "cd "
            + self.work_dir
            + ";radmc3d "
            + cmd_image
            + " incl "
            + str(incl)
            + " phi "
            + str(phi)
            + " npix "
            + str(self.nx)
            + " sizepc "
            + str(self.nx * self.dx / const.pc)
            + cmd_lambda
            + " doppcatch "
            + writepop
            + cmd_tail
        )
        if self.output_binary:
            img_ext = "bout"
        else:
            img_ext = "out"
        pop_ext = "dat"
        print(cmd)
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
        )
        output = (p.stdout.read()).decode("ascii")
        self.output = output
        fo = open(self.work_dir + "output.txt", "w")
        fo.write(output)
        fo.close()
        self.image_new_name = image_new_name
        if image_new_name is not None:
            cmd = "mv " + self.work_dir + "image." + img_ext + " " + image_new_name
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
        )
        print(cmd)
        print(p.stdout.read())
        if levelpop_new_name is not None:
            cmd = (
                "mv "
                + self.work_dir
                + "levelpop_"
                + self.species_main
                + "."
                + pop_ext
                + " "
                + levelpop_new_name
            )
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
        )
        print(cmd)
        print(p.stdout.read())
        return

    def readimage(self):
        # go to run directory
        img = Images(self.image_new_name)
        img.radmcData = copy.deepcopy(self)
        return img

    def _write_grid_binp(self, dir_out, nx, ny, nz, dx, xmin, ymin, zmin):
        """dx=dy=dz, cgs units, dir_out/"""
        fn_grid = dir_out + "amr_grid.binp"
        xi_arr = np.linspace(xmin, nx * dx + xmin, nx + 1)
        yi_arr = np.linspace(ymin, ny * dx + ymin, ny + 1)
        zi_arr = np.linspace(zmin, nz * dx + zmin, nz + 1)
        amr_header = (1, 0, 0, 0, 1, 1, 1, nx, ny, nz)
        fout = open(fn_grid, "wb")
        amr_header_b = struct.pack("q" * len(amr_header), *amr_header)
        xi_arr_b = struct.pack("d" * len(xi_arr), *xi_arr)
        yi_arr_b = struct.pack("d" * len(yi_arr), *yi_arr)
        zi_arr_b = struct.pack("d" * len(zi_arr), *zi_arr)
        fout.write(amr_header_b)
        fout.write(xi_arr_b)
        fout.write(yi_arr_b)
        fout.write(zi_arr_b)
        fout.close()
        return

    def _write_binp(self, arr0, nx, ny, nz, fn, ndim=1):
        arr = np.zeros(nx * ny * nz * ndim)
        i = 0
        if ndim > 1:
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        for idim in range(ndim):
                            arr[i] = arr0[ix, iy, iz, idim]
                            i = i + 1
        else:
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        arr[i] = arr0[ix, iy, iz]
                        i = i + 1
        header = (1, 8, nx * ny * nz)
        header_b = struct.pack("q" * len(header), *header)
        arr_b = struct.pack("d" * len(arr), *arr)
        fout1 = open(fn, "wb")
        fout1.write(header_b)
        fout1.write(arr_b)
        fout1.close()
        return


class Images:
    """Read RADMC-3D image output. Functions include plotting image and sed."""

    def __init__(self, fn_image, lam0=2600.75763346):
        """Read image from RADMC3D output.
        lam0: wavelength at line center in microns, used to calculate velocity v.
        Default: 2600.757: CO J1-0 line center."""

        ext = fn_image.split(".")[-1]
        if ext == "out":
            self.binary_output = False
        elif ext == "bout":
            self.binary_output = True
        else:
            raise RuntimeError("Extension {} not recognized.".format(ext))
        # read image output
        image_all = []
        if self.binary_output:
            fimage = open(fn_image, "rb")
            _ = fimage.read(8)  # iformat
            precis = 8  # precision
            real = "d"
            lines = fimage.read(8 * 3)
            lines = struct.unpack("q" * 3, lines)
            nximg = lines[0]
            nyimg = lines[1]
            nlam = lines[2]
            lines = fimage.read(precis * 2)
            lines = struct.unpack(real * 2, lines)
            pixsize_x = lines[0]
            pixsize_y = lines[1]
            lines = fimage.read(precis * nlam)
            lines = struct.unpack(real * nlam, lines)
            lam_arr = np.array(lines)
            for i in range(nlam):
                image = np.zeros((nximg, nyimg))
                for iy in range(nyimg):
                    for ix in range(nximg):
                        fi = fimage.read(precis)
                        fi = struct.unpack(real, fi)[0]
                        image[ix, iy] = fi
                # image1 = np.swapaxes(image, 0, 1)
                image_all.append(image)
            fimage.close()
        else:
            fimage = open(fn_image, "r")
            _ = fimage.readline()  # format 1
            ldim = fimage.readline().split()
            nximg = int(ldim[0])
            nyimg = int(ldim[0])
            # print nx, ny
            nlam = int(fimage.readline())
            # print nlam
            lpix = fimage.readline().split()
            pixsize_x = float(lpix[0])
            pixsize_y = float(lpix[1])
            # print pixsize_x/const.pc, pixsize_y/const.pc
            lam_arr = np.zeros(nlam)
            for i in range(nlam):
                llam = fimage.readline().split()
                lam_arr[i] = float(llam[0])
            # print lam_arr
            for i in range(nlam):
                _ = fimage.readline()  # empty line
                image = np.zeros((nximg, nyimg))
                for iy in range(nyimg):
                    for ix in range(nximg):
                        fi = float(fimage.readline())
                        image[ix, iy] = fi
                # image1 = np.swapaxes(image, 0, 1)
                image_all.append(image)
            fimage.close()
        # assign to class
        self.fn = fn_image
        self.nlam = nlam
        self.lam_micron = lam_arr
        self.lam_cm = lam_arr * 1.0e-4
        self.nx = nximg
        self.ny = nyimg
        self.dx = pixsize_x
        self.dy = pixsize_y
        self.image = image_all
        self.nu = const.c / self.lam_cm  # frequency in Hz
        # line center wavelength and frequency
        self.lam0_micron = lam0
        self.lam0_cm = lam0 * 1.0e-4
        self.nu0 = const.c / self.lam0_cm
        # velocity relative to the line center, in km/s
        self.vel_kms = -const.c * (self.nu - self.nu0) / self.nu0 / const.km
        # translate images in to Antena temperature T_A units (Kevin)
        self.image_TA = []
        for i in range(nlam):
            im = self.image[i]
            nui = self.nu[i]
            fac_TA = const.c**2 / (2.0 * const.kb * nui**2)
            self.image_TA.append(im * fac_TA)
        # get the average image
        self.image_avg = [im.mean() for im in self.image]
        self.image_TA_avg = [im.mean() for im in self.image_TA]
        return

    def clip(self, ixmin, ixmax, iymin, iymax):
        self.nx = ixmax - ixmin
        self.ny = iymax - iymin
        image1 = []
        image_TA1 = []
        for im in self.image:
            im1 = im[ixmin:ixmax, iymin:iymax]
            image1.append(im1)
        self.image = image1
        for im in self.image_TA:
            im1 = im[ixmin:ixmax, iymin:iymax]
            image_TA1.append(im1)
        self.image_TA = image_TA1
        self.image_avg = [im.mean() for im in self.image]
        self.image_TA_avg = [im.mean() for im in self.image_TA]
        return

    def write_fits(
        self,
        fn_fits="image.fits",
        clip_range=None,
        comment="",
        axis_name1="x",
        axis_name2="y",
        extent=None,
    ):
        if clip_range is not None:
            self.clip(clip_range[0], clip_range[1], clip_range[2], clip_range[3])
        if extent is None:
            xmin = -self.dx * self.nx / 2.0 / const.pc
            xmax = self.dx * self.nx / 2.0 / const.pc
            ymin = -self.dy * self.ny / 2.0 / const.pc
            ymax = self.dy * self.ny / 2.0 / const.pc
        else:
            xmin = extent[0]
            xmax = extent[1]
            ymin = extent[0]
            ymax = extent[1]
        data = np.array(self.image_TA).swapaxes(1, 2)
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header.set("author", "Munan Gong")
        primary_hdu.header.set("comment", comment)
        primary_hdu.header.set(axis_name1.upper() + "MIN", xmin)
        primary_hdu.header.set(axis_name1.upper() + "MAX", xmax)
        primary_hdu.header.set("D" + axis_name1.upper(), self.dx / const.pc)
        primary_hdu.header.set(axis_name2.upper() + "MIN", ymin)
        primary_hdu.header.set(axis_name2.upper() + "MAX", ymax)
        primary_hdu.header.set("D" + axis_name2.upper(), self.dy / const.pc)
        primary_hdu.header.set("VMIN", self.vel_kms[0])
        primary_hdu.header.set("VMAX", self.vel_kms[-1])
        primary_hdu.header.set("DV", self.vel_kms[1] - self.vel_kms[0])
        hdu_data = fits.ImageHDU(data)
        hdu_data.header.set("EXTNAME", "TA")
        for hdui in [primary_hdu, hdu_data]:
            hdui.header.set("CDELT1", self.dx / const.pc)
            hdui.header.set("CTYPE1", axis_name1)
            hdui.header.set("CUNIT1", "pc")
            hdui.header.set("CRVAL1", xmin)
            hdui.header.set("CRPIX1", 0)
            hdui.header.set("CDELT2", self.dy / const.pc)
            hdui.header.set("CTYPE2", axis_name2)
            hdui.header.set("CUNIT2", "pc")
            hdui.header.set("CRVAL2", ymin)
            hdui.header.set("CRPIX2", 0)
            hdui.header.set("CDELT3", self.vel_kms[1] - self.vel_kms[0])
            hdui.header.set("CTYPE3", "v")
            hdui.header.set("CUNIT3", "km/s")
            hdui.header.set("CRVAL3", self.vel_kms[0])
            hdui.header.set("CRPIX3", 0)
        hdul = fits.HDUList([primary_hdu, hdu_data])
        hdul.writeto(fn_fits)
        return

    def plotimage(self, ax, TA=True, ilam=0, **keys):
        """If TA=True, show in Antena temperature Kelvin units. Otherwise,
        show in orignal units of erg/s/cm^2/Hz/str"""
        if TA:
            imgi = self.image_TA[ilam]
            clabel = r"$T_A (\mathrm{K})$"
        else:
            imgi = self.image[ilam]
            clabel = r"$I_\nu (\mathrm{erg/s/cm^2/Hz/ster})$"
        cax = ax.imshow(
            np.swapaxes(imgi, 0, 1) + 1e-100,
            origin="lower",
            extent=[0, self.nx * self.dx / const.pc, 0, self.ny * self.dy / const.pc],
            cmap="magma",
            **keys,
        )
        cbar = ax.figure.colorbar(cax)
        cbar.solids.set_edgecolor("face")
        cbar.ax.set_ylabel(clabel)
        ax.set_xlabel("$x/\mathrm{pc}$")
        ax.set_ylabel("$y/\mathrm{pc}$")
        return

    def plotsed(self, ax, **keys):
        """plot velocity in km/s versus the average antena temprature over the
        whole image."""
        ax.plot(self.vel_kms, self.image_TA_avg, **keys)
        ax.set_xlabel(r"$v(\mathrm{km/s})$")
        ax.set_ylabel(r"$\langle T_A \rangle (\mathrm{K})$")
        return

    def getavgWCO(self, Tdect=0.6):
        """Calculate the WCO averaged over the whole image.
        Tdect: the detection limit in Kevin. Defaut 0.6 (K)."""
        # calculate the moment 0 image
        img = self.getimageWCO(Tdect=Tdect)
        return img.mean()

    def getimageWCO(self, Tdect=0.6, dv=0.07):
        """Return image in K km/s."""
        # get more finely spaced velocity array
        # spacing in velocity array
        vel0 = self.vel_kms
        vel_arr = np.arange(vel0[0], vel0[-1], dv)
        imageWCO = np.zeros((self.nx, self.ny))
        dv0 = np.diff(vel0)[0]
        # if small velocity channel, interpolation
        if dv <= dv0:
            # loop over each pixel
            for ix in range(self.nx):
                for iy in range(self.ny):
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    TA_dect = np.interp(vel_arr, vel0, TA_arr)
                    indx = TA_dect < Tdect
                    TA_dect[indx] = 0
                    # integration
                    imageWCO[ix, iy] = TA_dect.sum() * dv
        # if large velocity channel, need to calculate flux in each channel
        else:
            dv_fine = min(dv / 10.0, dv0)
            vel_arr_fine = np.arange(vel0[0], vel0[-1], dv_fine)
            for ix in range(self.nx):
                for iy in range(self.ny):
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    TA_dect_fine = np.interp(vel_arr_fine, vel0, TA_arr)
                    TA_dect = np.zeros(len(vel_arr))
                    for i in range(len(vel_arr)):
                        if i < len(vel_arr) - 1:
                            indx_channel = (vel_arr_fine >= vel_arr[i]) & (
                                vel_arr_fine < vel_arr[i + 1]
                            )
                        else:
                            indx_channel = vel_arr_fine >= vel_arr[i]
                        TA_dect[i] = TA_dect_fine[indx_channel].sum() * dv_fine / dv
                    indx = TA_dect < Tdect
                    TA_dect[indx] = 0
                    # integration
                    imageWCO[ix, iy] = TA_dect.sum() * dv
        return imageWCO

    def getTexc(self, Tdect=0.0, dv=None):
        Tmax_axis = self.getTpeak(Tdect=Tdect, dv=dv)
        Texc = 5.5 / np.log(1 + 5.5 / (Tmax_axis))
        return Texc

    def getTpeak(self, Tdect=0.0, dv=None):
        if dv is not None:
            # nv = self.nlam
            # else:
            vel0 = self.vel_kms
            dv0 = np.diff(vel0)[0]
            vel_arr = np.arange(vel0[0], vel0[-1], dv)
            # nv = len(vel_arr)
        Tmax = np.zeros((self.nx, self.ny))
        for ix in range(self.nx):
            for iy in range(self.ny):
                if dv is None:
                    ITA = np.array([m[ix, iy] for m in self.image_TA])
                elif dv < dv0:
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    ITA = np.interp(vel_arr, vel0, TA_arr)
                else:
                    vel_arr = np.arange(vel0[0], vel0[-1], dv)
                    dv_fine = min(dv / 10.0, dv0)
                    vel_arr_fine = np.arange(vel0[0], vel0[-1], dv_fine)
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    TA_dect_fine = np.interp(vel_arr_fine, vel0, TA_arr)
                    ITA = np.zeros(len(vel_arr))
                    for i in range(len(vel_arr)):
                        if i < len(vel_arr) - 1:
                            indx_channel = (vel_arr_fine >= vel_arr[i]) & (
                                vel_arr_fine < vel_arr[i + 1]
                            )
                        else:
                            indx_channel = vel_arr_fine >= vel_arr[i]
                        ITA[i] = TA_dect_fine[indx_channel].sum() * dv_fine / dv
                indx_nodect = ITA < Tdect
                ITA[indx_nodect] = 0.0
                Tmax[ix, iy] = ITA.max()
        return Tmax

    def plotimageWCO(
        self,
        ax,
        Tdect=0.6,
        extent=None,
        cax=None,
        orientation="vertical",
        cbticks=None,
        **keys,
    ):
        img = self.getimageWCO(Tdect=Tdect)
        if extent is None:
            extent = [
                0,
                self.nx * self.dx / const.pc / 1e3,
                0,
                self.ny * self.dy / const.pc / 1e3,
            ]
        cax1 = ax.imshow(
            np.swapaxes(img, 0, 1) + 1e-100, origin="lower", extent=extent, **keys
        )
        if cbticks is not None:
            cbar = ax.figure.colorbar(
                cax1, cax=cax, ticks=cbticks, orientation=orientation
            )
        else:
            cbar = ax.figure.colorbar(cax1, cax=cax, orientation=orientation)
        cbar.solids.set_edgecolor("face")
        cbar.set_label(r"$\mathrm{W_{CO}}(\mathrm{K\cdot km/s})$", fontsize=25)
        ax.set_xlabel("$x(\mathrm{kpc})$", fontsize=20)
        ax.set_ylabel("$y(\mathrm{kpc})$", fontsize=20)
        return

    def plotimageXCO(
        self,
        ax,
        imgNH2,
        Tdect=0.6,
        WCOmin=0.6,
        extent=None,
        cbticks=None,
        cax=None,
        orientation="vertical",
        **keys,
    ):
        """WCOmin: below which assign XCO=0 (non-detection.)"""
        img = self.getimageWCO(Tdect=Tdect)
        imgXCO = imgNH2 / 2.0e20 / img
        indx_WCO = img < WCOmin
        imgXCO[indx_WCO] = 0
        if extent is None:
            extent = [
                0,
                self.nx * self.dx / const.pc / 1e3,
                0,
                self.ny * self.dy / const.pc / 1e3,
            ]
        cax1 = ax.imshow(
            np.swapaxes(imgXCO, 0, 1), origin="lower", extent=extent, **keys
        )
        if cbticks is not None:
            cbar = ax.figure.colorbar(
                cax1, cax=cax, ticks=cbticks, orientation=orientation
            )
        else:
            cbar = ax.figure.colorbar(cax1, cax=cax, orientation=orientation)
        cbar.solids.set_edgecolor("face")
        cbar.set_label(
            r"$\mathrm{X_\mathrm{CO}}(\mathrm{2\times 10^{20}cm^{-2}K^{-1} km^{-1}s})$",
            fontsize=25,
        )
        ax.set_xlabel("$x(\mathrm{kpc})$", fontsize=20)
        ax.set_ylabel("$y(\mathrm{kpc})$", fontsize=20)
        return

    def getimage_velocities(self, Tdect=0.0, dv=None):
        """img_v_avg_kms: intensty weighted mean velocity, in unit of km/s.
        img_sigma_v: intensity weighted veolcity dispersion, in unit of km/s."""
        img_v_avg_kms = np.zeros((self.nx, self.ny))
        img_sigma_v_kms = np.zeros((self.nx, self.ny))
        vel0 = self.vel_kms
        dv0 = np.diff(vel0)[0]
        for ix in range(self.nx):
            for iy in range(self.ny):
                if dv is None:
                    ITA = np.array([m[ix, iy] for m in self.image_TA])
                    vel_arr = self.vel_kms
                elif dv < dv0:
                    # if small velocity channel, interpolation
                    vel_arr = np.arange(vel0[0], vel0[-1], dv)
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    ITA = np.interp(vel_arr, vel0, TA_arr)
                else:
                    # if large velocity channel,
                    # need to calculate flux in each channel
                    vel_arr = np.arange(vel0[0], vel0[-1], dv)
                    dv_fine = min(dv / 10.0, dv0)
                    vel_arr_fine = np.arange(vel0[0], vel0[-1], dv_fine)
                    TA_arr = np.array(
                        [
                            (self.image_TA[self.nlam - ilam - 1])[ix, iy]
                            for ilam in np.arange(self.nlam)
                        ]
                    )
                    TA_dect_fine = np.interp(vel_arr_fine, vel0, TA_arr)
                    ITA = np.zeros(len(vel_arr))
                    for i in range(len(vel_arr)):
                        if i < len(vel_arr) - 1:
                            indx_channel = (vel_arr_fine >= vel_arr[i]) & (
                                vel_arr_fine < vel_arr[i + 1]
                            )
                        else:
                            indx_channel = vel_arr_fine >= vel_arr[i]
                        ITA[i] = TA_dect_fine[indx_channel].sum() * dv_fine / dv
                indx_nodect = ITA < Tdect
                ITA[indx_nodect] = 0.0
                if ITA.sum() == 0.0:
                    img_v_avg_kms[ix, iy] = 0.0
                    img_sigma_v_kms[ix, iy] = 0.0
                else:
                    v_avg = np.average(vel_arr, weights=ITA)
                    vsq_avg = np.average(vel_arr**2, weights=ITA)
                    sigma_v = np.sqrt(vsq_avg - v_avg**2)
                    img_v_avg_kms[ix, iy] = v_avg
                    img_sigma_v_kms[ix, iy] = sigma_v
        return img_v_avg_kms, img_sigma_v_kms

    def rescale(self, factor, circular_beam=False, stamp=None):
        """Rescale image to resolution = orignal_resolution*factor.
        This is to see the effect of a larger observational beam."""
        # check whether nx and ny is dividable by factor
        if (self.nx % factor) != 0 or (self.ny % factor) != 0:
            print("ERROR: nx or ny not dividable by factor")
            return
        # make a copy of the class
        img_copy = copy.deepcopy(self)
        # rescale image
        img_copy.nx = int(self.nx / factor)
        img_copy.ny = int(self.ny / factor)
        img_copy.dx = self.dx * factor
        img_copy.dy = self.dy * factor
        img_copy.image = []
        img_copy.image_TA = []
        for i in range(self.nlam):
            im = self.image[i]
            im0 = im.reshape([img_copy.nx, factor, img_copy.ny, factor]).mean(3).mean(1)
            if circular_beam:
                im1 = mcb.stamp_avg(im0, stamp)
            else:
                im1 = im0
            img_copy.image.append(im1)
            nui = self.nu[i]
            fac_TA = const.c**2 / (2.0 * const.kb * nui**2)
            img_copy.image_TA.append(im1 * fac_TA)
        img_copy.image_avg = [im.mean() for im in img_copy.image]
        img_copy.image_TA_avg = [im.mean() for im in img_copy.image_TA]
        return img_copy


class LevelPop:
    def __init__(self, fn, dims=None):
        """dims = [ny, nx, nz], if not given, then 1D array, assuming ASCII
        output"""
        self.fn = fn
        f = open(fn, "r")
        _ = f.readline()  # skip first line
        self.ncells = int(f.readline())
        self.nlevels = int(f.readline())
        self.levels = [int(x) for x in f.readline().split()]
        if self.nlevels != len(self.levels):
            raise RuntimeError("File: " + fn + ", nlevels != len(levels)")
        self.levelpop = {}
        for ilevel in self.levels:
            if dims is None:
                self.levelpop[ilevel] = np.zeros(self.ncells)
            else:
                if dims[0] * dims[1] * dims[2] != self.ncells or len(dims) != 3:
                    raise RuntimeError(
                        "File: " + fn + ", dims is not consistent with ncells."
                    )
                self.levelpop[ilevel] = np.zeros(dims)
        if dims is None:
            for i in range(self.ncells):
                li = [float(x) for x in f.readline().split()]
                for ii, ilevel in zip(range(self.nlevels), self.levels):
                    (self.levelpop[ilevel])[i] = li[ii]
        else:
            for iz in range(dims[2]):
                for iy in range(dims[1]):
                    for ix in range(dims[0]):
                        li = [float(x) for x in f.readline().split()]
                        for ii, ilevel in zip(range(self.nlevels), self.levels):
                            (self.levelpop[ilevel])[ix, iy, iz] = li[ii]
        f.close()
        return

    def getTexc(self, Jl=1, Ju=2, Tul=5.5, gl=1, gu=3):
        """Return excitation temperature.
        Note that in RADMC, Ju is defined to start from 1.
        """
        nl = self.levelpop[Jl]
        nu = self.levelpop[Ju] + 1e-100
        Texc = Tul / np.log((nl / float(gl)) / (nu / float(gu)))
        return Texc


def bplanck(temp, nu):
    """Black body planck function B_nu(T), copied from radmc source code.
                   2 h nu^3 / c^2
       B_nu(T)  = ------------------    [ erg / cm^2 s ster Hz ]
                  exp(h nu / kT) - 1

    ARGUMENTS:
       temp  [K]             = Temperature
       nu    [Hz]            = Frequency
    """
    small_ = 1.0e-200
    # large_ = 709.78271
    if temp < small_:
        return 0.0
    xx = 4.7989e-11 * nu / temp
    ret = 1.47455e-47 * nu * nu * nu / (np.exp(xx) - 1.0) + small_
    return ret


def read_line_freq(fn_mol):
    """Read frequencies for molecular lines.
    input:
        fn_mol: filename for molecular data
    output:
        line_nu0: frequencies of transitions in Hz
    """
    f = open(fn_mol)
    lines = f.readlines()
    nlevel = int(lines[5].strip())
    ntrans = int(lines[8 + nlevel].strip())
    line_nu0 = np.zeros(ntrans)
    for i in np.arange(ntrans):
        indx = i + 10 + nlevel
        line_list = (lines[indx].strip()).split()
        line_nu0[i] = float(line_list[4]) * 1.0e9
    return line_nu0


def get_wavelength_micron(line_nu0, iline, widthkms, vkms, linenlam):
    nu0 = line_nu0[iline - 1]
    camera_nu = np.zeros(linenlam)
    for inu in np.arange(linenlam):
        camera_nu[inu] = nu0 * (
            1.0
            - vkms * const.km / const.c
            - (2.0 * inu / (linenlam - 1.0) - 1.0) * widthkms * const.km / const.c
        )
    wavelength_micron = 1.0e4 * const.c / camera_nu
    return wavelength_micron
