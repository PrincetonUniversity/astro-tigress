import yt
from yt.data_objects.level_sets.api import *
import yt.units as u
from yt import add_field
from yt import YTQuantity
from yt.utilities.physical_constants import \
    mh, \
    me, \
    sigma_thompson, \
    clight, \
    kboltz, \
    G
import numpy as np
import copy
import math
from yt.funcs import mylog
mylog.setLevel(50)

global Zdg, GxHe, GxC, GxO, GxSi
Gzdg = 1.
GxHe = 0.1
GxC = 1.6e-4 * Gzdg
GxO = 3.2e-4 * Gzdg
GxSi = 1.7e-6 * Gzdg

class coolftn(object):
    def __init__(self,fname='../module/coolftn.txt'):
        cdf=np.loadtxt(fname)
        self.cool=cdf[:, 0]
        self.heat=cdf[:, 1]
        self.temp=cdf[:, 2]
        self.T1=cdf[:, 3]
        self.Tmin=self.T1.min()
        self.Tmax=self.T1.max()
        self.dT=np.log10(self.T1[1]/self.T1[0])
        self.nT=len(self.T1)

    def get_Tidx(self,T):
        if type(T)==np.ndarray:
            Tidx=np.log10(T/self.Tmin)/self.dT
            Tidx[np.where(T<self.Tmin)]=0
            Tidx[np.where(T>=self.Tmax)]=self.nT-2
            return Tidx.astype(int)
        else:
            if T < self.Tmin: return 0
            if T >= self.Tmax: return self.nT-2
            Tidx=np.log10(T/self.Tmin)/self.dT
            return int(Tidx)

    def get_temp(self,T1):
        T1idx=self.get_Tidx(T1)
        Ti=self.temp[T1idx]
        Tip1=self.temp[T1idx+1]
        T1i=self.T1[T1idx]
        T1ip1=self.T1[T1idx+1]
        T=Ti+(Tip1-Ti)*(T1-T1i)/(T1ip1-T1i)

        return T

    def get_cool(self,T1):
        T1idx=self.get_Tidx(T1)
        Li=self.cool[T1idx]
        Lip1=self.cool[T1idx+1]
        T1i=self.T1[T1idx]
        T1ip1=self.T1[T1idx+1]
        L=Li+(Lip1-Li)*(T1-T1i)/(T1ip1-T1i)

        return L
    
    def get_heat(self,T1):
        T1idx=self.get_Tidx(T1)
        Gi=self.heat[T1idx]
        Gip1=self.heat[T1idx+1]
        T1i=self.T1[T1idx]
        T1ip1=self.T1[T1idx+1]
        G=Gi+(Gip1-Gi)*(T1-T1i)/(T1ip1-T1i)

        return G

# basic qunatities with renormalization
def _ndensity(field, data):
        return data["gas","density"]/(1.4271*mh)

def _ram_pok_z(field,data):
        return data["gas","density"]*data["gas","velocity_z"]**2/kboltz

# thermodynamics quantities
def _pok(field, data):
        return data["gas","pressure"]/kboltz

def _cs(field, data):
        return np.sqrt(data["gas","pressure"]/data["gas","density"])

def _T1(field, data):
        return data["gas","pressure"]/data["gas","nH"]/kboltz

def _mu(field, data):
    cf=coolftn()
    T1=data["gas","T1"].d
    temp=cf.get_temp(T1)
    return temp/T1

def _temperature(field,data):
    return data["gas","T1"]*data["gas","mu_cgk"]

# rotation
Omega=YTQuantity(28,"km/s/kpc")
def _dvelocity(field,data):
        return data["gas","velocity_y"]+data["gas","x"]*Omega

def _dvelocity_mag(field,data):
        return np.sqrt(data["gas","velocity_x"]**2+data["gas","dvelocity_y"]**2+data["gas","velocity_z"]**2)

def _dkinetic_energy(field,data):
    return 0.5*data['gas','dvelocity_magnitude']**2*data['gas','density']

# magnetic fields
def _mag_pok(field,data):
        return data["gas","magnetic_pressure"]/kboltz

#chemistry
def _Si(field, data):
    return GxSi - data["Si+"]

def _C(field, data):
    return GxC - data["CHx"] - data["CO"] - data["C+"] - data["HCO+"]

def _H(field, data):
    return ( 1. - (data["H2"]*2 + 3*data["H3+"] +
                        2*data["H2+"] + data["H+"] + data["HCO+"] +
                        data["CHx"] + data["OHx"] )  )

def _O(field, data):
    return GxO - data["OHx"] - data["CO"] - data["HCO+"]

def _He(field, data):
    return GxHe - data["He+"]

def _CO2Ctot(field, data):
    return data["CO"]/GxC

def _C2Ctot(field, data):
    return data["C"]/GxC

def _Cp2Ctot(field, data):
    return data["C+"]/GxC

def _tH2(field, data):
    return data["H2"]*2

def _electron(field, data):
    return  (data["He+"] + data["C+"] + data["HCO+"] + 
                            data["H+"] + data["H3+"] + data["H2+"])

def _temperature_chem(field, data):
    kb = 1.381e-16
    return yt.units.K * data["E"] / (1.5 * kb *  
                    (1. - data["H2"] + data["e"] + GxHe))

def _pressure(field, data):
    ds = data.ds
    pressure_unit = ds.mass_unit / (ds.length_unit * (ds.time_unit)**2)
    return data["press"]*pressure_unit

def _density(field, data):
    return data["rho"]


import matplotlib.pyplot as plt

def add_yt_fields(ds,chemistry=True, cooling=True,mhd=False,rotation=False,
                  zdg=1.):
    Gzdg = zdg
    GxHe = 0.1
    GxC = 1.6e-4 * Gzdg
    GxO = 3.2e-4 * Gzdg
    GxSi = 1.7e-6 * Gzdg
    ds.add_field(("gas","nH"),function=_ndensity, \
      units='cm**(-3)',display_name=r'$n_{\rm H}$', sampling_type="cell")
    ds.add_field(("gas","ram_pok_z"),function=_ram_pok_z, \
      units='K*cm**(-3)',display_name=r'$P_{\rm turb}/k_{\rm B}$',
      sampling_type="cell")
    ds.add_field(("gas", "density"), function=_density, display_name="density",
            sampling_type="cell")
    ds.add_field(("gas", "pressure"), function=_pressure, units="g/(cm* s**2)",
            display_name="pressure", sampling_type="cell")
    if chemistry:
        ds.add_field( ("gas", "Si"), function=_Si, 
                units="", display_name=r"${\rm Si}$", sampling_type="cell")
        ds.add_field( ("gas", "C"), function=_C, 
                units="", display_name=r"${\rm C}$", sampling_type="cell")
        ds.add_field( ("gas", "H"), function=_H, 
                units="", display_name=r"${\rm H}$", sampling_type="cell")
        ds.add_field( ("gas", "O"), function=_O, 
                units="", display_name=r"${\rm O}$", sampling_type="cell")
        ds.add_field( ("gas", "He"), function=_He, 
                units="", display_name=r"${\rm He}$", sampling_type="cell")
        ds.add_field( ("gas", "CO2Ctot"), function=_CO2Ctot, 
                units="", display_name=r"${\rm CO}/{\rm C}_{\rm tot}$",
                sampling_type="cell")
        ds.add_field( ("gas", "C2Ctot"), function=_C2Ctot, 
                units="", display_name=r"${\rm C}/{\rm C}_{\rm tot}$",
                sampling_type="cell")
        ds.add_field( ("gas", "C+2Ctot"), function=_Cp2Ctot, 
                units="", display_name=r"${\rm C^+}/{\rm C}_{\rm tot}$", 
                sampling_type="cell")
        ds.add_field( ("gas", "2H2"), function=_tH2, 
                units="", display_name=r"${\rm 2H_2}$", sampling_type="cell")
        ds.add_field( ("gas", "e"), function=_electron, 
                units="", display_name=r"${\rm e^{-}}$", sampling_type="cell")
        ds.add_field( ("gas", "temperature_chem"), function=_temperature_chem, 
                units="K", display_name=r"$T$", sampling_type="cell")
    if cooling:
        ds.add_field(("gas","pok"),function=_pok, \
          units='K*cm**(-3)',display_name=r'$P/k_{\rm B}$',
          sampling_type="cell")
        ds.add_field(("gas","cs"),function=_cs, \
          units='km*s**(-1)',display_name=r'$c_s$', sampling_type="cell")
        ds.add_field(("gas","T1"),function=_T1, \
          units='K',display_name=r'$T_1$', sampling_type="cell")
        ds.add_field(("gas","mu_cgk"),function=_mu, \
          units='',display_name=r'$\mu_{\rm cgk}$',force_override=True,
          sampling_type="cell")
        ds.add_field(("gas","temperature_init"),function=_temperature, \
          units='K',display_name=r'$T_{\rm init}$',force_override=True,
          sampling_type="cell")
    if rotation:
        ds.add_field(("gas","dvelocity_y"),function=_dvelocity, \
          units='km/s',display_name=r'$\delta v_y$',force_override=True,
          sampling_type="cell")
        ds.add_field(("gas","dvelocity_magnitude"),function=_dvelocity_mag, \
          units='km/s',display_name=r'$v$',force_override=True,
          sampling_type="cell")
        ds.add_field(("gas","dkinetic_energy"),function=_dkinetic_energy, \
          units='erg/cm**3',display_name=r'$E_k$',force_override=True,
          sampling_type="cell")
    if mhd:
        ds.add_field(("gas","mag_pok"),function=_mag_pok, \
          units='K*cm**(-3)',display_name=r'$P_{\rm mag}/k_{\rm B}$',
          sampling_type="cell")

def set_aux(model='solar'):
    aux={}
    aux['nH']=dict(label=r'$n_H\;[{\rm cm}^{-3}]$', \
        unit='cm**(-3)', limits=(1.e-6,1.e6), \
        cmap=plt.cm.Spectral_r,clim=(2.e-5,2.e2), \
        cticks=(1.e-4,1.e-2,1,1.e2), \
        n_bins=128, log=True)
    aux['pok']=dict(label=r'$P/k_B\;[{\rm K}\,{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=plt.cm.gnuplot2,clim=(10,5.e5), \
        n_bins=128, log=True)
    aux['temperature']=dict(label=r'$T\;[{\rm K}]$', \
        unit='K', limits=(1.e0,1.e9), \
        cmap=plt.cm.RdYlBu_r, \
        clim=(10,1.e8), \
        n_bins=128, log=True)
    aux['surface_density']=dict( \
        label=r'$\Sigma [{\rm M}_{\odot} {\rm pc}^{-2}]$', \
        cmap=plt.cm.pink_r,clim=(0.1,100),log=True)
    aux['dvelocity_magnitude']=dict(label=r'$v [{\rm km/s}]$', \
        unit='km/s', limits=(0.1,1.e4), \
        cmap=plt.cm.jet,clim=(1,1000), \
        n_bins=128, log=True)
    aux['velocity_z']=dict(label=r'$v_z [{\rm km/s}]$', \
        unit='km/s', limits=(-1500,1500), \
        cmap=plt.cm.RdBu_r,clim=(-200,200), \
        cticks=(-100,0,100), \
        n_bins=256, log=False)
    aux['magnetic_field_strength']=dict(label=r'$B [\mu{\rm G}]$', \
        unit='gauss', \
        cmap=plt.cm.viridis,clim=(0.01,10),factor=1.e6, \
        n_bins=128, log=True)
    aux['mag_pok']=dict(label=r'$P_{\rm mag}/k_B [{\rm K}{\rm cm}^{-3}]$',\
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        n_bins=128, log=True)
    aux['ram_pok_z']=dict(\
        label=r'$P_{\rm turb}/k_B [{\rm K}{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        n_bins=128, log=True)
    aux['plasma_beta']=dict(label=r'$\beta$', limits=(1.e-4,1.e16), \
        n_bins=256, log=True)

    if model is 'starburst':
        aux['nH']['clim']=(2.e-5,2.e3)
        aux['pok']['clim']=(10,5.e6)
        aux['surface_density']['clim']=(1,1000)
        aux['magnetic_field_strength']['clim']=(0.1,100)

    if model is 'multi_SN':
        aux['nH']['clim']=(2.e-5,2.e2)
        aux['pok']['clim']=(50,1.e5)
    return aux

def check_aux(fields):
    aux=set_aux()
    for f in fields:
        if f not in aux:
            print("auxiliary information for %s is missing",f)
            print(aux[f])


unit_base_pp={"length_unit": (1.0,"pc"), 
           "time_unit": (1.0,"s*pc/km"), 
           "mass_unit": (2.38858753789e-24,"g/cm**3*pc**3")}
unit_base={"length_unit": (1.0,"pc"), 
           "time_unit": (1.0,"s*pc/km"), 
           "mass_unit": (2.38858753789e-24,"g/cm**3*pc**3"),
           "velocity_unit": (1.0,"km/s"),
           "magnetic_unit": (5.4786746797e-07,"gauss")}

def load_athena4p2(filename):
    ds = yt.load(filename, units_override=unit_base)
    add_yt_fields(ds, chemistry=False)
    return ds

def load_athenapp(filename):
    ds = yt.load(filename, units_override=unit_base_pp)
    add_yt_fields(ds)
    return ds

def save_clumps(clumps_list, filename, fields_list=None):
    """The list of clumps must belong to the same base object."""
    base_object = clumps_list[0].data.base_object
    if fields_list == None:
        fl = flist=[ 
                ('gas', 'density'),
                ('gas', 'magnetic_field_x'),
                ('gas', 'magnetic_field_y'),
                ('gas', 'magnetic_field_z'),
                ('gas', 'pressure'),
                ('gas', 'velocity_x'),
                ('gas', 'velocity_y'),
                ('gas', 'velocity_z'),
                ('gravitational_potential'),
                ('gas', 'phi'),
                ('cell_mass'),
                ('cell_volume')
                ]
    else:
        fl = copy.deepcopy(fields_list)
    conditional_output = []
    #generate fields list
    for c in clumps_list:
        contour_id = c.data.conditionals[0].split('contours')[1].split("'")[0]
        conditional_i = c.data.conditionals[0]
        fl.append(  ("index", "contours"+contour_id) )
        conditional_output.append( ("contours"+contour_id, conditional_i)  )
    #save to file
    base_object.save_as_dataset(fields=fl, filename=filename)
    fo = open(filename+".cond", "w")
    for cond in conditional_output:
        fo.write("{:s},{:s}\n".format(cond[0], cond[1]))
    fo.close()
    return

def read_clumps(data_set, fn_cond=None):
    if fn_cond == None:
        fn_cond = data_set.filename_template + ".cond"
    clumps = []
    fo = open(fn_cond)
    for line in fo.readlines():
        ls = line.split(',')
        cli = ls[0]
        condi = ls[1].split('\n')[0]
        print(condi)
        cut = data_set.cut_region(data_set.all_data(), 
        [condi])
        c1 = Clump(cut, ('grid', 'phi') )
        clumps.append(c1)
    return clumps 

def get_clump_properties(cl):
    vol=cl['cell_volume'][0]
    totvol=cl.quantities.total_quantity('cell_volume').in_units('pc**3')
    radius = (totvol*3.*math.pi/4.)**(1./3.)
    mass=cl.quantities.total_quantity('cell_mass').in_units('Msun')
    nden=cl.quantities.weighted_average_quantity('density','cell_mass')
    vx_c = cl.quantities.weighted_average_quantity("velocity_x", "cell_mass").in_units("km/s")
    vy_c= cl.quantities.weighted_average_quantity("velocity_y", "cell_mass").in_units("km/s")
    vz_c = cl.quantities.weighted_average_quantity("velocity_z", "cell_mass").in_units("km/s")  
    vc = [vx_c, vy_c, vz_c]
    Ekin_cells = 0.5*cl["cell_mass"]*( (cl["velocity_x"]-vx_c)**2 + (cl["velocity_y"] - vy_c)**2 +
    (cl["velocity_z"] -vz_c)**2 )
    Ekin = Ekin_cells.sum()
    Emag=cl.quantities.total_quantity('magnetic_energy')*vol.in_cgs()
    Eth=cl.quantities.total_quantity('pressure')*vol.in_cgs()*1.5
    Egrav=-( cl["gravitational_potential"].max()-cl["gravitational_potential"].min() )*(u.km/u.s)**2 * mass
    B=yt.YTArray([cl.quantities.weighted_average_quantity('magnetic_field_x','cell_mass'),
                  cl.quantities.weighted_average_quantity('magnetic_field_y','cell_mass'),
                  cl.quantities.weighted_average_quantity('magnetic_field_z','cell_mass')])
    Eratio = -((Ekin+Emag+Eth)/Egrav)
    print('==========================================')
    print('Volume:',totvol)
    print('Effective radius:', radius)
    print('Mass:',mass)
    print('center of mass velocity', vc)
    print('Mass-weighted mean Density:',nden/u.mh)
    print('Mass-weighted mean B (micro G):',B.in_units('uG'))
    print('Ekin:',Ekin.in_units('erg'))
    print('Emag:',Emag.in_units('erg'))
    print('Eth:',Eth.in_units('erg'))
    print('Egrav:',Egrav.in_units('erg'))
    print('Eratio = (Ekin+Emag+Eth)/(-Egrav):', Eratio.in_cgs())
    print('==========================================')
    return {"volume": totvol, "radius": radius, "mtot": mass, "vc": vc,
            "density_mass_weighted": nden, "B_mass_weigthed": B, "Ekin": Ekin, "Emag": Emag, "Eth": Eth,
            "Egrav": Egrav, "Eratio": Eratio}

def get_covering_grid(ds, left_edge=None, dims=None):
    """get the covering grid of the whole domain."""
    if left_edge is None:
        left_edge = ds.domain_left_edge
    if dims is None:
        dims = ds.domain_dimensions
    return ds.covering_grid(level=0, left_edge=left_edge,
            dims=dims)

def get_extent(grid, axis, unit):
    """unit: use yt.units"""
    left_edge = grid.left_edge
    right_edge = grid.right_edge
    if axis == 2:
        xmin = (left_edge[0]/unit).value
        xmax = (right_edge[0]/unit).value
        ymin = (left_edge[1]/unit).value
        ymax = (right_edge[1]/unit).value
    elif axis == 1:
        xmin = (left_edge[0]/unit).value
        xmax = (right_edge[0]/unit).value
        ymin = (left_edge[2]/unit).value
        ymax = (right_edge[2]/unit).value
    elif axis == 0:
        xmin = (left_edge[1]/unit).value
        xmax = (right_edge[1]/unit).value
        ymin = (left_edge[2]/unit).value
        ymax = (right_edge[2]/unit).value
    extent = [float(xmin), float(xmax), float(ymin), float(ymax)]
    return extent

def get_surfd(grid, axis, unit=None):
    """unit: None: in cm^-2, otherwise, use yt units"""
    dx = grid["dx"][0, 0, 0]
    surfd = grid["density"].sum(axis=axis)*dx
    if unit is None:
        surfd_cgs = surfd.in_units("g/cm**2").to_ndarray()
        surfd_ret = surfd_cgs / (1.4 * const.mh)
    else:
        surfd_ret = surfd.in_units(unit).to_ndarray()
    return surfd_ret


