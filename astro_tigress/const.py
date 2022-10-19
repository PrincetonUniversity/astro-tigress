"""
Basic physical constants file. Units are mostly CGS.
"""

import math
import numpy as np
import scipy.integrate
import scipy.misc

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# CGS PHYSICAL CONSTANTS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

c = 2.99792458e10       # speed of light CGS
h = 6.6260755e-27       # Planck's constant CGS
g = 6.67259e-8          # Grav const CGS
kb = 1.380658e-16       # Boltzmann's const CGS
a = 7.56591e-15         # Radiation constant CGS
sb = 5.67051e-5         # sigma (stefan-boltzmann const) CGS
qe =  4.803206e-10      # Charge of electron CGS
ev =  1.60217733e-12    # Electron volt CGS
na =  6.0221367e23      # Avagadro's Number
me =  9.1093897e-28     # electron mass CGS
mp =  1.6726231e-24     # proton mass CGS
mn = 1.674929e-24       # neutron mass CGS
mh = 1.673534e-24       # hydrogen mass CGS
amu =  1.6605402e-24    # atomic mass unit CGS

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# ASTRONOMICAL CONSTANTS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

# GENERAL
au = 1.496e13           # astronomical unit CGS
pc = 3.0857e18          # parsec CGS
yr = 3.155815e7         # sidereal year CGS
Myr = 1.0e6 * yr        # million years CGS
ms = 1.98900e+33        # solar mass CGS
mj = 1.8986e30          # jupiter mass CGS
aj = 5.202545 * au           # jupiter semi-major axis CGS
rs = 6.9599e10          # sun's radius CGS
ls = 3.839e33           # sun's luminosity CGS
mm = 7.35000e+25        # moon mass CGS
mer = 5.97400e+27       # earth mass CGS
rer = 6.378e8           # earth's radius CGS
mmsn_er = 1700.0          # g/cm^2, MMSN at 1AU
medd = 3.60271e+34      # Eddington mass CGS
km = 1.0e5              # km in cm

# RADIO SPECIFIC
jy = 1.e-23                  # Jansky CGS
restfreq_hi = 1420405751.786 # 21cm transition (Hz)
restfreq_co = 115271201800.  # CO J=1-0 (Hz)
cm2perkkms_hi = 1.823e18     # HI column per intensity (thin)

# OTHER WAVELENGTHS
ksun = 3.28             # abs K mag of sun (Binney & Merrifield)

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# GEOMETRY
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

srdeg = 3283.0          # Degrees squared in a steradian
dtor =  0.0174532925199 # Degrees per radian
pi = 3.14159265359      # Pi

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# TELESCOPE STUFF
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

def hpbw(lam_cm=None, freq_hz=None, diam_m=None, units=True):
    """
    Return the half-power beam width for a telescope.
    """
    if (lam_cm == None and freq_hz == None) or \
            diam_m == None:
        hpbw = None
    else:
        if lam_cm == None:
            lam_cm = c/freq_hz
        hpbw = 1.2 * lam_cm / (diam_m*1e2)

    if units == True:
        return {"units":"sr",
                "value":hpbw}
    else:
        return hpbw

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# PLANCK FUNCTION CALCULATIONS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

def bb_inu(temp_k=None, lam_cm=None, units=True):
    """
    Planck function. Returns f_nu given T and lambda.
    """
    nu = c / lam_cm
    fnu = 2.0*h*nu**3/c**2 /(np.exp(h*nu/(kb*temp_k))-1.0)

    if units == True:
        return {"units":"erg/s/cm^2/sr/Hz",
                "value":fnu}
    else:
        return fnu

def bb_ilam(temp_k=None, lam_cm=None, units=True):
    """
    Planck function. Returns f_lambda given T and lambda.
    """
    flam = 2.0*h*c**2/lam_cm**5 / (np.exp(h*c/(lam_cm*kb*temp_k))-1.0)

    if units == True:
        return {"units":"erg/s/cm^2/sr/cm",
                "value":flam}
    else:
        return flam

### NOT FINISHED - FIX ###

def peak_inu(temp_k=None, units=True):
    """
    Wavelength for peak of F_nu given T. NOT FUNCTIONAL YET.
    """
    peak_lam = 2.821439372/kb/temp_k/h
    peak_nu = c / peak_lam

    if units == True:
        return {"units":"cm",
                "value":peak_lam}
    else:
        return peak_lam

### NOT FINISHED - FIX ###

def peak_ilam(temp_k=None, units=True):
    """
    Wavelength for peak of F_lambda given T. NOT FUNCTIONAL YET.
    """
    peak_lam = h*c/kb/temp_k/4.965114231
    if units == True:
        return {"units":"cm",
                "value":peak_lam}
    else:
        return peak_lam

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# TRANSLATE HMS/DMS STRINGS BAND AND FORTH FROM DEC. DEG
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

# CASA's QA tool should do this but sucks

def hms(val):
    """
    Decimal value to 3-element tuple of hour, min, sec. Just calls dms/15
    """
    return dms(val/15.0)

def dms(val):
    """
    Decimal value to 3-element tuple of deg, min, sec.
    """
    if type(val) == type(np.array([1,2,3])):
        val_abs = np.abs(val)
        d = np.floor(val_abs)
        m = np.floor((val_abs - d)*60.0)
        s = ((val_abs - d - m/60.0))*3600.0
        d *= +2.0*(val >= 0.0) - 1.0
    else:
        val_abs = abs(val)
        d = math.floor(val_abs)
        m = math.floor((val_abs - d)*60.0)
        s = ((val_abs - d - m/60.0))*3600.0
        d *= +2.0*(val >= 0.0) - 1.0
    return (d, m, s)

def ten(vec, hms=False):
    """
    Convert 3-element vector to decimal degrees. Use HMS for hours.
    """
    if type(val) == type(np.array([1,2,3])):
        dec = (np.abs(vec[0])+vec[1]/60.+vec[2]/3600.)
        dec *= +2.0*(vec[0] >= 0.0) - 1.0
    else:
        dec = (abs(vec[0])+vec[1]/60.+vec[2]/3600.)
        dec *= +2.0*(vec[0] >= 0.0) - 1.0
    if hms == True:
        return 15.0*dec
    else:
        return dec

def string_sixty(inp_str,hms=False):
    """
    Convert an HMS/DMS style string to a numeric vector. Again, QA
    should do this.
    """
    if hms==True:
        pass
    else:
        pass

def sixty_string(inp_val,hms=False,colons=False):
    """
    Convert a numeric vector to a string. Again, QA should do this.
    """
    if hms==True:
        out_str = "%02dh%02dm%05.3fs" % (inp_val[0], inp_val[1], inp_val[2])
    elif colons==True:
        out_str = "%+02d:%02d:%05.3f" % (inp_val[0], inp_val[1], inp_val[2])
    else:
        out_str = "%+02dd%02dm%05.3fs" % (inp_val[0], inp_val[1], inp_val[2])
    return out_str

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# LAPLACE COEFFICIENT AND ITS DIRIVATIVES CALCULATIONS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

def binte(phi, n, s, alpha):
    """the integration of Laplace coefficient"""
    return np.cos(n*phi) / (1 - 2*alpha*np.cos(phi) + alpha**2)**s

def laplace_b(n, s, alpha):
    """the Laplace coefficient b^n_s(alpha)"""
    return scipy.integrate.quad(binte, 0., 2.*math.pi, (n, s, alpha)
            )[0] / math.pi
def dlaplace_b_dalpha(n, s, alpha):
    def deri(alpha0, n0, s0):
        return laplace_b(n0, s0, alpha0)
    return scipy.misc.derivative(deri, alpha, dx=1.e-5, args=(n, s))

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# KEPLERIAN ORBIT COEFFICIENTS IN CGS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
def nkep(a, mstar=ms):
    """mean motion of orbit at semi-major axis a in CGS"""
    return math.sqrt(g * mstar / float(a)**3)
