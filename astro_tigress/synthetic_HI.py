import numpy as np
import tqdm

# save synthetic HI data to fits
from astropy.io import fits


def los_to_HI_axis_proj(dens, temp, vel, vchannel, deltas=1.0, los_axis=1):
    """
    inputs:
        dens: number density of hydrogen in units of 1/cm^3
        temp: temperature in units of K
        vel: line-of-sight velocity in units of km/s
        vchannel: velocity channel in km/s
    parameters:
        deltas: length of line segments in units of pc
        memlim: memory limit in GB
        los_axis: 0 -- z, 1 -- y, 2 -- x
    outputs: a dictionary
        TB: the brightness temperature
        tau: optical depth
    """

    # Nv=len(vchannel)

    Tlos = temp
    vlos = vel
    nlos = dens

    Tspin = Tlos

    ds = deltas * 3.085677581467192e18

    v_L = 0.21394414 * np.sqrt(Tlos)  # in units of km/s

    TB = []
    tau_v = []
    for vch in tqdm.tqdm(vchannel):
        phi_v = (
            0.00019827867
            / v_L
            * np.exp(-((1.6651092223153954 * (vch - vlos) / v_L) ** 2))
        )  # time
        kappa_v = 2.6137475e-15 * nlos / Tspin * phi_v  # area/volume = 1/length
        tau_los = kappa_v * ds  # dimensionless

        tau_cumul = tau_los.cumsum(axis=los_axis)

        # same unit with Tspin
        TB.append(
            np.nansum(
                Tspin * (1 - np.exp(-tau_los)) * np.exp(-tau_cumul), axis=los_axis
            )
        )
        # dimensionless
        tau_v.append(np.nansum(kappa_v * ds, axis=los_axis))
    return np.array(TB), np.array(tau_v)


def k10h(T2):
    """
    input: T/100K
    output: collisional excitation rate of hydrogen in c.g.s.
    """
    k10_1 = 1.19e-10 * T2 ** (0.74 - 0.20 * np.log(T2))
    k10_2 = 2.24e-10 * T2**0.207 * np.exp(-0.876 / T2)
    k10 = k10_1
    idx = k10_2 > k10_1
    k10[idx] = k10_2[idx]
    return k10


def k10e(T2):
    """
    input: T/100K
    output: collisional excitation rate of electrion in c.g.s.
    """
    Temp = T2 * 1.0e2
    k10 = -9.607 + 0.5 * np.log10(Temp) * np.exp(-((np.log10(Temp)) ** 4.5) / 1.8e3)
    return 10.0**k10


def create_fits(ytds):
    tMyr = ytds.current_time.to("Myr").value
    le = ytds.domain_left_edge.to("pc").value
    re = ytds.domain_right_edge.to("pc").value
    dx = (ytds.domain_width / ytds.domain_dimensions).to("pc").value
    hdr = fits.Header()
    hdr["time"] = float(tMyr)
    hdr["xmin"] = (le[0], "pc")
    hdr["xmax"] = (re[0], "pc")
    hdr["ymin"] = (le[1], "pc")
    hdr["ymax"] = (re[1], "pc")
    hdr["zmin"] = (le[2], "pc")
    hdr["zmax"] = (re[2], "pc")
    hdr["dx"] = (dx[0], "pc")
    hdr["dy"] = (dx[1], "pc")
    hdr["dz"] = (dx[2], "pc")
    hdu = fits.PrimaryHDU(header=hdr)
    return hdu


def add_header_for_glue(hdu, hdr, axis="xyz"):
    for i, ax in enumerate(axis):
        hdu.header["CDELT{}".format(i + 1)] = hdr["d{}".format(ax)]
        hdu.header["CTYPE{}".format(i + 1)] = ax
        hdu.header["CUNIT{}".format(i + 1)] = hdr.comments["d{}".format(ax)]
        hdu.header["CRVAL{}".format(i + 1)] = hdr["{}min".format(ax)]
        hdu.header["CRPIX{}".format(i + 1)] = (
            hdr["{}max".format(ax)] + hdr["{}min".format(ax)]
        )
    return


def save_to_fits(ytds, vchannel, TB, tau, axis="xyv", fitsname=None):
    """Function warpper for fits writing

    Parameters
    ----------
    ytds: yt.Dataset
        yt Dataset
    vchannel: array
        velocity channel array
    TB: array
        PPV data cube for brightness temperature
    tau: array
        PPV data cube for optical depth
    axis: str
        str for the PPV coordniate names
    fitsname: str
        if not None, fits file is created
    """
    hdul = fits.HDUList()
    hdu = create_fits(ytds)

    hdu.header["vmin"] = (vchannel.min(), "km/s")
    hdu.header["vmax"] = (vchannel.max(), "km/s")
    hdu.header["dv"] = (vchannel[1] - vchannel[0], "km/s")

    hdul.append(hdu)
    for fdata, label in zip([TB, tau], ["TB", "tau"]):
        hdul.append(fits.ImageHDU(name=label, data=fdata))

    hdr = hdu.header
    for hdu in hdul:
        add_header_for_glue(hdu, hdr, axis=axis)

    if fitsname is not None:
        hdul.writeto(fitsname, overwrite=True)
        print("FITS file named {} is created".format(fitsname))
    return hdul
