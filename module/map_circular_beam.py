import numpy as np
import math
from scipy.integrate import quad

def get_gauss_stamp(n):
    nx = n
    ny = n
    nxc = nx/2
    nyc = ny/2
    ix = np.zeros((nx, ny))
    iy = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            ix[i, j] = j-nyc
            iy[i, j] = i-nxc
    def gauss(x):
        sigmax = 0.5
        ret = 1./(np.sqrt(2*math.pi)*sigmax) * np.exp(-0.5*(x/sigmax)**2) 
        return ret
    stamp = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            x1 = ix[i, j] - 0.5
            x2 = x1 + 1.
            y1 = iy[i, j] - 0.5
            y2 = y1 + 1.
            inte = quad(gauss, x1, x2)[0] * quad(gauss, y1, y2)[0]
            stamp[i, j] = inte
    return stamp

#average data using a stamp
def stamp_avg(data, stamp):
    nxd, nyd = data.shape
    nxs, nys = stamp.shape
    ret = np.zeros(data.shape)
    nxc = nxs/2
    nyc = nys/2
    ix = np.zeros((nxs, nys))
    iy = np.zeros((nxs, nys))
    for i in range(nxs):
        for j in range(nys):
            ix[i, j] = i-nxc
            iy[i, j] = j-nyc
    for i in range(nxd):
        for j in range(nyd):
            for istamp in range(nxs):
                for jstamp in range(nys):
                    iret = i + ix[istamp, jstamp]
                    jret = j + iy[istamp, jstamp]
                    if (iret >= 0) and (iret < nxd
                        ) and (jret >= 0) and (jret < nyd):
                        ret[iret, jret] += data[i, j]*stamp[istamp, jstamp]
    return ret

def map_circular_beam(data, nstamp=9):
    stamp = get_gauss_stamp(nstamp)
    return stamp_avg(data, stamp)

