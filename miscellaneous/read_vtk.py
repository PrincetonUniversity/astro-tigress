"""
Read athena vtk file
"""

from __future__ import print_function

import os
import os.path as osp
import glob
import struct
import numpy as np
import xarray as xr

from units import Units


def read_vtk(filename, id0_only=False):
    """Convenience wrapper function to read Athena vtk output file
    using AthenaDataSet class.

    Parameters
    ----------
    filename : str
        Name of the file to open, including extension
    id0_only : bool
        Flag to enforce to read vtk file in id0 directory only.
        Default value is False.

    Returns
    -------
    ds : AthenaDataSet
    """

    return AthenaDataSet(filename, id0_only=id0_only)


class AthenaDataSet(object):
    def __init__(self, filename, id0_only=False, units=Units(), dfi=None):
        """Class to read athena vtk file.

        Parameters
        ----------
        filename : string
            Name of the file to open, including extension
        id0_only : bool
            Flag to enforce to read vtk file in id0 directory only.
            Default value is False.
        units : Units
            pyathena Units object (used for reading derived fields)
        dfi : dict
            Dictionary containing derived fields info
        """

        if not osp.exists(filename):
            raise IOError(("File does not exist: {0:s}".format(filename)))

        dirname, problem_id, num, suffix, ext, mpi_mode, nonzero_id = _parse_filename(
            filename
        )

        if id0_only:
            mpi_mode = False

        self.dirname = dirname
        self.problem_id = problem_id
        self.num = int(num)
        self.suffix = suffix
        self.ext = ext
        self.mpi_mode = mpi_mode
        self.fnames = [filename]
        self.u = units
        self.dfi = dfi
        if dfi is not None:
            self.derived_field_list = list(dfi.keys())
        else:
            self.derived_field_list = None

        # Find all vtk file names and add to flist
        if mpi_mode:
            if self.suffix is None:
                fname_pattern = osp.join(
                    dirname, "id*/{0:s}-id*.{1:s}.{2:s}".format(problem_id, num, ext)
                )
            else:
                fname_pattern = osp.join(
                    dirname,
                    "id*/{0:s}-id*.{1:s}.{2:s}.{3:s}".format(
                        problem_id, num, suffix, ext
                    ),
                )
            fnames = glob.glob(fname_pattern)
            self.fnames += fnames

        if nonzero_id:
            from collections import OrderedDict

            self.fnames = list(OrderedDict.fromkeys(self.fnames))

        self.grid = self._set_grid()
        self.domain = self._set_domain()
        self.set_region()

        # Need separte field_map for different grids
        if self.domain["all_grid_equal"]:
            self._field_map = _set_field_map(self.grid[0])
            for g in self.grid:
                g["field_map"] = self._field_map
        else:
            for g in self.grid:
                g["field_map"] = _set_field_map(g)
            self._field_map = self.grid[0]["field_map"]

        self.field_list = list(self._field_map.keys())

    def get_cc_pos(self):
        """Compute cell center positions

        Returns
        -------
        xc : numpy array
            Unique cell-centered cooridnates
        """

        xc = dict()
        for axis, le, re, dx in zip(
            ("x", "y", "z"), self.region["gle"], self.region["gre"], self.domain["dx"]
        ):
            xc[axis] = np.arange(le + 0.5 * dx, re + 0.5 * dx, dx)

        return xc

    def get_cc_ijk(s, x1, x2, x3):
        """Compute closest cell-center integer indices in which a particle resides

        Parameters
        ----------
        x1,x2,x3 : array of floats
            Particle position
        """
        domain = s.domain
        le1, le2, le3 = domain["le"]
        dx1, dx2, dx3 = domain["dx"]

        return (
            np.floor((x1 - le1) / dx1).astype(int),
            np.floor((x2 - le2) / dx2).astype(int),
            np.floor((x3 - le3) / dx3).astype(int),
        )

    def set_region(self, le=None, re=None):
        """Set region and find overlapping grids."""

        if le is None:
            le = self.domain["le"]
        if re is None:
            re = self.domain["re"]

        le = np.array(le)
        re = np.array(re)

        if (re < le).any():
            raise ValueError("Check left/right edge.")

        # Find all overlapping grids and their edges
        gle_all = []  # grid left edge
        gre_all = []  # grid right edge
        gidx = []  # grid indices that belongs to this region
        # print(self.grid,len(self.grid))
        for i, g in enumerate(self.grid):
            if (g["re"] >= le).all() and (g["le"] <= re).all():
                gidx.append(i)
                gle_all.append(g["le"])
                gre_all.append(g["re"])

        gidx = np.array(gidx)
        if len(gidx) == 0:
            raise ValueError(
                "Check left/right edges:",
                le,
                re,
                " Domain left/right edges are ",
                self.domain["le"],
                self.domain["re"],
            )

        gle_all = np.array(gle_all)
        gre_all = np.array(gre_all)

        # Find unique grid left/right edge coordinates
        gleu = [np.unique(gle_all[:, i]) for i in range(3)]
        greu = [np.unique(gre_all[:, i]) for i in range(3)]

        # Min/Max of gleu/greu
        gle = np.array([gle.min() for gle in gleu])
        gre = np.array([gre.max() for gre in greu])

        # Number of grids in each direction
        NGrid = np.array([len(gleu_) for gleu_ in gleu])

        # Number of cells per grid
        Nxg = np.concatenate(
            [
                (greu_ - gleu_) / dx
                for greu_, gleu_, dx in zip(greu, gleu, self.domain["dx"])
            ]
        )
        Nxg = np.array(np.array_split(Nxg, NGrid.cumsum()[:-1]), dtype=object).tolist()

        # Since floating point arithmetic may result in incorrect results, need
        # to round to the nearest integer
        for i, Nxg_ in enumerate(Nxg):
            Nxg[i] = np.rint(Nxg_).astype(int)

        # Number of cells in region
        Nxr = np.empty(len(Nxg), dtype=int)
        for i, Nxg_ in enumerate(Nxg):
            Nxr[i] = np.sum(Nxg_)

        # print(gidx,NGrid)
        assert len(gidx) == NGrid.prod(), print(
            "Unexpected error: Number of grids {0:d} != ".format(len(gidx))
            + "number of unique edges {0:d}.".format(NGrid.prod())
        )

        self.region = dict(
            le=le,
            re=re,
            gidx=gidx,
            gleu=gleu,
            greu=greu,
            gle=gle,
            gre=gre,
            NGrid=NGrid,
            Nxg=Nxg,
            Nxr=Nxr,
        )

    def get_slice(self, axis, field="density", pos="c", method="nearest"):
        """Read slice of fields.

        Parameters
        ----------
        axis : str
            Axis to slice along. 'x' or 'y' or 'z'
        field : (list of) str
            The name of the field(s) to be read.
        pos : float or str
            Slice through If 'c' or 'center', get a slice through the domain
            center. Default value is 'c'.
        method : str


        Returns
        -------
        slc : xarray dataset
            An xarray dataset containing slices.
        """

        axis_idx = dict(x=0, y=1, z=2)

        if pos is None:
            pos = "c"

        field = np.atleast_1d(field)
        axis = np.atleast_1d(axis)

        for ax in axis:
            le = np.copy(self.domain["le"])
            re = np.copy(self.domain["re"])
            if pos in ["c", "center"]:
                pos = self.domain["center"][axis_idx[ax]]

            # Let's make sure le < re always and truncation error does not cause
            # problem, although performance can be slowed down a bit.
            le[axis_idx[ax]] = pos - 0.5 * self.domain["dx"][axis_idx[ax]]
            re[axis_idx[ax]] = pos + 0.5 * self.domain["dx"][axis_idx[ax]]
            dat = self.get_field(field, le, re, as_xarray=True)

            slc = dat.sel(method="nearest", **{ax: pos})

        return slc

    def get_field(self, field="density", le=None, re=None, as_xarray=True):
        """Read 3d fields data.

        Parameters
        ----------
        field : (list of) string
            The name of the field(s) to be read.
        le : sequence of floats
           Left edge. Default value is the domain left edge.
        re : sequence of floats
           Right edge. Default value is the domain right edge.
        as_xarray : bool
           If True, returns results as an xarray Dataset. If False, returns a
           dictionary containing numpy arrays. Default value is True.

        Returns
        -------
        dat : xarray dataset
            An xarray dataset containing fields.
        """

        field = np.atleast_1d(field)

        # Derived field list
        dflist = set(field) - set(self.field_list)

        if not bool(dflist):
            # dflist is an empty set, we can read all fields directly from vtk
            # file
            return self._get_field(field, le, re, as_xarray)

        # If we are here, need to read all union of all input fields and those
        # required to calculate derived fields

        # Let's first make sure that we have all info about dflist
        if not dflist.issubset(set(self.dfi.keys())):
            tmp = []
            for f in dflist:
                if f not in self.dfi.keys():
                    tmp.append(f)

            raise KeyError("Unrecognized field name(s):", tmp)

        # Field names that are in the vtk file
        flist = set(field) - dflist

        # Fields that need to be read to calculate derived field
        flist_dep = set()
        for f in dflist:
            flist_dep = flist_dep | set(self.dfi[f]["field_dep"])

        # Fields names to be dropped later
        fdrop_list = flist_dep - flist

        # Need to adjust names for vector fields
        for f in list(fdrop_list):
            if self._field_map[f]["nvar"] > 1:
                for i in range(self._field_map[f]["nvar"]):
                    fdrop_list.add(f + str(i + 1))
                fdrop_list.remove(f)

        field = list(flist_dep | flist)
        dat = self._get_field(field, le, re, as_xarray)

        # Calculate derived fields
        for f in dflist:
            dat[f] = self.dfi[f]["func"](dat, self.u)

        # Drop fields that are not requested
        if as_xarray:
            dat = dat.drop(list(fdrop_list))
            dat.attrs["dfi"] = self.dfi
        else:
            for f in fdrop_list:
                del dat[f]

        return dat.squeeze()

    def _get_field(self, field="density", le=None, re=None, as_xarray=True):
        field = np.atleast_1d(field)

        # Create region
        self.set_region(le=le, re=re)
        arr = self._get_array(field)

        # Works only for 3d data
        if as_xarray:
            # Cell center positions
            x = dict()
            for axis, le, re, dx in zip(
                ("x", "y", "z"),
                self.region["gle"],
                self.region["gre"],
                self.domain["dx"],
            ):
                # May not result in correct number of elements due to truncation error
                # x[axis] = np.arange(le + 0.5*dx, re + 0.5*dx, dx)
                x[axis] = np.arange(le + 0.5 * dx, re + 0.25 * dx, dx)

            dat = dict()
            for k, v in arr.items():
                if len(v.shape) > self.domain["ndim"]:
                    for i in range(v.shape[-1]):
                        dat[k + str(i + 1)] = (("z", "y", "x"), v[..., i])
                else:
                    dat[k] = (("z", "y", "x"), v)

            attrs = dict()
            for k, v in self.domain.items():
                attrs[k] = v
            return xr.Dataset(dat, coords=x, attrs=attrs)
        else:
            if len(field) == 1:
                return arr[field[0]]
            else:
                # Return a dictionary of numpy arrays
                return arr

    def _get_array(self, field):
        arr = dict()
        for f in field:
            arr[f] = self._set_array(f)

        # Read from individual grids and copy to data
        le = self.region["gle"]
        dx = self.domain["dx"]
        for i in self.region["gidx"]:
            g = self.grid[i]
            il = (np.rint((g["le"] - le) / dx)).astype(int)
            iu = il + g["Nx"]
            slc = tuple([slice(l_, u_) for l_, u_ in zip(il[::-1], iu[::-1])])

            for f in field:
                arr[f][slc] = self._read_array(g, f)

        return arr

    def _read_array(self, grid, field):
        if field in grid["data"]:
            return grid["data"][field]
        elif field in self.field_list:
            fm = grid["field_map"][field]
            if "tarinfo" in grid:
                fp = self.tarfile.extractfile(grid["tarinfo"])
            else:
                fp = open(grid["filename"], "rb")
            fp.seek(fm["offset"])
            fp.readline()  # skip header
            if fm["read_table"]:
                fp.readline()
            grid["data"][field] = np.asarray(
                struct.unpack(">" + fm["ndata"] * fm["dtype"], fp.read(fm["dsize"]))
            )
            fp.close()

            if fm["nvar"] == 1:
                shape = np.flipud(grid["Nx"])
            else:
                shape = (*np.flipud(grid["Nx"]), fm["nvar"])

            grid["data"][field].shape = shape

            return grid["data"][field]

    def _set_array(self, field):
        dtype = self._field_map[field]["dtype"]
        nvar = self._field_map[field]["nvar"]
        Nxr = self.region["Nxr"]
        if "face_centered_B" in field:
            Nxr[int(field[-1]) - 1] += 1
        if nvar == 1:
            shape = np.flipud(Nxr)
        else:
            shape = (*np.flipud(Nxr), nvar)

        return np.empty(shape, dtype=dtype)

    def _set_domain(self):
        domain = dict()
        grid = self.grid
        ngrid = len(grid)

        # Grid left/right edges
        gle = np.empty((ngrid, 3), dtype="float32")
        gre = np.empty((ngrid, 3), dtype="float32")
        dx = np.empty((ngrid, 3), dtype="float32")
        Nx = np.ones_like(dx, dtype="int")

        for i, g in enumerate(grid):
            gle[i, :] = g["le"]
            gre[i, :] = g["re"]
            Nx[i, :] = g["Nx"]
            dx[i, :] = g["dx"]

        # Check if all grids have the equal size
        if (Nx[0] == Nx).all():
            domain["all_grid_equal"] = True
        else:
            domain["all_grid_equal"] = False

        # Set domain
        le = gle.min(axis=0)
        re = gre.max(axis=0)
        domain["ngrid"] = ngrid
        domain["le"] = le
        domain["re"] = re
        domain["dx"] = dx[0, :]
        domain["Lx"] = re - le
        domain["center"] = 0.5 * (le + re)
        domain["Nx"] = np.round(domain["Lx"] / domain["dx"]).astype("int")
        domain["ndim"] = 3  # always 3d

        # file = open(self.fnames[0], 'rb')
        # tmpgrid = dict()
        # tmpgrid['time'] = None
        # while tmpgrid['time'] is None:
        #     line = file.readline()
        #     _vtk_parse_line(line, tmpgrid)
        # file.close()
        # domain['time'] = tmpgrid['time']

        domain["time"] = grid[0]["time"]

        return domain

    def _set_grid(self):
        grid = []
        # Record filename and data_offset
        for i, fname in enumerate(self.fnames):
            file = open(fname, "rb")
            g = dict()
            g["data"] = dict()
            g["filename"] = fname
            g["read_field"] = None
            g["read_type"] = None

            while g["read_field"] is None:
                g["data_offset"] = file.tell()
                line = file.readline()
                _vtk_parse_line(line, g)

            file.close()
            g["Nx"] -= 1
            g["Nx"][g["Nx"] == 0] = 1
            g["dx"][g["Nx"] == 1] = 1.0

            # Right edge
            g["re"] = g["le"] + g["Nx"] * g["dx"]
            grid.append(g)

        return grid


def _parse_filename(filename):
    """Break up a filename into its component
    to check the extension and extract the output number.

    Parameters
    ----------
    filename : string
        Name of the file, including extension

    Returns
    -------
    tuple containing dirname, problem_id, output number, extension, mpi flag, nonzero_id flag

    Examples
    --------
    >>> _parse_filename('/basedir/id0/problem_id.0000.vtk')
    ('/basedir', 'problem_id', '0000', 'vtk', True, False)

    >>> _parse_filename('/basedir/id10/problem_id-id10.0000.d1.vtk')
    ('/basedir', 'problem_id', '0000', 'vtk', True, True)

    >>> _parse_filename('/basedir/problem_id.0000.vtk')
    ('/basedir', 'problem_id', '0000', 'vtk', False, False)
    """

    sep = osp.sep
    dirname = osp.dirname(filename)
    nonzero_id = False
    # Check if dirname ends with id0
    dirname_last = dirname.split(sep)[-1]
    if dirname_last.startswith("id") and dirname_last[2:].isdigit():
        dirname = sep.join(dirname.split(sep)[:-1])
        mpi_mode = True
    else:
        mpi_mode = False

    base = os.path.basename(filename)
    base_split = base.split(".")
    if len(base_split) == 3:
        problem_id = ".".join(base_split[:-2])
        num = base_split[-2]
        suffix = None
        ext = base_split[-1]
    else:
        try:
            inum = -3
            # test = int(base_split[inum])
            # If dirname is idXX where XX>0, (2d vtk slices)
            # need to remove idXX string from the problem_id
            suffix = base_split[-2]
        except ValueError:
            inum = -2
            suffix = None

        if mpi_mode and int(dirname_last[2:]) != 0:
            problem_id = ".".join(base_split[:inum])
            problem_id = problem_id.replace("-" + dirname_last, "")
            nonzero_id = True
        else:
            problem_id = ".".join(base_split[:inum])

        num = base_split[inum]
        ext = base_split[-1]

    return dirname, problem_id, num, suffix, ext, mpi_mode, nonzero_id


def _set_field_map(grid):
    fp = open(grid["filename"], "rb")
    fp.seek(0, 2)
    eof = fp.tell()
    offset = grid["data_offset"]
    fp.seek(offset)

    field_map = dict()
    if "Nx" in grid:
        Nx = grid["Nx"]

    while offset < eof:
        line = fp.readline()
        sp = line.strip().split()
        field = sp[1].decode("utf-8")
        field_map[field] = dict()
        field_map[field]["read_table"] = False
        field_map[field]["title"] = line
        if b"SCALARS" in line:
            tmp = fp.readline()
            field_map[field]["read_table"] = True
            field_map[field]["nvar"] = 1
        elif b"VECTORS" in line:
            field_map[field]["nvar"] = 3
        else:
            raise TypeError(sp[0] + " is unknown type.")

        field_map[field]["offset"] = offset
        field_map[field]["ndata"] = field_map[field]["nvar"] * grid["ncells"]
        if field == "face_centered_B1":
            field_map[field]["ndata"] = (Nx[0] + 1) * Nx[1] * Nx[2]
        elif field == "face_centered_B2":
            field_map[field]["ndata"] = Nx[0] * (Nx[1] + 1) * Nx[2]
        elif field == "face_centered_B3":
            field_map[field]["ndata"] = Nx[0] * Nx[1] * (Nx[2] + 1)

        if sp[2] == b"int":
            dtype = "i"
        elif sp[2] == b"float":
            dtype = "f"
        elif sp[2] == b"double":
            dtype = "d"

        field_map[field]["dtype"] = dtype
        field_map[field]["dsize"] = field_map[field]["ndata"] * struct.calcsize(dtype)
        fp.seek(field_map[field]["dsize"], 1)
        offset = fp.tell()
        tmp = fp.readline()
        if len(tmp) > 1:
            fp.seek(offset)
        else:
            offset = fp.tell()

    return field_map


def _vtk_parse_line(line, grid):
    sp = line.strip().split()

    if b"vtk" in sp:
        grid["vtk_version"] = sp[-1]
    elif b"time=" in sp:
        time_index = sp.index(b"time=")
        grid["time"] = float(sp[time_index + 1].rstrip(b","))
        if b"level" in sp:
            grid["level"] = int(sp[time_index + 3].rstrip(b","))
        if b"domain" in sp:
            grid["domain"] = int(sp[time_index + 5].rstrip(b","))
        if sp[0] == b"PRIMITIVE":
            grid["prim_var_type"] = True
    elif b"DIMENSIONS" in sp:
        grid["Nx"] = np.array(sp[-3:]).astype("int")
    elif b"ORIGIN" in sp:  # left_edge
        grid["le"] = np.array(sp[-3:]).astype("float64")
    elif b"SPACING" in sp:
        grid["dx"] = np.array(sp[-3:]).astype("float64")
    elif b"CELL_DATA" in sp:
        grid["ncells"] = int(sp[-1])
    elif b"SCALARS" in sp:
        grid["read_field"] = sp[1]
        grid["read_type"] = "scalar"
    elif b"VECTORS" in sp:
        grid["read_field"] = sp[1]
        grid["read_type"] = "vector"
    elif b"POINTS" in sp:
        grid["npoint"] = eval(sp[1])
