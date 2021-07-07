import numpy as np
import re
import struct
from pathlib import Path
from collections import namedtuple
# Many of these are simplified in Python 3.9:
from typing import Union
from typing import Optional
from typing import Any
from typing import Sequence
from typing import Literal
from typing import Dict
from typing import List
# import numpy.typing as npt


# -----------------------------------------------------------------------------

def loadtxt_iter(txtfile: Union[str, Path],
                 ncols: int,
                 fromiter_count: int = -1,
                 delimiter: Optional[str] = None,
                 skiprows: int = 0,
                 dtype: Any = np.float64,  # npt.DTypeLike = np.float64,
                 usecols: Optional[Sequence[int]] = None,
                 comment: str = '#'
                 ) -> Any:  # npt.ArrayLike:
    """Reads a simply formatted text file using Numpy's `fromiter` function.
    This function should perform faster than the `loadtxt` function.

    Parameters
    ----------
    txtfile
        Path to text file
    ncols
        Number of columns (not detected automatically)
    fromiter_count
        If passed, it might make numpy's fromiter more efficient (see doc)
    delimiter
        Passed to `split(delimiter=)` in every line of the text file.
        `None` means any number of white spaces
    skiprows
        Skip `skiprows` initial lines
    dtype

    Returns
    -------
    Numpy array with the data

    Notes
    -----
    Based on N. Schlomer function at:

    https://stackoverflow.com/questions/18259393/numpy-loading-csv-too-slow-compared-to-matlab

    and J. Kington at

    https://stackoverflow.com/questions/8956832/python-out-of-memory-on-large-csv-file-numpy
    """
    def iter_func():
        with open(txtfile, 'r') as infile:
            for _ in range(skiprows):
                next(infile)

            for line in infile:
                if line.startswith(comment):
                    continue
                # Might not be necessary to strip characters at start and end:
                line = line.strip().split(delimiter)
                # line = line.split(delimiter)
                # As a general solution we can also use regex:
                # re.split(" +", line)
                for item in line:
                    yield dtype(item)

    data = np.fromiter(iter_func(), dtype=dtype,
                       count=fromiter_count).flatten()

    if usecols is None:
        COLS = np.arange(ncols)
    else:
        COLS = usecols
    data = data.reshape((-1, ncols))[:, COLS]

    return data


# -----------------------------------------------------------------------------


class FieldData(object):
    """
    Class to extract the field data from an OOMMF output file (omf, ohf, ehf)
    which is ordered in a mesh grid. If the grid is regular, coordinates are
    generated in this class.

    Methods
    -------
    read_header
    generate_field
    generate_coordinates
    """

    def __init__(self, input_file: Union[str, Path]):

        self.input_file = input_file

        # Values from header of file:
        self.dx = self.dy = self.dz = 0.0
        self.xmin = self.ymin = self.zmin = 0.0
        self.xmax = self.ymax = self.zmax = 0.0
        self.valuedim = 0
        self.xbase = self.ybase = self.zbase = 0.0
        self.nx = self.ny = self.nz = 0
        self.data_format = ''

        self.read_header()

    def read_header(self) -> None:
        """
        Read header from file and store values as class variables. Some of the
        original names are simplified:

            [c]stepsize -> [c]x
            [c]nodes    -> n[c]

        with [c] = x, y, z
        """

        # Generate a single string with the whole header up to the line where
        # numerical Data starts
        # If not binary file, open as binary anyways:
        _file = open(self.input_file, 'rb')

        line = _file.readline().decode()
        data = ''
        while not line.startswith('# Begin: Data'):
            data += line
            line = _file.readline().decode()
        data += line

        HeaderDict = {'xstepsize': 'dx',  'ystepsize': 'dy', 'zstepsize': 'dz',
                      # 'xbase': 'xbase',  'ybase': 'ybase', 'zbase': 'zbase',
                      'xmin': 'xmin', 'ymin': 'ymin', 'zmin': 'zmin',
                      'xmax': 'xmax', 'ymax': 'ymax', 'zmax': 'zmax',
                      'valuedim': 'valuedim',
                      'xnodes': 'nx', 'ynodes': 'ny', 'znodes': 'nz'}

        # data_dict: Dict[str, Any] = {}

        # Regex search the attributes. Stepsizes are specified as dx, dy, dz
        for k in HeaderDict.keys():

            num_val = re.search(r'(?<={}: )[0-9\-\.e]+'.format(k), data)

            # .................................................................
            # Only OVF 2.0 file format allows arbitrary number of dimensions
            # of the field data. Otherwise just set it to default 3
            if k == 'valuedim':
                if num_val is None:
                    setattr(self, HeaderDict[k], 3)
                else:
                    setattr(self, HeaderDict[k], int(num_val.group(0)))

            # .................................................................
            # Rectangular grids specify the number of nodes/mesh-points in
            # every direction
            elif k == 'xnodes' or k == 'ynodes' or k == 'znodes':
                if num_val is None:
                    setattr(self, HeaderDict[k], None)
                else:
                    setattr(self, HeaderDict[k], int(num_val.group(0)))
            # .................................................................

            else:
                if num_val is None:
                    raise Exception(f'Unable to find value for {k}')
                else:
                    setattr(self, HeaderDict[k], float(num_val.group(0)))

        # Add data format from last line
        datafmt = re.search(r'(?<=Begin: Data ).+', data)
        if datafmt is None:
            raise Exception('Unable to find data format')
        else:
            setattr(self, 'data_format', datafmt.group(0))

        # Compute number of elements in each direction
        # Assuming if nx is None, then ny and nz also are None
        if self.nx is None:
            for c in ['x', 'y', 'z']:
                # Compute nx = (xmax - xmin) / dx :
                diff = (getattr(self, f'{c}max') - getattr(self, f'{c}min'))
                setattr(self, f'n{c}', round(diff / getattr(self, f'd{c}')))

        # Obtain binary data type from header and check the dtype
        # to use it in Numpy's methods
        # Based on: https://github.com/deparkes/OOMMFTools/blob/master/oommftools/core/oommfdecode.py
        if self.data_format == 'Binary 4':
            flag = _file.read(4)
            if struct.unpack('>f', flag)[0] == 1234567.0:
                # print(struct.unpack('>f', flag)[0])
                self._dtype = '>f4'
                self._dtype_st = '>f'
            elif struct.unpack('<f', flag)[0] == 1234567.0:
                self._dtype = '<f4'
                self._dtype_st = '<f'
            else:
                raise Exception('Cannot get binary 4 data dtype')

        # In this case we could use _dtype = '>f8' or '<f8' but this does not
        # work with struct.unpack, which requires a double(?) >d or <d
        elif self.data_format == 'Binary 8':
            flag = _file.read(8)
            if struct.unpack('>d', flag)[0] == 123456789012345.0:
                self._dtype = '>f8'
                self._dtype_st = '>d'
            elif struct.unpack('<d', flag)[0] == 123456789012345.0:
                self._dtype = '<f8'
                self._dtype_st = '<d'
            else:
                raise Exception('Cannot get binary 8 data dtype')

        # Check mesh type to calculate the base positions ---------------------

        # Get the first 3 values if using binary data
        first_num_data: List[float]
        if self.data_format.startswith('Binary'):
            first_num_data = []
            for i in range(3):
                # the data using the number of binary bits
                num_data = _file.read(int(self.data_format[-1]))
                # Decode the data using the dtype without binary bits number(?)
                first_num_data.append(struct.unpack(self._dtype_st,
                                                    num_data)[0])

        # else read the next line after: Begin: Data with the coords and spins
        else:
            # this should guess data is separated by any number of white spaces
            first_num_data = [float(v) for v in
                              _file.readline().decode()[1:].split()]

        meshtype = re.search('(?<=meshtype: )[a-z]+', data)
        if meshtype is None:
            raise Exception('Cannot find meshtype')
        else:
            self.meshtype = meshtype.group(0)

        for i, k in enumerate(['xbase', 'ybase', 'zbase']):
            if self.meshtype == 'irregular':
                setattr(self, k, first_num_data[i])
            elif self.meshtype == 'rectangular':
                num_val = re.search(r'(?<={}: )[0-9\-\.e]+'.format(k), data)
                if num_val is None:
                    raise Exception('Cannot set base mesh values')
                else:
                    setattr(self, k, float(num_val.group(0)))
        # ---------------------------------------------------------------------

        _file.close()

    def _generate_data(self) -> Any:  # npt.ArrayLike:
        """
        If the data is in binary format, we decode the information using
        Numpy's `fromfile` function. Otherwise load the text file
        using `loadtxt_iter`
        """

        if self.data_format is None:
            raise Exception('File data format not identified')

        if self.data_format == 'Binary 4' or self.data_format == 'Binary 8':

            with open(self.input_file, 'rb') as _file:
                # First read the initial comments until the data begins
                line = _file.readline().decode()
                while not line.startswith('# Begin: Data'):
                    line = _file.readline().decode()

                # Load the the binary data using the dtype decoded when
                # reading the header
                # TODO: check if self._dtype exists
                # We discarded the comments at the beginning of the file,
                # however we still have the comment finishing the data:
                # "# End: Data" which we need to discard
                data = np.fromfile(_file, dtype=self._dtype)

                # According to OOMMF the data starts with an endian flag thus
                # it is the first element read in the data variable, which
                # which depends on the binary type
                endianflag = data[0]
                # Count the total number of mesh sites, this is the total
                # number of field elements read in the binary file
                # We are assuming we only have field data: fx, fy, fz per line
                n_sites = self.nx * self.ny * self.nz
                # Discard the first element (flag) and discard the final data,
                # which is the final comment ending the data file (the binary
                # decoding of numpy transforms this comment in non-sense nums)
                # Finally reshape to have columns: fx fy fz
                if self.meshtype == 'irregular':
                    # Number of dimensions of data: coordinates + field
                    field_dim = 3 + self.valuedim
                    data = data[1:field_dim * n_sites + 1].reshape(-1, field_dim)
                    # only get the last cols with field data (non coordinates)
                    data = data[:, 3:]
                elif self.meshtype == 'rectangular':
                    field_dim = self.valuedim
                    data = data[1:field_dim * n_sites + 1].reshape(-1, field_dim)

        else:
            if self.meshtype == 'irregular':
                count = 6 * self.nx * self.ny * self.nz
                data = loadtxt_iter(self.input_file, ncols=6, comment='#',
                                    usecols=[3, 4, 5], fromiter_count=count)
            elif self.meshtype == 'rectangular':
                count = 3 * self.nx * self.ny * self.nz
                data = loadtxt_iter(self.input_file, ncols=3, comment='#',
                                    fromiter_count=count)

        return data

    def generate_field(self, normalise_field: bool = False) -> None:
        """
        Read the field data: any scalar or vector field assuming the data from
                             the field is always Nx3 in size

        normalise_field     :: Creates the self.nfield variable with the field
                               normalised per mesh site
        """
        self.field = self._generate_data()
        self.field_norm = np.sqrt(np.sum(self.field ** 2, axis=1))
        self.field_norm[self.field_norm == 0.0] = 0.0

        if normalise_field:
            self.nfield = np.copy(self.field)
            self.nfield_x = self.nfield[:, 0]
            self.nfield_y = self.nfield[:, 1]
            self.nfield_z = self.nfield[:, 2]

            _filter = self.field_norm != 0.0
            self.nfield_x[_filter] /= self.field_norm[_filter]
            self.nfield_y[_filter] /= self.field_norm[_filter]
            self.nfield_z[_filter] /= self.field_norm[_filter]

    def generate_coordinates(self) -> None:
        """
        Create the self.x, self.y, self.z arrays with the coordinates of the
        mesh sites. Unique values of the coordinates are stored in the
        correpsonding xs, yz, zs arrays
        """
        xs, ys, zs = (np.arange(float(self.nx)),
                      np.arange(float(self.ny)),
                      np.arange(float(self.nz))
                      )
        xs *= self.dx
        xs += self.xbase
        ys *= self.dy
        ys += self.ybase
        zs *= self.dz
        zs += self.zbase

        xs = np.tile(np.tile(xs, self.ny), self.nz)
        ys = np.tile(np.repeat(ys, self.nx), self.nz)
        zs = np.repeat(zs, self.nx * self.ny)

        # int(len(xs), len(ys), len(zs))

        # self.coordinates = np.ravel(np.column_stack((xs, ys, zs))) * 1e9
        self.coordinates = np.column_stack((xs, ys, zs)) * 1e9
        self.x, self.y, self.z = (self.coordinates[:, 0],
                                  self.coordinates[:, 1],
                                  self.coordinates[:, 2])

        self.xs = np.unique(self.x)
        self.ys = np.unique(self.y)
        self.zs = np.unique(self.z)

        # ---------------------------------------------------------------------
        # Generate the grid of vertices of the mesh of cuboids, useful for
        # saving VTKs
        self.grid = (np.linspace(self.xmin, self.xmax, self.nx + 1),
                     np.linspace(self.ymin, self.ymax, self.ny + 1),
                     np.linspace(self.zmin, self.zmax, self.nz + 1)
                     )


# -----------------------------------------------------------------------------


class MagnetisationData(FieldData):
    """
    Class to extract the magnetisation data from an OOMMF omf file. The
    magnetisation field is normalised when the generate_field() method is
    called.
    This class includes a method to compute the skyrmion number
    """

    def __init__(self, input_file: Union[str, Path]):

        super(MagnetisationData, self).__init__(input_file)

    def generate_field(self) -> None:
        """
        Compute the magnetisation field data from the given input file

        The data is stored in the self.field variable, and individual
        components in the self.field_i variables, with i in {x,y,z}

        The field norm self.field_norm is also stored in the alias self.Ms
        """
        self.field = self._generate_data()
        self.field_norm = np.sqrt(np.sum(self.field ** 2, axis=1))
        self.field_norm[self.field_norm == 0.0] = 0.0
        # Alternative name
        self.Ms = self.field_norm

        self.field_x, self.field_y, self.field_z = (self.field[:, 0],
                                                    self.field[:, 1],
                                                    self.field[:, 2])

        _filter = self.field_norm != 0.0
        self.field_x[_filter] /= self.field_norm[_filter]
        self.field_y[_filter] /= self.field_norm[_filter]
        self.field_z[_filter] /= self.field_norm[_filter]

        # Generate alternative names for the field data:
        self.mx = self.field_x
        self.my = self.field_y
        self.mz = self.field_z

    __PlaneOptions = Literal['xy', 'xz', 'yz']
    __MethodOptions = Literal['finite_differences', 'spin_lattice']

    def compute_sk_number(self,
                          index: int = 0,
                          plane: __PlaneOptions = 'xy',
                          method: __MethodOptions = 'finite_differences'
                          ) -> float:
        r"""

        Compute the skyrmion number S for a two dimensional layer. To do this,
        we convert the self.m array with the magnetization, into a (nx, ny, 3)
        matrix (a grid), using a 2D slice of the mesh. The slice lies in the
        specified plane and an index integer. For example, plane='xy' and
        index=0 means a slice at the z=self.zs[0] coordinate

        This function creates a self.sk_number array with the sk number density
        per mesh site and returns the total sum

        Calculation: there are two methods:

        * FINITE DIFFERENCES: defined as
                                _
                     1         /       dm     dm
             S   =  ---  *    /   m .  --  X  --   dx dy
                    4 PI   _ /         dx     dy

        A finite difference discretisation of the continuum magnetisation
        field, using central differences, and a simple midpoint rule
        for the 2D integration leads to:

            S = (1 / 4 * PI) *  (  M_i \cdot ( M_{i+1} \times M_{j+1} )
                                 + M_i \cdot ( M_{i-1} \times M_{j-1} )
                                 - M_i \cdot ( M_{i-1} \times M_{j+1} )
                                 - M_i \cdot ( M_{i+1} \times M_{j-1} )
                                 ) / 4

        * SPIN LATTICE: Kim and Mulkers [IOP SciNotes1(2020) 025211] algorithm
        to compute the topological charge of a finite differences grid using
        Berg and Luscher spin lattice approach, and averaged using 4 spin
        triangles. The charge for a triangle made of spins q_i, q_j, q_k is
        defined as q_ijk, thus the total charge is
                        __
            S =  1     \   q_ijk
                ---    /__
                4 PI  triangles

            tan( 1 q_ijk )                M_i * (M_j X M_k)
               ( -       )  =  --------------------------------------
               ( 2       )      1 + M_i * M_j + M_i * M_k + M_j * M_k

        The triangles need to be weighted by 1.0 when the spins are next to an
        empty site or at the boundary [IOP SciNotes1(2020) 025211], and by 0.5
        when they are within the sample.

        Parameters
        ----------
        plane
            'xy', 'xz' or 'yz'
        index
            any integer from 0 up to len(xs) or len(ys) or len(zs) depending on
            the slice plane
        method
            'finite_differences' or 'spin_lattice'
        """

        # Spin data in a grid, dimensions are: z, y, x, 3
        spin_grid = self.field.reshape(-1, self.ny, self.nx, 3)
        # Get the specified slice along the specified dimension
        if plane == 'xy':
            # y-direction in axis=1; x-direction in axis=2  transforms to
            # -----> y-direction in axis=0; x-direction in axis=1
            spin_grid = spin_grid[index, :, :, :]
        elif plane == 'xz':
            # z-direction in axis=0; x-direction in axis=2
            # -----> z-direction in axis=0; x-direction in axis=1
            spin_grid = spin_grid[:, index, :, :]
        elif plane == 'yz':
            # z-direction in axis=0 y-direction in axis=1
            # -----> z-direction in axis=0; y-direction in axis=1
            spin_grid = spin_grid[:, :, index, :]
        else:
            raise Exception('Specify a valid plane')

        # 2nd argument are how many zeroes (before,after) we pad in each axis
        # (we keep 3-spin-components, so we don't pad anything at axis=2)
        spin_pad = np.pad(spin_grid, ((1, 1), (1, 1), (0, 0)),
                          mode='constant', constant_values=0.0)
        # Same as doing:
        # spin_pad = np.zeros((nx +2, ny + 2, 3))
        # spin_pad[1:-1, 1:-1, :] = spin_grid

        # TODO: adding support for PBCs is straightforward if we copy the 1st
        # column to the N+1 column, and the Nth column to the 0th column
        # in the corresponding direction of PBC

        # Skyrmion number density
        self.sk_number = np.zeros((spin_pad.shape[0] - 2,
                                   spin_pad.shape[1] - 2))
        # A view of the self.sk_number array
        sk_num = self.sk_number[:, :]

        # Here we vectorise the cross products using the padded matrix to
        # obtain neighbours (which are zero spin components) at the boundary of
        # the sample
        # We obtain a grid with the sk number density
        if method == 'finite_differences':
            # Here we still use central difference at the boundary spins
            # Notice: if we use a 'xy' plane the
            #         x direction is in axis=1 and y direction in axis=0
            sk_num_cross = (np.cross(spin_pad[1:-1, 2:, :],   # s(i+1,j)
                                     spin_pad[2:, 1:-1, :],   # s(i,j+1)
                                     axis=2) +
                            np.cross(spin_pad[1:-1, :-2, :],  # s(i-1,j)
                                     spin_pad[:-2, 1:-1, :],  # s(i,j-1)
                                     axis=2) -
                            np.cross(spin_pad[1:-1, :-2, :],  # s(i-1,j)
                                     spin_pad[2:, 1:-1, :],   # s(i,j+1)
                                     axis=2) -
                            np.cross(spin_pad[1:-1, 2:, :],   # s(i+1,j)
                                     spin_pad[:-2, 1:-1, :],  # s(i,j-1)
                                     axis=2)
                            )

            # The dot product of every site with the cross product between
            # their neighbours that was already computed above
            # We save this quantity to the self.sk_number method

            # self.sk_number = np.sum(spin_grid * sk_num,
            #                         axis=2) / (16 * np.pi)
            np.einsum('ijk,ijk->ij', spin_grid, sk_num_cross, out=sk_num)
            np.multiply(0.25, sk_num, out=sk_num)

            # -----------------------------------------------------------------
            # The following approach uses backward or forward difference
            # for spins at the boundary

            # # s(i+1,j) - s(i-1,j)
            # fdiff_x = spin_pad[2:, 1:-1, :] - spin_pad[:-2, 1:-1, :]
            # # Forward difference for spins at the boundary +x
            # ngbs_x = np.linalg.norm(spin_pad[2:, 1:-1, :], axis=2)
            # fdiff_x[ngbs_x < 1e-6] += spin_grid[ngbs_x < 1e-6]
            # # Backward difference for spins at the boundary -x (re-use array)
            # ngbs_x[:] = np.linalg.norm(spin_pad[:-2, 1:-1, :], axis=2)
            # fdiff_x[ngbs_x < 1e-6] -= spin_grid[ngbs_x < 1e-6]
            #
            # # s(i+1,j) - s(i-1,j)
            # fdiff_y = spin_pad[1:-1, 2:, :] - spin_pad[1:-1, :-2, :]
            # # Forward difference for spins at the boundary +y
            # ngbs_y = np.linalg.norm(spin_pad[1:-1, 2:, :], axis=2)
            # fdiff_y[ngbs_y < 1e-6] += spin_grid[ngbs_y < 1e-6]
            # # Backward difference for spins at the boundary -y (re-use array)
            # ngbs_y[:] = np.linalg.norm(spin_pad[1:-1, :-2, :], axis=2)
            # fdiff_y[ngbs_y < 1e-6] -= spin_grid[ngbs_y < 1e-6]

            # np.einsum('ijk,ijk->ij', spin_grid,
            #           np.cross(fdiff_x, fdiff_y, axis=2), out=sk_num),
            # np.multiply(0.25, sk_num, out=sk_num)

        elif method == 'spin_lattice':
            # This list contains 2-tuples where each element of every tuple
            # is a slice (or Numpy slice np.s_) to obtain the components of
            # the spin in one of the triangle vertices, e.g.
            #   triangles[0][0] are the spin components of s(i+1, j) and
            #   triangles[0][1] are the spin components of s(i, j+1)
            #
            #              s(i+1,j)            s(i,j+1)        # right triangle
            triangles = [(np.s_[1:-1, 2:, :], np.s_[2:, 1:-1, :]),
                         # s(i-1,j)            s(i,j-1)        # bottom left
                         (np.s_[1:-1, :-2, :], np.s_[:-2, 1:-1, :]),
                         # s(i,j+1)           s(i-1,j)         # left triangle
                         (np.s_[2:, 1:-1, :], np.s_[1:-1, :-2, :]),
                         # s(i,j-1)            s(i+1,j)        # bottom right
                         (np.s_[:-2, 1:-1, :], np.s_[1:-1, 2:, :])]
            # These are weights associated to spins in the opposite corner of
            # the triangles to avoid double counting
            weight_ngbs = [np.s_[2:, 2:, :],    # s(i+1,j+1)
                           np.s_[:-2, :-2, :],  # s(i-1,j-1)
                           np.s_[2:, :-2, :],   # s(i-1,j+1)
                           np.s_[:-2, 2:, :]]   # s(i+1,j-1)

            for n, ngbs in enumerate(triangles):   # bottom right

                # Compute weights of digonally opposite neighbour: 1 or 1/2
                weights = np.linalg.norm(spin_pad[weight_ngbs[n]], axis=2)
                weights = np.where(weights < 1e-6, 1, 0.5)

                denom = 1 + (np.einsum('ijk,ijk->ij',               # s_i * s_j
                                       spin_pad[1:-1, 1:-1, :],
                                       spin_pad[ngbs[0]]) +
                             np.einsum('ijk,ijk->ij',               # s_i * s_k
                                       spin_pad[1:-1, 1:-1, :],
                                       spin_pad[ngbs[1]]) +
                             np.einsum('ijk,ijk->ij',               # s_j * s_k
                                       spin_pad[ngbs[0]],
                                       spin_pad[ngbs[1]])
                             )

                # s_j X s_k
                triangle_charge = np.cross(spin_pad[ngbs[0]],  # s_j
                                           spin_pad[ngbs[1]],  # s_k
                                           axis=2)
                # s_i * (s_j X s_k)
                triangle_charge = np.einsum('ijk,ijk->ij',
                                            spin_pad[1:-1, 1:-1, :],
                                            triangle_charge)
                np.arctan2(triangle_charge, denom, out=triangle_charge)
                np.multiply(2.0, triangle_charge, out=triangle_charge)
                # Multiply the weights
                np.multiply(weights, triangle_charge, out=triangle_charge)

                sk_num += triangle_charge
        else:
            raise Exception('Specify valid method for the calculation of Q_sk')

        # Total sk number (integral)
        return np.sum(self.sk_number.flatten()) / (4. * np.pi)

# -----------------------------------------------------------------------------


class OOMMFODTReader(object):

    """
    Class to read an ODT file from an OOMMF simulation output
    Columns from the ODT table can be obtained calling an element of this class
    from a  valid column name, e.g.
        data = OOMMFODTRead(my_odt_file)
        total_energy = data['Oxs_CGEvolve::Total energy']
    """

    def __init__(self, input_file: Union[str, Path]):

        self.input_file = input_file
        self.read_header()
        # Load the ODT file numerical data
        self.data = np.loadtxt(self.input_file)

    def read_header(self):
        f = open(self.input_file)

        # Read the third line which has the header names
        i = 0
        for i in range(4):
            line = f.readline()
        f.close()
        # Remove the starting "# Columns:" string
        line = line[11:]

        # Separate the header strings:
        # re.split('}\s{|\s\s\s\s|\s{|\sOxs', l)
        header = re.findall(r'Oxs_[A-Za-z\s\:}]+(?=Oxs|{Oxs|\n)', line)
        # Assign the name and column number
        self.columns = {}
        for i, h in enumerate(header):
            h = h.strip()
            h = h.strip('}')
            self.columns[h] = i

    def __getitem__(self, column_name: str) -> Any:  # npt.ArrayLike:
        """
        Returns the corresponding column from the name when calling
        an element of this Class through []
        """
        if column_name not in self.columns.keys():
            raise Exception('Invalid column name: {}. \n'.format(column_name) +
                            'Options:\n' + '\n'.join(self.columns.keys())
                            )

        return self.data[:, self.columns[column_name]]
