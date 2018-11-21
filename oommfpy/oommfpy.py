import numpy as np
import re
import colorsys
import struct


# -----------------------------------------------------------------------------
# Utilities to generate a HSL colourmap from the magnetisation field data

def convert_to_RGB(hls_color):
    return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                        hls_color[1],
                                        hls_color[2]))


def generate_colours(field_data, colour_model='rgb'):
    """
    field_data      ::  (n, 3) array
    """
    hls = np.ones_like(field_data)
    hls[:, 0] = np.arctan2(field_data[:, 1],
                           field_data[:, 0]
                           )
    hls[:, 0][hls[:, 0] < 0] = hls[:, 0][hls[:, 0] < 0] + 2 * np.pi
    hls[:, 1] = 0.5 * (field_data[:, 2] + 1)

    if colour_model == 'rgb':
        rgbs = np.apply_along_axis(convert_to_RGB, 1, hls)
        return rgbs

    elif colour_model == 'hls':
        return hls

    else:
        raise Exception('Specify a valid colour model: rgb or hls')


# -----------------------------------------------------------------------------


class OOMMFData(object):
    """
    Class to extract the magnetisation field data from an OOMMF file with
    a regular mesh grid (coordinates are generated in this class)
    """

    def __init__(self, input_file):

        self.input_file = input_file
        self.read_header()

    def read_header(self):
        """
        Read header from file
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

        attrs = {'xstepsize': 'dx',  'ystepsize': 'dy', 'zstepsize': 'dz',
                 # 'xbase': 'xbase',  'ybase': 'ybase', 'zbase': 'zbase',
                 'xmin': 'xmin', 'ymin': 'ymin', 'zmin': 'zmin',
                 'xmax': 'xmax', 'ymax': 'ymax', 'zmax': 'zmax',
                 }

        # Regex search the attributes. Stepsizes are specified as dx, dy, dz
        for k in attrs.keys():
            num_val = float(re.search('(?<={}: )[0-9\-\.e]+'.format(k),
                            data).group(0))
            setattr(self, attrs[k], num_val)

        # Add data format from last line
        setattr(self, 'data_format', re.search('(?<=Begin: Data ).+',
                                               data).group(0))

        # Compute number of elements in each direction
        self.nx = int((self.xmax - self.xmin) / self.dx)
        self.ny = int((self.ymax - self.ymin) / self.dy)
        self.nz = int((self.zmax - self.zmin) / self.dz)

        # Obtain binary data type from header and check the dtype
        # to use it in Numpy's methods
        # Based on: https://github.com/deparkes/OOMMFTools/blob/master/oommftools/core/oommfdecode.py
        if self.data_format == 'Binary 4':
            flag = _file.read(4)
            if struct.unpack('>f', flag)[0] == 1234567.0:
                # print(struct.unpack('>f', flag)[0])
                self._dtype = '>f4'
            elif struct.unpack('<f', flag)[0] == 1234567.0:
                self._dtype = '<f4'
            else:
                raise Exception('Cannot get binary data dtype')

        elif self.data_format == 'Binary 8':
            flag = _file.read(8)
            if struct.unpack('>d', flag)[0] == 123456789012345.0:
                self._dtype = '>f8'
            elif struct.unpack('<d', flag)[0] == 123456789012345.0:
                self._dtype = '<f8'
            else:
                raise Exception('Cannot get binary data dtype')

        # Check mesh type to calculate the base positions ---------------------

        # Get the first 3 values if using binary data
        if self.data_format.startswith('Binary'):
            first_num_data = []
            for i in range(3):
                # the data using the number of binary bits
                num_data = _file.read(int(self.data_format[-1]))
                # unpack using dtype without the num of bits
                first_num_data.append(struct.unpack(self._dtype[:-1],
                                                    num_data)[0])
        # else read the next line after: Begin: Data with the coords and spins
        else:
            first_num_data = _file.readline()[1:].split(' ')

        self.meshtype = re.search('(?<=meshtype: )[a-z]+', data).group(0)

        for i, k in enumerate(['xbase', 'ybase', 'zbase']):
            if self.meshtype == 'irregular':
                setattr(self, k, float(first_num_data[i]))
            elif self.meshtype == 'rectangular':
                num_val = float(re.search('(?<={}: )[0-9\-\.e]+'.format(k),
                                data).group(0))
                setattr(self, k, num_val)
        # ---------------------------------------------------------------------

        _file.close()

    def read_m(self):
        """
        Read the magnetisation data

        If the data is in binary format, we decode the information using
        the Numpy's fromfile function. Otherwise just load the text file
        with Numpy's loadtxt
        """

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
                # Count the total number of spins, this is the total number
                # of magnetisation elements read in the binary file
                # We are assuming we only have spin data: mx, my, mz per line
                n_spins = self.nx * self.ny * self.nz
                # Discard the first element (flag) and discard the final data,
                # which is the final comment ending the data file (the binary
                # decoding of numpy transforms this comment in non-sense numbers)
                # Finally reshape to have columns: mx my mz
                if self.meshtype == 'irregular':
                    data = data[1:6 * n_spins + 1].reshape(-1, 6)
                    # only get the last 3 cols with spin data
                    data = data[:, 3:]
                elif self.meshtype == 'rectangular':
                    data = data[1:3 * n_spins + 1].reshape(-1, 3)

        else:
            if self.meshtype == 'irregular':
                data = np.loadtxt(self.input_file, usecols=[3, 4, 5])
            elif self.meshtype == 'rectangular':
                data = np.loadtxt(self.input_file)

        self.Ms = np.sqrt(np.sum(data ** 2, axis=1))
        self.Ms[self.Ms == 0.0] = 0.0
        self.mx, self.my, self.mz = (data[:, 0], data[:, 1], data[:, 2])
        self.mx[self.Ms != 0.0] /= self.Ms[self.Ms != 0.0]
        self.my[self.Ms != 0.0] /= self.Ms[self.Ms != 0.0]
        self.mz[self.Ms != 0.0] /= self.Ms[self.Ms != 0.0]

    def set_coordinates(self):
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


# -----------------------------------------------------------------------------


class OOMMFODTRead(object):

    """
    Class to read an ODT file from an OOMMF simulation output
    Columns from the ODT table can be obtained calling an element of this class
    from a  valid column name, e.g.
        data = OOMMFODTRead(my_odt_file)
        total_energy = data['Oxs_CGEvolve::Total energy']
    """

    def __init__(self, input_file):

        self.input_file = input_file
        self.read_header()
        # Load the ODT file numerical data
        self.data = np.loadtxt(self.input_file)

    def read_header(self):
        f = open(self.input_file)

        # Read the third line which has the header names
        i = 0
        for i in range(4):
            l = f.readline()
        f.close()
        # Remove the starting "# Columns:" string
        l = l[11:]

        # Separate the header strings:
        # re.split('}\s{|\s\s\s\s|\s{|\sOxs', l)
        header = re.findall(r'Oxs_[A-Za-z\s\:}]+(?=Oxs|{Oxs|\n)', l)
        # Assign the name and column number
        self.columns = {}
        for i, h in enumerate(header):
            h = h.strip()
            h = h.strip('}')
            self.columns[h] = i

    def __getitem__(self, column_name):
        """
        Returns the correspondign column from the name when calling
        an element of this Class through []
        """
        if column_name not in self.columns.keys():
            raise Exception('Invalid column name: {}. \n'.format(column_name) +
                            'Options:\n' + '\n'.join(self.columns.keys())
                            )

        return self.data[:, self.columns[column_name]]
