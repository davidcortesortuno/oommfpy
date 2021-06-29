[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2611194.svg)](https://doi.org/10.5281/zenodo.2611194)

```
   .+------+-------+-------+-------+-------+-------+-------+
 .'      .'      .'      .'      .'      .'      .'      .'|
+---+--+'------+'------+'------+'------+'------+'------+'  -
|      |       |       |       |       |       |       |   |
|   O  +   O   +   M   +   M   +   F   +   P   +   Y   +   +
|      |       |       |       |       |       |       | .'
+------+'------+'------+'------+-------+-------+-------+'
```

# OOMMFPy

A very minimal and simple Python library to read and extract data from OOMMF
magnetisation files `omf`, which are also used in MuMax3. In addition to this
library we provide tools to plot `omf` files and convert them to `vtk` files.

Highlights:

- Read `omf` files in any format
- Can also read `ovf` files and MuMax3 files
- Painless conversion of the data in an `omf` file into Numpy arrays for data
  analysis
- Fast calculation (using Numpy) of the skyrmion number in a slice of the
  system in any plane orientation (`xy`, `xz`, `yz`)
- Fast reading of `omf` files in binary format (using Numpy's `fromfile`)
- Minimal and super fast tool to convert `omf` files to VTK format
- Plot functions
- Early support for Paraview plugin: read `omf` files directly!

## Install

The easiest is to use `pip` or `poetry` to install the package from
[PyPI](https://pypi.org/project/oommfpy)

    pip install oommfpy

The Github address can also be directly used to install the package via `pip`

    pip install git+https://github.com/davidcortesortuno/oommfpy

Alternatively, a `setup.py` file is provided to install this library

    git clone https://github.com/davidcortesortuno/oommfpy
    cd oommfpy
    pip install ./

If successful, the `plot_omf` and `omf2vtk` tools are installed in the
corresponding `bin` directory and can be called from the command line.

A C library is built with the installation process, thus the setup file tries
to install Cython if is not present in the system.

### Paraview plugin

A first version of a reader for Paraview is added in this last version. For now
the installation is a bit of a hack:

- After installing the `oommfpy` library, locate the `oommfpy` folder from
  the`site-packages` directory

- Download the latest version of Paraview with Python > 3.8 support

- Copy the `oommfpy` directory into the Paraview Python `site-packages` folder.
  For example, for Paraview 5.9.0 installed in the `home` folder:

  ```
  cp -r oommfpy $HOME/ParaView-5.9.0-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages/
  ```

- Open Paraview and go to `Tools -> Manage Plugins -> Load New` and select the
  Python file in the `tools/` folder of `oommfpy` (you can clone the
  repository)

- Now you can open any `omf` file without converting to VTK!

## Documentation

For now check the `doc/ipynb` folder which contains a tutorial with basic
functionality. To load a file with a magnetisation field, which is found more
commonly in simulations, use the `MagnetisationData` class. To load any field,
such as the dipolar field, use the `FieldData` class.

Scripts to convert `omf` to VTK can be called directly as, for example,

```
omf2vtk -i omfs/my_oommf_output.omf -o test.vtk
```

The input path can also be a directory or a path with a wildcard, *e.g.*
`omfs/*.omf`. This method assumes the files in the path come from the same
simulation as the tool loads the mesh from the first file in the path and then
only updates the magnetisation fields.

Similar options are provided for the `plot_omf` function. Use the `--help` for
details.

## TODO

- [ ] More tests
- [ ] Add pyproject.toml file to avoid manual installation of Cython in setup.py
- [ ] More options to plotting library
- [ ] Print `z` coordinate when computing sk number
- [ ] Allow Periodic boundaries for the skyrmion number calculation
- [ ] Add typing check
- [ ] Support for multiple OS

# Citation

If you find this library useful, please cite this repository as:

```
@Misc{Cortes2019,
  author       = {David Cort{\'e}s-Ortu{\~n}o},
  title        = {OOMMFPy},
  howpublished = {Zenodo doi:10.5281/zenodo.2611194. Github: https://github.com/davidcortesortuno/oommfpy},
  year         = {2019},
  doi          = {10.5281/zenodo.2611194},
  url          = {https://doi.org/10.5281/zenodo.2611194},
}
```
