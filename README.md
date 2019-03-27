[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2611194.svg)](https://doi.org/10.5281/zenodo.2611194)

# OOMMFPy

A very minimal and simple Python library to read and extract data from OOMMF
magnetisation files `omf`. In addition to this library we provide tools to plot
`omf` files and convert them to `vtk` files.

Highlights:

- Read `omf` files in any format
- Painless convertion of the data in an `omf` file into Numpy arrays for data
  analysis
- Fast calculation (using Numpy) of the skyrmion number in a slice of the
  system in any plane orientation (`xy`, `xz`, `yz`)
- Minimal tool to convert `omf` files to VTK
- Plot functions

## Install

A `setup.py` file is provided to install this library using `pip`

    git clone https://github.com/davidcortesortuno/oommfpy
    cd oommfpy
    pip install ./

If successful, the `plot_omf` and `omf2vtk` tools are installed in the
corresponding `bin` directory and can be called from the command line.

## Documentation

For now check the `doc/ipynb` folder which contains a tutorial with basic
functionality

Scripts to convert `omf` to VTK can be called directly as, for example,

```
omf2vtk -i omfs/my_oommf_output.omf -o test.vtk
```

Similarly with the `plot_omf` function.

## TODO

- [ ] More tests
- [ ] More documentation and function docstrings
- [ ] More options to plotting library
- [ ] Point data at cell centres of VTK file
- [ ] Print `z` coordinate when computing sk number
- [ ] Allow Periodic boundaries for the skyrmion number calculation

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
