# OOMMFPy

A very minimal and simple Python library to read and extract data from OOMMF
magnetisation files `omf`. In addition to this library we provide tools to plot
`omf` files and convert them to `vtk` files.

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
- [ ] More documentation
- [ ] More options to plotting library
- [ ] Point data at cell centres of VTK file
- [ ] Print `z` coordinate when computing sk number
- [ ] Allow Periodic boundaries

# Citation

Cite this repository as:

```
```
