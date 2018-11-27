from .. import OOMMFData
import numpy as np
import pyvtk
import click


# -----------------------------------------------------------------------------


@click.command()
@click.option('-i', '--input_omf_file', type=str,
              help='Path to OMF file', required=True)
@click.option('-o', '--output_vtk_file', type=str,
              help='Output VTK file name', required=True)
@click.option('-of', '--output_format', type=str,
              help='Output VTK file name')
def omf2vtk(input_omf_file,
            output_vtk_file,
            output_format='ascii'
            ):
    """
    """

    data = OOMMFData(input_omf_file)
    data.read_m()
    data.set_coordinates()

    grid = (np.linspace(data.xmin * 1e9, data.xmax * 1e9, data.nx + 1),
            np.linspace(data.ymin * 1e9, data.ymax * 1e9, data.ny + 1),
            np.linspace(data.zmin * 1e9, data.zmax * 1e9, data.nz + 1))

    structure = pyvtk.RectilinearGrid(* grid)
    vtk_data = pyvtk.VtkData(structure, "")

    # Save the magnetisation
    vtk_data.cell_data.append(pyvtk.Vectors(np.column_stack((data.mx,
                                                             data.my,
                                                             data.mz)),
                              "m"))

    # Save Ms as scalar field
    vtk_data.cell_data.append(pyvtk.Scalars(data.Ms, "Ms"))

    # Save to VTK file with specified output filename
    vtk_data.tofile(output_vtk_file, output_format)
