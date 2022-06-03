from .. import MagnetisationData
from . import clib
# import numpy as np
import click
import timeit
from pathlib import Path
import numpy as np


# -----------------------------------------------------------------------------

def omf2vtk(input_omf_file,
            output_vtk_file=None,
            output_format='binary',
            spatial_scale=1e9
            ):
    """
    Convert a given input_omf_file into an Image Data VTK file with numerical
    data in raw binary format (not base64 encoded and not compressed).
    Magnetisation (direction and magnitude) values are stored as cell values
    """
    # TODO: spatial_scale should be replaced by valuemultiplier from the header

    # Currently the C library accepts only binary format
    if output_format == 'ascii':
        raise Exception('Only raw binary format is supported in this version')

    data = MagnetisationData(input_omf_file)
    data.generate_field()
    data.generate_coordinates()

    # If no name is supplied use the same name of the input file
    # Add the vti extension to the output file name
    if not output_vtk_file:
        output_vtk_file = Path(input_omf_file).with_suffix('.vti')
    else:
        output_vtk_file = Path(output_vtk_file).with_suffix('.vti')

    # Save VTK in binary format with Ms and m data
    clib.WriteVTK_ImageData_XML_C(data.xmin * spatial_scale,
                                  data.ymin * spatial_scale,
                                  data.zmin * spatial_scale,
                                  data.dx * spatial_scale,
                                  data.dy * spatial_scale,
                                  data.dz * spatial_scale,
                                  data.field.reshape(-1),
                                  data.field_norm,
                                  data.nx, data.ny, data.nz,
                                  str(output_vtk_file))


def omf2vtk_batch(input_omfs,
                  output_format='binary',
                  spatial_scale=1e9,
                  vtk_format='xml'
                  ):
    """
    Convert all `omf` files (or with a different extension) into VTK files
    according to the input path, which can be a directory or a path with a
    wildcard expression. Progress of the conversion is printed on every
    iteration.
    IMPORTANT: this function assumes the OMF files in the folder come from
    the SAME simulation, since the mesh grid is not updated. The loop only
    updates the magnetisation field.

    input_omfs    :: a directory, or path with a wildcard, e.g. 'omfs/m_*.ovf'
    output_format :: only binary
    vtk_format    :: xml or legacy
    """

    # Currently the C library accepts only binary format
    if output_format == 'ascii':
        raise Exception('Only binary format is supported in this version')

    input_path = Path(input_omfs)
    if '*' in input_path.name:
        # Use the wildcard pattern from the path name
        _files = input_path.parent.glob(input_path.name)
    else:
        _files = input_path.parent.glob('*.omf')

    # Convert to a list (we might change this if there are too many files)
    # Try sorting (can be done properly in the future)
    _files = sorted(list(_files))

    if not _files:
        raise Exception('Could not find any files! Check the input path')

    data = MagnetisationData(_files[0])
    data.generate_field()

    if vtk_format == 'xml':
        data.generate_coordinates()
    elif vtk_format == 'legacy':
        data.generate_coordinates(compute_vertex_grid=True)
    else:
        raise Exception('Specify a correct file version')

    r_base = np.array([data.xmin, data.ymin, data.zmin]) * spatial_scale
    dr = np.array([data.dx, data.dy, data.dz]) * spatial_scale

    # Save VTK in binary format with Ms and m data
    len_total = len(str(len(_files)))
    print('Converting files in:', input_path.parent, '\n')
    for i, FILE in enumerate(_files):

        start = timeit.default_timer()
        data.input_file = FILE
        data.generate_field()

        if vtk_format == 'xml':
            clib.WriteVTK_ImageData_XML_C(*r_base, *dr,
                                          data.field.reshape(-1),
                                          data.field_norm,
                                          data.nx, data.ny, data.nz,
                                          str(FILE.with_suffix('.vti')))
        elif vtk_format == 'legacy':
            clib.WriteVTK_RectilinearGrid_C(
                data.vertex_grid[0] * spatial_scale,
                data.vertex_grid[1] * spatial_scale,
                data.vertex_grid[2] * spatial_scale,
                data.field.reshape(-1),
                data.field_norm,
                data.nx, data.ny, data.nz,
                str(FILE.with_suffix('.vtk')))

        end = timeit.default_timer()

        print(f'{str(i + 1).zfill(len_total)}/{len(_files)}:',
              FILE.name,
              f'({end - start:.1e} s)')

    print('\nFinished!')


# Command line interface ------------------------------------------------------

@click.command()
@click.argument('input_path', type=click.Path(), nargs=-1)
@click.option('-o', '--output_vtk_file', type=click.STRING, default=None,
              help='Output VTK file name (try using .vti suffix)',
              required=False)
@click.option('-of', '--output_format', type=str, default='binary',
              help='Format: ascii or binary. binary only in this version')
@click.option('-ss', '--spatial_scale', type=float, default=1e9,
              help='Scale for the spatial coordinates, e.g.: 1e9')
@click.option('-vf', '--vtk_format', type=str, default='xml',
              help='Format for the VTK file: xml or legacy. xml is modern.' +
              ' For a single file only xml is used, this option is skipped')
def omf2vtk_cli(input_path, output_vtk_file, output_format,
                spatial_scale, vtk_format):

    for IP in input_path:

        IP = Path(IP)

        if IP.is_dir() or ('*' in IP.name):
            omf2vtk_batch(IP, output_format=output_format,
                          spatial_scale=spatial_scale, vtk_format=vtk_format)

        elif IP.is_file():
            omf2vtk(IP, output_vtk_file, output_format=output_format,
                    spatial_scale=spatial_scale)

        else:
            raise Exception(f'No valid input path: {IP}')


if __name__ == '__main__':
    omf2vtk_cli()
