from .. import MagnetisationData
from . import clib
# import numpy as np
import click
import glob
import os
import timeit


# -----------------------------------------------------------------------------

def omf2vtk(input_omf_file,
            output_vtk_file=None,
            output_format='binary',
            spatial_scale=1e9
            ):
    """
    Convert a given input_omf_file into a VTK file in binary format
    Magnetisation (direction and magnitude) values are stored as cell values
    """

    # Currently the C library accepts only binary format
    if output_format == 'ascii':
        raise Exception('Only binary format is supported in this version')

    data = MagnetisationData(input_omf_file)
    data.generate_field()
    data.generate_coordinates()

    # If no name is supplied use the same name of the input file
    if not output_vtk_file:
        output_vtk_file = input_omf_file[:-3] + '.vtk'
    # Be sure to add the vtk extension to the output file name
    if not output_vtk_file.endswith('.vtk'):
        output_vtk_file += '.vtk'

    # Save VTK in binary format with Ms and m data
    clib.WriteVTK_RectilinearGrid_C(data.grid[0] * spatial_scale,
                                    data.grid[1] * spatial_scale,
                                    data.grid[2] * spatial_scale,
                                    data.field.reshape(-1),
                                    data.field_norm,
                                    data.nx, data.ny, data.nz,
                                    output_vtk_file)


def omf2vtk_batch(input_omfs,
                  output_format='binary',
                  spatial_scale=1e9
                  ):
    """
    Convert all OMF files (or with a different extension) into VTK files
    according to the input path, which can be a directory or a path with a
    wildcard expression. Progress of the conversion is printed on every
    iteration.
    IMPORTANT: this function assumes the OMF files in the folder come from
    the SAME simulation, since the mesh grid is not updated. The loop only
    updates the magnetisation field.

    input_omfs    :: a directory, or path with a wildcard, e.g. 'omfs/m_*.ovf'
    output_format :: only binary
    """

    # Currently the C library accepts only binary format
    if output_format == 'ascii':
        raise Exception('Only binary format is supported in this version')

    if '*' in input_omfs:
        _files = glob.glob(input_omfs)
    else:
        if not input_omfs.endswith('/'):
            input_omfs += '/'
        _files = glob.glob(input_omfs + '*.omf')

    if not _files:
        raise Exception('Could not find any files! Check the input path')

    # Try sorting (can be done properly in the future)
    _files = sorted(_files)

    data = MagnetisationData(_files[0])
    data.generate_field()
    data.generate_coordinates()

    # Save VTK in binary format with Ms and m data
    len_total = len(str(len(_files)))
    print('Converting files in:', os.path.dirname(_files[0]), '\n')
    for i, FILE in enumerate(_files):

        start = timeit.default_timer()
        data.input_file = FILE
        data.generate_field()

        clib.WriteVTK_RectilinearGrid_C(data.grid[0] * spatial_scale,
                                        data.grid[1] * spatial_scale,
                                        data.grid[2] * spatial_scale,
                                        data.field.reshape(-1),
                                        data.field_norm,
                                        data.nx, data.ny, data.nz,
                                        FILE[:-3] + 'vtk')
        end = timeit.default_timer()

        print(f'{str(i + 1).zfill(len_total)}/{len(_files)}:',
              os.path.basename(FILE),
              f'({end - start:.1e} s)')
    print('\nFinished!')


# Command line interface ------------------------------------------------------

@click.command()
@click.option('-i', '--input_path', type=str,
              help='Path to OMF file, directory or wildcard', required=True)
@click.option('-o', '--output_vtk_file', type=click.STRING, default=None,
              help='Output VTK file name', required=False)
@click.option('-of', '--output_format', type=str, default='binary',
              help='Output VTK file name')
@click.option('-ss', '--spatial_scale', type=float, default=1e9,
              help='Output VTK file name')
def omf2vtk_cli(input_path, output_vtk_file, output_format, spatial_scale):

    if os.path.isdir(input_path) or ('*' in input_path):
        omf2vtk_batch(input_path, output_format=output_format,
                      spatial_scale=spatial_scale)

    elif os.path.isfile(input_path):
        omf2vtk(input_path, output_vtk_file, output_format=output_format,
                spatial_scale=spatial_scale)

    else:
        raise Exception('No valid input path')


if __name__ == '__main__':
    omf2vtk_cli()
