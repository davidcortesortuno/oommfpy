import oommfpy.tools as ot
import oommfpy as op
import numpy as np
import os
import subprocess
import shutil
import timeit
import glob
import pyvtk
import pytest
this_dir = os.path.dirname(os.path.abspath(__file__))


def generate_omfs():

    """
    Simulate a skyrmion reducing the mesh discretisation
    using the OOMMF script: isolated_sk_DMI_Cnv.mif
    It is necessary to have OOMMF installed with the latest DMI
    modules.

    For the speed test we generate a mesh of 200 x 200 x 1 cells
    """

    OMF_DIR = os.path.join(this_dir, 'vtk_writer_omfs/')
    if os.path.exists(OMF_DIR):
        shutil.rmtree(OMF_DIR)
    os.makedirs(OMF_DIR)

    for n in [200]:
        SIM_NAME = 'vtk_writer_omfs/isolated_sk_Cnv_n_{:03d}'.format(n)
        SCRIPT = os.path.join(this_dir, 'isolated_sk_DMI_Cnv.mif')

        # if glob.glob(SIM_NAME + '*.omf'):
        #     continue

        job = ('oommf boxsi -threads 2 '
               '-parameters "NX {0} '
               'BASENAME {1}" '
               '"{2}"'.format(n, SIM_NAME, SCRIPT)
               )
        print(job)

        subprocess.call(job, shell=True)


def test_vtk_writer():
    nx, ny, nz = 3, 3, 1
    m = np.array([-1., 0., 1.,
                  0, 0.5, 1,
                  1, 0, 0,
                  -1, 1, 1,
                  0, 1, 1,
                  0, 1, 1,
                  -1, 0.75, 0,
                  1, 0, 0,
                  1, 0, 0
                  ])

    Ms = np.ones(9) * 5.
    Ms[3:5] = 8.
    Ms[6:9] = 15.

    xmin, xmax, ymin, ymax, zmin, zmax = 0., 3., 0., 3, 0., 1.
    grid = (np.linspace(xmin, xmax, nx + 1),
            np.linspace(ymin, ymax, ny + 1),
            np.linspace(zmin, zmax, nz + 1)
            )

    # Write VTK file in legacy format
    vtk_fname = 'test_c_vtk_writer_rect.vtk'
    ot.clib.WriteVTK_RectilinearGrid_C(grid[0], grid[1], grid[2],
                                       m, Ms,
                                       nx, ny, nz,
                                       vtk_fname)

    # Use different base coordinates:
    vtk_fname = 'test_c_vtk_writer_imdata.vti'
    ot.clib.WriteVTK_ImageData_XML_C(-1.5, -1.0, -2.0,
                                     1.0, 1.0, 1.0,
                                     m, Ms,
                                     nx, ny, nz,
                                     vtk_fname)

    # TODO:
    # Reading the VTK file and transorming the binary data to analyse it


@pytest.mark.skip(reason="Needs update: create a large VTK file")
def test_vtk_writer_speed():
    """
    Comparison of C backend with the old pyvtk interface

    For this we use a 200 x 200 x 1 mesh and convert the OMF file into a VTK
    file using both the C library and the PyVTK library. We run the same
    conversion for each method during 20 times and print a mean.

    """

    # C backend ---------------------------------------------------------------
    _file = glob.glob('./vtk_writer_omfs/isolated_sk_Cnv_n_200-Oxs*.omf')[0]
    data = op.MagnetisationData(_file)
    data.generate_field()
    data.generate_coordinates()
    output_vtk_file = 'vtk_writer_omfs/isolated_sk_Cnv_n_200.vtk'

    C_timings = []
    for i in range(20):
        start = timeit.default_timer()

        ot.clib.WriteVTK_RectilinearGrid_C(data.grid[0],
                                           data.grid[1],
                                           data.grid[2],
                                           data.field.reshape(-1),
                                           data.field_norm,
                                           data.nx, data.ny, data.nz,
                                           output_vtk_file)

        end = timeit.default_timer()
        C_timings.append(end - start)
    print('C writer (best of 20): ', np.mean(np.array(C_timings)))

    # PyVTK -------------------------------------------------------------------
    output_vtk_file = 'vtk_writer_omfs/isolated_sk_Cnv_n_200_pyvtk.vtk'

    pyvtk_timings = []
    for i in range(20):
        start = timeit.default_timer()

        structure = pyvtk.RectilinearGrid(* data.grid)
        vtk_data = pyvtk.VtkData(structure, "")

        # Save the magnetisation
        vtk_data.cell_data.append(pyvtk.Vectors(data.field, "m"))

        # Save Ms as scalar field
        vtk_data.cell_data.append(pyvtk.Scalars(data.field_norm, "Ms"))

        # Save to VTK file with specified output filename
        vtk_data.tofile(output_vtk_file, 'binary')

        end = timeit.default_timer()
        pyvtk_timings.append(end - start)
    print('PyVTK writer (best of 20): ', np.mean(np.array(pyvtk_timings)))

    # -------------------------------------------------------------------------


if __name__ == "__main__":

    # Generate the OMF file with the 200 x 200 x 1 mesh
    generate_omfs()
    # Compare the conversion of the OMF file into VTK using the C library
    # with a conversion using PyVTK
    test_vtk_writer_speed()

    # Test rect grid and image data VTK writers
    test_vtk_writer()
