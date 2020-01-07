import oommfpy.tools as ot
import numpy as np

nx, ny, nz = 3, 3, 1
m = np.array([-1., 0., 1.,
              0, 0, 1,
              1, 0, 0,
              -1, 1, 1,
              0, 1, 1,
              0, 1, 1,
              -1, 0, 0,
              1, 0, 0,
              1, 0, 0
              ])

Ms = np.ones(9) * 5.

xmin, xmax, ymin, ymax, zmin, zmax = 0., 3., 0., 3, 0., 1.
grid = (np.linspace(xmin, xmax, nx + 1),
        np.linspace(ymin, ymax, ny + 1),
        np.linspace(zmin, zmax, nz + 1)
        )

# Write VTK file
vtk_fname = 'test_c_vtk_writer.vtk'
ot.clib.WriteVTK_RectilinearGrid_C(grid[0], grid[1], grid[2],
                                   m, Ms,
                                   nx, ny, nz,
                                   vtk_fname)
