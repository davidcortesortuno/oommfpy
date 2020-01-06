#ifndef __CLIB__
#define __CLIB__

void WriteVTK_RectilinearGrid(double * gridx, double * gridy, double * gridz,
                              double * m,
                              double * Ms,
                              int nx, int ny, int nz,
                              char * fname
                              );

#endif
