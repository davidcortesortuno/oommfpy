#ifndef __CLIB__
#define __CLIB__

void WriteVTK_RectilinearGrid(double * gridx, double * gridy, double * gridz,
                              double * m,
                              double * Ms,
                              int nx, int ny, int nz,
                              char * fname
                              );

void WriteVTK_ImageData_XML(double r0x, double r0y, double r0z,
                            double dx, double dy, double dz,
                            double * m,
                            double * Ms,
                            int nx, int ny, int nz,
                            char * fname
                            );

#endif
