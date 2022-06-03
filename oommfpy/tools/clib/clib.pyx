#cython: language_level=3

# -----------------------------------------------------------------------------

cdef extern from "clib.h":

    void WriteVTK_RectilinearGrid(double * gridx, double * gridy, double * gridz,
                                  double * m, double * Ms,
                                  int nx, int ny, int nz,
                                  char * fname
                                  )

    void WriteVTK_ImageData_XML(double r0x, double r0y, double r0z,
                                double dx, double dy, double dz,
                                double * m,
                                double * Ms,
                                int nx, int ny, int nz,
                                char * fname
                                )

# -----------------------------------------------------------------------------

def WriteVTK_RectilinearGrid_C(double [:] gridx,
                               double [:] gridy,
                               double [:] gridz,
                               double [:] m,
                               double [:] Ms,
                               nx, ny, nz,
                               fname
                               ):
    cdef bytes py_bytes = fname.encode()
    cdef char* c_string = py_bytes

    # fname_ = fname.encode('utf-8')
    WriteVTK_RectilinearGrid(&gridx[0], &gridy[0], &gridz[0],
                             &m[0], &Ms[0], nx, ny, nz,
                             &c_string[0])

def WriteVTK_ImageData_XML_C(r0x, r0y, r0z, dx, dy, dz,
                             double [:] m,
                             double [:] Ms,
                             nx, ny, nz,
                             fname
                             ):
    cdef bytes py_bytes = fname.encode()
    cdef char* c_string = py_bytes

    # fname_ = fname.encode('utf-8')
    WriteVTK_ImageData_XML(r0x, r0y, r0z, dx, dy, dz,
                           &m[0], &Ms[0], nx, ny, nz,
                           &c_string[0])
