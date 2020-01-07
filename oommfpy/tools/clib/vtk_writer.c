#include "clib.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef unsigned char uchar;
// Reverse byte to BigEndian representation (reverse byte order?) from the
// LittleEndian repr given in Linux (in other architectures this might not be
// necessary)
double DoubleSwap( double f )
{
   union
   {
      double f;
      // byte b[4];
      uchar b[sizeof(double)];
   } dat1, dat2;

   dat1.f = f;
   dat2.b[0] = dat1.b[7];
   dat2.b[1] = dat1.b[6];
   dat2.b[2] = dat1.b[5];
   dat2.b[3] = dat1.b[4];
   dat2.b[4] = dat1.b[3];
   dat2.b[5] = dat1.b[2];
   dat2.b[6] = dat1.b[1];
   dat2.b[7] = dat1.b[0];
   return dat2.f;
}

/* Write array of doubles with size n into the file given by the pointer fptr
 * Values are written with the BigEndian order
 */
void WriteDouble(double * arr, int n, FILE * fptr) {
    for(int i = 0; i < n; i++) {
        double f = DoubleSwap(arr[i]);
        fwrite(&f, sizeof(f), 1, fptr);
    }
}

void WriteMagData(double * m,   // 3 * nx*ny*nz array
                  double * Ms,  // nx*ny*nz array
                  int n,
                  FILE * fptr,
                  char * header
                  ) {

    // Magnetisation **********************************************************

    sprintf(header,"\nCELL_DATA %d\n", n);
    fprintf(fptr, "%s", header);

    sprintf(header,"SCALARS Ms double 1\n");
    fprintf(fptr, "%s", header);
    sprintf(header,"LOOKUP_TABLE default\n");
    fprintf(fptr, "%s", header);
    for(int i = 0; i < n; i++) {
        // printf("%f\n", Ms[i]);
        double d = DoubleSwap(Ms[i]);
        fwrite(&d, sizeof(d), 1, fptr);
    }

    // Spin directions ********************************************************

    sprintf(header,"\nVECTORS m double\n");
    fprintf(fptr, "%s", header);
    for(int i = 0; i < n; i++) {
        double dx = DoubleSwap(m[3 * i    ]);
        double dy = DoubleSwap(m[3 * i + 1]);
        double dz = DoubleSwap(m[3 * i + 2]);
        fwrite(&dx, sizeof(dx), 1, fptr);
        fwrite(&dy, sizeof(dy), 1, fptr);
        fwrite(&dz, sizeof(dz), 1, fptr);
    }
}

// If the mesh is made up of nx * ny * nz cells, it has
// (nx + 1) * (ny + 1) * (nz + 1) vertices.
void WriteVTK_RectilinearGrid(double * gridx, double * gridy, double * gridz,
                              double * m,   // 3 * nx*ny*nz array
                              double * Ms,  // nx*ny*nz array
                              int nx, int ny, int nz,
                              char * fname
                              ) {

    char header[1024];
    FILE * fptr;
    // Create binary file
    fptr = fopen(fname, "wb");

    if(fptr == NULL) {
        printf("Error!");
        exit(1);
    }

    sprintf(header,"# vtk DataFile Version 2.0\n");
    sprintf(header + strlen(header),"OOMMFPY VTK Data\n");
    sprintf(header + strlen(header),"BINARY\n");
    sprintf(header + strlen(header),"DATASET %s\n","RECTILINEAR_GRID");
    fprintf(fptr, "%s", header);

    // COORDINATES ------------------------------------------------------------

    sprintf(header,"DIMENSIONS %d %d %d\n", nx + 1, ny + 1, nz + 1);
    fprintf(fptr, "%s", header);

    sprintf(header,"X_COORDINATES %d double\n", nx + 1);
    fprintf(fptr, "%s", header);
    WriteDouble(gridx, nx + 1, fptr);

    sprintf(header,"\nY_COORDINATES %d double\n", ny + 1);
    fprintf(fptr, "%s", header);
    WriteDouble(gridy, ny + 1, fptr);

    sprintf(header,"\nZ_COORDINATES %d double\n", nz + 1);
    fprintf(fptr, "%s", header);
    WriteDouble(gridz, nz + 1, fptr);

    // DATA -------------------------------------------------------------------

    int n_cell_data = nx * ny * nz;
    WriteMagData(m, Ms, n_cell_data, fptr, header);

    // ------------------------------------------------------------------------

    fclose(fptr);
}
