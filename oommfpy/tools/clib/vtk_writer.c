#include "clib.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef unsigned char uchar;
// Reverse byte order to BigEndian representation from the LittleEndian repr
// (eg in Linux) (in other architectures this might not be necessary)
double DoubleSwap_L2B(double f, unsigned short endn)
{
   // If big endian, do not do anything
   if (!endn) {return f;}

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

// For big endian double, swap bytes to big endian repr
// This function also converts from double to float
float DoubleToFloat_B2L(double db, unsigned short endn)
{

   union
   {
      float f;
      // byte b[4];
      uchar b[sizeof(float)];
   } dat1, dat2;
   dat1.f = db;  // implicit cast double to float

   // If little endian, do not swap
   if (endn) {return dat1.f;}

   dat2.b[0] = dat1.b[3];
   dat2.b[1] = dat1.b[2];
   dat2.b[2] = dat1.b[1];
   dat2.b[3] = dat1.b[0];
   return dat2.f;
}

/* Write array of doubles with size n into the file given by the pointer fptr
 * Values are written with the BigEndian order
 */
void WriteDouble(double * arr, int n, FILE * fptr, unsigned short endn) {
    for(int i = 0; i < n; i++) {
        double f = DoubleSwap_L2B(arr[i], endn);
        fwrite(&f, sizeof(f), 1, fptr);
    }
}

void WriteMagData(double * m,   // 3 * nx*ny*nz array
                  double * Ms,  // nx*ny*nz array
                  int n,
                  FILE * fptr,
                  char * header,
                  unsigned short endn
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
        double d = DoubleSwap_L2B(Ms[i], endn);
        fwrite(&d, sizeof(d), 1, fptr);
    }

    // Spin directions ********************************************************

    sprintf(header,"\nVECTORS m double\n");
    fprintf(fptr, "%s", header);
    for(int i = 0; i < n; i++) {
        double dx = DoubleSwap_L2B(m[3 * i    ], endn);
        double dy = DoubleSwap_L2B(m[3 * i + 1], endn);
        double dz = DoubleSwap_L2B(m[3 * i + 2], endn);
        fwrite(&dx, sizeof(dx), 1, fptr);
        fwrite(&dy, sizeof(dy), 1, fptr);
        fwrite(&dz, sizeof(dz), 1, fptr);
    }
}

// If the mesh is made up of nx * ny * nz cells, it has
// (nx + 1) * (ny + 1) * (nz + 1) vertices.
// This writer uses the VTK legacy format: .vtk (old)
void WriteVTK_RectilinearGrid(double * gridx, double * gridy, double * gridz,
                              double * m,   // 3 * nx*ny*nz array
                              double * Ms,  // nx*ny*nz array
                              int nx, int ny, int nz,
                              char * fname
                              ) {
    int x = 1;
    // Little endian: true
    unsigned short ENDIANNESS = *((char*)&x) == 1;

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
    WriteDouble(gridx, nx + 1, fptr, ENDIANNESS);

    sprintf(header,"\nY_COORDINATES %d double\n", ny + 1);
    fprintf(fptr, "%s", header);
    WriteDouble(gridy, ny + 1, fptr, ENDIANNESS);

    sprintf(header,"\nZ_COORDINATES %d double\n", nz + 1);
    fprintf(fptr, "%s", header);
    WriteDouble(gridz, nz + 1, fptr, ENDIANNESS);

    // DATA -------------------------------------------------------------------

    int n_cell_data = nx * ny * nz;
    WriteMagData(m, Ms, n_cell_data, fptr, header, ENDIANNESS);

    // ------------------------------------------------------------------------

    fclose(fptr);
}

// ----------------------------------------------------------------------------
// A more modern writer using XML format and using ImageData (smaller size of
// file taking advantage of the mesh regularity)
// https://vtk.org/Wiki/VTK_XML_Formats#The_VTKFile_Element

// ImageData is the most simple grid: regular cells in x, y and z positions
// Save in .vti format. This writer uses the more modern XML format.
// Data is saved as floats (more efficient although approximated)
void WriteVTK_ImageData_XML(double r0x, double r0y, double r0z,  // Origin
                            double dx, double dy, double dz,     // Spacings: dx, dy, dz
                            double * m,   // 3 * nx*ny*nz array
                            double * Ms,  // nx*ny*nz array
                            int nx, int ny, int nz,
                            char * fname
                            ) {

    int x = 1;
    unsigned short ENDIANNESS = *((char*)&x) == 1;

    char header[1024];
    FILE * fptr;
    // Create binary file
    fptr = fopen(fname, "wb");

    if(fptr == NULL) {
        printf("Error!");
        exit(1);
    }

    // HEADER .................................................................

    sprintf(header                 , "<?xml version=\"1.0\"?>\n");
    sprintf(header + strlen(header), "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    sprintf(header + strlen(header), "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%f %f %f\" Spacing=\"%f %f %f\">\n",
                                     nx,  ny, nz, r0x, r0y, r0z, dx, dy, dz);
    sprintf(header + strlen(header), "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx,  ny, nz);
    fprintf(fptr, "%s", header);

    // DATA ...................................................................
    // Will be appended, see: https://vtk.org/Wiki/VTK_XML_Formats
    //                        https://kitware.github.io/vtk-examples/site/VTKFileFormats/
    // If the data is NOT appended, it has to be base-64 encoded (we would need
    // a library or a header file for this)

    // Here we compute the lengths for the spin and the Ms data as int and byte
    union
    {
        unsigned int integer;
        unsigned char byte[4];
    } m_data_len, Ms_data_len;
    m_data_len.integer = 3 * nx * ny * nz * sizeof(float);  // float -> 4 bytes
    Ms_data_len.integer = nx * ny * nz * sizeof(float);

    // Start the data .........................................................
    sprintf(header, "      <CellData Vectors=\"m\" Scalars=\"Ms\">\n");

    sprintf(header + strlen(header),
            "        <DataArray type=\"Float32\" Name=\"m\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>\n");

    // offset bytes include the 4 bytes from the size of m after "_"
    sprintf(header + strlen(header),
            "        <DataArray type=\"Float32\" Name=\"Ms\" format=\"appended\" offset=\"%d\"/>\n", 4 + m_data_len.integer);

    fprintf(fptr, "%s", header);

    // ........................................................................

    sprintf(header, "      </CellData>\n");
    sprintf(header + strlen(header), "    </Piece>\n");
    sprintf(header + strlen(header), "  </ImageData>\n");
    fprintf(fptr, "%s", header);

    // ........................................................................
    // Appended data (not base64 encoded -> raw)
    // The format is
    //          NNNN<data>NNNN<data>NNNN<data>
    //          ^         ^         ^
    //          1         2         3
    //
    // where each "NNNN" is an unsigned 32-bit integer (4 bytes), and <data>
    // consists of a number of bytes equal to the preceding NNNN value.  The
    // corresponding DataArray elements must have format="appended" and offset
    // attributes equal to the following:
    //
    // 1.) offset="0"
    // 2.) offset="(4+NNNN1)"
    // 3.) offset="(4+NNNN1+4+NNNN2)"
    // ...

    // Start with "_"
    sprintf(header, "  <AppendedData encoding=\"raw\">\n    _");
    fprintf(fptr, "%s", header);

    // m data: (write the length as 4 bytes int first, then data)
    fwrite(&m_data_len.byte[0], sizeof(m_data_len.byte), 1, fptr);
    // printf("m Bytes: %u Size: %d\n", m_data_len.byte[0], m_data_len.integer);
    // Can't do a single line write since we have to cast double into float :
    // fwrite(&m[0], sizeof(float), 3 * nx * ny * nz, fptr);
    // So make a loop: (take into account endianness)
    for(int i = 0; i < (nx * ny * nz); i++) {
        float mx = DoubleToFloat_B2L(m[3 * i    ], ENDIANNESS);  // implicit cast double to float
        float my = DoubleToFloat_B2L(m[3 * i + 1], ENDIANNESS);
        float mz = DoubleToFloat_B2L(m[3 * i + 2], ENDIANNESS);
        fwrite(&mx, sizeof(mx), 1, fptr);
        fwrite(&my, sizeof(my), 1, fptr);
        fwrite(&mz, sizeof(mz), 1, fptr);
        // printf("mx %f my %f mz %f\n", mx, my, mz);
    }

    // Ms data:
    fwrite(&Ms_data_len.byte[0], sizeof(Ms_data_len.byte), 1, fptr);
    // printf("Ms Bytes: %u Size: %d\n", Ms_data_len.byte[0], Ms_data_len.integer);
    // fwrite(&Ms[0], sizeof(float), nx * ny * nz, fptr);
    for(int i = 0; i < (nx * ny * nz); i++) {
        float M = DoubleToFloat_B2L(Ms[i], ENDIANNESS);
        fwrite(&M, sizeof(M), 1, fptr);
    }

    sprintf(header, "\n  </AppendedData>\n");
    fprintf(fptr, "%s", header);

    // ........................................................................
    sprintf(header, "</VTKFile>");
    fprintf(fptr, "%s", header);

    // ------------------------------------------------------------------------

    fclose(fptr);
}
