#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <arpa/inet.h>
#include <stdint.h>

// #define _FILE "../../../../test/test_rect_txt.omf"
#define _FILE "../../../../test/test_rect_b4.omf"


bool StartsWith(const char *a, const char *b)
{
   if(strncmp(a, b, strlen(b)) == 0) return 1;
   return 0;
}

bool isBigEndian(){
   int number = 1;
   return (*(char*)&number != 1);
}

int read_ovf(void) {

    FILE *fp = fopen(_FILE, "rb");
    if(fp == NULL) {
        perror("Unable to open file!");
        return(-1);
    }

    // Chink size to read the first lines
    char chunk[128] = "0";

    while(!StartsWith(chunk, "# Begin: Data")) {
        if(fgets(chunk, 128, fp)!=NULL) {
        /* writing content to stdout */
        puts(chunk);
        }
    }
    // fgets(chunk, 128, fp);
    // puts(chunk);


    float a;
    char buffer[4];
    char buffer_new[4];
    int st;
    st = fread(buffer, 4, 1, fp);

    // // First chunk of data is in Bigendian format so we must reverse the order
    // // of the bytes
    // // My Linux machine is LittleEndian
    // //    echo -n I | od -to2 | head -n1 | cut -f2 -d" " | cut -c6
    // // which results in 1 (big endian -> 0) 
    // for(int i=0; i < sizeof(float); i++) {
    //     buffer_new[sizeof(float) - i - 1] = buffer[i];
    // }
    // memcpy(&a, buffer_new, sizeof(float));

    // Unsigned 32 bit integer (4 bytes)
    uint32_t unint;
    // Convert char to unint
    memcpy(&unint, buffer, sizeof(unint));
    // Network Byte Order to host -> we get an unisgned 4 byte int
    uint32_t ai = ntohl(unint);
    // Convert to float (4 bytes)
    memcpy(&a, &ai, sizeof(float));

    // TODO: if we use 8 bytes, we have to use a double (64 bits) and
    // uint64_t for ntohl
    // TODO: Wrap this reading process in a function

    // sscanf(buffer, "%f", &a);

    printf("Buff %s\n", buffer_new);
    printf("%f\n", a);
    printf("Size t %d\n", st);

    // double vx, vy, vz;
    // int n_data = 9;
    // for(int i=0; i < n_data; i++) {
    //     fscanf(fp, "%lf %lf %lf", &vx, &vy, &vz);
    //     printf("%lf %lf %lf \n", vx, vy, vz);
    // }

    fclose(fp);

    return 0;
}


int main(void) {

    read_ovf();

    return 0;
}
