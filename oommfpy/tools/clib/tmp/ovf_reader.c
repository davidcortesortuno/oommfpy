#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <arpa/inet.h>  // Windows might need winsock2.h
#include <stdint.h>

// https://stackoverflow.com/questions/3022552/is-there-any-standard-htonl-like-function-for-64-bits-integers-in-c
#if __BIG_ENDIAN__
// In this case it is not necessary to swap byte order
# define ntohll(x) (x)
#else
// If the system is LittleEndian (e.g. Linux) then define ntohll for 64 bits
# define ntohll(x) (((uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))
#endif

// #define _FILE "../../../../test/test_rect_txt.omf"
#define _FILE "../../../../test/test_rect_b8.omf"

typedef enum {
    TXT,
    BINARY4,
    BINARY8,
} data_type;


bool StartsWith(const char *a, const char *b)
{
   if(strncmp(a, b, strlen(b)) == 0) return 1;
   return 0;
}

// Check if machine is Bigendian
// bool isBigEndian(){
//    int number = 1;
//    return (*(char*)&number != 1);
// }

void read_binary8_N_times (FILE * fp, double * decoded_arr, int N) {

    char buffer[8];
    int st;
    // Unsigned 32 bit integer (4 bytes)
    uint64_t unint;
    
    for (int i = 0; i < N; ++i) {
        st = fread(buffer, 8, 1, fp);
        // Convert char to unint
        memcpy(&unint, buffer, sizeof(unint));
        // Network Byte Order to host -> we get an unsigned 8 byte int
        unint = ntohll(unint);
        // Convert to doouble (8 bytes)
        memcpy(&decoded_arr[i], &unint, sizeof(double));
    }
}

void read_binary4_N_times (FILE * fp, float * decoded_arr, int N) {

    char buffer[4];
    int st;
    // Unsigned 32 bit integer (4 bytes)
    uint32_t unint;
    
    for (int i = 0; i < N; ++i) {
        st = fread(buffer, 4, 1, fp);
        // Convert char to unint
        memcpy(&unint, buffer, sizeof(unint));
        // Network Byte Order to host -> we get an unsigned 4 byte int
        unint = ntohl(unint);
        // Convert to float (4 bytes)
        memcpy(&decoded_arr[i], &unint, sizeof(float));
    }
}


float read_binary4_chunk (FILE * fp) {

    float f;
    char buffer[4];
    int st;
    st = fread(buffer, 4, 1, fp);

    // Unsigned 32 bit integer (4 bytes)
    uint32_t unint;
    // Convert char to unint
    memcpy(&unint, buffer, sizeof(unint));
    // Network Byte Order to host -> we get an unsigned 4 byte int
    unint = ntohl(unint);
    // Convert to float (4 bytes)
    memcpy(&f, &unint, sizeof(float));

    return f;
}

double read_binary8_chunk (FILE * fp) {

    double d;
    char buffer[8];
    int st;
    uint64_t unint;
    st = fread(buffer, 1, 8, fp);
    memcpy(&unint, buffer, sizeof(unint));
    unint = ntohll(unint);
    memcpy(&d, &unint, sizeof(double));
    return d;
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

    double f;
    // The first chunk of data should print the signature 1234567 for binary4
    f = read_binary8_chunk(fp);
    printf("%lf\n", f);

    // char buffer[4];
    // char buffer_new[4];
    // int st;
    // st = fread(buffer, 4, 1, fp);
    // // First chunk of data is in Bigendian format so we must reverse the order
    // // of the bytes
    // // My Linux machine is LittleEndian
    // //    echo -n I | od -to2 | head -n1 | cut -f2 -d" " | cut -c6
    // // which results in 1 (big endian -> 0) 
    // for(int i=0; i < sizeof(float); i++) {
    //     buffer_new[sizeof(float) - i - 1] = buffer[i];
    // }
    // memcpy(&f, buffer_new, sizeof(float));

    // TODO: if we use 8 bytes, we have to use a double (64 bits) and
    // uint64_t for ntohl
    // TODO: Wrap this reading process in a function


    // Now we get the rest of the data
    // for(int i = 0; i < 6; i++) {
    //     f = read_binary4_chunk(fp);
    //     printf("%f\n", f);
    // }
    
     double data[27];

     read_binary8_N_times(fp, data, 27);
     for (int i = 0; i < 27; ++i) {
         printf("%lf\n", data[i]);
     }

    fclose(fp);

    return 0;
}


int main(void) {

    read_ovf();

    return 0;
}
