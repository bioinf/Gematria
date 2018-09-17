#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pthread.h"

#include "_structs.c"
#include "_reader.c"
#include "_decomposition.c"
#include "_sorter.c"
#include "_makegms.c"

int main(int argc, char *argv[])
{
    read_length = atoi(argv[2]);
    unsigned int threads = atoi(argv[3]);
    /* --------------------------------------------------------------------- */

    makegms(argv[1], threads);

    num blocks = (num)( ceil(seq->size/8) );
    unsigned char * GMS = (unsigned char *)malloc(sizeof(unsigned char)*blocks);
    num bytes_counter = 0;
    unsigned char current_byte = 0;
    
    for (num k = 0; k != seq->size; ++k) {
        current_byte += pow(2,7-k%8)*(seq->counts[k] > 1 ? 0 : 1);
        if ( k%8 == 7 ){
            *(GMS + bytes_counter) = current_byte;
            current_byte = 0;
            bytes_counter++;
        }
    }
    if (seq->size%8 != 0)
        *(GMS + bytes_counter) = current_byte;
    
    FILE * pFile = fopen ("./track.bin", "wb");
    fwrite (GMS , sizeof(unsigned char), blocks, pFile);
    fclose (pFile);
    
    /* --------------------------------------------------------------------- */
    free(GMS);
    free(seq->counts);
    free(seq);
    return 0;
}

/*
 Helpers:
 gcc -pthread main.c -lm -o exec
 valgrind --leak-check=full --show-leak-kinds=all \
 ./exec ../../test/example.fa 20 1
 */
