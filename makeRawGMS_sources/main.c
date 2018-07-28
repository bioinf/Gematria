#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "tsort.h"
#include "shared.h"

int main(int argc, char *argv[])
{
    unsigned long long int seq_lng = 0;
    unsigned long long int k = 0;
    unsigned long long int number_of_reads = 0;
    unsigned long long int number_of_byte_blocks = 0;
    unsigned long long int bytes_counter = 0;
    unsigned char current_byte = 0;

    unsigned long long int *ptr_positions = NULL;
    unsigned int *ptr_counts = NULL;
    unsigned short int *ptr_GMS_track = NULL;
    unsigned char *ptr_GMS_track_byte = NULL;

    int num_threads = 1;

    int flag = 0;
    time_t rawtime;
    struct tm * timeinfo;
    FILE * pFile;

    if ( argc != 5 ) {
        printf("Wrong number of input arguments: \n");
        printf("  1st argument - path to a dataset; \n");
        printf("  2nd argument - read length. \n");
        printf("  3rd argument - number of threads. \n");
        printf("  4th argument - working directory. \n");
        return -1;
    }
    else {
        read_lng = atoi(argv[2]);
        num_threads = atoi(argv[3]);
    }

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "Raw GMS-track builder was started at: %s", asctime (timeinfo) );

    // Counting nucleotides
    flag = count_seq_lng ( argv[1], &seq_lng );
    printf ("Total number of nucleotides is %llu \n", seq_lng);

    // Getting data from file
    ptr_DNA_seq = (char*)malloc(sizeof(char) * seq_lng);
    flag = load_sequence_FASTA ( argv[1], ptr_DNA_seq );
    printf ("The whole sequence was loaded. \n");

    // Filling auxiliary array
    number_of_reads = seq_lng - read_lng + 1;
    ptr_positions = (unsigned long long int*)malloc(sizeof(unsigned long long int)*number_of_reads);
    for (k = 0; k < number_of_reads; k++)
        *(ptr_positions + k) = k;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "Sorting was started at: %s", asctime (timeinfo) );

    // Sorting of k-mers
    if (num_threads == 1) {
        qsort(ptr_positions, number_of_reads, sizeof(unsigned long long int), (int(*) (const void *, const void *))read_cmpfunc);
    }
    else {
        tsort(ptr_positions, number_of_reads, sizeof(unsigned long long int), (int(*) (const void *, const void *))read_cmpfunc,num_threads);
    }

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "Counting of repeated reads was started at: %s", asctime (timeinfo) );

    // Counting repeated reads
    ptr_counts = (unsigned int*)malloc(sizeof(unsigned int)*number_of_reads);
    flag = count_reads(ptr_positions, ptr_counts, number_of_reads);

    free(ptr_DNA_seq);
    free(ptr_positions);

    // Building mappability track
    ptr_GMS_track = (unsigned short int *)malloc(sizeof(unsigned short int )*number_of_reads);
    for (k = 0; k < number_of_reads; k++)
        if ( *(ptr_counts + k) > 1 )
            *(ptr_GMS_track + k) = 0;
        else
            *(ptr_GMS_track + k) = 1;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    // Writing to binary file
    number_of_byte_blocks = (unsigned long long int)( ceil(number_of_reads/8) );
    ptr_GMS_track_byte = (unsigned char *)malloc(sizeof(unsigned char)*number_of_byte_blocks);
    bytes_counter = 0;
    current_byte = 0;
    for (k = 0; k < number_of_reads; k++){
        current_byte += pow(2,7-k%8)*( *(ptr_GMS_track + k) );
        if ( k%8 == 7 ){
            *(ptr_GMS_track_byte + bytes_counter) = current_byte;
            current_byte = 0;
            bytes_counter++;
        }
    }
    if (number_of_reads%8 != 0)
        *(ptr_GMS_track_byte + bytes_counter) = current_byte;

    pFile = fopen (strcat(argv[4],"/Cache/GMS_track.bin"), "wb");
    fwrite (ptr_GMS_track_byte , sizeof(unsigned char), number_of_byte_blocks, pFile);
    fclose (pFile);

    printf ( "Raw GMS-track builder has finished its task at: %s", asctime (timeinfo) );

    return 0;
}
