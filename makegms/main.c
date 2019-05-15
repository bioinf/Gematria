#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pthread.h"

#include "common/types.c"
#include "common/reader.c"

#include "bloom/filter_helpers.c"
#include "bloom/make_gms_bloom.c"

#include "qsort/decomposition.c"
#include "qsort/sorter.c"
#include "qsort/make_gms_qsort.c"

#include "common/make_gms_wrapper.c"

int main(int argc, char *argv[])
{
    read_length = 100;
    if (argc > 2) read_length = atoi(argv[2]);

    unsigned threads = 1;
    if (argc > 3) threads = atoi(argv[3]);

    unsigned quality = 0;
    if (argc > 4) quality = atoi(argv[4]);

    /* --------------------------------------------------------------------- */

	makeGmsWrapper(argv[1], threads, quality);
	
	big true_x = 0;
	big false_x = 0;
	
    for (big k = 0; k != seq->reads; ++k) {
		if (seq->counts[k] > 1) true_x += 1; else false_x += 1;
    	// printf("%u", seq->counts[k] > 1 ? 0 : 1);
    }

	printf("T: %llu\n", true_x);
	printf("F: %llu\n", false_x);

    /* --------------------------------------------------------------------- */

    if (seq) free(seq->counts);
    if (seq) free(seq);

    return 0;
}

/*
 Helpers:
 gcc src/main.c -std=c99 -m64 -O3 -lpthread -lm -o exec 

 Debug:
 valgrind --leak-check=full --show-leak-kinds=all \
 ./exec ../../test/example.fa 20 1
*/
