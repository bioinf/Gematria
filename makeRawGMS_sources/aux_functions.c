#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int load_sequence_FASTA( char *fname , char *ptr_DNA_seq )
{
    FILE * fp;
	char * line = NULL;
	size_t len = 0;
	size_t new_substring;
	int k = 0;
	unsigned int substring_lng = 0;
	unsigned long long int nuc_counter = 0;

	fp = fopen(fname, "r");
	if (fp == NULL)
		exit(EXIT_FAILURE);

	while ((new_substring = getline(&line, &len, fp)) != -1) {
		if (line[new_substring - 1] == '\n') {
            line[new_substring - 1] = '\0';
			substring_lng = new_substring-1; // Delete \n symbol
        }
        else
			substring_lng = new_substring;

		if (line[0] != '>') // Not name of read
			for (k = 0; k < substring_lng; k++) {
                *(ptr_DNA_seq + nuc_counter) = line[k];
                nuc_counter++;
            }
        else {
			printf("Reading nucleotides from chr. %s \n", line);
        }
	}

	fclose(fp);

	return 0;
}

int count_seq_lng ( char *fname, unsigned long long int *ptr_seq_lng )
{
    FILE * fp;
	char * line = NULL;
	size_t len = 0;
	size_t new_substring;

    unsigned int substring_lng = 0;
	unsigned long long int seq_lng = 0;

	fp = fopen(fname, "r");
	if (fp == NULL)
		exit(EXIT_FAILURE);

	while ((new_substring = getline(&line, &len, fp)) != -1) {
		/* Delete trailing newline */
		if (line[new_substring - 1] == '\n'){
            line[new_substring - 1] = '\0';
			substring_lng = new_substring-1; // Delete \n symbol
		}
        else
			substring_lng = new_substring;

		if (line[0] != '>') // Not name of read
			seq_lng += (unsigned long long int)substring_lng;
        else
			printf("Counting number of nucleotides in chr. %s \n", line);
	}

	*ptr_seq_lng = seq_lng;

	fclose(fp);

	return 0;
}

int read_cmpfunc(const unsigned long long int *a, const unsigned long long int *b)
{
    extern char *ptr_DNA_seq;
    extern unsigned int read_lng;

    int k = 0;

    char *ptr_read_1 = ptr_DNA_seq + *a;
    char *ptr_read_2 = ptr_DNA_seq + *b;

    for (k = 0; k < read_lng; k++){
        if (*ptr_read_1 > *ptr_read_2)
            return +1;
        if (*ptr_read_1 < *ptr_read_2)
            return -1;
        ptr_read_1++;
        ptr_read_2++;
    }
    return 0;
}

int fst_read_cmpfunc(const unsigned long long int a, const unsigned long long int b)
{
    extern char *ptr_DNA_seq;
    extern unsigned int read_lng;

    int k = 0;

    char *ptr_read_1 = ptr_DNA_seq + a;
    char *ptr_read_2 = ptr_DNA_seq + b;

    for (k = 0; k < read_lng; k++){
        if (*ptr_read_1 > *ptr_read_2)
            return +1;
        if (*ptr_read_1 < *ptr_read_2)
            return -1;
        ptr_read_1++;
        ptr_read_2++;
    }
    return 0;
}

int count_reads(unsigned long long int *ptr_positions, unsigned int *ptr_counts, unsigned long long int number_of_reads)
{
//    extern char *ptr_DNA_seq;
    unsigned long long int k = 0;
    unsigned long long int start_pos = 0;
    unsigned long long int current_pos = 0;
    unsigned int counts = 0;

//    printf("nor = %llu \n", number_of_reads);
    for (k = 0; k < number_of_reads-1; k++){
//        printf ("step = %llu \n", k);
        counts++;
        if ( fst_read_cmpfunc( *(ptr_positions+k), *(ptr_positions+k+1) ) != 0 ){
//            printf("%llu - ", k);
            for (current_pos = k + 1; current_pos-- > start_pos; )
//                printf("%llu \n", current_pos);
                *(ptr_counts + *(ptr_positions + current_pos )) = counts;
//                printf("!!! \n");
            start_pos = k + 1;
            counts = 0;
//            printf(">>> \n");
        }
//        printf("aaaa \n");
    }
//    printf("xxx \n");
//    k++;
    counts++;
    for (current_pos = k + 1; current_pos-- > start_pos; )
        *(ptr_counts + *(ptr_positions + current_pos)) = counts;

    return 0;
}
