char *ptr_DNA_seq = NULL;
unsigned int read_lng = 0;

int load_sequence_FASTA( char *fname, char *ptr_DNA_seq );
int count_seq_lng( char *fname, unsigned long long int *ptr_seq_lng );
int read_cmpfunc( const unsigned long long int *a, const unsigned long long int *b );
int fst_read_cmpfunc( const unsigned long long int a, const unsigned long long int b );
int count_reads( unsigned long long int *ptr_positions, unsigned int *ptr_counts, unsigned long long int number_of_reads );

