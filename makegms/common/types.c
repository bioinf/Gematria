typedef unsigned int num;
typedef unsigned long long int big;

typedef struct {
    num * positions;
    big pointer;
    big size;
} ThreadData;

typedef struct {
    big * positions;
    big pointer;
    big size;
} ThreadDataBig;

typedef struct {
    unsigned hash_count;
    big length;
    char * second;
    char * first;
    big * hashes;
} Filter;

typedef struct {
    char * sequence;
	unsigned short * counts;
    big length;
    big reads;
} Genome;

unsigned large_genome = 0;
unsigned read_length = 100;
Genome * seq = NULL;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
