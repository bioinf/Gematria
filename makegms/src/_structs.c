#define num unsigned long long int

typedef struct {
    char * sequence;
    num length;
    unsigned short * counts;
    num size;
} Genome;

typedef struct {
    num * positions;
    num pointer;
    num size;
} ThreadData;

unsigned read_length = 100;
Genome * seq = NULL;
// pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
