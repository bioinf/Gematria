int compare(const num * a, const num * b)
{
    extern Genome * seq;
    extern unsigned read_length;
    
    unsigned int k = 0;
    char * read_a = seq->sequence + *a;
    char * read_b = seq->sequence + *b;
    
    for (k = 0; k < read_length; k++) {
        if (*read_a > *read_b) return +1;
        if (*read_a < *read_b) return -1;
        read_a++;
        read_b++;
    }
    
    return 0;
}

int compare_big(const big * a, const big * b)
{
    extern Genome * seq;
    extern unsigned read_length;
    
    unsigned int k = 0;
    char * read_a = seq->sequence + *a;
    char * read_b = seq->sequence + *b;
    
    for (k = 0; k < read_length; k++) {
        if (*read_a > *read_b) return +1;
        if (*read_a < *read_b) return -1;
        read_a++;
        read_b++;
    }
    
    return 0;
}


void set_counts(num last, num k, num * positions, unsigned cnt)
{
    extern Genome * seq;

    for (num i = k + 1; i-- > last; ) {
        num pt = *(positions + i);
        if (pt < seq->reads) seq->counts[pt] = cnt;
    }
}

void set_counts_big(big last, big k, big * positions, unsigned cnt)
{
    extern Genome * seq;
	
    for (big i = k + 1; i-- > last; ) {
        big pt = *(positions + i);
        if (pt < seq->reads) seq->counts[pt] = cnt;
    }
}


void * sorter(void * data)
{
    extern Genome * seq;
    extern unsigned read_length;

    ThreadData * block = (ThreadData*) data;
    qsort(block->positions, block->pointer, sizeof(num),
          (int(*) (const void *, const void *)) compare);

    num k = 0;
    num last = 0;
    unsigned cnt = 0;

    if (block->pointer > 0) {
        for (k = 0; k < block->pointer - 1; k++) {
            cnt++;
            if (compare( &block->positions[k], &block->positions[k+1] ) != 0) {
                set_counts(last, k, block->positions, cnt);
                cnt  = 0;
                last = k + 1;
            }
        }
        cnt++;
        set_counts(last, k, block->positions, cnt);
    }

    pthread_exit(0);
}

void * sorter_big(void * data)
{
    extern Genome * seq;
    extern unsigned read_length;

    ThreadDataBig * block = (ThreadDataBig*) data;
    qsort(block->positions, block->pointer, sizeof(big),
          (int(*) (const void *, const void *)) compare_big);
    
    big k = 0;
    big last = 0;
    unsigned cnt = 0;
    
    if (block->pointer > 0) {
        for (k = 0; k < block->pointer - 1; k++) {
            cnt++;
            if (compare_big( &block->positions[k], &block->positions[k+1] ) != 0) {
                set_counts_big(last, k, block->positions, cnt);
                cnt  = 0;
                last = k + 1;
            }
        }
        cnt++;
        set_counts_big(last, k, block->positions, cnt);
    }

    pthread_exit(0);
}