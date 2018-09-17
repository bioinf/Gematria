int compare(const num * a, const num * b)
{
    extern Genome * seq;
    extern unsigned read_length;
    
    int k = 0;
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

void setCounts(num last, num k, num * positions, unsigned cnt){
    extern Genome * seq;
    for (num i = k + 1; i-- > last; ) {
        num pt = *(positions + i);
        if (pt < seq->size) seq->counts[pt] = cnt;
    }
}

void * sorter(void * data){
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
                setCounts(last, k, block->positions, cnt);
                cnt  = 0;
                last = k + 1;
            }
        }
        cnt++;
        setCounts(last, k, block->positions, cnt);
    }

    // pthread_mutex_lock(&mutex);
    // print debug
    // pthread_mutex_unlock(&mutex);

    pthread_exit(0);
}
