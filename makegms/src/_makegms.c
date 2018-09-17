void makegms(char * src, unsigned int threads)
{
    extern Genome * seq;
    extern unsigned read_length;

    readGenome(src);
    seq->size = seq->length/2 - read_length + 1;
    
    if (threads > seq->length/10) threads = seq->length/10;
    if (threads == 0) threads = 1;

    ThreadData * parts = decomposition(seq, read_length, threads);
    pthread_t * thread = (pthread_t*) malloc(threads * sizeof(pthread_t));

    for (unsigned t = 0; t < threads; t++){
        pthread_create( &(thread[t]), NULL, sorter, &parts[t]);
    }
    for (unsigned t = 0; t < threads; t++){
        pthread_join(thread[t], NULL);
    }

    for (unsigned t = 0; t < threads; t++){
        free(parts[t].positions);
    }

    free(parts);
    free(thread);
    free(seq->sequence);
}
