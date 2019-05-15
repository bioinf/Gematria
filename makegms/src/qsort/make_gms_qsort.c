void makeGmsQsort(unsigned threads)
{
    printf("\033[32m[+] Using an approach based on QuickSort\033[0m\n");

    extern Genome * seq;
    extern unsigned read_length;
	extern unsigned large_genome;

    if (threads > seq->length/100) threads = seq->length/100;
    if (threads < 1) threads = 1;
	
    pthread_t * thread = (pthread_t*) malloc(threads * sizeof(pthread_t));

	if (large_genome) {
		ThreadDataBig * parts_big = decomposition_big(seq, read_length, threads);

	    for (unsigned t = 0; t < threads; t++){
	        pthread_create( &(thread[t]), NULL, sorter_big, &parts_big[t]);
	    }
	    for (unsigned t = 0; t < threads; t++){
	        pthread_join(thread[t], NULL);
	    }
	    for (unsigned t = 0; t < threads; t++){
	        free(parts_big[t].positions);
	    }

	    free(parts_big);
	} else {
		
		ThreadData * parts = decomposition(seq, read_length, threads);

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
	}

    free(thread);
    free(seq->sequence);
}
