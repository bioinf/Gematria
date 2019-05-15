void makeGmsBloom(unsigned quality)
{
    printf("\033[32m[+] Using a BloomFilter-based approach\033[0m\n");

    extern Genome * seq;
    extern unsigned read_length;

    Filter * bloom = (Filter *) malloc(sizeof(Filter));
	
	if (quality < 3) quality = 3;
	if (quality > 9) quality = 9;

    bloom->hash_count = round(2.77259 * quality); // Log[2] * 4 = 2.77259
    bloom->hashes = (big *) malloc(sizeof(big) * bloom->hash_count);

    bloom->length = seq->reads * quality;
    bloom->first = malloc(bloom->length);
    bloom->second = malloc(bloom->length);

    memset(bloom->first, 0, bloom->length);
    memset(bloom->second, 0, bloom->length);
	
	// Write hashes to filter
    for (big k = 0; k != seq->length - read_length; k++) {
        bloom_insert(bloom, seq, k);
    }
	
	// Export hashes from filter
    for (big k = 0; k != seq->reads; ++k) {
		seq->counts[k] = bloom_exist(bloom, seq, k) > 0 ? 2 : 1;
    }

    free(bloom->first);
    free(bloom->second);
    free(bloom->hashes);
    free(bloom);
}
