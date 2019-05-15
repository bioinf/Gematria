void setbit(void * field, big index)
{
    char * byte = NULL;
    byte = (char *) field + index / 8;
    *byte = *byte | (1 << (index % 8));
}

unsigned checkbit(const void * field, big index)
{
    char * byte = NULL;
    byte = (char *) field + index / 8;
    return (*byte >> (index % 8)) & 1;
}

void make_hashes(Filter * bloom, const Genome * seq, const big i)
{
    extern unsigned read_length;
	
	// MurmurHash64A
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;

    const uint64_t * data = (const uint64_t *) &seq->sequence[i];
    const uint64_t * end = data + (read_length/8);

    for (unsigned seed = 0; seed < bloom->hash_count; seed++) {
	    bloom->hashes[seed] = (uint64_t) seed ^ (read_length * m);
	}

    while(data != end)
    {
        uint64_t k = *data++;
        k *= m;
        k ^= k >> r;
        k *= m;

	    for (unsigned seed = 0; seed < bloom->hash_count; seed++) {
		    bloom->hashes[seed] ^= k;
		    bloom->hashes[seed] *= m;
		}
    }

    const unsigned char * data2 = (const unsigned char*) data;

    for (unsigned seed = 0; seed < bloom->hash_count; seed++) {
	    switch(read_length & 7)
	    {
	        case 7: bloom->hashes[seed] ^= (uint64_t) data2[6] << 48;
	        case 6: bloom->hashes[seed] ^= (uint64_t) data2[5] << 40;
	        case 5: bloom->hashes[seed] ^= (uint64_t) data2[4] << 32;
	        case 4: bloom->hashes[seed] ^= (uint64_t) data2[3] << 24;
	        case 3: bloom->hashes[seed] ^= (uint64_t) data2[2] << 16;
	        case 2: bloom->hashes[seed] ^= (uint64_t) data2[1] << 8;
	        case 1: bloom->hashes[seed] ^= (uint64_t) data2[0];
	        bloom->hashes[seed] *= m;
	    };

	    bloom->hashes[seed] ^= bloom->hashes[seed] >> r;
	    bloom->hashes[seed] *= m;
	    bloom->hashes[seed] ^= bloom->hashes[seed] >> r;
		bloom->hashes[seed] = bloom->hashes[seed] % (bloom->length * 8); 
    }
}

unsigned bloom_exist(Filter * bloom, const Genome * seq, const num i)
{
    make_hashes(bloom, seq, i);

    unsigned h = 0;
    unsigned exist = 1;
    for (h = 0; h < bloom->hash_count; h++) {
        if (!checkbit(bloom->second, bloom->hashes[h])) exist = 0;
    }
    
    return exist;
}

void bloom_insert(Filter * bloom, Genome * seq, big i)
{
    make_hashes(bloom, seq, i);

    unsigned h = 0;
    unsigned exist = 1;
    for (h = 0; h < bloom->hash_count; h++) {
        if (!checkbit(bloom->first, bloom->hashes[h])) exist = 0;
        setbit(bloom->first, bloom->hashes[h]);
    }

    if (exist) {
        for (h = 0; h < bloom->hash_count; h++) {
            setbit(bloom->second, bloom->hashes[h]);
        }
    }
}
