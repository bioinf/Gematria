unsigned int MurmurHash2(char * key, unsigned len)
{
    const unsigned m = 0x5bd1e995;
    const int r = 24;
    unsigned h = 1 ^ len;
    const unsigned char * data = (const unsigned char *) key;

    unsigned k;
    while (len >= 4) {
        k  = data[0];
        k |= data[1] << 8;
        k |= data[2] << 16;
        k |= data[3] << 24;
        k *= m;
        k ^= k >> r;
        k *= m;
        h *= m;
        h ^= k;
        data += 4;
        len -= 4;
    }

    switch (len) {
        case 3: h ^= data[2] << 16;
        case 2: h ^= data[1] << 8;
        case 1: h ^= data[0]; h *= m;
    };

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
}

void append(ThreadData * item, big k)
{
    if (item->pointer >= item->size) {
        item->size += 50;
        item->positions = (num*) realloc (item->positions, sizeof(num) * item->size);
    }
    item->positions[item->pointer++] = k;
}

void append_big(ThreadDataBig * item, big k)
{
    if (item->pointer >= item->size) {
        item->size += 50;
        item->positions = (big*) realloc (item->positions, sizeof(big) * item->size);
    }
    item->positions[item->pointer++] = k;
}


ThreadData * decomposition(Genome * seq, unsigned read_length, unsigned threads)
{
    ThreadData * data = (ThreadData*) malloc(threads * sizeof(ThreadData));

    for (unsigned t = 0; t < threads; t++){
        data[t].pointer   = 0;
        data[t].size      = seq->length/threads;
        data[t].positions = (num*) malloc (sizeof(num) * data[t].size);
    }

    big reads = seq->length - read_length + 1;
    for (big k = 0; k < reads; k++) {
        append(&data[ MurmurHash2(seq->sequence + k, read_length) % threads ], k);
    }

    return data;
}

ThreadDataBig * decomposition_big(Genome * seq, unsigned read_length, unsigned threads)
{
    ThreadDataBig * data = (ThreadDataBig*) malloc(threads * sizeof(ThreadDataBig));

    for (unsigned t = 0; t < threads; t++){
        data[t].pointer   = 0;
        data[t].size      = seq->length/threads;
        data[t].positions = (big*) malloc (sizeof(big) * data[t].size);
    }

    big reads = seq->length - read_length + 1;
    for (big k = 0; k < reads; k++) {
        append_big(&data[ MurmurHash2(seq->sequence + k, read_length) % threads ], k);
    }

    return data;
}
