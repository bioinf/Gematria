// ./exec genome.fa 100 16 5 -> Qsort / Bloom  (auto, вepends on size)
// ./exec genome.fa 100 0  5 -> ! Bloom        (force)
// ./exec genome.fa 100 16 0 -> ! Qsort        (force)

void makeGmsWrapper(char * src, unsigned threads, unsigned quality)
{
    extern Genome * seq;
    extern unsigned read_length;
	extern unsigned large_genome;
	
    seq = reader(src);
    seq->reads = seq->length/2 - read_length + 1;
	seq->counts = (unsigned short *) malloc (sizeof(unsigned short) * seq->reads);
	
    // unsigned int max value = 4,294,967,295
    // genome max length = 2147483647
	if (seq->length > 2147483647) large_genome = 1;

    if (threads == 0 && quality > 0) {
		makeGmsBloom(quality);
	
	} else if (quality == 0 && threads > 0) {
		makeGmsQsort(threads);
	
	} else {
        if (large_genome) {
			makeGmsBloom(quality);

        } else {
			if (threads > 1 || quality > 5) {
				makeGmsQsort(threads);

			} else {
				makeGmsBloom(quality);
			}
        }
	}
}
