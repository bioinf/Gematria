char complement(char ch)
{
    if (ch == 'A') return 'T';
    if (ch == 'T') return 'A';
    if (ch == 'G') return 'C';
    if (ch == 'C') return 'G';
    return ch;
}

void readGenome(char * src)
{
    extern Genome * seq;
    extern unsigned read_length;

    FILE * f = fopen(src, "r");
    if (f == NULL) exit(EXIT_FAILURE);

    fseek(f, 0, SEEK_END);
    
    unsigned long long bytes = (unsigned long long) ftell(f);

    seq = (Genome*) malloc(sizeof(Genome));
    seq->sequence = (char *) malloc(sizeof(char) * bytes * 2);

    fseek(f, 0, SEEK_SET);
    
    char ch;
    char state = 'N';
    num i = 0;
    
    while ((ch = fgetc(f)) != EOF) {
        if (state == 'N' && ch == '>') {
            state = 'T'; // is title
        } else if (ch == '\n') {
            state = 'N'; // is newline
        } else if (state != 'T') {
            seq->sequence[i++] = toupper(ch);
        }
    }
    fclose(f);
    
    num r = i;
    while (i > 0) seq->sequence[r++] = complement(seq->sequence[--i]);

    seq->sequence[r] = '\0';
    seq->length = r;
    seq->counts = (unsigned short *) malloc (sizeof(unsigned short) * r);
    seq->size = seq->length/2 - read_length + 1;
}
