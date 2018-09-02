class Write():
    wig = "fixedStep chrom={0} start={1} step=1 span={2}\n{3}\n"
    bed = "{0}\t{1}\t{2}\t.\t{3}\t.\t{1}\t{2}\t{4}\n"
    stp = -1

    def __init__(self, file, fasta):
        output, self.ext = file
        init = getattr(self, 'init_' + self.ext)
        init(output + '.' + self.ext, fasta)

    # ----------------------------------------------------------------------- #
    def init_wig(self, f, fasta):
        self.h = open(f, 'w')

    def init_bed(self, f, fasta):
        self.h = open(f, 'w')

    def init_bw(self, f, fasta):
        import pyBigWig as bw
        self.h = bw.open(f, 'w')
        self.h.addHeader([(chr, lng) for chr, lng in fasta])

    # ----------------------------------------------------------------------- #
    def _wig(self, chr, pos, span, val):
        if span == self.stp:
            return self.h.write(str(val) + "\n")
        self.stp = span
        self.h.write(self.wig.format(chr, pos, span, val))

    def _bed(self, chr, pos, span, val):
        color = '{0},{0},255'.format(str(int(255 - round(255 * val/100))))
        self.h.write(self.bed.format(chr, pos-1, pos-1+span, int(val), color))

    def _bw(self, chr, pos, span, val):
        self.h.addEntries(chr, [pos-1], values=[val], span=span, step=1)

    # ----------------------------------------------------------------------- #
    def add(self, chr, gms):
        repeats = 1
        repeats_init = 0

        for i in range(len(gms) - 1):
            if gms[i] == gms[i+1]:
                repeats += 1
                continue
            add = getattr(self, '_' + self.ext)
            add(chr, repeats_init + 1, repeats, gms[i])
            repeats = 1
            repeats_init = i + 1
