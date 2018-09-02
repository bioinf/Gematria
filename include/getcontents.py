def getcontents(filename):
    chr = ['', 0]
    fasta = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            if line == "":
                continue
            if line[0] == '>':
                if chr[1] > 0:
                    fasta.append(chr)
                chr = [line[1:], 0]
            else:
                chr[1] += len(line)
        fasta.append(chr)

    return fasta
