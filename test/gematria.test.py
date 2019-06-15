#!/usr/bin/env python3
import os

root = os.path.dirname(os.path.abspath(__file__)) + '/../'
init = "{root}gematria.py {root}test/example.fa".format(root=root)
coli = "{root}gematria.py {root}test/ecoli.fa".format(root=root)


def run(cmd):
    print('-' * 79)
    os.system((' ').join(cmd).format(root=root))

run(["{root}gematria.py --help"])

run([init, "-l 11 -o wig,bigwig,bed,tdf,bigbed"])
run([init, "-l 12 -o bigwig,bigbed"])
run([init, "-l 20 -o bed"])

run(["rm -f _tests* && echo \"+ _tests* removed\""])

run([coli, "-l 12 -o tdf -t 2"])
run([coli, "-l 16 -o tdf -m soft"])

run(["rm -f _coli* && echo \"+ _coli* removed\""])
