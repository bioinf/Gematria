#!/usr/bin/env python3
import os

root = os.path.dirname(os.path.abspath(__file__)) + '/../'
init = "{root}gematria.py {root}test/example.fa".format(root=root)
coli = "{root}gematria.py {root}test/ecoli.fa".format(root=root)


def run(cmd):
    print('-' * 79)
    os.system((' ').join(cmd).format(root=root))

run(["{root}gematria.py --help"])

run([init, "-l 11 -o {root}_tests.output.1 -f wig,bigwig,bed,tdf,bigbed"])
run([init, "-l 12 -o {root}_tests.output.2 -f bigwig,bigbed"])
run([init, "-l 20 -o {root}_tests.output.3 -f bed"])

run(["rm -f _tests* && echo \"+ _tests* removed\""])

run([coli, "-l 12 -o {root}_coli-L12 -f tdf"])
run([coli, "-l 16 -o {root}_coli-L16 -f tdf"])

run(["rm -f _coli* && echo \"+ _coli* removed\""])
