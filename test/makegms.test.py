#!/usr/bin/env python3
import sys
import os

import makegms
import numpy as np


class Test():
    def __init__(self, desc, seq, read=10, quality=3):
        sys.stdout.write("\033[32m- {0}\033[0m".format(desc))
        self.write(seq)
        track = makegms.run(self.fp.name, read=read, quality=quality)
        self.validate(track, seq, read)
        os.remove(self.fp.name)

    def write(self, seq):
        txt = '\n'.join([seq[i:i+100] for i in range(0, len(seq), 100)])
        self.fp = open('./temp.fasta', 'w+')
        self.fp.write(">TMP\n{0}\n".format(txt))
        self.fp.close()

    def validate(self, track, seq, read):
        valid = 0
        for i in range(len(seq) - read + 1):
            if (seq[i:i+read] in seq) is not track[i]:
                valid += 1

        col = "31" if (valid != len(track)) else "32"
        tpl = " \033[1;" + col + "m:: {0}% valid: {1}/{2}\033[0m "
        tpl += "\033[32mSequence length: {3}bp\033[0m\n"

        pc = int(valid/len(track) * 100)
        sys.stdout.write(tpl.format(pc, valid, len(track), len(seq)))


class Unique(type):
    _seed = 1

    def __new__(cls, size, block=10, **kwargs):
        def rebase(num):
            seq = []
            while num:
                num, r = divmod(num, 3)
                seq.append(['A', 'T', 'G'][r])
            return 'C' * (block - len(seq)) + ''.join(reversed(seq))

        seq = ''
        while (len(seq) < size):
            seq += rebase(Unique._seed)
            Unique._seed += 1

        return seq


# --------------------------------------------------------------------------- #

print("\n- Repeats only")

Test('Track #1 | Bloom ', 'T' * 1000)
Test('Track #2 | Bloom ', 'AT' * 50000)
Test('Track #3 | Bloom ', 'GGGGC' * 100000)

Test('Track #1 | Qsort ', 'T' * 1000, 10, 0)
Test('Track #2 | Qsort ', 'AT' * 50000, 10, 0)
Test('Track #3 | Qsort ', 'GGGGC' * 100000, 10, 0)

print("\n- No repeats")

Test('Track #4 | Bloom ', Unique(100))
Test('Track #5 | Bloom ', Unique(1000))
Test('Track #6 | Bloom ', Unique(10000))

Test('Track #4 | Qsort ', Unique(100), 10, 0)
Test('Track #5 | Qsort ', Unique(1000), 10, 0)
Test('Track #6 | Qsort ', Unique(10000), 10, 0)

print("\n- Mixed reads")

Test('Track #7 | Bloom ', ('A' * 500) + Unique(500) + ('GC' * 250))
Test('Track #8 | Bloom ', Unique(500) + ('T' * 10000) + Unique(50))
Test('Track #9 | Bloom ', Unique(100) + ('AAAT' * 500))

Test('Track #7 | Qsort ', ('A' * 500) + Unique(500) + ('GC' * 250), 10, 0)
Test('Track #8 | Qsort ', Unique(500) + ('T' * 10000) + Unique(50), 10, 0)
Test('Track #9 | Qsort ', Unique(100) + ('AAAT' * 500), 10, 0)
