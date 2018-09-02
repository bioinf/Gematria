#!/usr/bin/env python3
import makegms
import numpy as np

from include.argparse import *
from include.getcontents import *
from include.write import *

# --------------------------------------------------------------------------- #

fasta = getcontents(app.argx['input'])
track = makegms.run(app.argx['input'],
                    read=app.argx['length'],
                    threads=int(app.argx['threads']))

h = Write([file, ext], fasta)
i = 0

for chr, lng in fasta:
    reads = lng - app.argx['length'] + 1
    subseq = track[i:i+reads]
    i += lng

    unique = np.convolve(subseq, kernel)
    left = np.concatenate(([0] * mdist, unique))[:reads]
    right = np.concatenate((unique, [0] * mdist))[-reads:]

    final = subseq + (np.ones(reads) - subseq) * (left + right) / 2
    gms = np.convolve(final, [100/app.argx['length']] * app.argx['length'])

    h.add(chr, np.append(np.round(gms), [-1]))
