#!/usr/bin/env python3
import os

import makegms
import numpy as np

from include.argparse import *
from include.getcontents import *
from include.write import *

# --------------------------------------------------------------------------- #
app.log('Gematria: Initialization')
igvtools, bed2bigbed = check_exe(__file__)

app.log('Get contents of fasta file: ' + app.argx['input'])
fasta = getcontents(app.argx['input'])

app.log('Making raw GMS-track')
track = makegms.run(app.argx['input'],
                    read=app.argx['length'],
                    threads=int(app.argx['threads']))

fs = {}
for k in outputs:
    if k in ['wig', 'bw', 'bed']:
        fs[k] = Write(outputs[k], k, fasta)

if 'tdf' in outputs and 'wig' not in fs:
    outputs['wig'] = '.temporary.wig'
    fs['wig'] = Write(outputs['wig'], 'wig', fasta)

if 'bigbed' in outputs and 'bed' not in fs:
    outputs['bed'] = '.temporary.bed'
    fs['bed'] = Write(outputs['bed'], 'bed', fasta)

app.log('Splitting raw GMS-track to chromosomes')
i = 0
for chr, lng in fasta:
    app.log(chr, ' >  Chr: ')
    reads = lng - app.argx['length'] + 1
    subseq = track[i:i+reads]
    i += lng

    unique = np.convolve(subseq, kernel)
    left = np.concatenate(([0] * mdist, unique))[:reads]
    right = np.concatenate((unique, [0] * mdist))[-reads:]

    final = subseq + (np.ones(reads) - subseq) * (left + right) / 2
    gms = np.convolve(final, [100/app.argx['length']] * app.argx['length'])

    for key in fs:
        fs[key].add(chr, np.append(np.round(gms), [-1]))

for key in fs:
    fs[key].h.close()

# --------------------------------------------------------------------------- #
sizes = False
if 'tdf' in outputs or 'bigbed' in outputs:
    sizes = '.temporary.chrom.sizes'
    data = ["{0}\t{1}".format(chr.split(' ')[0], lng) for chr, lng in fasta]
    with open(sizes, 'w+') as output:
        output.write('\n'.join(data))

if 'tdf' in outputs:
    app.log('Converting wig to tdf')

    os.system((' ').join([
      "java -Djava.awt.headless=true -Xmx1500m",
      "-jar {app}".format(app=igvtools),
      "toTDF {wig} {tdf}".format(wig=fs['wig'].h.name, tdf=outputs['tdf']),
      sizes, "> /dev/null 2>&1"
    ]))
    
    os.remove('igv.log')
    if fs['wig'].h.name[0] == '.':
        os.remove(fs['wig'].h.name)

if 'bigbed' in outputs:
    app.log('Converting bed to bigbed')

    os.system((' ').join([
      "sort -k1,1 -k2,2n",
      "{bed} > {bed}.sorted".format(bed=outputs['bed'])
    ]))

    os.system((' ').join([
      bed2bigbed, "-type=bed9",
      "{bed}.sorted".format(bed=outputs['bed']),
      sizes, outputs['bigbed'],
      "> /dev/null 2>&1"
    ]))

    if fs['bed'].h.name[0] == '.':
        os.remove(fs['bed'].h.name)
        os.remove(fs['bed'].h.name + '.sorted')

if sizes:
    os.remove(sizes)
    

app.success_log('Done')
