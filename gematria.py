#!/usr/bin/env python3
import os
import time

import makegms
import numpy as np

started = time.time()

from include.argparse import *
from include.getcontents import *
from include.write import *

# --------------------------------------------------------------------------- #
app.intro()
igvtools, bed2bigbed = check_exe(__file__)

begin = time.time()
app.log('Get contents of fasta file: ' + app.argx['input'])
fasta = getcontents(app.argx['input'])
app.success_log('File loaded: {0:.2f}sec.'.format(time.time()-begin))


begin = time.time()
app.log('Making raw GMS-track')
track = makegms.run(app.argx['input'],
                    read=app.argx['length'],
                    threads=int(app.argx['threads']))
app.success_log('GMS-track is created: {0:.2f}sec.'.format(time.time()-begin))


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
    begin = time.time()
    app.echo(' #  Chr: ' + chr, 'green')
    
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

    app.echo(' [{0:.2f}sec.]\n'.format(time.time()-begin), 'green_bold')

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
    begin = time.time()
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
    
    app.success_log('Done: {0:.2f}sec.'.format(time.time()-begin))

if 'bigbed' in outputs:
    begin = time.time()
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

    app.success_log('Done: {0:.2f}sec.'.format(time.time()-begin))

if sizes:
    os.remove(sizes)
    
app.echo('\nGematria has finished.\nElapsed time: ', 'white_bold')
app.echo('{0:.2f}sec.\n'.format(time.time()-started), 'white_bold')
