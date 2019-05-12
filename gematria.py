#!/usr/bin/env python3
import os
import time

import makegms
import numpy as np


from include.argparse import *
from include.write import *


started = time.time()

app.intro()
igvtools, bed2bigbed = check_exe(__file__)

if len(outputs) == 0:
    app.exit('The specified export formats cannot be generated.')

# --------------------------------------------------------------------------- #

begin = time.time()
app.log('Get contents of fasta file: ' + app.argx['input'])
fasta = app.fasta()
app.success_log('File loaded: {0:.2f}sec.'.format(time.time() - begin))

# --------------------------------------------------------------------------- #

fs = {}
for k in outputs:
    if k in ['wig', 'bw', 'bed']:
        fs[k] = Write(outputs[k], k)

if 'tdf' in outputs and 'wig' not in fs:
    outputs['wig'] = '.__temporary.wig'
    fs['wig'] = Write(outputs['wig'], 'wig')

if 'bigbed' in outputs and 'bed' not in fs:
    outputs['bed'] = '.__temporary.bed'
    fs['bed'] = Write(outputs['bed'], 'bed')

if 'bw' in outputs:
    fs['bw'].h.addHeader([(chr, size) for chr, size, name in fasta])

# --------------------------------------------------------------------------- #

begin = time.time()
app.log('Making raw GMS-track')
track = makegms.run(app.argx['input'],
                    read=app.argx['length'],
                    quality=int(app.argx['quality']),
                    threads=int(app.argx['threads']))

app.success_log('GMS-track is created: {0:.2f}sec.'.format(time.time()-begin))

# --------------------------------------------------------------------------- #

app.log('Splitting raw GMS-track to chromosomes')
mask = 100 * np.ones(app.argx['length']) / app.argx['length']
index = 0
for chr, lng, name in fasta:
    begin = time.time()
    app.echo(' -  Chr: ' + name + '\r', 'green')
    
    reads = lng - app.argx['length'] + 1
    subseq = track[index:index+reads]
    index += lng

    if app.argx['reads'][0] == 'S':
        final = subseq
    else:
        unique = np.convolve(subseq, kernel)
        left = np.concatenate((np.zeros(mdist), unique))[:reads]
        right = np.concatenate((unique, np.zeros(mdist)))[-reads:]
        final = subseq + (np.ones(reads) - subseq) * (left + right) / 2

    gms = np.append(np.round(np.convolve(final, mask)), [-1])
    gms = np.insert(gms, 0, -1)
    zer = np.nonzero(np.convolve(gms, np.array([-1, 1]))[2:-1])[0]
    for key in fs:
        fs[key].add(chr, gms, zer)

    app.echo('[+] Chr: ' + name, 'green')
    app.echo(' [{0:.2f}sec.]\n'.format(time.time()-begin), 'green_bold')

for key in fs:
    fs[key].h.close()


# --------------------------------------------------------------------------- #

size_ = False
if 'tdf' in outputs or 'bigbed' in outputs:
    size_ = '.__temporary.chrom.sizes'
    data = ["{0}\t{1}".format(chr, size) for chr, size, name in fasta]
    with open(size_, 'w+') as output:
        output.write('\n'.join(data))

if 'tdf' in outputs:
    begin = time.time()
    app.log('Converting wig to tdf')

    os.system((' ').join([
      "java -Djava.awt.headless=true -Xmx1500m",
      "-jar {app}".format(app=igvtools),
      "toTDF {wig} {tdf}".format(wig=fs['wig'].h.name, tdf=outputs['tdf']),
      size_, "> /dev/null 2>&1"
    ]))
    
    os.remove('igv.log')
    if fs['wig'].h.name[0:3] == '.__':
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
      size_, outputs['bigbed'],
      "> /dev/null 2>&1"
    ]))

    if fs['bed'].h.name[0:3] == '.__':
        os.remove(fs['bed'].h.name)

    os.remove(fs['bed'].h.name + '.sorted')

    app.success_log('Done: {0:.2f}sec.'.format(time.time()-begin))

if size_:
    os.remove(size_)

app.echo('\nGematria has finished. Results:\n', 'white_bold')
for f in outputs:
    app.echo('[#] ' + outputs[f] + '\n', 'white')

app.echo('\nElapsed time: ', 'white_bold')
app.echo('{0:.2f}sec.\n'.format(time.time()-started), 'white_bold')

# --------------------------------------------------------------------------- #

if app._debug:
    time.sleep(1)
    inf = 'cat /proc/{pid}/status | grep "VmHWM" | xargs'

    m = os.popen(inf.format(pid=os.getpid()))
    memory = int(m.read().split(' ')[1])/1024
    m.close()

    fsize = os.path.getsize(app.argx['input'])/1024/1024

    h = open(app._debug, 'a+')
    h.write(('  ').join([
      'Length: {0: <3}',
      'Filesize [MB]: {1: <5.2f}',
      'Memory [MB]: {2: <8.0f}',
      'Read Type: {3}',
      'Filename: {4}\n'
    ]).format(app.argx['length'], fsize, memory, 
              app.argx['reads'], app.argx['input']))
    h.close()
# --------------------------------------------------------------------------- #
