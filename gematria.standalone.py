#!/usr/bin/env python3
import makegms
import numpy as np

# --------------------------------------------------------------------------- #
# filename: include/argparse.py


import os
import sys

# --------------------------------------------------------------------------- #
# filename: include/app.py


import sys


class App():
    def __init__(self, init="", args=[], demo=[]):
        self.init = init
        self.demo = [[sys.argv[0], e] for e in demo]
        self.parse(args)

    def parse(self, args):
        names = {}
        self.args = []
        self.argx = {}
        for arg in args:
            name = arg[1].replace('-', '')
            names[arg[0]] = name
            names[arg[1]] = name
            self.argx[name] = None
            self.args.append(["{0}, {1} ".format(arg[0], arg[1]), arg[2]])
            [self.args.append(['', line]) for line in arg[3:]]

        argv = [names[i] if i in names else i for i in sys.argv]

        if len(argv) <= 1:
            self.exit('Kio okazas? (Specify arguments. Please)')

        if argv[1] == 'help':
            self.exit()

        first = -1 * len(argv) % 2
        for i in range(0, len(argv) - 1, 2):
            k, v = [argv[i + first], argv[i + first + 1]]
            if k in self.argx:
                self.argx[k] = v

    def echo(self, text, color=""):
        sys.stderr.write({
          '': "{0}",
          'red': "\033[31m{0}\033[0m",
          'green': "\033[32m{0}\033[0m",
          'white': "\033[37m{0}\033[0m",
          'red_bold': "\033[1;31m{0}\033[0m",
          'green_bold': "\033[1;32m{0}\033[0m"
        }[color].format(text))

    def params(self, items):
        space = max([len(name) for name, desc in items]) + 1
        for name, desc in items:
            name += ' ' * (space - len(name))
            self.echo('  ' + name, 'green_bold')
            self.echo(desc + '\n', 'green')
        self.echo('\n')

    def default(self, name, value):
        if name not in self.argx:
            self.argx[name] = value
        if self.argx[name] is None:
            self.argx[name] = value

    def exit(self, cause=False):
        if cause is not False:
            self.echo('Error:\n', 'red')
            self.echo('  ' + cause + '\n\n', 'red_bold')

        if self.init is not "":
            self.echo('Usage:\n', 'white')
            self.params([[sys.argv[0], self.init]])

        if len(self.args) > 0:
            self.echo('Optional arguments:\n', 'white')
            self.params(self.args)

        if len(self.demo) > 0:
            self.echo('Examples:\n', 'white')
            self.params(self.demo)

        sys.exit(1 if cause else 0)

init = '[fasta file] [Optional arguments]'
args = [
  ['-i', '--input', 'Path to `genome.fasta` file'],
  ['-l', '--length', 'Read length. Default: 100'],
  ['-t', '--threads', 'Number of threads. Default: auto'],
  ['-o', '--output', 'Output filename with an extension wig, bw or bed'],
  ['-r', '--reads', 'Reads type parameters in the following format:',
                    'S - for single-end reads',
                    'N:mu:sigma - for Normal distribution of insertion size',
                    'U:min:max - for Uniform distribution of insertion size'],
  ['-h', '--help', 'Show this help']]
demo = [
  'input.fasta',
  'input.fasta -l 35 -t 4 -o result.bed -r U:0:80',
  '-i input.fasta -l 50',
  '-i input.fasta -r N:40:20']

app = App(init, args, demo)

# --------------------------------------------------------------------------- #
app.default('input', sys.argv[1])
if not os.path.isfile(app.argx['input']):
    app.exit('Input file not found: ' + app.argx['input'])

app.default('output', app.argx['input'] + '.wig')
file, ext = os.path.splitext(app.argx['output'])
ext = ext[1:]
if ext not in ['wig', 'bw', 'bed']:
    ext = 'wig'
if ext == 'bw':
    try:
        __import__('imp').find_module('pyBigWig')
    except ImportError:
        stop('Python module pyBigWig not found')
    import pyBigWig as bw

app.default('reads', 'S')
try:
    # Single-end reads
    if app.argx['reads'][0] == 'S':
        mdist, kernel = [0, [1]]

    # Normal distribution of insertions size
    if app.argx['reads'][0] == 'N':
        from math import sqrt, pi, exp
        mu, s = map(int, app.argx['reads'][2:].split(':'))

        def k(i):
            return exp(-(float(i) - mu)**2 / (2 * s**2)) / (sqrt(2*pi)*s)

        mdist = mu - 3 * s
        ker = [k(i) for i in range(mdist, mu + 3 * s + 1)]
        mdist, kernel = [mdist if mdist > 0 else 0, ker]

    # Uniform distribution of insertion size
    # U:min:max
    if app.argx['reads'][0] == 'U':
        v_min, v_max = map(int, app.argx['reads'][2:].split(':'))
        ker = [1/(v_max - v_min)] * (v_max - v_min)
        mdist, kernel = [v_min, ker]

except:
    app.exit('Unable to parse reads type parameters')

app.default('threads', os.sysconf('SC_NPROCESSORS_ONLN'))
app.default('length', 100)
app.argx['length'] = int(app.argx['length'])
# --------------------------------------------------------------------------- #
# filename: include/getcontents.py


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
# --------------------------------------------------------------------------- #
# filename: include/write.py


class Write():
    wig = "fixedStep chrom={0} start={1} step=1 span={2}\n{3}\n"
    bed = "{0}\t{1}\t{2}\t.\t{3}\t.\t{1}\t{2}\t{4}\n"
    stp = ['', -1]

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
        if [span, chr] == self.stp:
            return self.h.write(str(val) + "\n")
        self.stp = [span, chr]
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
