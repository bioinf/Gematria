#!/usr/bin/env python3
import os

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

    def log(self, text, prefix='    '):
        self.echo(prefix + text + '\n', 'green')

    def error_log(self, text, prefix='[-] '):
        self.echo(prefix, 'red_bold')
        self.echo(text + '\n', 'red')

    def success_log(self, text, prefix='[+] '):
        self.echo(prefix, 'green_bold')
        self.echo(text + '\n', 'green')

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
  ['-o', '--output', 'Output filenames without extension'],
  ['-f', '--formats', 'Comma separated output formats',
                      'Acceptable formats: wig, bigwig, bed, tdf, bigbed'],
  ['-r', '--reads', 'Reads type parameters in the following format:',
                    'S - for single-end reads',
                    'N:mu:sigma - for Normal distribution of insertion size',
                    'U:min:max - for Uniform distribution of insertion size'],
  ['-h', '--help', 'Show this help']]
demo = [
  'input.fasta',
  'input.fasta -l 35 -t 4 -o result -f bed,tdf -r U:0:80',
  '-i input.fasta -l 50',
  '-i input.fasta -r N:40:20']

app = App(init, args, demo)

# --------------------------------------------------------------------------- #
app.default('input', sys.argv[1])
if not os.path.isfile(app.argx['input']):
    app.exit('Input file not found: ' + app.argx['input'])

app.default('output', app.argx['input'])
app.default('formats', 'wig')

outputs = {}
for ext in app.argx['formats'].split(','):
    ext = ext.replace('bigwig', 'bw')
    if ext not in ['wig', 'bed', 'tdf', 'bw', 'bigbed']:
        continue
    outputs[ext] = app.argx['output'] + '.' + ext

if 'wig' in outputs:
    try:
        import pyBigWig as bw
    except ImportError:
        stop('Python module pyBigWig not found')

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
def download(src, dst, iexec=False):
    import urllib.request

    app.log('Downloading: ' + src)
    res = urllib.request.urlopen(src)
    with open(dst, 'wb') as output:
        output.write(res.read())

    if os.path.isfile(dst):
        app.success_log('Download complete: ' + dst)
        if iexec:
            import stat
            os.chmod(dst, os.stat(dst).st_mode | stat.S_IEXEC)
        return True

    app.error_log('Download error: ' + dst)
    return False


def check_exe(root):
    root = os.path.dirname(os.path.abspath(root))

    ucsc = "http://hgdownload.cse.ucsc.edu/admin/exe/"
    bed2bigbed = "{root}/exe/bedToBigBed".format(root=root)

    if sys.platform == "darwin":
        src = ucsc + "macOSX.x86_64/bedToBigBed"
        bed2bigbed += "_darwin"

    if sys.platform == "linux":
        src = ucsc + "linux.x86_64/bedToBigBed"
        bed2bigbed += "_linux"

    if 'bigbed' in outputs:
        if os.path.isfile(bed2bigbed) and os.access(bed2bigbed, os.X_OK):
            app.success_log('Executable `bedToBigBed` found')
        else:
            app.error_log('Executable `bedToBigBed` not found')
            if not download(src, bed2bigbed, True):
                del(outputs['bigbed'])

    repo = "https://github.com/evgeny-bakin/GeMaTrIA/"
    igvtools = "{root}/exe/igvtools.jar".format(root=root)
    src = repo + "raw/master/exe/igvtools.jar"

    if 'tdf' in outputs:
        if os.path.isfile(igvtools):
            app.success_log('Executable `igvtools` found')
        else:
            app.error_log('Executable `igvtools` not found')
            if not download(src, igvtools):
                del(outputs['tdf'])

    return [igvtools, bed2bigbed]
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

    def __init__(self, name, ext, fasta):
        self.ext = ext
        init = getattr(self, 'init_' + self.ext)
        init(name, fasta)

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
        chr = chr.split(' ')[0]
        self.stp = [span, chr]
        self.h.write(self.wig.format(chr, pos, span, val))

    def _bed(self, chr, pos, span, val):
        chr = chr.split(' ')[0]
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
