#!/usr/bin/env python3
import os
import time

import makegms
import numpy as np

# --------------------------------------------------------------------------- #
# filename: include/argparse.py


import sys
import os

# --------------------------------------------------------------------------- #
# filename: include/app.py


import sys


class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)


class App():
    def __init__(self, init="", args=[], demo=[]):
        self._debug = False
        self.version = 'v1.0 (cli) (built: 16.06.2019)'
        self.stderr = Unbuffered(sys.stderr)
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

        if argv[1] == 'version':
            self.echo('Gematria {}\n'.format(self.version), 'white_bold')
            sys.exit()

        if argv[1] == 'help':
            self.exit()

        first = -1 * len(argv) % 2
        for i in range(0, len(argv) - 1, 2):
            k, v = [argv[i + first], argv[i + first + 1]]
            if k in self.argx:
                self.argx[k] = v

    def intro(self):
        self.echo('Gematria\nCommand executed:\n', 'white_bold')
        self.echo('{app} {genome}'.format(app=sys.argv[0], genome=self.argx['input']), 'white')
        if self.argx['quality'] == 0:
            show = ['length', 'output', 'threads', 'paired']
        else:
            show = ['length', 'output', 'paired', 'lowmem']
        for nm in show:
            self.echo(" --{0} {1}".format(nm, self.argx[nm]), 'white')
            
        self.echo('\n\n')
        
    def fasta(self):
        chr = ['', 0]
        _fasta = []
    
        with open(self.argx['input'], 'r') as f:
            for line in f:
                line = line.replace('\n', '')
                if line == "":
                    continue
                if line[0] == '>':
                    if chr[1] > 0:
                        _fasta.append(chr)
                    chr = [line[1:], 0]
                else:
                    chr[1] += len(line)
            _fasta.append(chr)
        
        sizes = [lng for chr, lng in _fasta]
        names = [chr for chr, lng in _fasta]
        short = [chr.split(' ')[0] for chr, lng in _fasta]

        lng, chr, name = zip(*sorted(zip(sizes, short, names), reverse=True))
        return list(zip(chr, lng, name))

    def echo(self, text, color=""):
        self.stderr.write({
          '': "{0}",
          'red': "\033[31m{0}\033[0m",
          'green': "\033[32m{0}\033[0m",
          'white': "\033[37m{0}\033[0m",
          'white_bold': "\033[1;37m{0}\033[0m",
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
  ['-l', '--length', 'Read length. Default: 100'],
  ['-t', '--threads', 'Number of threads. Default: auto'],
  ['-o', '--output', 'Output formats: [wig],bigwig,bed,tdf,bigbed,all'],
  ['-p', '--paired', 'To use paired end reads, specify the insert size',
                     'and the standard deviation (normal distribution model)',
                     'If you don’t know where to start, use -p 300,100'],
  ['-m', '--lowmem', 'Use the RAM-memory saving algorithm. Default: none',
                     '  -m hard / will be used ~ 5 × Genome Size',
                     '  -m soft / will be used ~ 7 × Genome Size',
                     'Specifying this parameter in the value hard or soft ',
                     'makes it impossible to use multithreading.'],
  ['-h', '--help', 'Show this help'],
  ['-v', '--version', 'Show version number']
]
demo = [
  'genome.fa',
  'genome.fa -l 10 -o bw,bed,tdf -t 8',
  'genome.fa --length 15 --paired 50,10 --lowmem hard --output all']

app = App(init, args, demo)

# --------------------------------------------------------------------------- #
app.default('input', sys.argv[1])
if not os.path.isfile(app.argx['input']):
    app.exit('Input file not found: ' + app.argx['input'])

app.default('output', 'wig')
basename = os.path.splitext(app.argx['input'])[0]

outputs = {}
exts = app.argx['output'].lower().replace('all', 'wig,bed,tdf,bw,bigbed')
for ext in exts.split(','):
    ext = ext.replace('bigwig', 'bw')
    if ext not in ['wig', 'bed', 'tdf', 'bw', 'bigbed']:
        continue
    outputs[ext] = basename + '.' + ext

if 'wig' in outputs:
    try:
        import pyBigWig as bw
    except ImportError:
        stop('Python module pyBigWig not found')

app.default('paired', 'S')

try:
    # Single-end reads
    if app.argx['paired'][0] == 'S':
        mdist, kernel = [0, [1]]

    # Normal distribution of insertions size
    else:
        from math import sqrt, pi, exp
        mu, s = map(int, app.argx['paired'].split(','))
        def k(i):
            return exp(-(float(i) - mu)**2 / (2 * s**2)) / (sqrt(2*pi)*s)
        mdist = mu - 3 * s
        ker = [k(i) for i in range(mdist, mu + 3 * s + 1)]
        mdist, kernel = [mdist if mdist > 0 else 0, ker]
    
except:
    app.exit('Unable to parse reads type parameters [--reads]')

app.default('threads', os.sysconf('SC_NPROCESSORS_ONLN'))
app.default('length', 100)
app.argx['length'] = int(app.argx['length'])

app.default('lowmem', 'none')
app.argx['quality'] = {'hard': 3, 'soft': 5}.get(app.argx['lowmem'], 0)
if app.argx['quality'] != 0:
    app.argx['threads'] = 0

# --------------------------------------------------------------------------- #

if '--debug' in sys.argv:
    app._debug = sys.argv[sys.argv.index('--debug') + 1]
    print(app._debug)

# --------------------------------------------------------------------------- #
def download(src, dst, iexec=False):
    import urllib.request

    app.log('Downloading: ' + src)
    os.system("curl -s {0} > {1}".format(src, dst))

    if os.path.isfile(dst):
        app.success_log('Download complete: ' + dst)
        if iexec:
            import stat
            os.chmod(dst, os.stat(dst).st_mode | stat.S_IEXEC)
        return True

    app.error_log('Download error: ' + dst)
    return False

# --------------------------------------------------------------------------- #
def check_exe(root):
    root = os.path.dirname(os.path.abspath(root))
    
    exec_dir = root + '/exe'
    if not os.path.exists(exec_dir):
        os.makedirs(exec_dir)

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

    repo = "https://raw.githubusercontent.com/evgeny-bakin/GeMaTrIA/"
    igvtools = "{root}/exe/igvtools.jar".format(root=root)
    src = repo + "master/exe/igvtools.jar"

    if 'tdf' in outputs:
        if os.path.isfile(igvtools):
            app.success_log('Executable `igvtools` found')
        else:
            app.error_log('Executable `igvtools` not found')
            if not download(src, igvtools):
                del(outputs['tdf'])

        java_exists = os.popen("java 2>&1")
        if java_exists.read().find('Usage') == -1:
            app.error_log("".join([
              'Install Java JDK if you want to export results as tdf file\n'
              '    https://www.oracle.com/technetwork/java/javase/downloads/index.html'
            ]))
            del(outputs['tdf'])
        java_exists.close()

    return [igvtools, bed2bigbed]
# --------------------------------------------------------------------------- #
# filename: include/write.py


class Write():
    wig = "fixedStep chrom={0} start={1} step=1 span={2}\n{3}\n"
    bed = "{0}\t{1}\t{2}\t.\t{3}\t.\t{1}\t{2}\t{4}\n"
    stp = ['', -1]

    def __init__(self, name, ext):
        self.ext = ext
        init = getattr(self, 'init_' + self.ext)
        init(name)

    # ----------------------------------------------------------------------- #
    def init_wig(self, f):
        self.h = open(f, 'w')

    def init_bed(self, f):
        self.h = open(f, 'w')

    def init_bw(self, f):
        import pyBigWig as bw
        self.h = bw.open(f, 'w')

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
    def add(self, chr, gms, zer):
        prev = 0
        for i in zer:
            add = getattr(self, '_' + self.ext)
            add(chr, int(prev + 1), int(i - prev + 1), float(gms[prev + 1]))
            prev = i + 1

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

    if app.argx['paired'] == 'S':
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

    os.remove(fs['bed'].h.name + '.sorted')

    app.success_log('Done: {0:.2f}sec.'.format(time.time()-begin))

if size_:
    os.remove(size_)

app.echo('\nGematria has finished. Results:\n', 'white_bold')
for f in outputs:
    if outputs[f][0:3] == '.__':
        os.remove(outputs[f])
    else:
        app.echo(outputs[f] + '\n', 'white')

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
