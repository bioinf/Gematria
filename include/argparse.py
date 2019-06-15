import sys
import os

from include.app import *

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
