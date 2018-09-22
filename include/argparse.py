import sys
import os

from include.app import *

init = '[fasta file] [Optional arguments]'
args = [
  ['-i', '--input', 'Path to `genome.fasta` file'],
  ['-l', '--length', 'Read length. Default: 100'],
  ['-t', '--threads', 'Number of threads. Default: auto'],
  ['-o', '--output', 'Output filenames without extension'],
  ['-f', '--formats', 'Comma separated output formats',
                      'Acceptable: wig, bigwig, bed, tdf, bigbed, all'],
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
exts = app.argx['formats'].lower().replace('all', 'wig,bed,tdf,bw,bigbed')
for ext in exts.split(','):
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
    os.system("curl -s {0} > {1}".format(src, dst))

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

    return [igvtools, bed2bigbed]
