# Gematria
> Genome Mappability Track Instant Analysis

## What is mappability?

Mappability is a genome-wide function showing if a region in a genome could be 
identified unambiguously using a read of a particular length. This function 
typically ranges between 0 and 100. If exact genome subsequence may be found 
in more than one location, then the mappability in that location of the genome 
is set to zero. Otherwise, if the subsequence is unique, mappabity is close to 
100.

Obviously, mappability depends on read length (the more read length the more 
unambiguous locus is), and insertion size for paired-end reads (growing of 
insertion size improves mappability as well).

## Quick Start

### Installation (for Python3)

```bash
pip3 install --user git+https://github.com/latur/Makegms
pip3 install numpy pyBigWig
curl https://raw.githubusercontent.com/evgeny-bakin/GeMaTrIA/master/gematria.standalone.py > gematria.py
chmod +x gematria.py
```

### Usage

```bash
./gematria.py [fasta file] [Optional arguments]
```

### Optional arguments

```bash
-i, --input    Path to `genome.fasta` file
-l, --length   Read length. Default: 100
-t, --threads  Number of threads. Default: auto
-o, --output   Output filenames without extension
-f, --formats  Comma separated output formats
               Acceptable: wig, bigwig, bed, tdf, bigbed, all
-r, --reads    Reads type parameters in the following format:
               S - for single-end reads
               N:mu:sigma - for Normal distribution of insertion size
               U:min:max - for Uniform distribution of insertion size
-h, --help     Show this help
```

### Examples:

```bash
# Run of Gematria for 5bp single-end reads 
# with a result saving in a bigWig, bed, tdf files:
./gematria.py ./test/example.fa -l 5 -o result -f bw,bed,tdf

# Calculate mappability for 7bp paired-end reads with 10..25bp insertion size:
# the result will be saved to file ./test/example.fa.wig
./gematria.py -i ./test/example.fa -l 7 -r U:10:25
```

## Installation without a superuser

if you are working on a shared server, and you do not have root access, you 
can locally install the python and the required extensions as follows:

```bash
# Python
git clone https://github.com/python/cpython.git
cd cpython && ./configure --with-pydebug && make -j

# Pip
wget https://bootstrap.pypa.io/get-pip.py
./cpython/python get-pip.py --user

# Packages
./cpython/python .local/bin/pip3.8 install numpy pyBigWig git+https://github.com/latur/Makegms --user

# Gematria
curl https://raw.githubusercontent.com/evgeny-bakin/GeMaTrIA/master/gematria.standalone.py > gematria.py
./cpython/python gematria.py
```

## Development

```bash
git clone https://github.com/evgeny-bakin/GeMaTrIA.git ./gematria
cd gematria
```

### Makegms

The `makegms` module generates a raw binary track where the `false` correspond to repeating regions, and the `true` correspond to unique regions in the genome for a given length of the read.

Build as a python module:

```bash
cd makegms/
rm -rf makegms.egg-* dist build
python3 setup.py build
python3 setup.py install

# Testing:
cd ../
python3 ./test/makegms.test.py
```

Example of module usage

```python
import makegms
track = makegms.run('input.file.fasta', read=100, threads=4)
```

Build as an independent application:

```bash
cd makegms/src/
gcc -pthread -std=c99 -m64 -O main.c -lm \
  -o ../../exe/makegms.exe

# Memory leak testing:
cd ../../
valgrind --leak-check=full --show-leak-kinds=all \
 ./exe/makegms.exe ./test/example.fa 10 1
```

### Gematria.py

Python wrapper for makegms. The script uses the contents of the folder `/include` and auxiliary applications from the `/exe` folder.

```bash
# Running the script on test data with various parameters
python3 ./test/gematria.test.py
```

### Gematria.standalone.py

If you find it convenient to compile all the code into one script (`gematria.standalone.py`) and use it, run the command:

```bash
python3 build.py
chmod +x gematria.standalone.py
```

## Our team

This tool was developed by joint efforts of Evgeny Bakin<sup>1</sup>, Natalia Zorina<sup>2</sup>, Elizaveta Skalon<sup>2</sup>, Michael Utman<sup>3</sup> and Vyacheslav Batunov<sup>4</sup> under general supervision of Alexander Predeus<sup>1</sup>. 

[1] Bioinformatics institute, Saint-Petersburg, Russia  
[2] Saint-Petersburg State University (SPbU), Russia  
[3] Saint-Petersburg Academic University (SPbAU), Russia  
[4] Saint-Petersburg State University of Aerospace Instrumentation (SUAI), Russia

## TBD
 
 [-] Improve speed of output files generation.  
 [-] Debug multithreading for large genomes.  
 [+] Make code more readable.  
 [-] Add more reallistic insertion size distributions.  
 [-] Add erroneous reads support.  
 [-] Add automatic generation of Circos plot.
