![Gematria Logo](https://raw.githubusercontent.com/latur/Gematria/master/assets/logo-center.jpg)

-----------------------------------------------------------

# Genome Mappability Track Instant Analysis

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
pip3 install --user git+https://github.com/bioinf/Gematria numpy pyBigWig
curl https://raw.githubusercontent.com/bioinf/Gematria/master/gematria.standalone.py > gematria.py
chmod +x gematria.py
```

### Usage

```bash
./gematria.py [fasta file] [Optional arguments]
```

### Optional arguments

```
-l, --length   Read length. Default: 100
-t, --threads  Number of threads. Default: auto
-o, --output   Output formats: [wig],bigwig,bed,tdf,bigbed,all
-p, --paired   To use paired end reads, specify the insert size
               and the standard deviation (normal distribution model)
               If you don’t know where to start, use -p 300,100
-m, --lowmem   Use the RAM-memory saving algorithm. Default: none
                 -m hard / will be used ~ 5 × Genome Size
                 -m soft / will be used ~ 7 × Genome Size
               Specifying this parameter in the value hard or soft 
               makes it impossible to use multithreading.
-h, --help     Show this help
-v, --version  Show version number
```

## Examples:

Run with default parameters

```bash
./gematria.py ./test/example.fa
```

Read size determination `-l 20`

```bash
./gematria.py ./test/example.fa -l 20
```

By default, the result is saved to the .wig file. 
You can specify additional export formats `-o all`

```bash
./gematria.py ./test/example.fa -o bigbed
./gematria.py ./test/example.fa -o bigwig,bed,tdf
./gematria.py ./test/example.fa -o all
```

The default mappability calculation is for single-end reads `-p S`. 

For paired reads, you can specify the insert size. 
For example, if the insert size is normally distributed 
with an average = 30 and sigma = 10, use the following parameter: `-p 30,10`

```bash
./gematria.py ./test/example.fa -p 30,10
```

Inside the program there are [two ways](https://github.com/latur/MakeGMS#algorithm-description) to calculate mappability track.
If you want to use exactly the [approach with Filters](https://github.com/latur/MakeGMS#bloom-based-method), specify `-m hard` or `-m hard`
This approach saves memory but does not have a multithreaded run.

```bash
./gematria.py ./test/example.fa -m hard
```

If you want to use the [QuickSort-approach](https://github.com/latur/MakeGMS#qsort-based-method), specify number of threads `-t 16`

```bash
./gematria.py ./test/example.fa -t 16
```


## Installation without a superuser

if you are working on a shared server, and you do not have root access, you 
can locally install the python and the required extensions as follows:

### Python3, pip, Packages

```bash
git clone https://github.com/python/cpython.git
cd cpython && git checkout origin/3.7

./configure --with-pydebug && make -j && cd ../
wget https://bootstrap.pypa.io/get-pip.py

./cpython/python get-pip.py --user
./cpython/python .local/bin/pip3.7 \
  install --user git+https://github.com/bioinf/Gematria numpy pyBigWig
```

### Gematria

```bash
curl https://raw.githubusercontent.com/bioinf/Gematria/master/gematria.standalone.py > gematria.py
./cpython/python gematria.py
```

## Development

```bash
git clone https://github.com/bioinf/Gematria.git ./gematria
cd gematria
```

### Quick tests

Run gematria with different parameters:

```bash
python3 ./test/gematria.test.py
```

### Gematria.py

The script uses the contents of the folder `/include` and auxiliary applications from the `/exe` folder.
If you find it convenient to compile all the code into one script (`gematria.standalone.py`) and use it, run the command:

```bash
python3 build.py && chmod +x gematria.standalone.py
```
