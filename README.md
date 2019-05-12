# Gematria
**Ge**nome **Ma**ppability **Tr**ack **I**nstant **A**nalysis

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
pip3 install --user git+https://github.com/latur/MakeGMS numpy pyBigWig
curl https://raw.githubusercontent.com/latur/Gematria/master/gematria.standalone.py > gematria.py
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
-q, --quality  Memory usage: quality * genome size. Default: 12
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

Run with default parameters

```bash
./gematria.py ./support/example.fa
```

Read size determination `-l 20`

```bash
./gematria.py ./support/example.fa -l 20
```

Define prefix for result files `-o prefix-string`

```bash
./gematria.py ./support/example.fa -o results-filename
```

By default, the result is saved to the .wig file. 
You can specify additional export formats `-f all`

```bash
./gematria.py ./support/example.fa -f bigbed
./gematria.py ./support/example.fa -f bigwig,bed,tdf
./gematria.py ./support/example.fa -f all
```

The default mappability calculation is for single-end reads `-r S`. 

For paired reads, you can specify the insert size. 
For example, if the insert size is normally distributed 
with an average = 30 and sigma = 10, use the following parameter: `-r N:30:10`

```bash
./gematria.py ./support/example.fa -r N:30:10
```

For uniform distribution use the parameter: `-r U:from:to`

```bash
./gematria.py ./support/example.fa -r U:10:20
```

Inside the program there are [two ways](https://github.com/latur/MakeGMS#algorithm-description) to calculate mappability track. 
You can determine the number of threads `-t 8` and the quality of filtering `-q 4`.
The method will be automatically selected depending on these parameters.

Do not specify the number of threads more than you have free processors. 
Be careful when specifying the quality of the filter. The amount of memory allocated to an application depends on this parameter.

```bash
./gematria.py ./support/example.fa -t 12 -q 3
```

If you want to use exactly the [approach with Filters](https://github.com/latur/MakeGMS#bloom-based-method), specify `-t 0`

```bash
./gematria.py ./support/example.fa -t 0 -q 4
```

If you want to use the [QuickSort-approach](https://github.com/latur/MakeGMS#qsort-based-method), specify `-q 0`

```bash
./gematria.py ./support/example.fa -t 24 -q 0
```


## Installation without a superuser

if you are working on a shared server, and you do not have root access, you 
can locally install the python and the required extensions as follows:

### Python3, pip, Packages

```bash
git clone https://github.com/python/cpython.git
cd cpython && git checkout origin/3.7

./configure --with-pydebug && make -j
wget https://bootstrap.pypa.io/get-pip.py

./cpython/python get-pip.py --user
./cpython/python .local/bin/pip3.7 install --user git+https://github.com/latur/MakeGMS numpy pyBigWig
```

### Gematria

```bash
curl https://raw.githubusercontent.com/latur/Gematria/master/gematria.standalone.py > gematria.py
./cpython/python gematria.py
```

## Development

```bash
git clone https://github.com/latur/Gematria.git ./gematria
cd gematria
```

### Quick tests

Run gematria with different parameters:

```bash
python3 ./support/gematria.test.py
```

### Gematria.py

The script uses the contents of the folder `/include` and auxiliary applications from the `/exe` folder.
If you find it convenient to compile all the code into one script (`gematria.standalone.py`) and use it, run the command:

```bash
python3 build.py
chmod +x gematria.standalone.py
```
