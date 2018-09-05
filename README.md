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

### Installation

```bash
pip install makegms numpy pyBigWig
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
-t, --threads  Number of threads. Default: autodetect
-o, --output   Output filename with an extension wig, bw or bed
-r, --reads    Reads type parameters in the following format:
               S - for single-end reads
               N:mu:sigma - for Normal distribution of insertion size
               U:min:max - for Uniform distribution of insertion size
-h, --help     Show this help
```

### Examples:

```bash
# Run of Gematria for 100bp single-end reads with a result saving in a wig format:
./gematria.py input.fasta

# .. Result saving in a bigWig file:
./gematria.py -i test/ecoli.fasta -o result.bw

# Calculate mappability for 35bp paired-end reads with 200..300bp insertion size:
./gematria.py input.fasta -l 35 -r U:200:300

```

## Installation without a superuser

If you are on a shared server, and you do not have root access, you can 
locally install the python and the required extensions as follows:

```bash
# Python
git clone https://github.com/python/cpython.git
cd cpython && ./configure --with-pydebug && make -j

# Pip
wget https://bootstrap.pypa.io/get-pip.py
./cpython/python get-pip.py --user

# Packages
./cpython/python .local/bin/pip3.8 install makegms numpy pyBigWig --user

# Gematria
curl https://raw.githubusercontent.com/evgeny-bakin/GeMaTrIA/master/gematria.standalone.py > gematria.py
./cpython/python gematria.py
```

## Requirements

Gematria interface part was developed via Python 3 (scripts gematria.py for developing and gematria.standalone.py for production) while the most computationally-consuming algorithms were implemented in C (`makegms`).
The following set of extra Python libraries are required (in brackets we give the versions of the libraries for which our tool was tested): `numpy`

---

[edit]

---

Out of the box it can produce output files with a mappability track in two formats: Wig and Bed (in the latter case a mappability will be represented with a color intensity). However a user may enhance capabilites of GeMaTrIA output generation in the following way:

* For bigWig: install Python [pyBigWig library](https://github.com/deeptools/pyBigWig)
* For bigBed: 
	- create directory _utils_ in GeMaTrIA main directory;
	- download there [bigBedToBed](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed) and [faSize](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize) utilities from UCSC official web site.
* For TDF: 
	- create directory _utils_ in GeMaTrIA main directory;
	- download there [faSize](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize) from UCSC official web site;
	- download there [igvtools](http://software.broadinstitute.org/software/igv/download) from Broadinstitute web site.

Thus, for the full set of options, your GeMaTrIA directory should look as follows:

GeMaTrIA<br/>
|<br/>
+--&nbsp; utils/<br/>
|&nbsp; &nbsp; &nbsp; |<br/>
|&nbsp; &nbsp; &nbsp; +-- bedToBigBed<br/>
|&nbsp; &nbsp; &nbsp; +-- faSize<br/>
|&nbsp; &nbsp; &nbsp; +-- igvtools<br/>
|&nbsp; &nbsp; &nbsp; +-- igvtools.jar<br/>
+--&nbsp; GMS\_aux\_lib.py<br/>
+--&nbsp; makeRawGMS<br/>
+--&nbsp; gematria.py<br/> 

## Installation

Just clone this GitHub project and download extra utilities if required (see section **Requirements**). Check, that file makeRawGMS has permissions to be executed (`chmod u+x makeRawGMS`).

## Tool options 

```shell
python3 gematria.py [-h] -i FILE -o Str [-l Int] [-r Str] [-f Str] [-t Int] [-m]                   
```

Arguments:
*  `-i, --input`    Name of FASTA-file with a genome.
*  `-o, --output`   Prefix of output files generated with GeMaTrIA.
*  `-l, --length`   Read length (default is 100).
*  `-r, --read`     Reads type parameters in the following format: S (for single-end reads) **or** N:mu:sigma (for Normal distribution of insertions size) **or** U:min:max (for Uniform distribution of insertion size). Default is a single-end mode.
*  `-f, --formats`  List of required output formats separated with commas (**without spaces**), e.g. format1,format2,... etc. Supported formats are Wig, bigWig, Bed, bigBed, TDF, ALL (for generation of track in all available formats). Default is bigWig.
*  `-m, --mat`      Save debug data in Matlab MAT format.
*  `-t, --threads`  Number of threads (**please, note, that multithreading for large genomes is still under debugging**).

## Usage examples

Easiest run of GeMaTrIA for 100bp single-end reads with a result saving in a bigWig file:
```shell
python3 gematria.py -i e_coli.fasta -o e_coli
```

Save data in Wig and TDF formats:
```shell
python3 gematria.py -i e_coli.fasta -o e_coli -f Wig,TDF
```

Calculate mappability for 100bp paired-end reads with 200..300bp insertion size:
```shell
python3 gematria.py -i e_coli.fasta -o e_coli -f Wig,TDF -r U:200:300
```
Calculate mappability for 200bp paired-end reads with 400..600bp insertion size:
```shell
python3 gematria.py -i e_coli.fasta -o e_coli -f Wig,TDF -l 200 -r U:400:600
```

Do all this, using 4 threads:
```shell
python3 gematria.py -i e_coli.fasta -o e_coli -f Wig,TDF -l 200 -r U:400:600 -t 4
```

## Our team

This tool was developed by joint efforts of Evgeny Bakin<sup>1</sup>, Natalia Zorina<sup>2</sup>, Elizaveta Skalon<sup>2</sup>, Michael Utman<sup>3</sup> and Vyacheslav Batunov<sup>4</sup> under general supervision of Alexander Predeus<sup>1</sup>. 

[1] Bioinformatics institute, Saint-Petersburg, Russia

[2] Saint-Petersburg State University (SPbU), Russia

[3] Saint-Petersburg Academic University (SPbAU), Russia

[4] Saint-Petersburg State University of Aerospace Instrumentation (SUAI), Russia

## TBD 
* Improve speed of output files generation.
* Debug multithreading for large genomes.
* Make code more readable.
* Add more reallistic insertion size distributions.
* Add erroneous reads support.
* Add automatic generation of Circos plot.
