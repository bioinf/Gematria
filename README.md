# GeMaTrIA
This project is devoted to development of fast and easy tool for Genome Mappability Track Instant Analysis (GeMaTrIA).

## What is mappability?

Mappability is a genome-wide function showing if a region in a genome could be identified unambiguously using a read of a particular length. This function typically ranges between 0 and 100. If exact genome subsequence may be found in more than one location, then the mappability in that location of the genome is set to zero. Otherwise, if the subsequence is unique, mappabity is close to 100.

Obviously, mappability depends on read length (the more read length the more unambiguous locus is), and insertion size for paired-end reads (growing of insertion size improves mappability as well).

## Requirements

GeMaTrIA interface part was developed via Python 3 (scripts _gematria.py_ and _GMS_aux_lib.py_) while the most computationally-consuming algorithms were implemented in C (_makeRawGMS_).

The following set of extra Python libraries are required (in brackets we give the versions of the libraries for which our tool was tested):
* numpy (1.10.4)
* biopython (1.68)
* scipy (0.17.1)

Out of the box it can produce output files with a mappability track in two formats: Wig and Bed (in this case mappability will be reflected with color intensity). However a user may enhance capabilites of GeMaTrIA output formats in the following way:

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

This 

## Tool options 

TBD

## Usage examples

TBD

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
