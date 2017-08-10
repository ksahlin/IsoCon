IsoCON
========

Supported on Linux / OSX with python 2.7, and 3 [![Build Status](https://travis-ci.org/ksahlin/BESST.svg?branch=master)](https://travis-ci.org/ksahlin/IsoCon)

Codecov:        https://img.shields.io/codecov/c/github/codecov/example-python.svg
Codecov branch:     https://img.shields.io/codecov/c/github/codecov/example-python/master.svg

IsoCon is a tool for deriving *finished transcript sequences* from *Iso-Seq* reads. Input is a set of full-length-non-chimeric reads and output is a set of predicted transcripts. IsoCon can be run as follows

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

predicted transcripts is found in file */path/to/output/final_candidates.fa*. Reads that could not be corrected or clustered (e.g., due to chimeric products or less than three reads sequenced from the transcript) are found in */path/to/output/not_converged.fa*. For more instructions see below.

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† (2017) "Title." Journal (bioRxiv for now?)(17), 2215-2222 [Link](here)

INSTALLATION
----------------

The preferred way to install IsoCon is with pythons package installer pip.

#### Using pip (recommended)
---------------------------------------------------

Type in terminal:
```
sudo pip install IsoCon
```

If you do not have *pip* for python, you can get *pip* as follows:
```
curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
python get-pip.py
``` 
With proper installation of **IsoCon**, you should be able to run:

```
IsoCon
```

to view user instructions.

This installation will install the dependencies automatically for you. Customized installation of latest master branch is also easy, see below.


#### Dependencies

* scipy
* [edlib](https://github.com/Martinsos/edlib "edlib's Homepage"), for installation see [link](https://github.com/Martinsos/edlib/tree/master/bindings/python#installation)
* [ssw](https://github.com/vishnubob/ssw "Python wrapper for SSW"), for installation see [link](https://github.com/vishnubob/ssw#installation)
* [networkx](https://networkx.github.io/)

These libraries can all be installed with `pip install`.


#### Source from GitHub

This assumes that you have the dependencies already installed.

```sh
git clone https://github.com/ksahlin/IsoCon.git
cd IsoCon
./IsoCon
```

Usage
-------

Running IsoCon
```
    IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

```
$ IsoCon --help
usage: Pipeline for obtaining non-redundant haplotype specific transcript isoforms using PacBio IsoSeq reads.
       [-h] [--version]
       {pipeline,remove_barcodes,get_candidates,stat_filter} ...

positional arguments:
  {pipeline,remove_barcodes,get_candidates,stat_filter}
                        help for subcommand
    pipeline            Run the entire pipeline, i.e., runs: 1.
                        remove_barcodes() if --barcodes is specified 2.
                        get_candidates() to find candidate transcripts out of
                        CCS full length reads 3. stat_filter() to filter
                        candidates and only output significant ones
    remove_barcodes     Remove barcodes from sequences.
    get_candidates      Get candidate transcripts with minimizer approach.
    stat_filter         Get candidate transcripts with minimizer approach.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

Example datasets
-----------------

We ran IsoCon on the publicly availabe datasets from pacbio [MCF7](link) and [Alzheimer](link) using the following command

```
    IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --minimizer_search_depth 100
```

IsoCon predicted X and Y transcripts and had a runtime of X and Y, for MCF7 and Alzheimer datasets respectively, on a 64 core machine with a peak memory usage of Z Gb. The predicted transcripts as well as the reads that could not be clustered or corrected are found [here](link).



