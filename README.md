IsoCon
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

Type in terminal:
```
pip install IsoCon
```

With proper installation of **IsoCon**, you should be able to run:

```
    IsoCon
```

to view user instructions. This installation will install the dependencies automatically for you. Customized installation of latest master branch is also easy, see below.

#### Downloading source from GitHub

##### Dependencies

* scipy
* [edlib](https://github.com/Martinsos/edlib "edlib's Homepage"), for installation see [link](https://github.com/Martinsos/edlib/tree/master/bindings/python#installation)
* [ssw](https://github.com/vishnubob/ssw "Python wrapper for SSW"), for installation see [link](https://github.com/vishnubob/ssw#installation)
* [networkx](https://networkx.github.io/)

These libraries can all be installed with `pip install`. With the dependencies already installed.

```sh
git clone https://github.com/ksahlin/IsoCon.git
cd IsoCon
./IsoCon
```

Usage
-------

Running IsoCon end-to-end (preffered)
```
    IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```


Example datasets
-----------------

We ran IsoCon on the publicly availabe datasets from pacbio [MCF7](link) and [Alzheimer](link) using the following command

```
    IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --minimizer_search_depth 100
```

| Dataset | runtime  | peak memory | final_candidates | corr | not_corr | *TOFU* | *nr original CCS* | 
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| [MCF7](http://www.pacb.com/blog/data-release-human-mcf-7-transcriptome/) | 24h43m  | <1.9Gb  | 2169 | 18458 | 401885 | 518701 | here |
|[Alzheimer](http://www.pacb.com/blog/data-release-alzheimer-brain-isoform-sequencing-iso-seq-dataset/)| time | Content Cell  | Content Cell  |

IsoCon predicted 2169 and Y transcripts and had a runtime of 24h43m and Y, for MCF7 and Alzheimer datasets respectively, on a 64 core machine with a peak memory usage of Z Gb. IsoCons output are found [here](link).

Manual BLAT of 20 sequences from the "corrected_but_not_converged" predictions to human show an alignment identity increase from 94-98% of the ccs reads up to 99.5-99.9% for the corrected reads.


Detailed usage
----------------

IsoCOn has three different modes `pipeline` (end to end), `get_candidates` only correct and cluster reads (no statistical testing), as well as `stat_filter` (only statistical testing of a set of candidate transcripts).

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

#### Commands

```
    $ ./run3CO pipeline --help
usage: Pipeline for obtaining non-redundant haplotype specific transcript isoforms using PacBio IsoSeq reads. pipeline
       [-h] -fl_reads FL_READS -outfolder OUTFOLDER [--barcodes BARCODES]
       [--develop_mode] [--minimizer_search_depth MINIMIZER_SEARCH_DEPTH]
       [--single_core] [--min_candidate_support MIN_CANDIDATE_SUPPORT]
       [--ignore_ends_len IGNORE_ENDS_LEN] [--cleanup]
       [--prefilter_candidates]

optional arguments:
  -h, --help            show this help message and exit
  --barcodes BARCODES   If targeted approach. A fasta file with barcodes to
                        search for in ends of transcripts.
  --develop_mode        This will print more information abount workflow and
                        provide plots of similarity network etc.
  --minimizer_search_depth MINIMIZER_SEARCH_DEPTH
                        Maximum number of pairwise alignments in search matrix
                        to find minimizer. [default =2**32]
  --single_core         Force working on single core.
  --min_candidate_support MIN_CANDIDATE_SUPPORT
                        Required minimum number of reads converged to the same
                        sequence to be included in statistical test.
  --ignore_ends_len IGNORE_ENDS_LEN
                        Number of bp to ignore in ends. If two candidates are
                        identical expept in ends of this size, they are
                        collapses and the longest common substing is chosen to
                        represent them. In statistical test step, minimizers
                        are found based on ignoring the ends of this size.
                        Also indels in ends will not be tested. [default
                        ignore_ends_len=15].
  --cleanup             Remove everything except logfile.txt,
                        candidates_converged.fa and final_candidates.fa in
                        output folder.
  --prefilter_candidates
                        Filter candidates if they are not consensus over at
                        least one base pair, this can reduce runtime without
                        significant loss in true candidates.
```


