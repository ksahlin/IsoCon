IsoCon
========

Supported on Linux / OSX with python 2.7, and 3 [![Build Status](https://travis-ci.org/ksahlin/BESST.svg?branch=master)](https://travis-ci.org/ksahlin/IsoCon)

Codecov:        https://img.shields.io/codecov/c/github/codecov/example-python.svg
Codecov branch:     https://img.shields.io/codecov/c/github/codecov/example-python/master.svg

IsoCon is a tool for deriving *finished transcript sequences* from *Iso-Seq* reads. Input is a set of full-length-non-chimeric reads and output is a set of predicted transcripts. IsoCon can be run as follows

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```

predicted transcripts is found in file */path/to/output/final_candidates.fa*. Reads that could not be corrected or clustered (e.g., due to chimeric products or less than three reads sequenced from the transcript) are found in */path/to/output/not_converged.fa*. For more instructions see below.

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† (2017) "Title." Journal (bioRxiv for now?)(17), 2215-2222 [Link](here)

<!-- Table of Contents
=================
    * [INSTALLATION](#installation)
        * [Using pip (recommended)](#using-pip-recommended)
        * [Downloading source from GitHub](#downloading-source-from-github)
            * [Dependencies](#dependencies)
    * [Usage](#usage)
    * [IsoCon on general Iso-Seq datasets](#isocon-on-general-iso-seq-datasets)
    * [Detailed usage](#detailed-usage)
        * [Commands](#commands)

TOC created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc) -->

INSTALLATION
----------------

The preferred way to install IsoCon is with pythons package installer pip.

#### Using pip (recommended)

Type in terminal `pip install IsoCon` . With proper installation of **IsoCon**, you should be able to issue the command ` IsoCon` to view user instructions. This installation will install the dependencies automatically for you. Customized installation of latest master branch is also easy, see below.

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

Running IsoCon end-to-end with base call quality values (preffered)
```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```

IsoCon takes two input files: (1) a fasta file of full length non-chimeric (flnc) CCS reads and (2) the bam file of CCS reads containing predicted base call quality values. IsoCon also supports taking only the fasta read file as input. Without the base call quality values, IsoCon will use an empirically estimated error model. The ability to take only the flnc fasta file as input is useful when the reads have been altered after the CCS base calling algorithm, \emph{e.g.}, from correction using illumina reads. We highly recommend supplying the quality values to IsoCon if CCS reads has not gone through any additional correction. IsoCon is available at \url{[https://github.com/ksahlin/IsoCon]} and distributed as a python package for Unix and OSX based systems.

Running IsoCon end-to-end without base call quality values (e.g., if CCS reads have been altered after they were constructed such as Illumina corrected)
```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

Output
-------

The final high quality transcripts are written to the file `final_candidates.fa` in the output folder. If there was only one read coming from a transcript, which is sufficiently different from other reads (exon difference), it will be output in the file `not_converged.fa`. This file may also contain other erroneous CCS reads such as chimeras.

<!-- IsoCon on general Iso-Seq datasets
-----------------------------------

IsoCon has been evaluated on targeted Iso-Seq data [cite], but the algorithm makes no explicit assumptions, or depends on the data being targeted. Therefore it also runs on any general Iso-Seq data. To demonstrate this We ran IsoCon on the publicly availabe datasets from pacbio [MCF7](link) and [Alzheimer](link) using the following command

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --neighbor_search_depth 100
```

| Dataset | runtime  | peak memory | final_candidates | corr | not_corr | *TOFU* | *nr original CCS* | 
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| [MCF7](http://www.pacb.com/blog/data-release-human-mcf-7-transcriptome/) | 24h43m  | <1.9Gb  | 2,169 | 18,458 | 401,885<sup>[1](#myfootnote1)</sup> |  55,770 | 518,701 |
|[Alzheimer](http://www.pacb.com/blog/data-release-alzheimer-brain-isoform-sequencing-iso-seq-dataset/)| time | Content Cell  | Content Cell  |

TODO1: Output three separate files: validated_transcripts.fa,  corrected_not_validated.fa, not_converged.fa 
TODO2: Show set intersection (allowing end differences between final_cand + corr to TOFU transcripts). Also do other analysis.

IsoCon predicted 2169 and Y transcripts and had a runtime of 24h43m and Y, for MCF7 and Alzheimer datasets respectively, on a 64 core machine with a peak memory usage of Z Gb. IsoCons output are found [here](link).

Manual BLAT of 20 sequences from the "corrected_but_not_converged" predictions to human show an alignment identity increase from 94-98% of the ccs reads up to 99.5-99.9% for the corrected reads.

<a name="footnote_not_converged">1</a>: Footnote content goes here -->

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
    get_candidates      Get candidate transcripts with nearest_neighbor approach.
    stat_filter         Get candidate transcripts with nearest_neighbor approach.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

#### Commands

```
    $ ./IsoCon pipeline --help
usage: Pipeline for obtaining non-redundant haplotype specific transcript isoforms using PacBio IsoSeq reads. pipeline
       [-h] -fl_reads FL_READS -outfolder OUTFOLDER [--barcodes BARCODES]
       [--develop_mode] [--neighbor_search_depth neighbor_search_depth]
       [--single_core] [--min_candidate_support MIN_CANDIDATE_SUPPORT]
       [--ignore_ends_len IGNORE_ENDS_LEN] [--cleanup]
       [--prefilter_candidates]

optional arguments:
  -h, --help            show this help message and exit
  --barcodes BARCODES   If targeted approach. A fasta file with barcodes to
                        search for in ends of transcripts.
  --develop_mode        This will print more information abount workflow and
                        provide plots of similarity network etc.
  --neighbor_search_depth neighbor_search_depth
                        Maximum number of pairwise alignments in search matrix
                        to find nearest_neighbor. [default =2**32]
  --single_core         Force working on single core.
  --min_candidate_support MIN_CANDIDATE_SUPPORT
                        Required minimum number of reads converged to the same
                        sequence to be included in statistical test.
  --ignore_ends_len IGNORE_ENDS_LEN
                        Number of bp to ignore in ends. If two candidates are
                        identical expept in ends of this size, they are
                        collapses and the longest common substing is chosen to
                        represent them. In statistical test step, nearest_neighbors
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


