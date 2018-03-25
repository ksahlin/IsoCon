IsoCon
========

IsoCon is distributed as a python package supported on Linux / OSX with python v>=2.7, and 3.4-3.6, 3.5-dev and 3.6-dev [![Build Status](https://travis-ci.org/ksahlin/IsoCon.svg?branch=master)](https://travis-ci.org/ksahlin/IsoCon)


IsoCon is a tool for deriving *finished transcripts* from *Iso-Seq* reads from *targeted* sequencing. Input is either a set of full-length-non-chimeric reads in fasta format and the CCS base call values as a bam file, or a fastq file of CCS reads and their quality values provided by the ccs caller. IsoCon can only run on targeted Iso-Seq datasets due to underlying alignment approach --- runtime does not scale for nontargeted sequencing. The output is a set of predicted transcripts. IsoCon can be run as follows

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```
or

```
IsoCon pipeline -fl_reads <flnc.fastq> -outfolder </path/to/output>
```

predicted transcripts are found in file **/path/to/output/final_candidates.fa**. Reads that could not be corrected or clustered are found in /path/to/output/not_converged.fa. For more instructions see below.


Table of Contents
=================

  * [Table of Contents](#table-of-contents)
  * [INSTALLATION](#installation)
    * [Using pip](#using-pip)
    * [Downloading source from GitHub](#downloading-source-from-github)
    * [Dependencies](#dependencies)
  * [USAGE](#usage)
    * [Pipline](#pipline)
      * [Output](#output)
    * [get_candidates](#get_candidates)
    * [stat_filter](#stat_filter)
    * [Parameters](#parameters)
  * [CREDITS](#credits)
  * [LICENCE](#licence)


INSTALLATION
----------------

The preferred way to install IsoCon is with pythons package installer pip.

### Using pip 

`pip` is pythons official package installer. This section assumes you have `python` (v2.7 or >=3.4) and a recent version of `pip` installed which should be included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 

With `python` and `pip` available, create a file `requirements.txt` with contents copied from [this file](https://github.com/ksahlin/IsoCon/blob/master/requirements.txt). Then, type in terminal 

```
pip install --requirement requirements.txt IsoCon
```

This should install IsoCon. With proper installation of **IsoCon**, you should be able to issue the command `IsoCon pipeline` to view user instructions. You should also be able to run IsoCon on this [small dataset](https://github.com/ksahlin/IsoCon/tree/master/test/data). Simply download the test dataset and run:

```
IsoCon pipeline -fl_reads [path/simulated_pacbio_reads.fa] -outfolder [output path]
```

`pip` will install the dependencies automatically for you. IsoCon has been built with python 2.7, 3.4-3.6 on Linux systems using [Travis](https://travis-ci.org/). For customized installation of latest master branch, see below.

### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as IsoCon has not been tested with earlier versions of these libraries. However, IsoCon may also work with earliear versions of these libaries.
* [edlib](https://github.com/Martinsos/edlib "edlib's Homepage"), for installation see [link](https://github.com/Martinsos/edlib/tree/master/bindings/python#installation) (>= v1.1.2)
* [networkx](https://networkx.github.io/) (>= v1.10)
* [parasail](https://github.com/jeffdaily/parasail-python)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/IsoCon.git
cd IsoCon
./IsoCon
```


USAGE
-------

IsoCon's algorithm consists of two main phases; the error correction step and the statistical testing step. IsoCon can run these two steps in one go using `IsoCon pipeline`, or it can run only the correction or statistical test steps using `IsoCon get_candidates` and `IsoCon stat_filter` respectively. The preffered and most tested way is to use the entire pipeline `IsoCon pipeline`, but the other two settings can come in handy for specific cases. For example, running only `IsoCon get_candidates` will give more sequences if one is not concerned about precision and will also be faster, while one might use only `IsoCon stat_filter` using different parameters for a set of already constructed candidates in order to prevent rerunning the error correction step.


### Pipeline

IsoCon takes two input files: (1) a fasta file of full length non-chimeric (flnc) CCS reads and (2) the bam file of CCS reads containing predicted base call quality values. The fasta file containing flnc can be obtained from PacBios Iso-Seq pipeline [ToFU](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki) and the bam file is the output of running the consensus caller algorthm [ccs](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md) on the Iso-Seq reads (ccs takes bam files so if you have bax files, convert them using [bax2bam](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md#input) ). IsoCon can then be run as

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> [--ccs </path/to/filename.ccs.bam>]
```

IsoCon also supports taking only the fasta read file as input. Without the base call quality values in `--ccs`, IsoCon will use an empirically estimated error model. The ability to take only the flnc fasta file as input is useful when the reads have been altered after the CCS base calling algorithm, \emph{e.g.}, from error correction using Illumina reads. **However, we highly recommend supplying the CCS quality values to IsoCon if CCS reads has not gone through any additional correction.** 

Simply omit the `--ccs` parameter if running IsoCon without base call quality values as
```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

#### Output

The final high quality transcripts are written to the file `final_candidates.fa` in the output folder. If there was only one or two reads coming from a transcript, which is sufficiently different from other reads (exon difference), it will be output in the file `not_converged.fa`. This file may contain other erroneous CCS reads such as chimeras. The output also contains a file `cluster_info.tsv` that shows for each read which candidate it was assigned to in `final_candidates.fa`.

### get_candidates

Runs only the error correction step. The output is the converged candidates in a fasta file.

```
IsoCon get_candidates -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

### stat_filter

Runs only the statistical filtering of candidates.

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> -candidates <candidate_transcripts.fa>  [--ccs </path/to/filename.ccs.bam>]
```
Observe that `candidate_transcripts.fa` does not have to come from IsoCon's error correction algorithm. For example, this could either be a set of already validated transcripts to which one would like to see if they occur in the CCS reads, or they could be Illumina (or in other ways) corrected CCS reads.


### Parameters

```
    $ IsoCon pipeline --help
usage: Pipeline for obtaining non-redundant haplotype specific transcript isoforms using PacBio IsoSeq reads. pipeline
       [-h] -fl_reads FL_READS -outfolder OUTFOLDER [--ccs CCS]
       [--nr_cores NR_CORES] [--verbose]
       [--neighbor_search_depth NEIGHBOR_SEARCH_DEPTH]
       [--min_exon_diff MIN_EXON_DIFF]
       [--min_candidate_support MIN_CANDIDATE_SUPPORT]
       [--p_value_threshold P_VALUE_THRESHOLD]
       [--min_test_ratio MIN_TEST_RATIO]
       [--max_phred_q_trusted MAX_PHRED_Q_TRUSTED]
       [--ignore_ends_len IGNORE_ENDS_LEN] [--cleanup]
       [--prefilter_candidates]

optional arguments:
  -h, --help            show this help message and exit
  --ccs CCS             BAM/SAM file with CCS sequence predictions.
  --nr_cores NR_CORES   Number of cores to use. [default = 16]
  --verbose             This will print more information abount workflow and
                        provide plots of similarity network etc.
  --neighbor_search_depth NEIGHBOR_SEARCH_DEPTH
                        Maximum number of pairwise alignments in search matrix
                        to find nearest_neighbor. [default =2**32]
  --min_exon_diff MIN_EXON_DIFF
                        Minimum consequtive base pair difference between two
                        neigborss in order to remove edge. If more than this
                        nr of consequtive base pair difference, its likely an
                        exon difference. [default =20]
  --min_candidate_support MIN_CANDIDATE_SUPPORT
                        Required minimum number of reads converged to the same
                        sequence to be included in statistical test. [default
                        2]
  --p_value_threshold P_VALUE_THRESHOLD
                        Threshold for statistical test, filter everythin below
                        this threshold . [default = 0.01]
  --min_test_ratio MIN_TEST_RATIO
                        Don't do tests where candidate c has more than
                        <min_test_ratio> reads assigned to itself compared to
                        the reference t, calculated as test_ratio = c/t,
                        because c will likely be highly significant [default =
                        5]
  --max_phred_q_trusted MAX_PHRED_Q_TRUSTED
                        Maximum PHRED quality score trusted (T), linerarly
                        remaps quality score interval [0,93] --> [0, T].
                        Quality scores may have some uncertainty since T is
                        estimated from a consensus caller algorithm.
  --ignore_ends_len IGNORE_ENDS_LEN
                        Number of bp to ignore in ends. If two candidates are
                        identical except in ends of this size, they are
                        collapsed and the longest common substing is chosen to
                        represent them. In statistical test step, the nearest
                        neighbors are found based on ignoring the ends of this
                        size. Also indels "hanging off" ends of this size will
                        not be tested. [default 15].
  --cleanup             Remove everything except logfile.txt,
                        candidates_converged.fa and final_candidates.fa in
                        output folder. [default = False]
  --prefilter_candidates
                        Filter candidates if they are not consensus over any
                        base pair in the candidate transcript formed from
                        them, this can reduce runtime without significant loss
                        in true candidates. [default = False]

required arguments:
  -fl_reads FL_READS    Fast<a/q> file pacbio Reads of Insert.
  -outfolder OUTFOLDER  Outfolder.
```

CREDITS
----------------

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† (2018) "IsoCon: Deciphering highly similar multigene family transcripts from Iso-Seq data", bioRxiv [Link](https://www.biorxiv.org/content/early/2018/01/10/246066).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/IsoCon/blob/master/LICENCE.txt).

