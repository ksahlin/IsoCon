isocon_nontargeted
========

isocon_nontargeted is distributed as a python package supported on Linux / OSX with python v>=2.7, and 3.4-3.6, 3.5-dev and 3.6-dev [![Build Status](https://travis-ci.org/ksahlin/isocon_nontargeted.svg?branch=master)](https://travis-ci.org/ksahlin/isocon_nontargeted)


isocon_nontargeted is a tool for deriving *finished transcript sequences* from *Iso-Seq* reads. Input is a set of full-length-non-chimeric reads in fasta format and the CCS base call values as a bam file. The output is a set of predicted transcripts. isocon_nontargeted can be run as follows

```
isocon_nontargeted pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```

predicted transcripts are found in file **/path/to/output/final_candidates.fa**. Reads that could not be corrected or clustered are found in /path/to/output/not_converged.fa. For more instructions see below.


Table of Contents
=================

   * [isocon_nontargeted](#isocon_nontargeted)
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

The preferred way to install isocon_nontargeted is with pythons package installer pip.

### Using pip 

`pip` is pythons official package installer. This section assumes you have `python` (v2.7 or >=3.4) and a recent version of `pip` installed which should be included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 

With `python` and `pip` available, create a file `requirements.txt` with contents copied from [this file](https://github.com/ksahlin/isocon_nontargeted/blob/master/requirements.txt). Then, type in terminal 

```
pip install --requirement requirements.txt isocon_nontargeted
```

This should install isocon_nontargeted. With proper installation of **isocon_nontargeted**, you should be able to issue the command `isocon_nontargeted pipeline` to view user instructions. You should also be able to run isocon_nontargeted on this [small dataset](https://github.com/ksahlin/isocon_nontargeted/tree/master/test/data). Simply download the test dataset and run:

```
isocon_nontargeted pipeline -fl_reads [path/simulated_pacbio_reads.fa] -outfolder [output path]
```

`pip` will install the dependencies automatically for you. isocon_nontargeted has been built with python 2.7, 3.4-3.6 on Linux systems using [Travis](https://travis-ci.org/). For customized installation of latest master branch, see below.

### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as isocon_nontargeted has not been tested with earlier versions of these libraries. However, isocon_nontargeted may also work with earliear versions of these libaries.
* [edlib](https://github.com/Martinsos/edlib "edlib's Homepage"), for installation see [link](https://github.com/Martinsos/edlib/tree/master/bindings/python#installation) (>= v1.1.2)
* [networkx](https://networkx.github.io/) (>= v1.10)
* [ssw](https://github.com/vishnubob/ssw "Python wrapper for SSW"), for installation see [link](https://github.com/vishnubob/ssw#installation)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* [BioPython](http://biopython.org/wiki/Download) (>= v1.66)


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/isocon_nontargeted.git
cd isocon_nontargeted
./isocon_nontargeted
```


USAGE
-------

isocon_nontargeted's algorithm consists of two main phases; the error correction step and the statistical testing step. isocon_nontargeted can run these two steps in one go using `isocon_nontargeted pipeline`, or it can run only the correction or statistical test steps using `isocon_nontargeted get_candidates` and `isocon_nontargeted stat_filter` respectively. The preffered and most tested way is to use the entire pipeline `isocon_nontargeted pipeline`, but the other two settings can come in handy for specific cases. For example, running only `isocon_nontargeted get_candidates` will give more sequences if one is not concerned about precision and will also be faster, while one might use only `isocon_nontargeted stat_filter` using different parameters for a set of already constructed candidates in order to prevent rerunning the error correction step.


### Pipline

isocon_nontargeted takes two input files: (1) a fasta file of full length non-chimeric (flnc) CCS reads and (2) the bam file of CCS reads containing predicted base call quality values. The fasta file containing flnc can be obtained from PacBios Iso-Seq pipeline [ToFU](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki) and the bam file is the output of running the consensus caller algorthm [ccs](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md) on the Iso-Seq reads (ccs takes bam files so if you have bax files, convert them using [bax2bam](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md#input) ). isocon_nontargeted can then be run as

```
isocon_nontargeted pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> [--ccs </path/to/filename.ccs.bam>]
```

isocon_nontargeted also supports taking only the fasta read file as input. Without the base call quality values in `--ccs`, isocon_nontargeted will use an empirically estimated error model. The ability to take only the flnc fasta file as input is useful when the reads have been altered after the CCS base calling algorithm, \emph{e.g.}, from error correction using Illumina reads. **However, we highly recommend supplying the CCS quality values to isocon_nontargeted if CCS reads has not gone through any additional correction.** 

Simply omit the `--ccs` parameter if running isocon_nontargeted without base call quality values as
```
isocon_nontargeted pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

#### Output

The final high quality transcripts are written to the file `final_candidates.fa` in the output folder. If there was only one or two reads coming from a transcript, which is sufficiently different from other reads (exon difference), it will be output in the file `not_converged.fa`. This file may contain other erroneous CCS reads such as chimeras. The output also contains a file `cluster_info.tsv` that shows for each read which candidate it was assigned to in `final_candidates.fa`.

### get_candidates

Runs only the error correction step. The output is the converged candidates in a fasta file.

```
isocon_nontargeted get_candidates -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

### stat_filter

Runs only the statistical filtering of candidates.

```
isocon_nontargeted pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> -candidates <candidate_transcripts.fa>  [--ccs </path/to/filename.ccs.bam>]
```
Observe that `candidate_transcripts.fa` does not have to come from isocon_nontargeted's error correction algorithm. For example, this could either be a set of already validated transcripts to which one would like to see if they occur in the CCS reads, or they could be Illumina (or in other ways) corrected CCS reads.


### Parameters

```
    $ isocon_nontargeted pipeline --help
usage: Pipeline for obtaining non-redundant haplotype specific transcript isoforms using PacBio IsoSeq reads. pipeline
       [-h] -fl_reads FL_READS -outfolder OUTFOLDER [--ccs CCS]
       [--verbose] [--neighbor_search_depth NEIGHBOR_SEARCH_DEPTH]
       [--min_exon_diff MIN_EXON_DIFF] [--p_value_threshold P_VALUE_THRESHOLD]
       [--nr_cores NR_CORES] [--max_phred_q_trusted MAX_PHRED_Q_TRUSTED]
       [--min_candidate_support MIN_CANDIDATE_SUPPORT]
       [--ignore_ends_len IGNORE_ENDS_LEN] [--cleanup]
       [--prefilter_candidates]

optional arguments:
  -h, --help            show this help message and exit
  --ccs CCS             BAM/SAM file with CCS sequence predictions.
  --verbose        This will print more information abount workflow and
                        provide plots of similarity network etc.
  --neighbor_search_depth NEIGHBOR_SEARCH_DEPTH
                        Maximum number of pairwise alignments in search matrix
                        to find nearest_neighbor. [default =2**32]
  --min_exon_diff MIN_EXON_DIFF
                        Minimum consequtive base pair difference between two
                        neigborss in order to remove edge. If more than this
                        nr of consequtive base pair difference, its likely an
                        exon difference. [default =20]
  --p_value_threshold P_VALUE_THRESHOLD
                        Threshold for statistical test, filter everythin below
                        this threshold . [default = 0.01]
  --nr_cores NR_CORES   Number of cores to use.
  --max_phred_q_trusted MAX_PHRED_Q_TRUSTED
                        Maximum PHRED quality score trusted (T), linerarly
                        remaps quality score interval [0,93] --> [0, T].
                        Quality scores may have some uncertainty since T is
                        estimated from a consensus caller algorithm.
  --min_candidate_support MIN_CANDIDATE_SUPPORT
                        Required minimum number of reads converged to the same
                        sequence to be included in statistical test. [default
                        2]
  --ignore_ends_len IGNORE_ENDS_LEN
                        Number of bp to ignore in ends. If two candidates are
                        identical except in ends of this size, they are
                        collapsed and the longest common substing is chosen to
                        represent them. In the statistical test step,
                        the nearest neighbors are found based on ignoring the ends
                        of this size. Also indels "hanging off" ends of this size will not be tested.
                        [default 15].
  --cleanup             Remove everything except logfile.txt,
                        candidates_converged.fa and final_candidates.fa in
                        output folder. [default = False]
  --prefilter_candidates
                        Filter candidates if they are not consensus over any base pair 
                        in the candidate transcript formed from them, this can reduce runtime
                        without significant loss in true candidates. [default = False]

required arguments:
  -fl_reads FL_READS    Fast<a/q> file pacbio Reads of Insert.
  -outfolder OUTFOLDER  Outfolder.
```

CREDITS
----------------

Please cite [1] when using isocon_nontargeted.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† (2018) "isocon_nontargeted: Deciphering highly similar multigene family transcripts from Iso-Seq data", bioRxiv [Link](http/link)

LICENCE
----------------


<!-- isocon_nontargeted on general Iso-Seq datasets
-----------------------------------

isocon_nontargeted has been evaluated on targeted Iso-Seq data [cite], but the algorithm makes no explicit assumptions, or depends on the data being targeted. Therefore it also runs on any general Iso-Seq data. To demonstrate this We ran isocon_nontargeted on the publicly availabe datasets from pacbio [MCF7](link) and [Alzheimer](link) using the following command

```
isocon_nontargeted pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --neighbor_search_depth 100
```

| Dataset | runtime  | peak memory | final_candidates | corr | not_corr | *TOFU* | *nr original CCS* | 
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| [MCF7](http://www.pacb.com/blog/data-release-human-mcf-7-transcriptome/) | 24h43m  | <1.9Gb  | 2,169 | 18,458 | 401,885<sup>[1](#myfootnote1)</sup> |  55,770 | 518,701 |
|[Alzheimer](http://www.pacb.com/blog/data-release-alzheimer-brain-isoform-sequencing-iso-seq-dataset/)| time | Content Cell  | Content Cell  |

TODO1: Output three separate files: validated_transcripts.fa,  corrected_not_validated.fa, not_converged.fa 
TODO2: Show set intersection (allowing end differences between final_cand + corr to TOFU transcripts). Also do other analysis.

isocon_nontargeted predicted 2169 and Y transcripts and had a runtime of 24h43m and Y, for MCF7 and Alzheimer datasets respectively, on a 64 core machine with a peak memory usage of Z Gb. isocon_nontargeteds output are found [here](link).

Manual BLAT of 20 sequences from the "corrected_but_not_converged" predictions to human show an alignment identity increase from 94-98% of the ccs reads up to 99.5-99.9% for the corrected reads.

<a name="footnote_not_converged">1</a>: Footnote content goes here -->


<!-- ### Entry points in isocon_nontargeted
isocon_nontargeted has three different modes `pipeline` (end to end), `get_candidates` only correct and cluster reads (no statistical testing), as well as `stat_filter` (only statistical testing of a set of candidate transcripts).

```
$ isocon_nontargeted --help
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
``` -->