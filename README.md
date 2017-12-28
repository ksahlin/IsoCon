IsoCon
========

IsoCon is distributed as a python package supported on Linux / OSX with python v>=2.7, and 3.4-3.6, 3.5-dev and 3.6-dev [![Build Status](https://travis-ci.org/ksahlin/IsoCon.svg?branch=master)](https://travis-ci.org/ksahlin/IsoCon)


IsoCon is a tool for deriving *finished transcript sequences* from *Iso-Seq* reads. Input is a set of full-length-non-chimeric reads in fasta format and the CCS base call values as a bam file. The output is a set of predicted transcripts. IsoCon can be run as follows

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```

predicted transcripts are found in file **/path/to/output/final_candidates.fa**. Reads that could not be corrected or clustered are found in /path/to/output/not_converged.fa. For more instructions see below.


Table of Contents
=================

   * [IsoCon](#isocon)
      * [INSTALLATION](#installation)
            * [Using pip](#using-pip)
            * [Downloading source from GitHub](#downloading-source-from-github)
               * [Dependencies](#dependencies)
      * [USAGE](#usage)
         * [Output](#output)
         * [Parameters](#parameters)
      * [CREDITS](#credits)
      * [LICENCE](#licence)

INSTALLATION
----------------

The preferred way to install IsoCon is with pythons package installer pip.

#### Using pip (12/27/17 -- under construction) 

`pip` is pythons official package installer. This section assumes you have `python` and a recent version of `pip` installed which should be included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 

Now, create a file `requirements.txt` with these [lines](https://github.com/ksahlin/IsoCon/blob/master/requirements.txt) and type in terminal `pip install --requirement requirements.txt IsoCon` . This should install IsoCOn. With proper installation of **IsoCon**, you should be able to issue the command `IsoCon pipeline` to view user instructions. `pip` will install the dependencies automatically for you. For customized installation of latest master branch, see below.

#### Downloading source from GitHub

##### Dependencies

Install the following dependencies, either with `pip install` or using instructions in links below. Versions in parenthesis are suggested as IsoCon has not been tested with earlier versions of these libraries, IsoCon may however work with earliear versions of these libaries too.
* [edlib](https://github.com/Martinsos/edlib "edlib's Homepage"), for installation see [link](https://github.com/Martinsos/edlib/tree/master/bindings/python#installation) (>= v1.1.2)
* [networkx](https://networkx.github.io/) (>= v1.10)
* [ssw](https://github.com/vishnubob/ssw "Python wrapper for SSW"), for installation see [link](https://github.com/vishnubob/ssw#installation)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* [BioPython](http://biopython.org/wiki/Download) (>= v1.66)


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/IsoCon.git
cd IsoCon
./IsoCon
```


USAGE
-------

IsoCon takes two input files: (1) a fasta file of full length non-chimeric (flnc) CCS reads and (2) the bam file of CCS reads containing predicted base call quality values. The fasta file containing flnc can be obtained from PacBios Iso-Seq pipeline [ToFU](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki) and the bam file is the output of running the consensus caller algorthm [ccs](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md) on the Iso-Seq reads (ccs takes bam files so if you have bax files, convert them using [bax2bam](https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md#input) ). IsoCon can then be run as

```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```

IsoCon also supports taking only the fasta read file as input. Without the base call quality values, IsoCon will use an empirically estimated error model. The ability to take only the flnc fasta file as input is useful when the reads have been altered after the CCS base calling algorithm, \emph{e.g.}, from **error correction using Illumina reads**. However, we highly recommend supplying the CCS quality values to IsoCon if CCS reads has not gone through any additional correction. 

Simply omit the `--ccs` parameter if running IsoCon without base call quality values as
```
IsoCon pipeline -fl_reads <flnc.fasta> -outfolder </path/to/output>
```

### Output

The final high quality transcripts are written to the file `final_candidates.fa` in the output folder. If there was only one or two reads coming from a transcript, which is sufficiently different from other reads (exon difference), it will be output in the file `not_converged.fa`. This file may contain other erroneous CCS reads such as chimeras.


### Parameters

```
    $ IsoCon pipeline --help
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

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† (2017) "IsoCon: Deciphering highly similar multigene family transcripts from Iso-Seq data", bioRxiv [Link](http/link)

LICENCE
----------------


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


<!-- ### Entry points in IsoCon
IsoCon has three different modes `pipeline` (end to end), `get_candidates` only correct and cluster reads (no statistical testing), as well as `stat_filter` (only statistical testing of a set of candidate transcripts).

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
``` -->