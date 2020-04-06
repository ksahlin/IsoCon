IsoCon
========

IsoCon is distributed as a python package supported on Linux / OSX with python v2.7, and version >=3.6 [![Build Status](https://travis-ci.org/ksahlin/IsoCon.svg?branch=master)](https://travis-ci.org/ksahlin/IsoCon)


IsoCon is a tool for reconstructing highly similar sequences present in a dataset of from long noisy reads. Its original use case was transcripts from highly similar gene copies ([paper here](https://www.nature.com/articles/s41467-018-06910-x)), however the methodology extends to any dataset where sequences spans the region(s) of interest end-to-end. IsoCon use examples: 

* Deriving *finished transcripts* from *Iso-Seq* or ONT reads from *targeted* sequencing of gene families using primers. 
* Deriving consensus sequence from several passes of long noisy reads (e.g., pacbio polymerase reads to CCS or ONT Rolling Circle Amplification to Concatemeric Consensus (R2C2)).
* Deriving viral strains from  reads (assuming the reads spans the viral sequence, e.g., as for HIV).
* Deriving consensus ribosomal RNA.
* Deriving consensus from any targeted amplicone based sequencing technique.


Simplest usage is an input file of fastq or fasta containing reads. IsoCon can be run as follows

```
IsoCon pipeline -fl_reads <reads.fastq> -outfolder </path/to/output>
```

or 

```
IsoCon pipeline -fl_reads <reads.fasta> -outfolder </path/to/output> --ccs </path/to/filename.ccs.bam>
```



predicted transcripts are found in file **/path/to/output/final_candidates.fa**. Reads that could not be corrected or clustered are found in /path/to/output/not_converged.fa. 

* Can IsoCon be run on nontargeted Iso-Seq datasets? [see here](https://github.com/ksahlin/IsoCon/issues/2). 
* How does my data set affect the runtime? [see here](https://github.com/ksahlin/IsoCon/issues/3) 

For more instructions see below.


Table of Contents
=================

  * [Table of Contents](#Table-of-Contents)
  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Using pip](#Using-pip)
    * [Downloading source from GitHub](#Downloading-source-from-GitHub)
    * [Dependencies](#Dependencies)
  * [USAGE](#USAGE)
    * [Pipline](#Pipline)
    * [Output](#Output)
    * [get_candidates](#get_candidates)
    * [stat_filter](#stat_filter)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)


INSTALLATION
----------------

### Using conda
Conda is the preferred way to install IsoCon.

1. Create and activate a new environment called IsoCon

```
conda create -n IsoCon python=3 pip 
conda activate IsoCon
```

2. Install IsoCon 

```
pip install IsoCon
```

3. You should now have 'IsoCon' installed; try it:

```
IsoCon --help
```

Upon start/login to your server/computer you need to activate the conda environment "IsoCon" to run IsoCon as:

```
conda activate IsoCon
```


### Using pip 

`pip` is pythons official package installer. This section assumes you have `python` (v2.7 or >=3.6) and a recent version of `pip` installed which should be included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 

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

Using quality values (fastq) is preferred over fasta as IsoCon uses the quality values for statistical analysis.


```
IsoCon pipeline -fl_reads <reads.fast[a/q]> -outfolder </path/to/output>
```


#### Output

The final high quality transcripts are written to the file `final_candidates.fa` in the output folder. If there was only one or two reads coming from a transcript, which is sufficiently different from other reads (exon difference), it will be output in the file `not_converged.fa`. This file may contain other erroneous reads such as chimeras. The output also contains a file `cluster_info.tsv` that shows for each read which candidate it was assigned to in `final_candidates.fa`.

### get_candidates

Runs only the error correction step. The output is the converged candidates in a fasta file.

```
IsoCon get_candidates -fl_reads <flnc.fast[a/q]> -outfolder </path/to/output>
```

### stat_filter

Runs only the statistical filtering of candidates.

```
IsoCon pipeline -fl_reads <flnc.fast[a/q]> -outfolder </path/to/output> -candidates <candidate_transcripts.fa> 
```
Observe that `candidate_transcripts.fa` does not have to come from IsoCon's error correction algorithm. For example, this could either be a set of already validated transcripts to which one would like to see if they occur in the reads, or they could be Illumina (or in other ways) corrected CCS reads.



CREDITS
----------------

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin*, Marta Tomaszkiewicz*, Kateryna D. Makova†, Paul Medvedev† Deciphering highly similar multigene family transcripts from iso-seq data with isocon. Nature Communications, 9(1):4601, 2018. [Link](https://www.nature.com/articles/s41467-018-06910-x).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/IsoCon/blob/master/LICENCE.txt).

