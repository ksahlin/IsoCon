# Change Log
Notable changes will be documented in this file.



## [Unreleased]
Changes appearing in the develop branch.



## [0.3.0] - 2018-02-02

Major updates to speed and code readability, minor bugfixes. Previous versions are deprecated. 

### Fixed
- Bugfix when —nearest_neighbor_depth set. Previous version would not explore less than specified number of sequences when assigning reads to candidates (statistical test step).

### Added
- Added parameter --min_test_ratio X (default 5). This parameter omits testing candidates c to t if c has X times more support. This will speed up the algorithm because it omits tests that will (most likely) be significant as c is dominant to t.
- Added CHANGELOG and GPL LICENSE.

### Changed
- IsoCon does not longer built the multialignment matrix (MAM) in the statistical testing step. We now get support, base probabilities etc. by (1) aligning only c to t to get positions where they differ, (2) Obtain base qualities over these positions by using (the already created) read alignments to c and t (a read is assignet to either c or t). This implementation therefor skips the realignment of all reads to the reference t as well as the creation of the MAM. This gives a speed-up of 5-20x for the statistical test (most speedup for longer noisier sequences).
- Improved speed (~2-3x) in multialignment function (Now used only in the correction step).
- Now IsoCon store read alignments to the candidate when entering the hypothesis test alignments to t to avoid realignment.
- significant re-write of code in some regions, such as the statisitcal test.



## [0.2.5.1] - 2017-01-07

### Fixed
- Bugfix in if-statement described in Supplementary Section A: "Estimating the probability of a sequencing error"[biorxiv-suppl].
- Bugfix in how to treat tiebreakers, described in Supplementary Section A: "implementation details" in [biorxiv-suppl]. 

### Added
- Testdata and instructions for running IsoCon on testdata.
- Automatic builds and testing with Travis.
- Installation through pip now possible.
- Parameter --verbose and also removed lot of prints to stdout.
- Added parameter --min_exon_diff to break alignment with this many consecutive '-' (an indel). Previously this was hardcoded.
- Made the mapping quality values upper bound T (described in Supplementary Section A: "Estimating the probability of a sequencing error"[biorxiv-suppl]) a parameter --max_phred_q_trusted instead of hardcoded value.

### Changed
- change parameter --single_core (a flag that was false by default and IsoCon would use all cores available) to a more flexible format --nr_cores where the used can specify how many cores.
- Changed terminology. All occurences of "minimizer" is changed to "nearest_neighbor" (or "neighbor" in parameters) to adapt for new notation of nearest neighbor graph instead of minimizer graph.


### Removed
- Removed option --barcodes as it serves no purpose anymore --- if we have barcodes they have been detected and reads have been split into batches in a upstream step.



## [0.2.4] - 2017-10-16
Version used for results in the bioRxiv preprint [biorxiv].





[biorxiv-suppl]: https://www.biorxiv.org/content/biorxiv/suppl/2018/01/10/246066.DC1/246066-1.pdf
[biorxiv]: https://www.biorxiv.org/content/early/2018/01/10/246066

[Unreleased]: https://github.com/ksahlin/IsoCon/tree/develop
[0.2.4]:   https://github.com/ksahlin/IsoCon/tree/85eb122cba106b7f747bd1dafaf7f232a16b52e9
[0.2.5.1]: https://github.com/ksahlin/IsoCon/tree/d4567f3f0d85a31aa406f20ef51558d4f8a8c4e0
[0.3.0]: https://github.com/ksahlin/IsoCon