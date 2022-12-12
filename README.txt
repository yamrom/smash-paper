SMASH sequencing mapping and binning code

see mumdex.com for later versions of general-purpose genomic mapping software

Instructions for using this package:

to compile C++ code, first edit shared.mk to point GCC_DIR
to a recent version of GCC
this code was tested with GCC 4.9.2

set environment variables to point to your reference genome fasta
and the smash code directories, in a manner similar to this:

export SMASH_CODE=/mnt/wigtop2/data/safe/yamrom/analysis/smash-copy-number/smash-paper-code
export SMASH_REF=/mnt/wigtop2/data/safe/yamrom/analysis/smash-copy-number/chrAll.fa

to generate a binary index and assorted other files for your genome, run
$SMASH_CODE/index_setup.sh

the genome chrAll.fa used in the SMASH paper is available for download
by ftp at wigserv4.cshl.edu with user smashpaper and password smashpaper

to generate bin definition files, choose a method and follow the
format of the sample bin files in the sample_bins subdirectory.
The columns of the files bins.txt gc.txt and bad.txt are described in binning.sh

possible bin definition methods include:
1. fixed number of loci per bin
2. fixed number of mappings per bin using all reference genome kmers
3. theoretical expectation using genome mappability lengths and SMASH
   fragment size distribution
4. theoretical expectation using restriction site information
5. empirically determined bins using a reference sample or samples

we believe that #5 with many samples offers the best result.
Sample bins included here have 50k, 100k and 500k bins using method #2
and kmer size 50 mapped with bowtie.

to map your data and extract mapped positions, run on your fastq input files:
$SMASH_CODE/smash_mapping.sh some_id fastq_r1.gz fastq_r2.gz

to bin your data, GC correct, segment and plot:
$SMASH_CODE/binning.sh some_id bins_dir

