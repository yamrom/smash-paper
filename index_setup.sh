#! /bin/bash

if [ "$SMASH_REF" = "" ] || [ ! -f $SMASH_REF ] ; then
    echo export SMASH_REF variable as fasta file path 1>&2
    exit 1
fi

if [ "$SMASH_CODE" = "" ] ; then
    echo export SMASH_CODE variable as code directory 1>&2
    exit 1
fi

if [ -e $SMASH_REF.bin ] ; then
    echo binary index directory $SMASH_REF.bin already exists - quitting 1>&2
    exit 1
fi

# Create suffix array binary files
$SMASH_CODE/mummer -verbose -rcref $SMASH_REF dummy

# Create mappability binary file
$SMASH_CODE/mummer -verbose -rcref -mappability $SMASH_REF $SMASH_REF.bin/map.bin

# Create fai file
samtools faidx $SMASH_REF

# Make chromosome sizes text file
n=0; cat $SMASH_REF.fai | awk '{print $1,$2}' | grep -v _ | while read chr size ; do echo -e "$chr\t$size\t$n" ; n=$((n+size)) ; done > $SMASH_REF.bin/chrom_sizes.txt

# Make sam header file
cat $SMASH_REF.fai | awk '{print $1,$2}' | while read chr size ; do echo -e "@SQ\tSN:$chr\tLN:$size" ; done > $SMASH_REF.bin/sam_header.txt

