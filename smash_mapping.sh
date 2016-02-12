#! /bin/bash

if [ "$SMASH_REF" = "" ] || [ ! -f $SMASH_REF ] ; then
    echo export SMASH_REF variable as fasta file path 1>&2
    exit 1
fi

if [ "$SMASH_CODE" = "" ] ; then
    echo export SMASH_CODE variable as code directory 1>&2
    exit 1
fi

# input consists of two space-separated lists of gzipped fastq files
id=$1
reads_1_gz="$2"
reads_2_gz="$3"

# map input data fastq.gz -> sam parts
$SMASH_CODE/mummer -verbose -rcref -qthreads 12 -nomap -samin -samout $SMASH_REF  <($SMASH_CODE/fastqs_to_sam <(zcat $reads_1_gz) <(zcat $reads_2_gz) 1)
mv mapout $id.mapout

# sam parts -> mappability tagged sam -> namesorted bam
$SMASH_CODE/mappability_tag $SMASH_REF <(cat $id.mapout/*.txt | head -n 100 | grep ^@ ; cat $id.mapout/*.txt | grep -v ^@ | perl -pe 's/^(\S+?)\/\S+\/\d+/\1/' ) | samtools view -Sbu - | samtools sort -n - $id.namesort

# namesorted bam -> good, unique mappings -> smash output
python $SMASH_CODE/smashMEM.py $id.namesort.bam 0 0 10000 4 > $id.smash.txt

# smash output -> map positions for major chromosomes only 
cat $id.smash.txt | awk '{print $4, $5}' | perl -ne 'print if /^chr(\d+|[XY]) \d+$/' > $id.positions.txt

