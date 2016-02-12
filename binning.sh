#! /bin/bash

if [ "$SMASH_REF" = "" ] || [ ! -f $SMASH_REF ] ; then
    echo export SMASH_REF variable as fasta file path 1>&2
    exit 1
fi

if [ "$SMASH_CODE" = "" ] ; then
    echo export SMASH_CODE variable as code directory 1>&2
    exit 1
fi

# sample identifier
id=$1

# list of input positions: chr pos
input_positions=$id.positions.txt

# directory with bin information files bins.txt bad.txt and gc.txt
bindir=$2

# bin information:
# chr start_chrpos start_abspos stop_chrpos bin_len n_maps_expected
bins=$bindir/bins.txt

# bad bins (optionally exists):
# 1-based bad bin indexes
bad=$bindir/bad.txt

# gc information
# chr start_chrpos stop_chrpos bin_len	int int int gc_content int int
# columns labeled int are unused by this codebase
gc=$bindir/gc.txt

# input chr, pos -> binned data
$SMASH_CODE/varbin.py $input_positions $bins varbin.txt $id.stats.txt $SMASH_REF.bin/chrom_sizes.txt > $id.varbin.out.txt

# binned data -> gc corrected binned data and segmented results amd plots
R CMD BATCH --no-restore --no-save '--args sample="'$id'" gc="'$gc'" bad="'$bad'"' $SMASH_CODE/cbs.r
mv cbs.r.Rout $id.cbs.r.Rout
mv varbin.txt $id.varbin.txt

exit 0
