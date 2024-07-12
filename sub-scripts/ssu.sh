#!/bin/bash

for sa in metabat/*depth.txt
do
sampleName=$(echo $sa | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

grep 'bac_16SrRNA' metaerg/$sampleName\/data/rRNA.tab.txt >> results/16S.txt
grep 'arc_16SrRNA' metaerg/$sampleName\/data/rRNA.tab.txt >> results/16S.txt
#grep 'euk_18SrRNA' metaerg/$sampleName\/data/rRNA.tab.txt >> results/18S.txt

done
