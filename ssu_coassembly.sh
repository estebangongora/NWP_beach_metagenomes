#!/bin/bash

for sa in metabat/coassembly/*depth.txt
do
sampleName=$(echo $sa | sed 's/_[^_]*$//g' | cut -f 3 -d '/')

grep 'bac_16SrRNA' metaerg/coassembly/$sampleName\/data/rRNA.tab.txt >> results/16S_coassembly.txt
grep 'arc_16SrRNA' metaerg/coassembly/$sampleName\/data/rRNA.tab.txt >> results/16S_coassembly.txt
#grep 'euk_18SrRNA' metaerg/coassembly/$sampleName\/data/rRNA.tab.txt >> results/18S_coassembly.txt

done
