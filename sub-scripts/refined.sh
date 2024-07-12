#!/bin/bash

mkdir -p refinem/coassembly/refined

for bam in bbmap_coassembly/*sorted.bam
do
sampleName=$(echo $bam | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

cp refinem/coassembly/$sampleName\/filtered_bins_taxon/*.fa refinem/coassembly/refined

done
