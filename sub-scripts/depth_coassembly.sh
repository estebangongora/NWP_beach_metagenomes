#!/bin/bash

mkdir -p metabat/coassembly

for ba in bbmap_coassembly/*.bam
do
sampleName=$(echo $ba | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

jgi_summarize_bam_contig_depths --outputDepth metabat/coassembly/$sampleName\_depth.txt $ba

done
