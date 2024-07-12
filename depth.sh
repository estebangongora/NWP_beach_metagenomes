#!/bin/bash

for ba in bbmap/*.bam
do
sampleName=$(echo $ba | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

jgi_summarize_bam_contig_depths --outputDepth metabat/$sampleName\_depth.txt bbmap/$sampleName\_sorted.bam

done
