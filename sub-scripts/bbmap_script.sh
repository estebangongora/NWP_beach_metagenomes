#!/bin/bash

#Edit the path to include the location of your BBMap folder

for fq in trimmed/*1P.fastq.gz
do
sampleName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

<path_to_bbmap_folder>/bbmap.sh in1=$fq in2=trimmed/$sampleName\_2P.fastq.gz ref=spades/spades_$sampleName\/contigs.fasta out=bbmap/$sampleName\.sam bamscript=bs.sh; sh bs.sh

done
