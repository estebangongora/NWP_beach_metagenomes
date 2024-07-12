#!/bin/bash

#Edit the path to include the location of your BBMap folder

mkdir -p bbmap_coassembly
count=0

for fq in trimmed/concatenated/*_1P.fastq.gz
do
sampleName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_')
beachName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_' | cut -f 3 -d '/')
count=$((count + 1))

<path_to_bbmap_folder>/bbmap.sh in1=$fq in2=$sampleName\_2P.fastq.gz ref=spades/coassembly/$beachName\/contigs.fasta build=$count path=bbmap_coassembly/ out=bbmap_coassembly/$beachName\.sam bamscript=bs.sh; sh bs.sh

<path_to_bbmap_folder>/pileup.sh in=bbmap_coassembly/$beachName\.sam out=bbmap_coassembly/$beachName\_stats.txt

done
