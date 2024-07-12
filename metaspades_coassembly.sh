#!/bin/bash

#Choose the number of threads that you want to use

for fq in trimmed/concatenated/*1P.fastq.gz
do
sampleName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 3 -d '/')
#beachName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_' | cut -f 3 -d '/')

metaspades.py -1 $fq -2 trimmed/concatenated/$sampleName\_2P.fastq.gz -t <number_of_threads> -m 500 -o spades/coassembly/$sampleName

done
