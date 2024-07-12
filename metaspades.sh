#!/bin/bash

#Choose the number of threads that you want to use

for fq in trimmed/*1P.fastq.gz
do
sampleName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

metaspades.py -1 $sampleName\_1P.fastq.gz -2 $sampleName\_2P.fastq.gz -t <number_of_threads> -m 300 -o spades/spades_$sampleName

done
