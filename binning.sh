#!/bin/bash

#Choose the number of threads that you want to use

for de in metabat/*.txt
do
sampleName=$(echo $de | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

metabat2 -i spades/spades_$sampleName\/contigs.fasta -o metabat/$sampleName -a metabat/$sampleName\_depth.txt -t <number_of_threads> -v

done
