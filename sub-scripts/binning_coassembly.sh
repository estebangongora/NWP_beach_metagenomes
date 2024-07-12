#!/bin/bash

#Choose the number of threads that you want to use

for de in metabat/coassembly/*_depth.txt
do
sampleName=$(echo $de | sed 's/_[^_]*$//g' | cut -f 3 -d '/')

metabat2 -i spades/coassembly/$sampleName\/contigs.fasta -o metabat/coassembly/$sampleName -a $de -t <number_of_threads> -v

done
