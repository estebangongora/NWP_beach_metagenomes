#!/bin/bash

#Choose the number of threads that you want to use

for de in metabat/*.txt
do
sampleName=$(echo $de | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

metaerg.pl --prefix $sampleName --outdir metaerg/$sampleName --locustag $sampleName --cpus <number_of_threads> --depth metabat/$sampleName\_depth.txt spades/spades_$sampleName\/contigs.fasta

done
