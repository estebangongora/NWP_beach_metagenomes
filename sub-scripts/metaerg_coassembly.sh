#!/bin/bash

#Choose the number of threads that you want to use

mkdir -p metaerg/coassembly

for de in metabat/coassembly/*.txt
do
sampleName=$(echo $de | sed 's/_[^_]*$//g' | cut -f 3 -d '/')

metaerg.pl --prefix $sampleName --outdir metaerg/coassembly/$sampleName --locustag $sampleName --cpus <number_of_threads> --depth $de spades/coassembly/$sampleName\/contigs.fasta

done
