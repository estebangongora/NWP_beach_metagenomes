#!/bin/bash

#Choose the number of threads that you want to use

mkdir -p phylogenomics/contig-database
phylo=phylogenomics/contig-database

mkdir -p phylogenomics/taxonomy

for fa in phylogenomics/original/*.fa
do
fileName=$(echo $fa | cut -f 3 -d '/' | cut -f 1,2 -d '.')
sampleName=$(echo $fa | cut -f 3 -d '/' | cut -f 1,2 -d '.' | sed 's/\./_/g')

anvi-gen-contigs-database -f $fa -o $phylo/$sampleName\.db -n $sampleName -T <number_of_threads>

anvi-run-hmms -c $phylo/$sampleName.db -T <number_of_threads>

anvi-run-scg-taxonomy -c $phylo/$sampleName.db -T <number_of_threads>

anvi-estimate-scg-taxonomy -c $phylo/$sampleName.db -o phylogenomics/taxonomy/$sampleName.txt -T <number_of_threads>

done
