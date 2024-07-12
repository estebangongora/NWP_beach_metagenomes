#!/bin/bash

#Choose the number of threads that you want to use

for fi in Rhodococcus/original/*.filtered.filtered.fa
do
fileName=$(echo $fi | cut -f 1,2 -d '.' | sed 's/\./_/g' | sed 's/$/.fa/g')

mv "$fi" $(echo "$fileName")

done

for fa in Rhodococcus/original/*.fa
do
sampleName=$(echo $fa | cut -f 3 -d '/' | cut -f 1 -d '.')

anvi-script-reformat-fasta $fa -o Rhodococcus/contig-fasta/$sampleName.fa --simplify-names --prefix $sampleName --seq-type NT

anvi-gen-contigs-database -f Rhodococcus/contig-fasta/$sampleName.fa -o Rhodococcus/contig-database/$sampleName\.db -n $sampleName -T <number_of_threads>

anvi-run-hmms -c Rhodococcus/contig-database/$sampleName.db -T <number_of_threads>

anvi-run-hmms -c Rhodococcus/contig-database/$sampleName.db -H Rhodococcus/CANTHYD/ --hmmer-output-dir Rhodococcus/CANTHYD_results/$sampleName --domain-hits-table -T <number_of_threads> --hmmer-program hmmsearch

anvi-run-scg-taxonomy -c Rhodococcus/contig-database/$sampleName.db -T <number_of_threads>

anvi-estimate-scg-taxonomy -c $phylo/$sampleName.db -o phylogenomics/taxonomy/$sampleName.txt -T <number_of_threads>

done
