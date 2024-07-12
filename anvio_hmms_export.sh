#!/bin/bash

mkdir -p Rhodococcus/CANTHYD_results/noise_cutoff/export

for db in Rhodococcus/contig-database/*.db
do
sampleName=$(echo $db | cut -f 3 -d '/' | cut -f 1 -d '.')

anvi-export-table $db --table hmm_hits -f 'gene_callers_id, source, gene_name, e_value' -o Rhodococcus/CANTHYD_results/noise_cutoff/export/$sampleName\_hmms.tsv

done
