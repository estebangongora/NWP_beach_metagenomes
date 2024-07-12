#!/bin/bash

for db in Rhodococcus/contig-database/*.db
do
sampleName=$(echo $db | cut -f 3 -d '/' | cut -f 1 -d '.')

anvi-import-functions -c $db -i Rhodococcus/CANTHYD_results/noise_cutoff/$sampleName/$sampleName\_hmms.tsv

done
