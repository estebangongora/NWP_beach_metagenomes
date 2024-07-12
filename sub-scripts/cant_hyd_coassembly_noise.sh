#!/bin/bash

#Edit the path to include the location of your CANT-HYD folder containing the HMMs
#Choose the number of threads that you want to use

mkdir -p cant-hyd_coassembly/noise
mkdir -p cant-hyd_coassembly/noise/metagenomes
mkdir -p cant-hyd_coassembly/noise/bins

for sa in metabat/coassembly/*depth.txt
do
sampleName=$(echo $sa | sed 's/_[^_]*$//g' | cut -f 3 -d '/')

hmmsearch -o cant-hyd_coassembly/noise/metagenomes/$sampleName\.out --tblout cant-hyd_coassembly/noise/metagenomes/$sampleName\.tblout --cut_nc --cpu <number_of_threads> <path_to_CANT-HYD_folder>/CANT-HYD.hmm metaerg/coassembly/$sampleName\/data/cds.faa

done

for bi in refinem/coassembly/refined/*.fa
do
binName=$(echo $bi | sed 's/\.[^.]*$//' | cut -f 4 -d '/' | cut -f 1,2 -d '.')

hmmsearch -o cant-hyd_coassembly/noise/bins/$binName\.out --tblout cant-hyd_coassembly/noise/bins/$binName\.tblout --cut_nc --cpu <number_of_threads> <path_to_CANT-HYD_folder>/CANT-HYD.hmm metaerg/bins/coassembly/$binName\/data/cds.faa

done



echo "Beach	Gene	Count	CDS" > results/cant-hyd_noise_metagenome_results.txt
echo "MAG	Beach	Gene	Count	CDS" > results/cant-hyd_noise_bin_results.txt

for sa in cant-hyd_coassembly/noise/metagenomes/*.tblout
do
sampleName=$(echo $sa | sed 's/\.[^.]*$//' | cut -f 4 -d '/')
cds=$(grep "CDS count" metaerg/coassembly/$sampleName/data/master.stats.txt | cut -f 2 -d '	')

while read gene
do
echo "$sampleName	$gene	$(grep -cw "$gene" $sa)	$cds"
done < cant-hyd_genes.txt >> results/cant-hyd_noise_metagenome_results.txt

done

for bi in cant-hyd_coassembly/noise/bins/*.tblout
do
binName=$(echo $bi | sed 's/\.[^.]*$//' | cut -f 4 -d '/')
beachName=$(echo $binName | cut -f 1 -d '.')
cds=$(grep "CDS count" metaerg/bins/coassembly/$binName/data/master.stats.txt | cut -f 2 -d '	')

while read gene
do
echo "$binName	$beachName	$gene	$(grep -cw "$gene" $bi)	$cds"
done < cant-hyd_genes.txt >> results/cant-hyd_noise_bin_results.txt

done
