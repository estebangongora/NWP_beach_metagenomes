#!/bin/bash

#Choose the number of threads that you want to use
#Edit the path to include the location of your RefineM folder with the required databases

mkdir -p refinem/coassembly

for bam in bbmap_coassembly/*sorted.bam
do
sampleName=$(echo $bam | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

refinem scaffold_stats spades/coassembly/$sampleName\/contigs.fasta bins_coassembly/$sampleName refinem/coassembly/$sampleName $bam -x fa -c <number_of_threads>

refinem outliers refinem/coassembly/$sampleName\/scaffold_stats.tsv refinem/coassembly/$sampleName\/outliers/

refinem filter_bins bins_coassembly/$sampleName refinem/coassembly/$sampleName\/outliers/outliers.tsv refinem/coassembly/$sampleName\/filtered_bins -x fa

refinem call_genes refinem/coassembly/$sampleName\/filtered_bins/ refinem/coassembly/$sampleName\/called -x fa -c <number_of_threads>

refinem taxon_profile refinem/coassembly/$sampleName\/called/ refinem/coassembly/$sampleName\/scaffold_stats.tsv <path_to_refineM_database_folder>/gtdb_r95_protein_db.2020-07-30.faa.dmnd <path_to_refineM_database_folder>/gtdb_r95_taxonomy.2020-07-30.tsv refinem/coassembly/$sampleName\/taxon_profile -c <number_of_threads>

refinem taxon_filter refinem/coassembly/$sampleName\/taxon_profile/ refinem/coassembly/$sampleName\/taxon_filter.tsv -c <number_of_threads>

refinem filter_bins refinem/coassembly/$sampleName\/filtered_bins/ refinem/coassembly/$sampleName\/taxon_filter.tsv refinem/coassembly/$sampleName\/filtered_bins_taxon -x fa

refinem ssu_erroneous refinem/coassembly/$sampleName\/filtered_bins_taxon refinem/coassembly/$sampleName\/taxon_profile <path_to_refineM_database_folder>/ssu <path_to_refineM_database_folder>/gtdb_r95_taxonomy.2020-07-30.tsv refinem/coassembly/$sampleName\/filtered_bins_ssu -x fa -c <number_of_threads>

done
