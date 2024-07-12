#!/bin/bash

echo "Beach	Domain	Phylum	Class	Order	Family	Genus	Species	Count	Count_pct	Abund	Abund_pct	Abund	Abund_pct" > results/metaerg_cds_taxonomy.tsv

for ba in bbmap_coassembly/*.bam
do
sampleName=$(echo $ba | sed 's/_[^_]*$//g' | cut -f 2 -d '/')

sed "s/^/$sampleName\t/" metaerg/coassembly/$sampleName/data/taxon.cds.profile.tab.txt | tail -n+2 >> results/metaerg_cds_taxonomy.tsv

sed -i "s/d__//" results/metaerg_cds_taxonomy.tsv
sed -i "s/;p__/\t/" results/metaerg_cds_taxonomy.tsv
sed -i "s/;c__/\t/" results/metaerg_cds_taxonomy.tsv
sed -i "s/;o__/\t/" results/metaerg_cds_taxonomy.tsv
sed -i "s/;f__/\t/" results/metaerg_cds_taxonomy.tsv
sed -i "s/;g__/\t/" results/metaerg_cds_taxonomy.tsv
sed -i "s/;s__/\t/" results/metaerg_cds_taxonomy.tsv

done
