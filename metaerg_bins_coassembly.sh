#!/bin/bash

mkdir -p metaerg/bins/coassembly

for bi in refinem/coassembly/refined/*.fa
do
binName=$(echo $bi | sed 's/\.[^.]*$//' | cut -f 4 -d '/' | cut -f 1,2 -d '.')
sampleName=$(echo $binName | sed 's/\.[^.]*$//')

#Step1, extracting the gff-format annotations for the contigs from the total metaerg dataset annotation:
fastaContig2Gff.pl -c $bi -g metaerg/coassembly/$sampleName\/data/all.gff > metaerg/bins/coassembly/$binName\.gff

#Step 2, generating the html reports for the extracted contig subset
output_reports.pl -g metaerg/bins/coassembly/$binName\.gff -f $bi -o metaerg/bins/coassembly/$binName

done
