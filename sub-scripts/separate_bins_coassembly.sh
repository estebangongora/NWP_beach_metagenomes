#!/bin/bash

cp -r metabat/coassembly bins_coassembly
cd bins_coassembly

for de in *.txt
do
sampleName=$(echo $de | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_')

mkdir $sampleName

find . -name "*.fa" | grep "$sampleName" | xargs mv -t "$sampleName"

rm $de

done
