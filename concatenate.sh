#!/bin/bash

cd trimmed/

pigz -dk *.gz
mkdir fastq
mv *.fastq fastq/
mkdir concatenated

for fq in fastq/*e_1P.fastq
do
sampleName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_')
beachName=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 1,2 -d '_' | cut -f 2 -d '/')

cat $sampleName\_e_1P.fastq $sampleName\_i_1P.fastq > concatenated/$beachName\_1P.fastq
cat $sampleName\_e_2P.fastq $sampleName\_i_2P.fastq > concatenated/$beachName\_2P.fastq

echo "$beachName forward" >> concatenation_results.txt
awk '{s++}END{print s/4}' $sampleName\_e_1P.fastq $sampleName\_i_1P.fastq >> concatenation_results.txt
awk '{s++}END{print s/4}' concatenated/$beachName\_1P.fastq >> concatenation_results.txt

echo "$beachName reverse" >> concatenation_results.txt
awk '{s++}END{print s/4}' $sampleName\_e_2P.fastq $sampleName\_i_2P.fastq >> concatenation_results.txt
awk '{s++}END{print s/4}' concatenated/$beachName\_2P.fastq >> concatenation_results.txt

done

cd concatenated/

pigz *.fastq

cd ../..

echo "done"
