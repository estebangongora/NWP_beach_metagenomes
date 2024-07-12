#!/bin/bash

#Fastq files must be named with the following format: <base_name>_<year>_1.fastq.gz for the forward reads and <base_name>_<year>_2.fastq.gz
#e.g., Tank_2019_1.fastq.gz, Ale_2018_2.fastq.gz, (also works with samples without a year like AB_HT_1.fastq.gz), etc.
#Edit the path to include the location of your trimommatic folder
#Choose the number of threads that you want to use
#Trimmed reads will be outputted to the trimmed folder

for fq in fastq/*1.fastq.gz
do
base=$(echo $fq | sed 's/_[^_]*$//g' | cut -f 2 -d '/')
sampleName=$(echo $base | sed -r 's/.{12}//g')

java -jar <path_to_trimmomatic_folder>/trimmomatic.jar PE -threads <number_of_threads> -basein $fq -baseout trimmed/$sampleName\.fastq.gz SLIDINGWINDOW:4:20 LEADING:10 TRAILING:10 ILLUMINACLIP:NexteraPE-Pe.fa:2:30:10

done