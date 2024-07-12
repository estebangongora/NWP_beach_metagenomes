#!/bin/bash

echo "Sample	Beach	Site	Site_long	DNA	Coverage	Domain	Phylum	Class	Order	Family	Genus	Species" > results/16S_taxonomy.txt

cat results/16S.txt | while read line
do

sample=$(echo $line | cut -f 1 -d '|')
beach=$(echo $line | cut -f 1 -d '|' | sed 's/_.*//g')
site=$(echo $line | cut -f 1 -d '|' | cut -f 1,2 -d '_' )
if [ $site = "AB_HT" ]
then site_long=$(echo "Assistance Bay - High tide")
elif [ $site = "AB_LT" ]
then	site_long=$(echo "Assistance Bay - Low tide")
elif [ $site = "Ale_1" ]
then	site_long=$(echo "Alert - 1")
elif [ $site = "Cam_Cam" ]
then	site_long=$(echo "Cambridge Bay")
elif [ $site = "Dump_2018" ]
then	site_long=$(echo "Dump Beach - 2018")
elif [ $site = "Dump_2019" ]
then	site_long=$(echo "Dump Beach - 2019")
elif [ $site = "Dyna_2018" ]
then	site_long=$(echo "Dynamite Beach - 2018")
elif [ $site = "Dyna_2019" ]
then	site_long=$(echo "Dynamite Beach - 2019")
elif [ $site = "Nani_A1" ]
then	site_long=$(echo "Nanisivik - A1")
elif [ $site = "Nani_B1" ]
then	site_long=$(echo "Nanisivik - B1")
elif [ $site = "Tank_2018" ]
then	site_long=$(echo "Tank Farm - 2018")
elif [ $site = "Tank_2019" ]
then	site_long=$(echo "Tank Farm - 2019")
elif [ $site = "Tupi_2018" ]
then	site_long=$(echo "Tupirvik - 2018")
elif [ $site = "Tupi_2019" ]
then	site_long=$(echo "Tupirvik - 2019")
else
	site_long=$(echo "remove")
fi
dna=$(echo $line | cut -f 1 -d '|' | cut -f 3 -d '_' )
if [ $dna = "i" ]
then
	dna_type=$(echo "iDNA")
else
	dna_type=$(echo "eDNA")
fi
taxonomy=$(echo $line | cut -f 8 -d ' ')
coverage=$(echo $line | grep -Eo 'cov_[0-9]+\.[0-9]+' | cut -f 2 -d '_')

echo -e $sample'\t'$beach'\t'$site'\t'$site_long'\t'$dna_type'\t'$coverage'\t'$taxonomy >> results/16S_taxonomy.txt

done

sed -i 's/\;/\t/g' results/16S_taxonomy.txt
