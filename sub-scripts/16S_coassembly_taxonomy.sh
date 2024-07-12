#!/bin/bash

echo "Sample	Beach	Site	Coverage	Domain	Phylum	Class	Order	Family	Genus	Species" > results/16S_coassembly_taxonomy.txt

cat results/16S_coassembly.txt | while read line
do

sample=$(echo $line | cut -f 1 -d '|')
beach=$(echo $line | cut -f 1 -d '|' | sed 's/_.*//g')
if [ $sample = "AB_HT" ]
then site=$(echo "Assistance Bay - High tide")
elif [ $sample = "AB_LT" ]
then	site=$(echo "Assistance Bay - Low tide")
elif [ $sample = "Ale_1" ]
then	site=$(echo "Alert - 1")
elif [ $sample = "Cam_Cam" ]
then	site=$(echo "Cambridge Bay")
elif [ $sample = "Dump_2018" ]
then	site=$(echo "Dump Beach - 2018")
elif [ $sample = "Dump_2019" ]
then	site=$(echo "Dump Beach - 2019")
elif [ $sample = "Dyna_2018" ]
then	site=$(echo "Dynamite Beach - 2018")
elif [ $sample = "Dyna_2019" ]
then	site=$(echo "Dynamite Beach - 2019")
elif [ $sample = "Nani_A1" ]
then	site=$(echo "Nanisivik - A1")
elif [ $sample = "Nani_B1" ]
then	site=$(echo "Nanisivik - B1")
elif [ $sample = "Tank_2018" ]
then	site=$(echo "Tank Farm - 2018")
elif [ $sample = "Tank_2019" ]
then	site=$(echo "Tank Farm - 2019")
elif [ $sample = "Tupi_2018" ]
then	site=$(echo "Tupirvik - 2018")
elif [ $sample = "Tupi_2019" ]
then	site=$(echo "Tupirvik - 2019")
else
	site=$(echo "remove")
fi
taxonomy=$(echo $line | cut -f 8 -d ' ')
coverage=$(echo $line | grep -Eo 'cov_[0-9]+\.[0-9]+' | cut -f 2 -d '_')

echo -e $sample'\t'$beach'\t'$site'\t'$coverage'\t'$taxonomy >> results/16S_coassembly_taxonomy.txt

done

sed -i 's/\;/\t/g' results/16S_coassembly_taxonomy.txt
