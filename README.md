# NWP_beach_metagenomes
 *This repository was last updated on July 12, 2024*
 
 This repository provides the accompanying scripts for the manuscript "Metagenomic survey reveals hydrocarbon biodegradation potential of Canadian high Arctic beaches" (Góngora et al., 2024) which is currently under review in the journal *Environmental Microbiome*.

- `baseline_metagenomes.txt` is the main script used to analyze the metagenomes with various Linux-based bioinformatic pipelines
  - The folder `sub-scripts` contains the sub-scripts that automate certain parts of the `baseline_metagenomes.txt` script
- `baseline.R` is the R script used to perform the statistical analyses and produce the manuscript figures after obtaining the various outputs from `baseline_metagenomes.txt` and its sub-scripts
- The folder `examples` contains example files that need to be manually created for certain parts of the scripts to work

The sequencing data required to reproduce the study was deposited in the NCBI Sequence Read Archive (SRA) repository under the BioProject accession number PRJNA1046404.