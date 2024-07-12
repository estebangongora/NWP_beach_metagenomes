rm(list = ls())

#Load required libraries
library(tidyverse)
library(phyloseq)
library(vegan)

######################## 16S analysis for eDNA vs iDNA #############################################

#General data preparation
##Load data frame
taxa_16s <- read.csv("16S_taxonomy.txt", sep = "\t", header = T, stringsAsFactors = T,
                     na.strings = "")

##Remove samples not used in e/iDNA analysis because there is no 16S data in their counterpart
taxa_16s <- subset(taxa_16s, Sample!="Dump_2019_e" & Sample!="AB_HT_i")
taxa_16s$Sample <- taxa_16s$Sample[ , drop=TRUE]

##Remove chloroplast and mitochondria entries
chloroplast_ei <- which(taxa_16s$Order == "Chloroplast")
###Check that the chloroplast_ei object is not empty before running the next line
####If it is empty, don't run the next line or the whole data frame will be wiped
taxa_16s <- taxa_16s[-chloroplast_ei, ]
mitochondria_ei <- which(taxa_16s$Family == "Mitochondria")
###Check that the mitochondria_ei object is not empty before running the next line
####If it is empty, don't run the next line or the whole data frame will be wiped
####Because we did not find any mitochondria hits, there is no need to run the next line
#taxa_16s <- taxa_16s[-mitochondria_ei, ]

##Change entries with an empty Phylum into an "Other" level
taxa_16s$Phylum <- as.character(taxa_16s$Phylum)
taxa_16s["Phylum"][is.na(taxa_16s["Phylum"])] <- "Other"
phyla <- taxa_16s$Phylum[!duplicated(taxa_16s$Phylum)]
phyla <- phyla[phyla != "Other"]
phyla <- c(phyla, "Other")
taxa_16s$Phylum <- factor(taxa_16s$Phylum, levels = phyla)

##Merge Bacteroidetes and Bacteroidota into Bacteroidota
levels(taxa_16s$Phylum)[levels(taxa_16s$Phylum)=="Bacteroidetes"] <- "Bacteroidota"

##Change Actinobacteria to Actinomycetota
levels(taxa_16s$Phylum)[levels(taxa_16s$Phylum)=="Actinobacteria"] <- "Actinomycetota"


#Barplots for e/iDNA data

##Concatenate repeated species entries by sample
taxa_conc_ieDNA <- taxa_16s %>%
  group_by(Sample, Beach, Site, DNA, Phylum) %>%
  summarise(Coverage = sum(Coverage))
taxa_conc_ieDNA <- taxa_conc_ieDNA %>% group_by(Sample) %>%
  mutate(Relative_coverage = ((Coverage / sum(Coverage))*100))

##Grouped taxa with relative abundance lower than 5% into an "Other" category
include16s <- subset(taxa_conc_ieDNA, Relative_coverage > 5)

exclude16s <- taxa_conc_ieDNA %>% anti_join(include16s)
phylum_rows = nrow(exclude16s)
for (i in 1:phylum_rows){
  exclude16s$Phylum = "Other"
}
exclude16s$Phylum <- as.factor(exclude16s$Phylum)

ssu <- rbind(include16s, exclude16s)

##Changed the labeling for the beaches and and rearranged them so that the Resolute Bay beaches are next to each other
levels(ssu$Beach) <- c("Assistance Bay", "Alert", "Cambridge Bay", "Dump beach", "Dynamite beach",
                       "Nanisivik", "Tank farm", "Tupirvik")
ssu$Beach <- factor(ssu$Beach, levels = c("Assistance Bay", "Dump beach", "Dynamite beach",
                                          "Tank farm", "Tupirvik", "Alert", "Cambridge Bay",
                                          "Nanisivik"))

##Removed the unused levels for the pooled phyla
###Also renamed the phyla to the current nomenclature and rearranged phyla alphabetically
ssu$Phylum <- droplevels(ssu$Phylum)
levels(ssu$Phylum) <- c("Pseudomonadota", "Bacteroidota", "Actinomycetota", "Planctomycetota",
                        "Acidobacteriota", "Other")
ssu$Phylum <- factor(ssu$Phylum, levels = c("Acidobacteriota", "Actinomycetota",
                                            "Bacteroidota", "Planctomycetota",
                                            "Pseudomonadota","Other"))

##Barplots for i/eDNA data (Fig. S1)
ieDNA_barplots <- ggplot(data = ssu, aes(x = Sample, y = Coverage, fill = Phylum)) +
  geom_bar(stat = "identity",
           position = "fill") +
  labs(fill = "Phylum") +
  facet_wrap(~ Beach,
             scales = "free_x",
             nrow = 1,
             strip.position = "bottom") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(vjust=3, color = "black", size = 14),
        axis.title.x = element_text(vjust=1, color = "black", size = 14),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15)) +
  scale_y_continuous(labels=scales::percent) +
  scale_x_discrete(labels=c("AB_LT_e" = "Intertidal - 2019 - eDNA", "AB_LT_i" = "Intertidal - 2019 - iDNA",
                            "Ale_1_e" = "2018 - eDNA", "Ale_1_i" = "2018 - iDNA",
                            "Cam_Cam_e" = "eDNA", "Cam_Cam_i" = "iDNA",  
                            "Dump_2018_e" = "2018 - eDNA", "Dump_2018_i" = "2018 - iDNA",
                            "Dyna_2018_e" = "2018 - eDNA", "Dyna_2018_i" = "2018 - iDNA",
                            "Dyna_2019_e" = "2019 - eDNA", "Dyna_2019_i" = "2019 - iDNA",
                            "Nani_A1_e" = "East - eDNA", "Nani_A1_i" = "East - iDNA",
                            "Nani_B1_e" = "West - eDNA", "Nani_B1_i" = "West - iDNA",
                            "Tank_2018_e" = "2018 - eDNA", "Tank_2018_i" = "2018 - iDNA",
                            "Tank_2019_e" = "2019 - eDNA", "Tank_2019_i" = "2019 - iDNA",
                            "Tupi_2018_e" = "2018 - eDNA", "Tupi_2018_i" = "2018 - iDNA",
                            "Tupi_2019_e" = "2019 - eDNA", "Tupi_2019_i" = "2019 - iDNA")) +
  ylab("Relative abudance") +
  scale_fill_viridis_d(direction = 1)
ieDNA_barplots


#Statistical comparison of eDNA vs iDNA

##Concatenate repeated species entries by sample
taxa_conc <- taxa_16s %>%
  group_by(Sample, Beach, Site, DNA, Domain, Phylum, Class, Order, Family, Genus) %>%
  summarise(Coverage = sum(Coverage))

##Add "OTU" IDs for unique taxonomic entries
taxa_conc <- taxa_conc %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus) %>%
  mutate(OTU = cur_group_id())
taxa_conc$OTU <- paste0("OTU", taxa_conc$OTU)

##Create "OTU" table as matrix for phyloseq
OTU_table <- subset(taxa_conc, select = c(Sample, OTU, Coverage))
OTU_table <- spread(OTU_table, Sample, Coverage, fill = 0)
OTU_matrix <- subset(OTU_table, select = - OTU)
OTU_matrix <- data.matrix(OTU_matrix)
rownames(OTU_matrix) <- OTU_table$OTU

##Create taxonomy table as matrix for phyloseq
taxa_table <- subset(taxa_conc, select = c(OTU, Domain, Phylum, Class, Order, Family, Genus))
taxa_table <- distinct(taxa_table, OTU)
taxa_matrix <- subset(taxa_table, select = -OTU)
taxa_matrix <- as.matrix(taxa_matrix)
rownames(taxa_matrix) <- taxa_table$OTU

##Create metadata table
metadata <- subset(taxa_16s, select = c('Beach', 'Site', 'Site_long','DNA'))
metadata <- distinct(metadata)
rownames(metadata) <- levels(taxa_16s$Sample)

##Create phyloseq object
OTU <- otu_table(OTU_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_matrix)
MET <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, MET)

##Alpha diversity
###Calculate Shannon diversity
shannon <- estimate_richness(ps, split = TRUE, measures = "Shannon")
shannon$DNA <- c(rep(c("eDNA", "iDNA"), 12))

###Test for differences in alpha diversity between DNA types
eDNA <- subset(shannon,  DNA == "eDNA", Shannon, drop = TRUE)
iDNA <- subset(shannon,  DNA == "iDNA", Shannon, drop = TRUE)

differences <- with(shannon, Shannon[DNA == "eDNA"] - Shannon[DNA == "iDNA"])

shapiro.test(differences)
t.test(eDNA, iDNA, paired = TRUE)

###Plot Shannon diversity by DNA type (Fig. S2)
boxplot(Shannon~DNA, data = shannon, xlab = "DNA type", ylab = "Shannon Diversity")

##Test differences in Bray-Curtis dissimilarities between DNA types

###Run PERMDISP and PERMANOVA comparing eDNA vs iDNA
####Convert to relative abundances
ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))
bray_dist <- distance(ps_ra, method = "bray")

####PERMDISP
dispersions <- betadisper(bray_dist, sample_data(ps_ra)$DNA)
permutest(dispersions)

####PERMANOVA
adonis2(bray_dist ~ sample_data(ps_ra)$DNA)

####Plot ordination (Fig. S3)
bray <- ordinate(ps_ra, "NMDS", "bray")

bray_plot <- plot_ordination(ps_ra, bray, type = "sample", color = "Site_long", shape = "DNA")
bray_plot +
  geom_point(size = 5) +
  labs(shape = "DNA Type", colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        strip.background = element_blank())

######################## 16S analysis for co-assembly ##############################################

#General data preparation
##Load data frame
taxa_16s_coass <- read.csv("16S_coassembly_taxonomy.txt", sep = "\t", header = T,
                           stringsAsFactors = T, na.strings = "")

##Remove chloroplast and mitochondria entries
chloroplast <- which(taxa_16s_coass$Order == "Chloroplast")
###Check that the chloroplast object is not empty before running the next line
####If it is empty, don't run the next line or the whole data frame will be wiped
taxa_16s_coass <- taxa_16s_coass[-chloroplast, ]
mitochondria <- which(taxa_16s_coass$Family == "Mitochondria")
###Check that the mitochondria object is not empty before running the next line
####If it is empty, don't run the next line or the whole data frame will be wiped
taxa_16s_coass <- taxa_16s_coass[-mitochondria, ]

##Change entries with an empty Phylum into an "Other" level
taxa_16s_coass$Phylum <- as.character(taxa_16s_coass$Phylum)
taxa_16s_coass["Phylum"][is.na(taxa_16s_coass["Phylum"])] <- "Other"
phyla_coass <- taxa_16s_coass$Phylum[!duplicated(taxa_16s_coass$Phylum)]
phyla_coass <- phyla_coass[phyla_coass != "Other"]
phyla_coass <- c(phyla_coass, "Other")
taxa_16s_coass$Phylum <- factor(taxa_16s_coass$Phylum, levels = phyla_coass)

###Merge Bacteroidetes and Bacteroidota into Bacteroidota
levels(taxa_16s_coass$Phylum)[levels(taxa_16s_coass$Phylum)=="Bacteroidetes"] <- "Bacteroidota"


#Barplots for co-assembly data
##Concatenate repeated species entries by sample
taxa_coass_conc <- taxa_16s_coass %>%
  group_by(Sample, Beach, Site, Phylum, Class) %>%
  summarise(Coverage = sum(Coverage))
taxa_coass_conc <- taxa_coass_conc %>% group_by(Sample) %>%
  mutate(Relative_coverage = ((Coverage / sum(Coverage))*100))

##Grouped taxa with relative abundance lower than 5% into an "Other" category
include16s_coass <- subset(taxa_coass_conc, Relative_coverage > 2.5)
exclude16s_coass <- taxa_coass_conc %>% anti_join(include16s_coass)

include16s_coass[, 'Taxa'] <- NA
include16s_coass$Taxa <- paste(include16s_coass$Phylum, include16s_coass$Class, sep = ", ")
include16s_coass <- include16s_coass[c(1:5, 8, 6:7)]

exclude16s_coass[, 'Taxa'] <- NA
exclude16s_coass <- exclude16s_coass[c(1:5, 8, 6:7)]
phylum_rows_coass = nrow(exclude16s_coass)
for (i in 1:phylum_rows_coass){
  exclude16s_coass$Phylum = "Other"
  exclude16s_coass$Class = "Other"
  exclude16s_coass$Taxa = "Other"
  }

ssu_coass <- rbind(include16s_coass, exclude16s_coass)

##Changed the labeling for the beaches and and rearranged them so that the Resolute Bay beaches are next to each other
levels(ssu_coass$Beach) <- c("AB", "Ale", "CB", "Dump", "Dynamite",
                       "Nanisivik", "Tank", "Tupirvik")
ssu_coass$Beach <- factor(ssu_coass$Beach, levels = c("AB", "Dump", "Dynamite",
                                          "Tank", "Tupirvik", "Ale", "CB",
                                          "Nanisivik"))

##Removed the unused levels for the pooled taxa
ssu_coass$Taxa <- as.factor(ssu_coass$Taxa)
ssu_coass$Taxa <- droplevels(ssu_coass$Taxa)
levels(ssu_coass$Taxa) <- c("Acidobacteriota, Holophagae", "Actinomycetota, Acidimicrobiia",
                            "Actinomycetota, Actinomycetia", "Bacteroidota, Bacteroidia",
                            "Thermodesulfobacteriota, Desulfobulbia", "Other",
                            "Pseudomonadota, Alphaproteobacteria", "Pseudomonadota, Gammaproteobacteria",
                            "Verrucomicrobiota, Verrucomicrobiia")
ssu_coass$Taxa <- factor(ssu_coass$Taxa, levels = c("Acidobacteriota, Holophagae",
                                               "Actinomycetota, Acidimicrobiia",
                                               "Actinomycetota, Actinomycetia",
                                               "Bacteroidota, Bacteroidia",
                                               "Thermodesulfobacteriota, Desulfobulbia",
                                               "Pseudomonadota, Alphaproteobacteria",
                                               "Pseudomonadota, Gammaproteobacteria",
                                               "Verrucomicrobiota, Verrucomicrobiia",
                                               "Other"))

##Barplots by sample (Fig. 2)
barplots_coass <- ggplot(data = ssu_coass, aes(x = Sample, y = Coverage, fill = Taxa)) +
  geom_bar(stat = "identity",
           position = "fill") +
  labs(fill = "Taxonomy") +
  facet_wrap(~ Beach,
             scales = "free_x",
             nrow = 1,
             strip.position = "bottom") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(vjust=1.5, color = "black", size = 14),
        axis.title.x = element_text(vjust=1, color = "black", size = 14),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, vjust = 2)) +
  scale_y_continuous(labels=scales::percent) +
  scale_x_discrete(labels=c("AB_HT" = "Supratidal - 2019", "AB_LT" = "Intertidal - 2019",
                            "Ale_1" = "2018", "Cam_Cam" = "2018",
                            "Dump_2018" = "2018", "Dump_2019" = "2019",
                            "Dyna_2018" = "2018", "Dyna_2019" = "2019",
                            "Nani_A1" = "East - 2018", "Nani_B1" = "West - 2018",
                            "Tank_2018" = "2018", "Tank_2019" = "2019",
                            "Tupi_2018" = "2018", "Tupi_2019" = "2019")) +
  ylab("Relative abundance") +
  scale_fill_viridis_d(option = "viridis", direction = 1)
barplots_coass


#For taxonomy based on MetaErg annotation
##Load dataframe
taxa_metaerg <- read.csv("metaerg_cds_taxonomy.tsv", sep = "\t", header = T, stringsAsFactors = T,
                         na.strings = "")

##Change entries with an empty Phylum into an "Other" level
taxa_metaerg$Phylum <- as.character(taxa_metaerg$Phylum)
taxa_metaerg["Phylum"][is.na(taxa_metaerg["Phylum"])] <- "Other"
phyla_metaerg <- taxa_metaerg$Phylum[!duplicated(taxa_metaerg$Phylum)]
phyla_metaerg <- phyla_metaerg[phyla_metaerg != "Other"]
phyla_metaerg <- c(phyla_metaerg, "Other")
taxa_metaerg$Phylum <- factor(taxa_metaerg$Phylum, levels = phyla_metaerg)

##Barplots for metaerg data
###Concatenate repeated species entries by sample
taxa_metaerg <- taxa_metaerg %>%
  group_by(Sample, Beach, Site, Phylum, Class) %>%
  summarise(Relative_abundance = sum(Relative_abundance))


###Grouped taxa with relative abundance lower than 5% into an "Other" category
include_metaerg <- subset(taxa_metaerg, Relative_abundance > 2.5)
exclude_metaerg <- taxa_metaerg %>% anti_join(include_metaerg)

include_metaerg[, 'Taxa'] <- NA
include_metaerg$Taxa <- paste(include_metaerg$Phylum, include_metaerg$Class, sep = ", ")
include_metaerg <- include_metaerg[c(1:5, 7, 6)]

exclude_metaerg[, 'Taxa'] <- NA
exclude_metaerg <- exclude_metaerg[c(1:5, 7, 6)]
phylum_rows_metaerg = nrow(exclude_metaerg)
for (i in 1:phylum_rows_coass){
  exclude_metaerg$Phylum = "Other"
  exclude_metaerg$Class = "Other"
  exclude_metaerg$Taxa = "Other"
}
exclude_metaerg <- exclude_metaerg %>%
  group_by(Sample, Beach, Site, Phylum, Class, Taxa) %>%
  summarise(Relative_abundance = sum(Relative_abundance))

ssu_metaerg <- rbind(include_metaerg, exclude_metaerg)

##Changed the labeling for the beaches and and rearranged them so that the Resolute Bay beaches are next to each other
levels(ssu_metaerg$Beach) <- c("AB", "Ale", "CB", "Dump", "Dynamite",
                             "Nanisivik", "Tank", "Tupirvik")
ssu_metaerg$Beach <- factor(ssu_metaerg$Beach, levels = c("AB", "Dump", "Dynamite",
                                                      "Tank", "Tupirvik", "Ale", "CB",
                                                      "Nanisivik"))

##Removed the unused levels for the pooled taxa
ssu_metaerg$Taxa <- as.factor(ssu_metaerg$Taxa)
ssu_metaerg$Taxa <- droplevels(ssu_metaerg$Taxa)
levels(ssu_metaerg$Taxa) <- c("Acidobacteriota, Thermoanaerobaculia",
                              "Acidobacteriota, Vicinamibacteria",
                              "Actinomycetota, Acidimicrobiia",
                              "Actinomycetota, Actinomycetia",
                              "Bacteroidota, Bacteroidia",
                              "Gemmatimonadota, Gemmatimonadia",
                              "Myxococcota, Polyangiia", "Other",
                              "Planctomycetota, Planctomycetia",
                              "Pseudomonadota, Alphaproteobacteria",
                              "Pseudomonadota, Gammaproteobacteria",
                              "Unclassified")
ssu_metaerg$Taxa <- factor(ssu_metaerg$Taxa,
                           levels = c("Acidobacteriota, Thermoanaerobaculia",
                                      "Acidobacteriota, Vicinamibacteria",
                                      "Actinomycetota, Acidimicrobiia",
                                      "Actinomycetota, Actinomycetia",
                                      "Bacteroidota, Bacteroidia",
                                      "Gemmatimonadota, Gemmatimonadia",
                                      "Myxococcota, Polyangiia",
                                      "Planctomycetota, Planctomycetia",
                                      "Pseudomonadota, Alphaproteobacteria",
                                      "Pseudomonadota, Gammaproteobacteria",
                                      "Other", "Unclassified"))

###Barplots by sample using MetaErg taxonomy (Fig. S4)
barplots_metaerg <- ggplot(data = ssu_metaerg, aes(x = Sample,
                                                   y = Relative_abundance,
                                                   fill = Taxa)) +
  geom_bar(stat = "identity",
           position = "fill") +
  labs(fill = "Taxonomy") +
  facet_wrap(~ Beach,
             scales = "free_x",
             nrow = 1,
             strip.position = "bottom") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(vjust=1.5, color = "black", size = 14),
        axis.title.x = element_text(vjust=1, color = "black", size = 14),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, vjust = 2)) +
  scale_y_continuous(labels=scales::percent) +
  scale_x_discrete(labels=c("AB_HT" = "Supratidal - 2019", "AB_LT" = "Intertidal - 2019",
                            "Ale_1" = "2018", "Cam_Cam" = "2018",
                            "Dump_2018" = "2018", "Dump_2019" = "2019",
                            "Dyna_2018" = "2018", "Dyna_2019" = "2019",
                            "Nani_A1" = "East - 2018", "Nani_B1" = "West - 2018",
                            "Tank_2018" = "2018", "Tank_2019" = "2019",
                            "Tupi_2018" = "2018", "Tupi_2019" = "2019")) +
  ylab("Relative abundance") +
  scale_fill_viridis_d(option = "viridis", direction = 1)
barplots_metaerg


#Bubble plots for top 20 genera
##Concatenate repeated species entries by sample
taxa_bubble_coass <- taxa_16s_coass %>%
  group_by(Sample, Beach, Beach_long, Site, Region, Year, Domain, Phylum, Class, Order, Family, Genus) %>%
  summarise(Coverage = sum(Coverage))

##Add relative abundance column
taxa_top20 <- taxa_bubble_coass %>% group_by(Sample) %>% mutate(Relative_coverage = ((Coverage / sum(Coverage))*100))

##Estimate overall average abundances of all genera and select top 20
coverage <- taxa_top20 %>% group_by(Genus) %>% summarise(average_coverage = mean(Relative_coverage), sd_coverage = sd(Relative_coverage))
top_20_genera <- top_n(coverage, 20, average_coverage)
for (i in 1:20){
  top_20_genera$samples[i] = sum(taxa_top20$Genus == top_20_genera$Genus[i], na.rm = TRUE)
}
###Table S3
write.csv(top_20_genera, "top20.csv", row.names = FALSE)

top_20 <- subset(taxa_top20, Genus %in% levels(droplevels(top_20_genera$Genus)))
top_20$Genus <- droplevels(top_20$Genus)
top_20$Phylum <- droplevels(top_20$Phylum)

##Changed the labeling for the beaches and and rearranged them so that the Resolute Bay beaches are next to each other
levels(top_20$Beach) <- c("AB", "Alert", "CB", "Dump", "Dynamite",
                          "Nanisivik", "Tank", "Tupirvik")
top_20$Beach <- factor(top_20$Beach, levels = c("AB", "Dump", "Dynamite",
                                                "Tank", "Tupirvik", "Alert", "CB",
                                                "Nanisivik"))

##Fixed taxonomy and rearranged the phyla so that they are in alphabetical order
levels(top_20$Phylum) <- c("Actinomycetota", "Pseudomonadota", "Bacteroidota")
top_20$Phylum <- factor(top_20$Phylum, levels = c("Actinomycetota", "Bacteroidota",
                                                  "Pseudomonadota"))

##Create bubble plot (Fig. S5)
bubble_plot_top20 = ggplot(top_20, aes(x=Sample, y=Genus)) +
  geom_point(aes(size = Relative_coverage, fill = Phylum), alpha = 1, shape = 21) +
  scale_size_continuous(limits = c(0,72),
                        range = c(1,10),
                        breaks = c(5,15,25,50),
                        name = "Relative abundance (%)") +
  facet_grid( Phylum ~ Beach,
              scales = "free",
              space = "free",
              switch = "x") +
  theme_bw() +
  theme(plot.margin = margin(t = 15,
                             r = 10,
                             b = 15,
                             l = 15),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 14),
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 16),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_rect(fill = NA, colour = "lightgray"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_blank(),
        panel.border = element_rect(colour = "lightgray"),
        axis.title.y = element_text(vjust=3.5, color = "black", size = 16),
        axis.title.x = element_text(vjust=1, color = "black", size = 16)) +
  scale_x_discrete(labels=c("AB_HT" = "Supratidal - 2019", "AB_LT" = "Intertidal - 2019",
                            "Ale_1" = "2018", "Cam_Cam" = "2018",
                            "Dump_2018" = "2018", "Dump_2019" = "2019",
                            "Dyna_2018" = "2018", "Dyna_2019" = "2019",
                            "Nani_A1" = "East - 2018", "Nani_B1" = "West - 2018",
                            "Tank_2018" = "2018", "Tank_2019" = "2019",
                            "Tupi_2018" = "2018", "Tupi_2019" = "2019")) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("#440154", "#fde725", "#21918c")) +
  guides(fill = guide_legend(override.aes = list(size=6, alpha = 1)))
bubble_plot_top20


#Statistical comparison of beaches

##Add "OTU" IDs for unique taxonomic entries
taxa_conc_coass <- taxa_bubble_coass %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus) %>%
  mutate(OTU = cur_group_id())
taxa_conc_coass$OTU <- paste0("OTU", taxa_conc_coass$OTU)

##Create "OTU" table as matrix for phyloseq
OTU_table_coass <- subset(taxa_conc_coass, select = c(Sample, OTU, Coverage))
OTU_table_coass <- spread(OTU_table_coass, Sample, Coverage, fill = 0)
OTU_matrix_coass <- subset(OTU_table_coass, select = - OTU)
OTU_matrix_coass <- data.matrix(OTU_matrix_coass)
rownames(OTU_matrix_coass) <- OTU_table_coass$OTU

##Create taxonomy table as matrix for phyloseq
taxa_table_coass <- subset(taxa_conc_coass, select = c(OTU, Domain, Phylum, Class, Order, Family, Genus))
taxa_table_coass <- distinct(taxa_table_coass, OTU)
taxa_matrix_coass <- subset(taxa_table_coass, select = -OTU)
taxa_matrix_coass <- as.matrix(taxa_matrix_coass)
rownames(taxa_matrix_coass) <- taxa_table_coass$OTU

##Create metadata table
metadata_coass <- subset(taxa_16s_coass, select = c('Beach', 'Beach_long', 'Site', 'Region', 'Year'))
metadata_coass <- distinct(metadata_coass)
rownames(metadata_coass) <- levels(taxa_16s_coass$Sample)
metadata_coass$Year <- as.factor(metadata_coass$Year)
metadata_coass$PAH <- c("Absent", "Absent", "Present", "Present", "Absent", "Absent",
                        "Absent", "Present", "Present", "Absent", "Present", "Absent",
                        "Absent", "Absent")

##Create phyloseq object
OTU_coass <- otu_table(OTU_matrix_coass, taxa_are_rows = TRUE)
TAX_coass <- tax_table(taxa_matrix_coass)
MET_coass <- sample_data(metadata_coass)
ps_coass <- phyloseq(OTU_coass, TAX_coass, MET_coass)

##Alpha diversity
###Calculate Shannon diversity
shannon_coass <- estimate_richness(ps_coass, split = TRUE, measures = "Shannon")
shannon_coass$Group <- c("HT", "LT", "drop", "drop", 2018, 2019, 2018,
                         2019, "A1", "B1", 2018, 2019, 2018, 2019)
shannon_coass$Region <- c("Resolute", "Resolute", "Alert", "Cambridge",
                          "Resolute", "Resolute", "Resolute", "Resolute",
                          "Nanisivik", "Nanisivik", "Resolute", "Resolute",
                          "Resolute", "Resolute")

###Test for differences in alpha diversity between years
y2018 <- subset(shannon_coass,  Group == 2018, Shannon, drop = TRUE)
y2019 <- subset(shannon_coass,  Group == 2019, Shannon, drop = TRUE)

differences_coass_year <- with(shannon_coass, Shannon[Group == 2018] - Shannon[Group == 2019])

shapiro.test(differences_coass_year)
wilcox.test(y2018, y2019, paired = TRUE)

####Plot Shannon diversity by year (Fig. S6)
years <- subset(shannon_coass, Group == 2018 | Group == 2019)
boxplot(Shannon~Group, data = years, xlab = "Group", ylab = "Shannon Diversity")

###Test for differences in alpha diversity between regions
shapiro.test(shannon_coass$Shannon)
lm_region=lm(data=shannon_coass, Shannon~Region)
car::ncvTest(lm_region)
model_region=aov(data=shannon_coass,Shannon~Region)
summary(model_region)

###Plot Shannon diversity by region (Fig. S8)
boxplot(Shannon~Region, data = shannon_coass, xlab = "Region", ylab = "Shannon Diversity")

##Test differences in Bray-Curtis dissimilarities between years and regions
###Run PERMDISP and PERMANOVA comparing years
ps_years <- subset_samples(ps_coass, Year!=0)

####Convert to relative abundances and calculate dissimilarities
ps_years_ra <- transform_sample_counts(ps_years, function(x) x/sum(x))
bray_dist_years <- distance(ps_years_ra, method = "bray")

####PERMDISP
dispersions_years <- betadisper(bray_dist_years, sample_data(ps_years_ra)$Year)
permutest(dispersions_years)

####PERMANOVA
adonis2(bray_dist_years ~ sample_data(ps_years_ra)$Year)

###Check that there are no differences within compared to between years
####Subset samples with paired years
otus <- OTU_matrix_coass[,c(5:8,11:14),drop=FALSE]

####Remove any unused OTUs
otus <- otus[rowSums(otus[])>0,] 
meta <- ps_years@sam_data

####Calculate distances
braycurtis_years <- vegdist(otus)

####Next line will output distances for all sample pairs. Distances for all 2018, 2019,
####and 2018-2019 pairs need to be extracted and manually sorted into year_diff.csv
meandist(braycurtis_years, meta$Site)

####Test for difference among sample pairs
year_diff <- read.csv("year_diff.csv", header = T, stringsAsFactors = T, na.strings = "")
shapiro.test(year_diff$Bray)
kruskal.test(data=year_diff, Bray~Year)

####Plot ordination (Fig. S7)
bray_years <- ordinate(ps_years_ra, "NMDS", "bray")

bray_plot_years <- plot_ordination(ps_years_ra, bray_years, type = "samples", color = "Beach_long", shape = "Year", label = "Site")
bray_plot_years +
  geom_point(size = 5) +
  labs(shape = "Year", colour = "Beach") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        strip.background = element_blank())

###Run PERMDISP and PERMANOVA comparing regions

####Convert to relative abundances
ps_coass_ra <- transform_sample_counts(ps_coass, function(x) x/sum(x))
bray_dist_coass <- distance(ps_coass_ra, method = "bray")

####PERMDISP
dispersions_coass <- betadisper(bray_dist_coass, sample_data(ps_coass_ra)$Region)
permutest(dispersions_coass)
#PERMANOVA will be overly conservative because the dispersion of Resolute is larger

####PERMANOVA
adonis2(bray_dist_coass ~ sample_data(ps_coass_ra)$Region)
pairwiseAdonis::pairwise.adonis(t(otu_table(ps_coass_ra)), metadata_coass$Region)

####Plot ordination separating samples by beach and region (Fig. S9)
bray_coass <- ordinate(ps_coass_ra, "NMDS", "bray")

bray_plot_coass <- plot_ordination(ps_coass_ra, bray_coass, type = "samples", color = "Beach_long", shape = "Region", label = "Site")
bray_plot_coass +
  geom_point(size = 5) +
  labs(shape = "Region", colour = "Beach") +
  scale_shape_manual(values = c(18, 17, 16, 15)) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        strip.background = element_blank())

###Run PERMDISP and PERMANOVA comparing samples with/without PAH genes

####PERMDISP
dispersions_coass_PAH <- betadisper(bray_dist_coass, sample_data(ps_coass_ra)$PAH)
permutest(dispersions_coass_PAH)
#Multivariate dispersions are homogeneous

####PERMANOVA
adonis2(bray_dist_coass ~ sample_data(ps_coass_ra)$PAH)


#Stacked barplots for novelty by taxonomic rank (Fig. 3b)
##The novelty.csv file was created using the GTDB-Tk information
novelty <- read.csv("novelty.csv", header = T, stringsAsFactors = T, na.strings = "")
novelty$Taxa <- factor(novelty$Taxa, levels = c("Species", "Genus", "Family", "Order"))

novelty_bars <- ggplot(data = novelty, aes(x = Percentage, y = Taxa, fill = Class)) +
  geom_bar(stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, colour = "black", size = 12),
        axis.text.y = element_text(hjust = 1.5, colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(vjust=1, color = "black", size = 12)) +
  scale_fill_manual(values = c("#1a8872", "white"))
novelty_bars

######################## Degradation genes by metagenome (CANT-HYD) ################################

#General data preparation

##Load data frame
###The cant-hyd_noise_metagenome_results.txt file needs to be manually edited to remove all entries with 0 counts
###Hits per million (HPM) must also be calculated
###It then must be saved as hpm_canthyd_metagenomes_noise.csv
###The cant-hyd_noise_metagenome_results.xlsx file in the GitHub repository gives an example of how this can be done
beaches <- read.csv("hpm_canthyd_metagenomes_noise.csv", header = T, stringsAsFactors = T, na.strings = "")

##Add data column with log-converted HPM values to make figure more readable
beaches$logHPM <- log10(beaches$HPM)

##Rearranged the Site levels so that the Resolute Bay beaches are next to each other and renamed them to make figure smaller
levels(beaches$Site) <- c("Alert", "AB", "CB", "Dump", "Dynamite", "Nanisivik", "Tank", "Tupirvik")
beaches$Site <- factor(beaches$Site, levels = c("AB", "Dump", "Dynamite", "Tank", "Tupirvik",
                                                "Alert", "CB", "Nanisivik"))

##Rearranged the Compound levels so that propane and toluene are next to alkane and MAH, respectively
beaches$Compound <- factor(beaches$Compound, levels = c("Alkane", "EB", "MAH", "Toluene", "PAH"))

##Renamed beaches to make figure more readable
beaches$Beach_name <- droplevels(beaches$Beach_name)
levels(beaches$Beach_name) <- c("2018", "2019", "East - 2018", "2018", "West - 2018", "2018",
                                "Supratidal - 2019", "Intertidal - 2019")

#Heatmap of hydrocarbon degradation gene counts (HPM) by sample (Fig. 4)
heatmap_log_HPM <- ggplot(beaches, aes(Beach_name, Gene, fill=logHPM)) +
  geom_tile() +
  labs(fill = "Hits per million genes (log)") +
  facet_grid(Pathway + Compound ~ Site,
             scales = "free",
             space = "free_y",
             switch = "x") +
  theme_bw() +
  theme(plot.margin = margin(t = 8,
                             r = 10,
                             b = -15,
                             l = 15),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "points"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        legend.text = element_text(vjust = 14.5, colour = "black", size = 14),
        legend.title = element_text(vjust = 1, hjust = 1, colour = "black", size = 16),
        legend.key.width= unit(1.5, 'cm'),
        legend.position = "bottom",
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background.x = element_rect(fill = NA, colour = "lightgray"),
        strip.background.y = element_rect(fill = NA, colour = "lightgray"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(angle = 0, size = 14),
        panel.border = element_rect(colour = "lightgray"),
        axis.title.y = element_text(vjust=3, colour = "black", size = 18),
        axis.title.x = element_blank()) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  scale_y_discrete(limits=rev)
heatmap_log_HPM

######################## Degradation genes by MAG (CANT-HYD) #######################################

#General data preparation

##Load data frame
###The cant-hyd_noise_bin_results.txt file needs to be manually edited to remove all genes with 0 counts
###Hits per million (HPM) must also be calculated
###It then must be saved as genes_canthyd_bins_noise.csv
###The genes_canthyd_bins_noise.xlsx file in the GitHub repository gives an example of how this can be done
bins <- read.csv("genes_canthyd_bins_noise.csv", header = T, stringsAsFactors = T, na.strings = "")

##Add data column with log-converted CPM values to make figure more readable
bins$logHPM <- log10(bins$HPM)

##Rearranged the Site levels so that the Resolute Bay beaches are next to each other
levels(bins$Site) <- c("Alert", "AB", "Cambridge Bay", "DB", "Dynamite Beach",
                       "Nanisivik", "Tank Farm", "Tupirvik")
bins$Site <- factor(bins$Site, levels = c("AB", "DB", "Dynamite Beach",
                                                "Tank Farm", "Tupirvik", "Alert", "Cambridge Bay",
                                                "Nanisivik"))

##Rearranged the Compound levels so that toluene is next to MAH
bins$Compound <- factor(bins$Compound, levels = c("Alkane", "EB", "MAH", "Toluene", "PAH"))


##Rearranged the Class levels so that the Unclassified samples are at the end
bins$Class <- factor(bins$Class, levels = c("Acidimicrobiia", "Actinomycetia",
                                            "Alphaproteobacteria",
                                            "Gammaproteobacteria",
                                            "Xanthomonadales", "Unclassified"))

##Rearranged the Gene levels so that they are ordered alphabetically
bins$Gene <- factor(bins$Gene, levels = c("LadA_alpha", "PrmC", "PrmA", "CYP153",
                                          "AlmA_GroupI", "AlkB", "MAH_beta",
                                          "MAH_alpha", "NdoC", "non_NdoB_type",
                                          "NdoB", "TmoE", "TmoA_BmoA", "AhyA",
                                          "EbdA"))

#Bubble plot (Fig. 5)
bubble_plot_bins_HPM = ggplot(bins, aes(x=MAG_ID, y=Gene)) +
  geom_point(aes(size = logHPM, fill = Class), alpha = 1, shape = 21) +
  scale_size_continuous(limits = c(1.9,3.5),
                        range = c(1,13),
                        breaks = c(2,2.5,3,3.5),
                        name = "Hits per million (log)") +
  facet_grid(Pathway + Compound ~ Site,
             scales = "free",
             space = "free",
             switch = "x") +
  theme_bw() +
  theme(plot.margin = margin(t = 15,
                             r = 10,
                             b = 15,
                             l = 15),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 14),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_rect(fill = NA, colour = "lightgray"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(angle = 0, size = 11.5),
        panel.border = element_rect(colour = "lightgray"),
        axis.title.y = element_text(vjust=3, color = "black", size = 14),
        axis.title.x = element_text(vjust=-0, color = "black", size = 14)) +
  xlab("MAG/Bin") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(size=5, alpha = 1), order = ))
bubble_plot_bins_HPM
