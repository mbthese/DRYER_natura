#libraries


library(metabaR)
library(dada2)

load("../resources/Metabarlist_natura_clean_ITS2_traits_alpha.Rdata") # ITS 

seqs<- natura_clean$motus$sequence
taxa <- assignTaxonomy(seqs, "../resources/sh_general_release_dynamic_16.10.2022.fasta") #from UNITE

write.csv(taxa, "../results/taxa.ITS.csv")