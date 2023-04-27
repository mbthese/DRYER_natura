reads <- read_delim("E:/DRYER_natura_data/16S/DRYER_natura_R1R2_good_demultiplexed_derepl_basicfilt_cl_agg.tab", escape_double = F, trim_ws = T, delim ="\t")

reads <- reads[,which(grepl('sample',colnames(reads))==T)] 
colnames(reads) <- gsub('sample:','',colnames(reads))
#reads2 <- t(reads2)
# reads3 <- t(reads3)
pcrs <- read.delim2("E:/DRYER_natura_data/16S/pcrs.txt")

list1<- pcrs$sample_id
length(list1)
list2 <- colnames(reads)
length(list2)
list2[!(list2 %in% list1)]
#dplyr::setdiff(list2, list1) #does not work for whatever reason

read2 <-as.data.frame(reads2)
dplyr::setdiff(pcrs$sample_id, colnames(read2))
pcrs <- dplyr::filter(pcrs, !sample_id %in% setdiff(as.vector(pcrs$sample_id),as.vector(colnames(reads2))))

reads3 <- colnames(reads2)
reads3 <- as.data.frame(reads3)
nrow(pcrs)#407
colnames(reads3) 
essai <- as.data.frame(reads3) 
essai <- essai %>% rename(sample_id = reads3 )

essai$col2 <- '5'#just to add another column and see the missing individuals
A <- anti_join(essai, pcrs, by = 'sample_id')
B <- setdiff(as.vector(essai$sample_id), as.vector(pcrs$sample_id))

#missing are:

DRYER_natura_NA_16S_control_PCR_PCR-tags_21
DRYER_natura_NA_16S_control_PCR_PCR-tags_22
DRYER_natura_NA_16S_control_extraction_EXT_5

DRYER_natura_NA_ITS2_control_extraction_EXT_5
DRYER_natura_NA_ITS2_control_PCR_PCR-tags_21


#removing missing from pcrs and file_samples
pcrs <- pcrs %>% filter(sample_id != "DRYER_natura_NA_16S_control_PCR_PCR-tags_21")
pcrs <- pcrs %>% filter(sample_id != "DRYER_natura_NA_16S_control_PCR_PCR-tags_22")
pcrs <- pcrs %>% filter(sample_id != "DRYER_natura_NA_16S_control_extraction_EXT_5")
pcrs <- pcrs %>% filter(sample_id != "DRYER_natura_NA_ITS2_control_extraction_EXT_5")
pcrs <- pcrs %>% filter(sample_id != "DRYER_natura_NA_ITS2_control_PCR_PCR-tags_21")

file_samples <- file_samples %>% filter(sample_id != "DRYER_natura_NA_16S_control_PCR_PCR-tags_21")
file_samples <- file_samples %>% filter(sample_id != "DRYER_natura_NA_16S_control_PCR_PCR-tags_22")
file_samples <- file_samples %>% filter(sample_id != "DRYER_natura_NA_16S_control_extraction_EXT_5")
file_samples <- file_samples %>% filter(sample_id != "DRYER_natura_NA_ITS2_control_extraction_EXT_5")
file_samples <- file_samples %>% filter(sample_id != "DRYER_natura_NA_ITS2_control_PCR_PCR-tags_21")


motus <- data.frame(sequence=reads$sequence)

rownames(motus) <- rownames(reads2)

reads2 <- t(reads2)
rownames(pcrs) <- rownames(reads2)

rownames(file_samples) <- file_samples$sample_id
reads2 <- gsub('sample:','',reads2)


metabarList <- metabarlist_generator(reads2,motus,pcrs,file_samples)

# Compute the number of reads per pcr
metabarList$pcrs$nb_reads <- rowSums(metabarList$reads)

# Compute the number of motus per pcr
metabarList$pcrs$nb_motus <- rowSums(metabarList$reads>0)

# Load requested package for plotting
library(ggplot2)
library(reshape2)

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- melt(metabarList$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the soil_euk$pcrs table

ggplot(metabarList$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")
