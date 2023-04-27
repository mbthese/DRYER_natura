library(readr)
library(dplyr)
#change ngsfilter format for sophie's pipeline


ngsfilter <- read.csv("~/Desktop/DRYER_natura_data/ngs_filter_final2.csv")

ngsfilter <-as.data.frame(sapply(ngsfilter, function(x) gsub("\"", "", x)))

ngsfilter <- ngsfilter %>%
  mutate(tag = paste0(barcode_forward, ":", barcode_reverse)) %>%
  mutate(extra_information = paste0("F"," ","@"," ",paste0("species=",Family,",", Genus,",", Species), sep=";")) %>%
  select(project, samples_name, tag, forward_primers_sequence, reverse_primers_sequence, extra_information)

#The forward primer creates copies of the 5’-3’ strand 
# whereas the reverse primer makes copies of the complementary (runs 3’-5’) strand.
#the workflow wants : The tags correspond to short and specific sequences added on the 5’ (so the forward) end of each primer to distinguish the different samples
#trop bizarre, i'll keep both



# Write data to txt file: tab separated values
# sep = "\t"
write.table(ngsfilter, file = "resources/DRYER_natura_ngsfilter.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
