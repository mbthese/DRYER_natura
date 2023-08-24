library(ggpubr)
library(ade4)
library(vegan)
library(scales)
library(ggpubr)
library(BiocManager)
library(TITAN2)
library(purrr)

load("../resources/Metabarlist_natura_clean_16S.Rdata")# 16S

TITcom <- natura_clean$reads # get my 16s community
TITcom <- TITcom[,colSums(TITcom != 0)>=3] # remove otus present in less than three samples, don't have 
TITcom <- TITcom[,colSums(TITcom)>100] # keep otus with more than 100 reads, don't have
#rownames(TITcom) <- as.numeric(gsub("_.+$","",rownames(TITcom))) # rename rows as numbers
TITenv <- natura_clean$samples # get leaf data of samples for which i have the community
names(natura_clean$samples)
TITenv <- natura_clean$samples %>% dplyr::select(LSWC, SLA, TLP, MajVLA, LT, SRL, RootShoot, StemDiameter, Height, LA_rwc, gmin)
#rownames(TITenv) <- rownames(data_leaf[data_leaf$number%in%rownames(TITcom),]) # rename rows as numbers
dim(TITcom);dim(TITenv) # check that dims are compatible
# [1] 132 764
# [1] 132  48

z <- TITcom # rename
quiet_titan <- quietly(titan)

x <- TITenv
n <- names(TITenv)


walk2(TITenv, names(TITenv), function(x, n) { # iterate titan for all traits of TITenv
  y <- z[!is.na(x),] # rm when having na in env
  y <- y[,apply(y, 2, function(x) sum(x>0))>3] # rm otus present less than 3 times in the data (alrdy done)
  x <- x[!is.na(x)]
  titan(x, y, nBoot = 500, ncpus = 3, messaging = FALSE) %>%
    saveRDS(paste("../results/titan/", n, "16s.RDS", sep="_")) # save the results in a .RDS file
})
