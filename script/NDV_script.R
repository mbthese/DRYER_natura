#https://rpubs.com/longroad/882321
# Adjusted from Stegen et al 2012 GEB


#libraries

library(tidyverse)
library(labdsv)
library(dplyr)
library(FactoMineR)
library(car)
library(missMDA)
library(metabaR)
library(ggplot2)
library(corrplot)
library(factoextra)
library(MASS)
library(cowplot)
library(nlme)
library(MuMIn)
library(multcomp)
library(ggpubr)
library(ade4)
library(hillR)
library(adespatial)
library(reldist)
library(bipartite)
library(vegan)
library(scales)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(BiocManager)
library(microbiome)  


# leaf ITS 


load("./resources/Metabarlist_natura_clean_ITS2_traits_alpha.Rdata")
leaf <- subset_metabarlist(natura_clean, table = "samples",
                           indices = natura_clean$samples$organ == "leaf")


# Prepare and calculate abundance beta-null deviation metric
comm.t=t(leaf$reads) #transpose reads
bbs.sp.site <- comm.t
patches=nrow(bbs.sp.site)
rand <- 3

# Create empty matrices/vectors, the dimensions of these matrices/vectors are determined based on the number of species (OTU) or sites (number of seedlings) in the data.
null.alphas <- matrix(NA, ncol(comm.t), rand)
null.alpha <- matrix(NA, ncol(comm.t), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)

# Abundance Beta-Null Deviation Calculation
bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) #normalisation of abundance data by dividing to max abundance value
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) # mean.alpha : mean abundance per site (seedlings)
gamma <- ncol(bbs.sp.site) # gamma diversity : The total number of species : number of columns in bbs.sp.site, because we transposed it at the beginning
obs_beta <- 1-mean.alpha/gamma #observed beta-diveristy : 1- mean alpha/gamma: the effective species turnover expressed as a proportion of total species diversity (H. Tuomisto 2014)
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma #difference with obs_beta?

#loop that iterates *rand* = 999 times to calculate the null deviation metric.
#Create a null distribution of the community data, maintain the total abundance per species and randomly reassign individuals to sites/seedlings 
for (randomize in 1:rand) {  
  null.dist = comm.t
  for (species in 1:ncol(null.dist)) {
    tot.abund = sum(null.dist[,species]) #for each species (OTU), calculates the total abundance by summing the values in the corresponding column of "null.dist". 
    null.dist[,species] = 0 # It then sets all values in that column to 0.
    for (individual in 1:tot.abund) {
      sampled.site = sample(c(1:nrow(bbs.sp.site)), 1) #randomly reassign individuals from the total abundance to sites/seedlings 
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }
  }
  ##Calculate null deviation metrics (alpha, beta, gamma) for null patches and store
  null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
  null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
  expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
  null.alpha <- mean(null.alphas[,randomize])
  null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
  bucket_bray <- as.matrix(vegdist(null.dist, "bray")) #dissmiliarity matrix using bray-curtis distances
  diag(bucket_bray) <- NA
  bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE) #Mean dissimilarity values for each species are calculated and stored.
} ## end randomize loop 

## Calculate beta-diversity for observed metacommunity
beta_comm_abund <- vegdist(comm.t, "bray")
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
diag(res_beta_comm_abund) <- NA #diagonale null comparaisson deux a deux?

#output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE) #on the matrix 10099*10099, aplpy mean on the columns, and output is a list 10099. 

# output abundance beta-null deviation
bray_abund_null_dev_fun <- beta_div_abund_stoch - mean(bucket_bray_res) #NDV : difference between the mean observed dissimilarity values and the mean dissimilarity values from the null distributions.

#store in R project
save(bray_abund_null_dev_fun,file="./results/bray_abund_null_dev_fun.RData")
load("./results/null/bray_abund_null_dev_fun.RData")

boxplot(bray_abund_null_dev_fun,outline=FALSE)

df_bray_fun_endo<-data.frame(bray_abund_null_dev_fun)
colnames(df_bray_fun_endo)<-c("NDV")

#root ITS

root <- subset_metabarlist(natura_clean, table = "samples",
                           indices = natura_clean$samples$organ == "root")

# Prepare and calculate abundance beta-null deviation metric
comm.t.root=t(root$reads) #transpose reads
bbs.sp.site.root <- comm.t.root
patches=nrow(bbs.sp.site.root)
rand <- 99

# Create empty matrices/vectors, the dimensions of these matrices/vectors are determined based on the number of species (OTU) or sites (number of seedlings) in the data.
null.alphas <- matrix(NA, ncol(comm.t.root), rand)
null.alpha <- matrix(NA, ncol(comm.t.root), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)

# Abundance Beta-Null Deviation Calculation
bbs.sp.site.root = ceiling(bbs.sp.site.root/max(bbs.sp.site.root)) #normalisation of abundance data by dividing to max abundance value
mean.alpha = sum(bbs.sp.site.root)/nrow(bbs.sp.site.root) # mean.alpha : mean abundance per site (seedlings)
gamma <- ncol(bbs.sp.site.root) # gamma diversity : The total number of species : number of columns in bbs.sp.site.root, because we transposed it at the beginning
obs_beta <- 1-mean.alpha/gamma #observed beta-diveristy : 1- mean alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site.root)/gamma #difference with obs_beta?

#loop that iterates *rand* = 999 times to calculate the null deviation metric.
#Create a null distribution of the community data, maintain the total abundance per species and randomly reassign individuals to sites/seedlings 
for (randomize in 1:rand) {  
  null.dist = comm.t.root
  for (species in 1:ncol(null.dist)) {
    tot.abund = sum(null.dist[,species]) #for each species (OTU), calculates the total abundance by summing the values in the corresponding column of "null.dist". 
    null.dist[,species] = 0 # It then sets all values in that column to 0.
    for (individual in 1:tot.abund) {
      sampled.site = sample(c(1:nrow(bbs.sp.site.root)), 1) #randomly reassign individuals from the total abundance to sites/seedlings 
      null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
    }
  }
  ##Calculate null deviation metrics (alpha, beta, gamma) for null patches and store
  null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
  null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
  expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
  null.alpha <- mean(null.alphas[,randomize])
  null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
  bucket_bray <- as.matrix(vegdist(null.dist, "bray")) #dissmiliarity matrix using bray-curtis distances
  diag(bucket_bray) <- NA
  bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE) #Mean dissimilarity values for each species are calculated and stored.
} ## end randomize loop 

## Calculate beta-diversity for observed metacommunity
beta_comm_abund <- vegdist(comm.t.root, "bray")
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
diag(res_beta_comm_abund) <- NA #diagonale null comparaisson deux a deux?

#output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)

# output abundance beta-null deviation
bray_abund_null_dev_fun_root <- beta_div_abund_stoch - mean(bucket_bray_res) #NDV : difference between the mean observed dissimilarity values and the mean dissimilarity values from the null distributions.

#store in R project
save(bray_abund_null_dev_fun_root,file="./results/bray_abund_null_dev_fun_root.RData")
load("./results/null/bray_abund_null_dev_fun_root.RData")

boxplot(bray_abund_null_dev_fun, bray_abund_null_dev_fun_root,outline=FALSE)
t.test(bray_abund_null_dev_fun,bray_abund_null_dev_fun_root, alternative = "two.sided", var.equal = FALSE)
#test the (null) hypothesis that two populations have equal means. It is named for its creator, Bernard Lewis Welch, is an adaptation of Student's t-test,[1] and is more reliable when the two samples have unequal variances and possibly unequal sample sizes.[2][3] 

#Test statistic (t-value): The t-value is -193.71. This value measures the difference between the means of the two groups relative to the variability within each group. In this case, the t-value indicates a substantial difference between the means of the two samples.

# Degrees of freedom (df): The degrees of freedom for this test are 15000. It represents the number of independent pieces of information available for estimating the population parameter. A higher degree of freedom allows for a more precise estimation.
# 
# p-value: The p-value is < 2.2e-16, which is a very small value. The p-value represents the probability of observing a test statistic as extreme as the one calculated, assuming the null hypothesis is true. In this case, the p-value is extremely small, suggesting strong evidence against the null hypothesis (the true difference in means being zero).
# 
# Alternative hypothesis: The alternative hypothesis is that the true difference in means is not equal to 0. The results support this alternative hypothesis.

df_bray_fun_endo_root<-data.frame(bray_abund_null_dev_fun_root)
colnames(df_bray_fun_endo_root)<-c("NDV")

#Plot

#ITS
# load("./results/null/bray_abund_null_dev_ITS_leaf_nontranspose.RData") #leaf
# load("./results/null/bray_abund_null_dev_ITS_root_nontranspose.RData") #root
load("./results/null/bray_abund_null_dev_ITS_leaf_nontranspose_hel.RData") #leaf
load("./results/null/bray_abund_null_dev_ITS_root_nontranspose_hel.RData") #root


ndv_ITS_leaf <-data.frame(bray_abund_null_dev_fun)
ndv_ITS_root <-data.frame(bray_abund_null_dev_fun_root)
colnames(ndv_ITS_leaf)<-c("NDV")
colnames(ndv_ITS_root)<-c("NDV")

boxplot(bray_abund_null_dev_fun, bray_abund_null_dev_fun_root,outline=FALSE)
t.test(bray_abund_null_dev_fun,bray_abund_null_dev_fun_root, alternative = "two.sided", var.equal = FALSE)

plot_ITS <- ggplot()+
  geom_violin(data=ndv_ITS_leaf,aes(x=factor(0),y=NDV),colour="black", fill="#009E73",draw_quantiles = c(0.25,0.5,0.75),trim=FALSE) +
  geom_violin(data=ndv_ITS_root,aes(x=factor(1),y=NDV),colour="black",fill="#E69F00",draw_quantiles = c(0.25,0.5,0.75),trim=FALSE)+
  scale_x_discrete(breaks=c(0,1),labels=c("Leaf","Root"))+
 # ggtitle("Fungi null beta-distribution")+
  annotate("text",x=1,y=0.7,label="T.test, p<2.2e-16", size = 6)+
  theme_minimal(base_size =25)+
  theme(axis.title.x = element_blank())+
  ylim(0.2, 0.8)
plot_ITS


ggsave("./results/plot_ITS_ndv_hell.jpeg", plot_ITS)

#16S

# load("./results/null/bray_abund_null_dev_16S_leaf_nontranspose.RData") #leaf
# load("./results/null/bray_abund_null_dev_16S_root_nontranspose.RData") #root
load("./results/null/bray_abund_null_dev_16S_leaf_nontranspose_hel.RData") #leaf
load("./results/null/bray_abund_null_dev_16S_root_nontranspose_hel.RData") #root


boxplot(bray_abund_null_dev_fun, bray_abund_null_dev_fun_root,outline=FALSE)
t.test(bray_abund_null_dev_fun,bray_abund_null_dev_fun_root, alternative = "two.sided", var.equal = FALSE)

ndv_16S_leaf <-data.frame(bray_abund_null_dev_fun)
ndv_16S_root <-data.frame(bray_abund_null_dev_fun_root)
colnames(ndv_16S_leaf)<-c("NDV")
colnames(ndv_16S_root)<-c("NDV")

plot_16S <- ggplot()+
  geom_violin(data=ndv_16S_leaf,aes(x=factor(0),y=NDV),colour="black", fill="#009E73",draw_quantiles = c(0.25,0.5,0.75),trim=FALSE) +
  geom_violin(data=ndv_16S_root,aes(x=factor(1),y=NDV),colour="black",fill="#E69F00",draw_quantiles = c(0.25,0.5,0.75),trim=FALSE)+
  scale_x_discrete(breaks=c(0,1),labels=c("Leaf","Root"))+
  #ggtitle("Bacteria null beta-distribution")+
  annotate("text",x=1,y=0.7,label="T.test, p-value=0.43", size = 6)+
  theme_minimal(base_size = 25)+
  theme(axis.title.x = element_blank())+
  ylim(0.2, 0.8)

plot_16S
ggsave("./results/plot_16S_ndv_nontranspose_boxplot_hell.jpeg", plot_16S)

NDV_final_plot <- ggarrange(plot_16S, plot_ITS, labels = c("Bacteria", "Fungi"))
NDV_final_plot

ggsave("./results/null/plot_final_ndv_hel.jpeg", NDV_final_plot, width = 10, height = 8 )
