## PREP EXCLUSIONS ## 

#IMPT: this is based on release 2.0.1 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)

#Import data ====
covs <- readRDS("abcd201_mdd_covs_allsubj.rds")
disease_antidep <- readRDS("abcd201_disease_antidep_allsubj.rds")

#Determine unrelated sample ====
tmp <- merge(covs, disease_antidep, by = "src_subject_id")

#remove neurologically unhealthy indiv before determining unrelatedness 
#if not neurologically unhealthy may end up getting retained and removed later on
id_keep <- tmp %>% 
  filter(neuro_status == 0) %>% 
  select(src_subject_id, rel_family_id) %>% 
  as.data.frame()

#select unrelated individuals following method in GS-Imaging (using rnorm)
set.seed(2021) 
rand <- id_keep %>%  
  group_by(rel_family_id) %>% 
  mutate(random = sample.int(n())) #assign random number to each individual within a family

max_rand <- rand %>% 
  group_by(rel_family_id) %>% 
  summarise(max = max(random)) 

select_indiv <- merge(rand, max_rand, by="rel_family_id") %>% 
  mutate(keep = ifelse(random == max, 1, 0)) %>% 
  filter(keep==1) 

#final ids (i.e., no neuro condition and unrelated
unrel_id <- select_indiv["src_subject_id"]

#Determine whiteonly sample in whole sample ====
#using kmeans clustering of whole sample genetic PCs 

#import PCs generated from whole sample PCA (done using plink on Eddie)
setwd("C:/Users/gladi/OneDrive/Desktop/ABCD/MDD-PRS/dataprep/")
pca <- read.delim("ABCD3.0_wholesample_genotyped_PCA.eigenvec", header=F, sep="") %>% .[,-1] #first col redundant
colnames(pca)[1] <- "src_subject_id"
colnames(pca)[2:16] <- rep(1:15)
colnames(pca)[2:16] <-paste0('PC_',colnames(pca)[2:16])

#fix the 2IDs that are messed up -> `NDAR_INVF3FYXH1G and NDARINVPWLFYWTX
pca$src_subject_id[grep("`NDAR_INVF3FYXH1G",pca$src_subject_id)] <- "NDAR_INVF3FYXH1G"
pca$src_subject_id[grep("NDARINVPWLFYWTX",pca$src_subject_id)] <- "NDAR_INVPWLFYWTX"

#do kmeans 
#follow steps: https://static-content.springer.com/esm/art%3A10.1038%2Fng.3768/MediaObjects/41588_2017_BFng3768_MOESM40_ESM.pdf
set.seed(2021)
km_pc1 <- kmeans(pca[2],4,nstart=10) 
km_pc2 <- kmeans(pca[3],4,nstart = 10)

pca_cluster <- data.frame(pca, km_pc1 = km_pc1$cluster, km_pc2 = km_pc2$cluster)
race_reported <- data.frame(src_subject_id = covs$src_subject_id, 
                            race = covs$race_ethnicity) #white, black, hispanic, asians,others in order
pca_race_cluster <- merge(pca_cluster, race_reported, by="src_subject_id", all.x=T)

#create intersection between PCs
#table output:$1 is row, $2 is col
table(pca_race.cluster$race, pca_race_cluster$km_pc1) 
table(pca_race_cluster$race, pca_race_cluster$km_pc2) 

#select white based on numbers - ie white have the largest sample
#km_pc white cluster will differ each time you run, so check before running the next  line
pca_race_cluster$intersection_white <- ifelse(pca_race_cluster$km_pc1 ==3 & pca_race_cluster$km_pc2 ==2,1,0) #actually PC1 seems sufficient
table(pca_race_cluster$intersection_white)

#white only - select those that report ethnicity as white + white cluster + parents report white
white_only <- pca_race_cluster %>% 
  filter(race == 1 & intersection_white == 1) %>% 
  select(src_subject_id)

white_id <- merge(white_only, covs[covs$parent_eth==1,]) %>% select(src_subject_id)

#determine included sample (aka white only, neuro healthy, unrelated)
id_keep <- merge(unrel_id, white_id, by = "src_subject_id")

## Save ====
saveRDS(id_keep, "abcd_noneuro_unrel_whiteonly_ids.rds")