## PREP ABCD MDD PHENO + COVS ##

#IMPT: this is based on release 2.0.1 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)

## MDD pheno ====
abcd_cbcls01 <- readRDS("abcd_cbcls01.rds") 

abcd_cbcl <- abcd_cbcls01 %>%
  filter(eventname == "baseline_year_1_arm_1") %>% 
  select(contains(c("src_subject_id","interview_age","gender","dsm5","internal","external","totprob","somatic"))) %>% #DSM scales
  select(contains(c("src_subject_id","interview_age","gender","_r")))  #took raw scores 
colnames(abcd_cbcl)[2:3] <- c("age","sex")

#check distribution 
abcd_cbcl %>%
  select(starts_with("cbcl")) %>%
  gather() %>%
  ggplot(aes(as.numeric(value)))+
  facet_wrap(~ key, scales = "free") +
  geom_histogram(fill = "deepskyblue4") #skewed as expected, majority at 0

## Covariates ====

#demographics ==
demo <- readRDS("acspsw03.rds") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% #take baseline info
  select(src_subject_id,interview_age,gender,race_ethnicity,rel_family_id) %>% 
  mutate_at(vars(gender:rel_family_id),as.factor)
colnames(demo)[2:3] <- c("age","sex")

#parents ethnicity (needed for subsetting whiteonly sample) ==
parents_eth <- readRDS("abcd_meim01.rds") %>% 
  select(src_subject_id, meim_ethnic_id_p) %>% 
  mutate_at(vars(meim_ethnic_id_p), as.factor)
colnames(parents_eth)[2] <- "parent_eth"
table(parents_eth$parent_eth) #white is coded as 1

#study site ==
abcd_lt01 <- readRDS("abcd_lt01.rds")
site_id <- abcd_lt01 %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  separate(site_id_l,sep=4,into=c(NA,"site.id")) %>%  
  mutate_at(vars(site.id),as.numeric) %>% 
  mutate_at(vars(site.id),as.factor) %>% .[,c(1,6)]

#education and income ==
pdem02 <- readRDS("pdem02.rds")
cols <- c('src_subject_id', 'demo_prnt_ed_v2$', 'demo_comb_income_v2$') 

edu_income <- pdem02 %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% #no baseline
  .[,grep(paste0(cols,collapse="|"),colnames(.))] %>% 
  mutate_at(vars(demo_prnt_ed_v2,demo_comb_income_v2), as.factor)
colnames(edu_income)[2:3] <- c("parent_edu","fam_income")

edu_income$parent_edu <- as.factor(gsub("777","NA",edu_income$parent_edu))
edu_income$fam_income <- as.factor(gsub("777","NA",edu_income$fam_income))
edu_income$fam_income <- as.factor(gsub("999","NA",edu_income$fam_income))
apply(edu_income, 2, function(x) sum(is.na(x)))

#genetic covariates ==

#batch info ==
batch_info <- read.delim("ABCD_release3.0_.batch_info.txt") %>% .[,2:4] %>% 
  mutate_at(vars(Axiom_Plate:BATCH),as.factor)
colnames(batch_info) <- c("src_subject_id","genotype_plate_no","batch_no")
#fix the 2IDs that are messed up -> `NDAR_INVF3FYXH1G and NDARINVPWLFYWTX
batch_info$src_subject_id[grep("`NDAR_INVF3FYXH1G",batch_info$src_subject_id)] <- "NDAR_INVF3FYXH1G"
batch_info$src_subject_id[grep("NDARINVPWLFYWTX",batch_info$src_subject_id)] <- "NDAR_INVPWLFYWTX"

#pca 
pca <- read.delim("ABCD3.0_imputed_whiteonlywparents_MAF0.005_PCA.eigenvec", header=F, sep="") %>% .[,-1] #%>% .[,1:16] #select top 15 PCs
colnames(pca)[1] <- "src_subject_id"
colnames(pca)[2:16] <- rep(1:15)
colnames(pca)[2:16] <-paste0('PC_',colnames(pca)[2:16])


#Compile mdd pheno and covariates together ====
abcd_mdd_covs <- full_join(abcd_cbcl, demo) %>% 
  full_join(., parents_eth) %>% 
  full_join(., site_id) %>%
  full_join(., edu_income) %>% 
  full_join(., batch_info) %>% 
  full_join(., pca)

saveRDS(abcd_mdd_covs, "abcd201_mdd_covs_allsubj.rds")
