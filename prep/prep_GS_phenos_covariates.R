## PREP GS PHENOTYPES EXCL IMAGING ##

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(broom);library(dplyr);library(stringr);library(plyr)

#Load necessary data ====
load("masterDB.RData") #should see totaldata in env
linkage_id <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")

#Covariates ====

#demographic variables at imaging assessment 
demo <- read_excel("STRADL_collated.xlsx",sheet="demographics") %>% 
  select(ID,AgeFaceToFace,Sex)
colnames(demo) <- c("ID","age","sex")

#education level 
basic_info <- read_excel("STRADL_collated.xlsx",sheet="basicinfo")
edu <- basic_info %>% select(ID,Education) %>% mutate_at(vars(Education),as.factor)
colnames(edu) <- c("ID","edu")

#income
gs_income <- totaldata %>% 
  select(id, income) %>% #(1 - less than 10,000, 2 - between 10,000 and 30,000, 3 - between 30,000 and 50,000, 4 - between 50,000 and 70,000, 5 - more than 70,000, 6 - prefer not to answer, 99 - not known)
  mutate_at(vars(income), as.factor) %>% 
  merge(., linkage_id[2:3], by.x="id", by.y="genotype_ID") %>% #select those w stradl only
  select(stradl_ID, income)
colnames(gs_income) <- c("ID","GS_income")

gs_income$GS_income <- as.factor(gsub("6","NA",gs_income$GS_income))
gs_income$GS_income <- as.factor(gsub("7","NA",gs_income$GS_income))
gs_income$GS_income <- as.factor(gsub("99","NA",gs_income$GS_income))

table(gs_income$GS_income)

#imaging related covariates 
img_cov <- read_csv("STRADL_Measures_FreeSurfer_Main.csv") %>% 
  select(id,site,batch,edited) %>% 
  mutate_at(vars(site:edited),as.factor)
colnames(img_cov)[1] <- "ID"

icv_std <- read_csv("STRADL_Measures_Standardised_ICV.csv") #std by site 
colnames(icv_std)[2] <- "ICV"

#genetic PCs
mds <- read.table("HM3mds.mds",header=T) %>% select(IID, C1:C15) 
mds_white <- mds[c(1:19994),] %>% 
  merge(linkage_id[2:3],., by.x = "genotype_ID", by.y = "IID")

#compile
covariates <- full_join(demo,edu) %>% 
  full_join(.,gs_income) %>% 
  full_join(.,img_cov) %>% 
  full_join(.,icv_std) %>% 
  full_join(., mds_white, by.x ="ID", by.y = "stradl_ID")

#mdd ====
qids <- read_excel("STRADL_collated.xlsx",sheet="QIDS") %>% select(ID,TotalQIDS,SeverityDepression)
scid <- read_excel("STRADL_collated.xlsx",sheet="SCID") %>% select(ID,DiagnosisGiven) #binary score of no history of depression (0) or lifetime episode of depression (1).
colnames(scid)[2] <- "scid_stradl"

load("STRADL.RData")
cidi <- x %>%  
  select(id,CIDI_MDD) %>% 
  filter(CIDI_MDD == 0 | CIDI_MDD == 1) %>% 
  merge(linkage.id[2:3],.,by.x="genotype_ID",by.y="id")  %>% 
  select(stradl_ID,CIDI_MDD) 
colnames(cidi) <- c("ID","cidi")

mdd <- full_join(qids, scid) %>% 
  mutate(lifetime_MDD = ifelse(cidi ==1 | scid_stradl ==1,1,
                               ifelse(cidi ==0 & scid_stradl==0,0,NA)))
table(mdd$lifetime_MDD)

#mdd-prs ====
prs_p0.001 <- read.table("Meta_predict_FullGS.S4.profile",header=T) %>% select(IID,SCORE_0.001)
colnames(prs_p0.001)[2] <-"p_0.001"

prs_p0.01 <- read.table("Meta_predict_FullGS.S5.profile",header=T)%>% select(IID,SCORE_0.01)
colnames(prs_p0.01)[2] <-"p_0.01"

prs_p0.05 <- read.table("Meta_predict_FullGS.S6.profile",header=T)%>% select(IID,SCORE_0.05)
colnames(prs_p0.05)[2] <-"p_0.05"

prs_p0.1 <- read.table("Meta_predict_FullGS.S7.profile",header=T)%>% select(IID,SCORE_0.1)
colnames(prs_p0.1)[2] <-"p_0.1"

prs_p0.2 <- read.table("Meta_predict_FullGS.S8.profile",header=T)%>% select(IID,SCORE_0.2)
colnames(prs_p0.2)[2] <-"p_0.2"

prs_p0.5 <- read.table("Meta_predict_FullGS.S9.profile",header=T)%>% select(IID,SCORE_0.5)
colnames(prs_p0.5)[2] <-"p_0.5"

prs_p1 <- read.table("Meta_predict_FullGS.S10.profile",header=T)%>% select(IID,SCORE_1)
colnames(prs_p1)[2] <-"p_1"

prs_fullsample <- full_join(prs_p0.001,prs_p0.01) %>% 
  full_join(.,prs_p0.05) %>% 
  full_join(.,prs_p0.1) %>% 
  full_join(.,prs_p0.2) %>% 
  full_join(.,prs_p0.5) %>% 
  full_join(.,prs_p1) %>% 
  merge(linkage_id[2:3],.,by.x="genotype_ID",by.y="IID") #select those with stradl ID only

#compile everything ==
mdd_covs_prs <- full_join(mdd, prs, by.x = "ID", by.y = "stradl_ID") %>% 
  full_join(., covariates, by = "ID")

saveRDS(mdd_covs_prs, "stradl_mdd_covs_prs_allsubj.rds")
