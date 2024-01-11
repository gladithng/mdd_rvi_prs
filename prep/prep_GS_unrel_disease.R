## PREP UNRELATED, DISEASE STATUS IN GS ##

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(broom);library(dplyr);library(stringr);library(readxl)

#Exclusions ====

#neurological conditions
basic_info <- read_excel("STRADL_collated.xlsx",sheet="basicinfo")
neuro <- basic_info[grep("epilepsy|stroke|parkinson's|aneurysm|seizure|haemorrhage|hematoma|cerebral palsy|brain tumour",
                         basic_info$IllnessDiagnosis,ignore.case = T),] 
loc <- grep("misdiagnosis|viral",neuro$IllnessDiagnosis,ignore.case = T)
neuro_clean <- neuro[-loc,]

#relatedness 
#stradl is family rich, hence remove duplicate familyid to get unrelated
#generate random number for each participant
#for duplicated family id, keep participant with higher random number 

#get famid
load("masterDB.RData") #should see totaldata in env
linkage_id <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")

full_id <- totaldata[, c("id", "famid")] %>% 
  merge(linkage_id[2:3], ., by.x = "genotype_ID", by.y = "id")

#select unrelated
set.seed(2021)
full_id <- full_id %>% 
  mutate(random = rnorm(n = nrow(.))) 

related <- full_id %>%  #to list ALL duplicated rows
  select(stradl_ID,famid,random) %>% 
  group_by(famid) %>% 
  filter(n()>1) 

unrelated <- full_id[full_id$stradl_ID %in% related$stradl_ID ==F,] 

max_rnorm_byfamid <- related %>% 
  group_by(famid) %>% 
  summarise(max=max(random)) %>% 
  ungroup()

select_indiv <- merge(related,max_rnorm_byfamid) %>% 
  mutate(keep = ifelse(random == max, 1,0)) %>% 
  filter(keep==1) 

unrel_keep <- full_join(select_indiv["stradl_ID"],unrelated["stradl_ID"])
colnames(unrel_keep) <- "ID"

#final ids to keep - no neuro and unrelated
ids_keep <- full_id %>% 
  .[.$stradl_ID %in% neuro_clean$ID == F,] %>% #remove neuro
  .[.$stradl_ID %in% unrel_keep$ID == T,] #remove related

saveRDS(ids_keep, "stradl_noneuro_unrel_ids.rds")

#Psych, antidepressants status ====
psych <- basic_info[grep("ADHD|bipolar|anxiety|depression|schizophrenia|autism|intellectual disability|substance", 
                         basic_info$IllnessDiagnosis,ignore.case = T),] 

antidep <- read_csv("STRADL_Medications_Aleks.csv") %>% 
  .[.$Meds_Antidepressant == 1,]

#compile disease status
disease_status <- basic_info %>% 
  mutate(neuro_status = ifelse(basic_info$ID %in% neuro_clean$ID ==T,1,0), #0 meaning no neuro problems 
         psych_status = ifelse(basic_info$ID %in% psych$ID ==T,1,0), #0 meaning no psych problems
         antidepressant_status = ifelse(basic_info$ID %in% antidep$ID ==T,1,0)) %>% #0 meaning no antidep use
  select(ID,neuro_status,psych_status, antidepressant_status) %>% 
  mutate_at(vars(neuro_status:antidepressant_status),as.factor)
apply(disease_status[-1], 2, function(x) table(x))

saveRDS(disease_status, "stradl_disease_antidep_allsubj.rds")


