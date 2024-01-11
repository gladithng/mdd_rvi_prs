## PREP ABCD DISEASE STATUS / EXCLUSIONS ## 

#IMPT: this is based on release 2.0.1 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)

#Disease status ====
disease <- readRDS("abcd_screen01.rds")

#select neurological and mental disorders
#reverse code if needed
tmp <- disease %>% 
  .[,grep("subject|interview|cpalsy|tumor|stroke|aneurysm|hemorrhage|hemotoma|
         medcond|epls|commondx|schiz|asd$|sud|intdiasb|psychdx",colnames(.))] %>% 
  mutate_at(vars(scrn_cpalsy:scrn_psychdx_other),as.numeric) %>% 
  mutate(scrn_cpalsy_reverse = ifelse(scrn_cpalsy == 0, 1, 0),
         scrn_tumor_reverse = ifelse(scrn_tumor == 0, 1, 0),
         scrn_stroke_reverse = ifelse(scrn_stroke == 0, 1, 0),
         scrn_aneurysm_reverse = ifelse(scrn_aneurysm == 0, 1, 0),
         scrn_hemorrhage_reverse = ifelse(scrn_hemorrhage == 0, 1, 0),
         scrn_hemotoma_reverse = ifelse(scrn_hemotoma == 0, 1, 0)) #reverse code questions so that 0 means no

disease_status <- tmp %>% 
  mutate(psych_condition = rowSums(tmp %>% select(scrn_commondx:scrn_psychdx_other), na.rm=T),#0 means no, 1 means yes
         neuro_condition = rowSums(tmp %>% select(scrn_cpalsy_reverse:scrn_hemotoma_reverse, scrn_epls), na.rm=T), 
         psych_status = ifelse(psych_condition == 0, 0, 1), #0 meaning no psych problems
         neuro_status = ifelse(neuro_condition == 0, 0, 1)) #0 meaning no neuro problems 


#antidepressant status ==
medi_dat <- readRDS("medsy01.rds")
img_dat <- readRDS("abcd_smrip201.rds") #to use as reference to match date 

#check which instance to use to match with imaging assessment
medi_dat$interview_date <- as.Date(medi_dat$interview_date, "%m/%d/%Y")
img_dat$interview_date <- as.Date(img_dat$interview_date, "%m/%d/%Y")
img_dat_tomerge <- img_dat[,c('src_subject_id','interview_date','interview_age')]
medi_img_dat <- merge(medi_dat, img_dat_tomerge, by = c('src_subject_id','interview_age','interview_date'))

#extract medication names 
medi_name <- medi_img_dat[,grep('src_subject_id|brought_medications|_rxnorm_p',colnames(medi_dat))]
medi_prescribed <- medi_name[,grep('_rxnorm_p',colnames(medi_name))]
medi_prescribed <- medi_prescribed[,!grepl('_otc_',colnames(medi_prescribed))] #questionnaire asked for another kind of medication (otc?) which you dont need

medi_name_all <- as.vector(medi_prescribed)
medi_name_all <- medi_name_all[!is.na(medi_name_all)]
medi_name_all <- medi_name_all[medi_name_all!='']
medi_name_all <- unique(medi_name_all)
medi_name_all <- tolower(medi_name_all) #change to lower case to suit ls.antidepressant

#load antidepressant list
ls_antidepres <- read.table('ls.antidepressants.txt',header=F,stringsAsFactors=F)
ls_antidepres <- ls_antidepres$V1
ls_anti_ABCD <- medi_name_all[grep(paste0(ls_antidepres,collapse = '|'),medi_name_all)] #getting list of antidepressants that is present in abcd cohort

detec_antidepre <- function(dat, ls_tofind){
  k = sum(tolower(dat) %in% ls_tofind,na.rm=T) #get the sum of # of antidepressants
  if(k==0){output = 0}else{output = 1}
  return(output)
}

medi_name$anti_count = apply(medi_name, 1, detec_antidepre, ls.anti.ABCD) #apply function over rows
table(medi_name$anti_count)

med_status <- medi_name %>% select(src_subject_id,anti_count)
colnames(med_status)[2] <- "antidepressant_status"

#Merge disease status and antidepressant status 
#controls for RVI analysis will be mentally and physically healthy (no serious condition) individuals wo antidep use
disease_antidep <- full_join(disease_status, med_status, by = "src_subject_id") 
  mutate_at(vars(psych_status:antidepressant_status),as.factor)

#save ====
saveRDS(disease_antidep, "abcd201_disease_antidep_allsubj.rds")
  