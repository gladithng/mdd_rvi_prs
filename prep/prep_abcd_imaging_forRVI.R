## PREP ABCD DATA##

#IMPT: this is based on release 2.0.1 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)

#Load necessary data 
id_keep <- readRDS("abcd_noneuro_unrel_whiteonly_ids.rds") #df; technically also took into consideration subj that have genetic data
disease_antidep <- readRDS("abcd201_disease_antidep_allsubj.rds") #df
covs <- readRDS("abcd201_mdd_covs_allsubj.rds")

#Select ids that passed processing ====
#follow recommended QC criteria by ABCD in 2.NDA 2.01 Release Notes Imaging Instrument
qc_pre <- readRDS("mriqcrp102.rds") %>% filter(eventname == "baseline_year_1_arm_1") #for satisfactory t1 raw images
qc_post <- readRDS("freesqc01.rds")  %>% filter(eventname == "baseline_year_1_arm_1") #for passed freesurfer QC
qc_criteria <- full_join(qc_pre,qc_post)

t1_ids_pass <- qc_criteria %>% 
  select(src_subject_id,iqc_t1_ok_ser,fsqc_qc) %>% 
  filter(iqc_t1_ok_ser >0) %>% #satisfactory t1 raw scans
  filter(fsqc_qc==1) #pass freesurfer qc
t1_ids_pass <- t1_ids_pass$src_subject_id

dti_ids_pass <- qc_criteria %>% 
  select(src_subject_id,iqc_dmri_ok_ser,iqc_t1_ok_ser,fsqc_qc) %>% 
  filter(iqc_dmri_ok_ser>0) %>% 
  filter(iqc_t1_ok_ser>0) %>% 
  filter(fsqc_qc==1) 
dti_ids_pass <- dti_ids_pass$src_subject_id

#Subcortical ====
abcd_smrip201 <- readRDS("abcd_smrip201.rds")

abcd_subcort <- abcd_smrip201 %>% 
  .[,grep("subject|smri_vol_scs_", colnames(.))] %>%  #select smri vol data
  .[,grep("subject|aa|amygdala|caudate|hpus|pallidum|putamen|tp|latventricles|cranial",colnames(.))] #select columns that overlap with enigma data

#sum up columns for left and right hemisphere and average
abcd_subcort_sum <- abcd_subcort %>% 
  mutate(Accumbens = as.character(rowSums(abcd_subcort %>% select(contains("aa")))/2),
         Amygdala = as.character(rowSums(abcd_subcort %>% select(contains("amygdala")))/2),
         Caudate = as.character(rowSums(abcd_subcort %>% select(contains("caudate")))/2),
         Hippocampus = as.character(rowSums(abcd_subcort %>% select(contains("hpus")))/2),
         Pallidum = as.character(rowSums(abcd_subcort %>% select(contains("pallidum")))/2),
         Putamen = as.character(rowSums(abcd_subcort %>% select(contains("putamen")))/2),
         Thalamus = as.character(rowSums(abcd_subcort %>% select(contains("tp")))/2),
         ICV = as.character(abcd_subcort$smri_vol_scs_intracranialv)) %>% 
  mutate_at(vars(Accumbens:ICV),as.numeric)
#glimpse(abcd_subcort_sum)

#select those that passed QC check and merge with control.status
abcd_subcort_final<- abcd_subcort_sum %>% 
  .[.$src_subject_id %in% id_keep$src_subject_id,] %>% #no neuro, unrel, whiteonly
  .[.$src_subject_id %in% t1_ids_pass,] %>%  #select indiv that passed QC
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#Cortical vol and SA (APARc ROI) ====

abcd_smrip101 <- readRDS(paste0(getwd(),"/abcd_smrip101.rds"))

#select cortical thickness - match to ENIGMA 
abcd_corthick <- abcd_smrip101 %>% 
  .[,grep("subject|thick_cdk", colnames(.))] %>% 
  .[,grep("subject|lh$|rh$", colnames(.))] #to choose hemispheric variables only 
colnames(abcd_corthick) <- sub("smri_thick_cdk_","",colnames(abcd_corthick))
#put lh and rh at the front instead 
colnames(abcd_corthick)[grep('lh$',colnames(abcd_corthick))]=paste0('lh.',gsub('lh$','',colnames(abcd_corthick)))[grep('lh$',colnames(abcd_corthick))]
colnames(abcd_corthick)[grep('rh$',colnames(abcd_corthick))]=paste0('rh.',gsub('rh$','',colnames(abcd_corthick)))[grep('rh$',colnames(abcd_corthick))]

#select cortical SA - match to ENIGMA 
abcd_corSA <- abcd_smrip101 %>% 
  .[,grep("subject|area_cdk", colnames(.))] %>% 
  .[,grep("subject|lh$|rh$", colnames(.))] #to choose hemispheric variables only 
colnames(abcd_corSA) <- sub("smri_area_cdk_","",colnames(abcd_corSA))
colnames(abcd_corSA)[grep('lh$',colnames(abcd_corSA))]=paste0('lh.',gsub('lh$','',colnames(abcd_corSA)))[grep('lh$',colnames(abcd_corSA))]
colnames(abcd_corSA)[grep('rh$',colnames(abcd_corSA))]=paste0('rh.',gsub('rh$','',colnames(abcd_corSA)))[grep('rh$',colnames(abcd_corSA))]

#select those that passed QC check 
abcd_corthick_final<- abcd_corthick %>% 
  .[.$src_subject_id %in% id_keep$src_subject_id,] %>% #no neuro, unrel, whiteonly
  .[.$src_subject_id %in% t1_ids_pass,] %>%  #select indiv that passed QC
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

abcd_corSA_final<- abcd_corSA %>% #para is skewed
  .[.$src_subject_id %in% id_keep$src_subject_id,] %>% #no neuro, unrel, whiteonly
  .[.$src_subject_id %in% t1_ids_pass,] %>%  #select indiv that passed QC
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#merge hemi for cortical measures 
avg_hemi_FUN <- function(dat){
  ls <- gsub("([a-z]{2}\\.)(.+)","\\2",colnames(dat %>% select(starts_with(c("lh","rh")))))
  names <- ls[!duplicated(ls)] #get unique names
  new <- as.data.frame(sapply(names,function(x) (rowSums(dat[,grep(paste0("\\b",x,"\\b"),colnames(dat))],na.rm=T)/2))) %>%
    cbind(dat[1],.,dat[,-grep("subject|^lh|^rh",colnames(dat))])
  return(new)}

abcd_corSA_final_avg <- avg_hemi_FUN(abcd_corSA_final)
abcd_corthick_final_avg <- avg_hemi_FUN(abcd_corthick_final)

#MD ====
abcd_dti_p101 <- readRDS("abcd_dti_p101.rds")
dti_def <- read_csv("abcd_dti_p101_definitions.csv") #subset to MD and FA only

#MD 
abcd_dti_MD<- abcd_dti_p101 %>% 
  .[,grep("subject|_dtimd_", colnames(.))] %>% 
  .[,grep('subject|fiberat',colnames(.))] %>% 
  .[,(colnames(.) %in% dti_def$ElementName[grep("subject|thalamic|cingulum|corticospinal|uncinate|fronto-occipital|foreceps|corpus|fasiculus|fornix|all", dti_def$ElementDescription, ignore.case = T)])]
#select columns that overlap with enigma data

#sum up columns for left and right hemisphere and average - match to ENIGMA for MD (only matched to half of the regions?)
abcd_dtiMD_sum <- abcd_dti_MD%>% 
  mutate(CC = as.character(abcd_dti_MD$dmri_dtimd_fiberat_cc),
         CGC = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_cgcrh, dmri_dtimd_fiberat_cgclh))/2),
         CGH = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_cghrh, dmri_dtimd_fiberat_cghlh))/2),
         CST = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_cstrh, dmri_dtimd_fiberat_cstlh))/2),
         FX = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_fxrh, dmri_dtimd_fiberat_fxlh))/2), 
         GCC = as.character(abcd_dti_MD$dmri_dtimd_fiberat_fmin), #genu map to foreceps minor
         SCC = as.character(abcd_dti_MD$dmri_dtimd_fiberat_fmaj), #splenium map to foreceps major
         ALIC = as.character(rowSums(abcd_dti_MD%>% select(dmri_dtimd_fiberat_atrrh, dmri_dtimd_fiberat_atrlh))/2), #map ALIC to anterior thalamic radiation
         IFO = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_iforh,dmri_dtimd_fiberat_ifolh))/2),
         SLF = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_slfrh, dmri_dtimd_fiberat_slflh))/2),
         UNC = as.character(rowSums(abcd_dti_MD %>% select(dmri_dtimd_fiberat_uncrh, dmri_dtimd_fiberat_unclh))/2),
         AvgMD = as.character(abcd_dti_MD$dmri_dtimd_fiberat_allfibers))
#glimpse(abcd_dtiMD_sum)

abcd_dtiMD_final <- abcd_dtiMD_sum %>% select(!grep("dmri_",colnames(.))) %>%
  .[.$src_subject_id %in% id_keep$src_subject_id,] %>% #no neuro, unrel, whiteonly
  .[.$src_subject_id %in% dti_ids_pass,] %>%  #select indiv that passed QC
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) %>% 
  mutate_at(vars(CC:AvgMD),as.numeric) 

#check dist
abcd_dtiMD_final %>%
  .[2:13] %>%
  gather() %>%
  ggplot(aes(as.numeric(value)))+
  facet_wrap(~ key, scales = "free") +
  geom_histogram(fill = "deepskyblue4") #very skewed

#remove outlier (for release2.01)
abcd_dtiMD_final <- filter(abcd_dtiMD_final,!is.na(AvgMD))
abcd_dtiMD_final[abs(scale(abcd_dtiMD_final$AvgMD))>5,]=NA #for all rows with avgMD >5sd, set NA to everything 
#check dist again -> looks fine

#FA ====
abcd_dti_FA<- abcd_dti_p101 %>% 
  .[,grep("subject|_dtifa_", colnames(.))] %>% 
  .[,grep('subject|fiberat',colnames(.))] %>% 
  .[,(colnames(.) %in% dti_def$ElementName[grep("subject|thalamic|cingulum|corticospinal|uncinate|fronto-occipital|foreceps|corpus|fasiculus|fornix|all", dti_def$ElementDescription, ignore.case = T)])]
#select columns that overlap with enigma data

#sum up columns for left and right hemisphere and average - match to ENIGMA for FA
abcd_dtiFA_sum <- abcd_dti_FA%>% 
  mutate(CC = as.character(abcd_dti_FA$dmri_dtifa_fiberat_cc),
         CGC = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_cgcrh, dmri_dtifa_fiberat_cgclh))/2),
         CGH = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_cghrh, dmri_dtifa_fiberat_cghlh))/2),
         CST = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_cstrh, dmri_dtifa_fiberat_cstlh))/2),
         FX = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_fxrh, dmri_dtifa_fiberat_fxlh))/2), 
         GCC = as.character(abcd_dti_FA$dmri_dtifa_fiberat_fmin), #genu map to foreceps minor
         SCC = as.character(abcd_dti_FA$dmri_dtifa_fiberat_fmaj), #splenium map to foreceps major
         ALIC = as.character(rowSums(abcd_dti_FA%>% select(dmri_dtifa_fiberat_atrrh, dmri_dtifa_fiberat_atrlh))/2), #map ALIC to anterior thalamic radiation
         IFO = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_iforh,dmri_dtifa_fiberat_ifolh))/2),
         SLF = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_slfrh, dmri_dtifa_fiberat_slflh))/2),
         UNC = as.character(rowSums(abcd_dti_FA %>% select(dmri_dtifa_fiberat_uncrh, dmri_dtifa_fiberat_unclh))/2),
         AvgFA = as.character(abcd_dti_FA$dmri_dtifa_fiberat_allfibers))
#glimpse(abcd_dtiFA_sum)

abcd_dtiFA_final <- abcd_dtiFA_sum %>% select(!grep("dmri_",colnames(.))) %>% 
  .[.$src_subject_id %in% id_keep$src_subject_id,] %>% #no neuro, unrel, whiteonly
  .[.$src_subject_id %in% dti_ids_pass,] %>%  #select indiv that passed QC
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) %>% 
  mutate_at(vars(CC:AvgFA),as.numeric) 

#check dist
# abcd.dtiFA.final %>%
#   .[2:13] %>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

#remove outliers (for release2.01)
abcd_dtiFA_final <- filter(abcd_dtiFA_final,!is.na(AvgFA))
abcd_dtiFA_final[abs(scale(abcd_dtiFA_final$AvgFA))>5,]=NA #for all rows with avgMD >5sd, set NA to everything 
#check dist again -> looks fine

#Compile all imaging measures ====
abcd_img_forRVI <- list(abcd_subcort = abcd_subcort_final, 
                        abcd_corthick = abcd_corthick_final_avg, 
                        abcd_corSA = abcd_corSA_final_avg,
                        abcd_dtiMD = abcd_dtiMD_final, 
                        abcd_dtiFA = abcd_dtiFA_final)

saveRDS(abcd_img_forRVI, "abcd201_img_forRVI.rds")

