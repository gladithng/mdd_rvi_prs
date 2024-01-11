## CALCULATE RVI FOR ABCD ##

#script here shown for abcd baseline, repeat for abcd 2-year and GS (not shown)

library(dplyr);library(tidyverse);library(broom);library(RVIpkg)

#Load data ====
#abcd
abcd_img <- readRDS("abcd201_img_forrvi_rds") #list
enigma <- readRDS("enigma_abcd201.rds")

#RVI subcort (without latvent) ====
abcd_dat <- abcd_img[["abcd_subcort"]] 
img <- abcd_dat %>% 
  select(Accumbens:Thalamus) %>% 
  .[,order(colnames(.))] %>% #to organise according to the order in dataframe
  cbind(.,abcd_dat %>% select(src_subject_id, age, sex, site_id, ICV, psych_status, antidepressant_status)) %>% 
  mutate_at(vars(Accumbens:Thalamus),as.numeric) #IMPT!!!*

eng_dat <- enigma[["eng_subcort"]] %>% 
  .[order(.$V1),]
ref_beta <- eng_dat$V2

control_id <- "psych_status == 0 & antidepressant_status == 0"
cov <- c("age","sex","site_id","ICV")

calc <- RVI_func(ID = "src_subject_id",
                 DXcontrol = control_id,
                 covariates = cov,
                 resp.range = c(1:7),
                 EP = ref_beta,
                 data = drop_na(img))

rvi_subcort_adult <- as.data.frame(calc[["RVI"]][1:2])
colnames(rvi_subcort_adult)[2] <- "rvi_subcort_adult" 

#corSA
abcd_dat <- abcd_img[["abcd_corSA"]] 
img <- abcd_dat %>% 
  select(bankssts:total) %>% 
  .[,order(colnames(.))] %>% #to organise according to the order in dataframe 
  cbind(.,abcd_dat %>% select(src_subject_id,age,sex,site_id,psych_status, antidepressant_status))  %>% 
  mutate_at(vars(bankssts:trvtm),as.numeric) #IMPT!!!* 

eng_dat <- enigma[["eng_corSA"]] %>% 
  .[order(.$V1),]
ref_beta <- eng_dat$V2

control_id <- "psych_status == 0 & antidepressant_status == 0"
cov <- c("age","sex","site_id")

calc <- RVI_func(ID="src_subject_id",
                 DXcontrol = control_id,
                 covariates = cov,
                 resp.range = c(1:35),
                 EP = ref_beta,
                 data = drop_na(img))

rvi_corSA_adult <- as.data.frame(calc[["RVI"]][1:2])
colnames(rvi_corSA_adult)[2] <- "rvi_corSA_adult"

#corthick
abcd_dat <- abcd_img[["abcd_corthick"]]
img <- abcd_dat %>% 
  select(bankssts:mean) %>% 
  .[,order(colnames(.))] %>% #to organise according to the order in dataframe 
  cbind(.,abcd_dat %>% select(src_subject_id,age,sex,site_id,psych_status, antidepressant_status)) %>% 
  mutate_at(vars(bankssts:trvtm),as.numeric) #IMPT!!!* 

eng_dat <- enigma[["eng_corthick"]] %>% 
  .[order(.$V1),]
ref_beta <- eng_dat$V2

control_id <- "psych.status ==0 & antidepressant.status ==0"
cov <- c("age","sex","site_id")

calc <- RVI_func(ID="src_subject_id",
                 DXcontrol = control_id,
                 covariates = cov,
                 resp.range = c(1:35),
                 EP = ref_beta,
                 data = drop_na(img))

rvi_corthick_adult <- as.data.frame(calc[["RVI"]][1:2])
colnames(rvi_corthick_adult)[2] <- "rvi_corthick_adult"

#dtiMD
abcd_dat <- abcd_img[["abcd_dtiMD"]]
img <- abcd_dat %>% 
  select(CC:AvgMD) %>%
  .[,order(colnames(.))] %>% #to organise according to the order in dataframe 
  cbind(.,abcd_dat %>% select(src_subject_id,age,sex,site_id,psych_status, antidepressant_status)) %>% 
  mutate_at(vars(ALIC:UNC),as.numeric) #IMPT!!!* 

eng_dat <- enigma[["eng_dtiMD"]] %>% 
  .[order(.$V1),]
ref_beta <- eng_dat$V2

control_id <- "psych.status ==0 & antidepressant.status ==0"
cov <- c("age","sex","site_id","I(age^2)","age*sex","I(age^2)*sex") #enigma controlled for age, sex, agexsex, age2, age2xsex and scansite 

calc <- RVI_func(ID="src_subject_id",
                 DXcontrol = control_id,
                 covariates = cov,
                 resp.range = c(1:12),
                 EP = ref_beta,
                 data = drop_na(img))

rvi_dtiMD_adult <- as.data.frame(calc[["RVI"]][1:2])
colnames(rvi_dtiMD_adult)[2] <- "rvi_dtiMD_adult"

#dtiFA
abcd_dat <- abcd_img[["abcd_dtiFA"]]
img <- abcd_dat %>% 
  select(CC:AvgFA) %>% 
  .[,order(colnames(.))] %>% #to organise according to the order in dataframe 
  cbind(.,abcd_dat %>% select(src_subject_id,age,sex,site_id,psych_status, antidepressant.status)) %>% 
  mutate_at(vars(ALIC:UNC),as.numeric) #IMPT!!!* 

eng_dat <- enigma[["eng_dtiFA"]] %>% 
  .[order(.$V1),]
ref_beta <- eng_dat$V2

control_id <- "psych.status ==0 & antidepressant.status ==0"
cov <- c("age","sex","site_id","I(age^2)","age*sex","I(age^2)*sex") #enigma controlled for age, sex, agexsex, age2, age2xsex and scansite 

calc <- RVI_func(ID="src_subject_id",
                 DXcontrol = control_id,
                 covariates = cov,
                 resp.range = c(1:12),
                 EP = ref_beta,
                 data = drop_na(img))

rvi_dtiFA_adult <- as.data.frame(calc[["RVI"]][1:2])
colnames(rvi_dtiFA_adult)[2] <- "rvi_dtiFA_adult"                                  


#Combine RVI scores----
rvi_combined <- full_join(rvi_subcort_adult, rvi_corthick_adult) %>% 
  full_join(., rvi_corSA_adult) %>% 
  full_join(., rvi_dtiMD_adult) %>% 
  full_join(., rvi_dtiFA_adult) 
colnames(rvi_combined)[1] <- "src_subject_id"

rvi_combined <- rvi_combined %>% 
  mutate(rvi_wholebrain = (rowSums(rvi_combined %>% select(rvi_subcort_adult,rvi_corthick_adult,rvi_corSA_adult,rvi_dtiFA_adult,rvi_dtiMD_adult)))/5)

#check dist
# rvi_combined %>%
#   select(-src_subject_id) %>% 
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

#Save RDS ----
saveRDS(rvi_combined,"abcd201_enigma_rvi_rds")

