## PREP IMAGING DATA FOR GS ##

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(dplyr);library(broom);library(pbapply);library(nlme)

#Load necessary data
id_keep <- readRDS("stradl_noneuro_unrel_ids.rds")
disease_antidep <- readRDS("stradl_disease_antidep_allsubj.rds") 
covs <- readRDS("stradl_mdd_covs_prs_allsubj.rds")

#Subcortical====
stradl_FS_main <- read_csv("STRADL_Measures_FreeSurfer_Main.csv")

stradl_subcort <- stradl_FS_main %>% .[grep("^id|scv.bilat",colnames(.))] #select bilateral 
stradl_subcort <- stradl_subcort[,1:ncol(stradl_subcort)-1] #drop ventraldc cuz no matching in enigma
colnames(stradl_subcort) <- gsub('scv.bilat.','',colnames(stradl_subcort))
colnames(stradl_subcort) <- c("ID","Accumbens","Amygdala","Caudate","Hippocampus","Pallidum","Putamen","Thalamus")

# stradl_subcort %>%
#   select(!starts_with(c("id","site")))%>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

stradl_subcort_final <- stradl_subcort %>% 
  .[.$ID %in% ids_keep$ID == T,] %>% #no neuro, unrelated
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#cortical SA ====
stradl_corSA <- stradl_FS_main %>% .[grep("^id|hem.lh.csa|hem.rh.csa|csa.bilat",colnames(.))]
colnames(stradl_corSA) <- gsub('csa.bilat.','',colnames(stradl_corSA))
colnames(stradl_corSA) <- gsub('hem.','',colnames(stradl_corSA))
colnames(stradl_corSA)[1] <- "ID"

avg_hemi_FUN <- function(dat){
  ls <- gsub("([a-z]{2}\\.)(.+)","\\2",colnames(dat %>% select(starts_with(c("lh","rh")))))
  names <- ls[!duplicated(ls)] #get unique names
  new <- as.data.frame(sapply(names,function(x) (rowSums(dat[,grep(paste0("\\b",x,"\\b"),colnames(dat))],na.rm=T)/2))) %>%
    cbind(dat[,-grep("subject|^lh|^rh",colnames(dat))],.)
  return(new)}

stradl_corSA_avg <- avg_hemi_FUN(stradl_corSA) #to get the avg SA for both hemisphere
colnames(stradl_corSA_avg)[36] <- "AvgCSA"

# stradl_corSA %>%
#   select(!starts_with(c("id","site")))%>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

stradl_corSA_final <- stradl_corSA_avg %>% 
  .[.$ID %in% ids_keep$ID == T,] %>% #no neuro, unrelated
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#cortical thickness ====
stradl_corthick <- stradl_FS_main %>% .[grep("^id|hem.lh.ct|hem.rh.ct|ct.bilat",colnames(.))]
colnames(stradl_corthick) <- gsub('ct.bilat.','',colnames(stradl_corthick))
colnames(stradl_corthick) <- gsub('hem.','',colnames(stradl_corthick))
colnames(stradl_corthick)[1] <- "ID"

stradl_corthick_avg <- avg_hemi_FUN(stradl_corthick) #to get the avg SA for both hemisphere
colnames(stradl_corthick_avg)[36] <- "AvgCT"

# stradl_corthick %>%
#   select(!starts_with(c("id","site")))%>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

stradl_corthick_final <- stradl_corthick_avg %>% 
  .[.$ID %in% ids_keep$ID == T,] %>% #no neuro, unrelated
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#dtiMD ====
stradl_dtiMD <- read_csv("STRADL_Measures_DTI_MD.csv")
stradl_dtiMD <- stradl_dtiMD %>% (function(x) x[,!grepl('\\-L$|\\-R$',colnames(x))])
colnames(stradl_dtiMD) <- gsub('AverageFA','AvgMD',colnames(stradl_dtiMD))
colnames(stradl_dtiMD)[1] <- "ID"

# stradl_dtiMD %>%
#   select(!starts_with(c("id","site")))%>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

stradl_dtiMD_final <- stradl_dtiMD_avg %>% 
  .[.$ID %in% ids_keep$ID == T,] %>% #no neuro, unrelated
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#dtiFA ==== 
stradl_dtiFA <- read_csv("STRADL_Measures_DTI_FA.csv")
stradl_dtiFA <- stradl_dtiFA %>% (function(x) x[,!grepl('\\-L$|\\-R$',colnames(x))])
colnames(stradl_dtiFA)=gsub('AverageFA','AvgFA',colnames(stradl_dtiFA))
colnames(stradl_dtiFA)[1] <- "ID"

# stradl_dtiFA %>%
#   select(!starts_with(c("id","site")))%>%
#   gather() %>%
#   ggplot(aes(as.numeric(value)))+
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram(fill = "deepskyblue4")

stradl_dtiFA_final <- stradl_dtiFA_avg %>% 
  .[.$ID %in% ids_keep$ID == T,] %>% #no neuro, unrelated
  left_join(.,disease_antidep) %>% 
  left_join(.,covs) 

#Compile all imaging measures ====
stradl_img_forRVI <- list(stradl_subcort = stradl_subcort_final, 
                        stradl_corthick = stradl_corthick_final_avg, 
                        stradl_corSA = stradl_corSA_final_avg,
                        stradl_dtiMD = stradl_dtiMD_final, 
                        stradl_dtiFA = stradl_dtiFA_final)

saveRDS(stradl_img_forRVI, "stradl_img_forRVI.rds")

#ENIGMA ====

#define functions
source("stradl_rename_enigma_FUN.R")

avg_hemi_ENIGMA_FUN <- function(dat){
  ls <- gsub("([a-z]{2}\\.)(.+)","\\2",dat$V1)
  names <- ls[!duplicated(ls)] %>% as.vector() #get unique names
  new <- pivot_wider(data,names_from=V1, values_from=V2)
  new.df <- as.data.frame(sapply(names,function(x) (rowSums(new[,grep(paste0("\\b",x,"\\b"),colnames(new))],na.rm=T)/2))) %>% 
    rownames_to_column()
  colnames(new.df) <- c("V1","V2")
  return(new.df)
}

#subcortical
enigma_subcort_adult <- read.delim("enigma_adult_subcort_cohend.txt",header=F,stringsAsFactors=F) 
enigma_subcort_adult <- enigma_subcort_adult[2:8,] #remove latven cuz no matching in stradl, remove ICV too cuz not "regional"?

#corSA
enigma_corSA_adult <- read.delim("enigma_adult_corSA_cohend.txt",header=F,stringsAsFactors=F) 
#enigma_corSA_adult$V2 <- gsub("?^'","-", enigma_corSA_adult$V2)
head(enigma_corSA_adult)

enigma_corSA_adult <- renameFUN(enigma_corSA_adult) %>% 
  select(annotations, V2) %>% 
  rename(., V1 = annotations) %>% 
  mutate_at(colnames(.)[2], as.numeric) %>% 
  drop_na(.) 

enigma_corSA_adult_avg <- avg_hemi_enigma_FUN(enigma_corSA_adult)

#corthick
enigma_corthick_adult <- read.delim("enigma_adult_corticalthick_cohend.txt",header=F,stringsAsFactors=F) 
enigma_corthick_adult$V2 <- gsub("?^'","-", enigma_corthick_adult$V2)
head(enigma_corthick_adult) 

enigma_corthick_adult <- renameFUN(enigma_corthick_adult) %>% 
  select(annotations, V2) %>% 
  rename(., V1 = annotations) %>% 
  mutate_at(colnames(.)[2], as.numeric) %>% 
  drop_na(.)

enigma_corthick_adult_avg <- avg_hemi_enigma_FUN(enigma_corthick_adult)

#dtiMD
enigma_MD_adult <- read.delim("enigma_dtiMD_adult.txt",header=F,stringsAsFactors=F) #txt file without colnames
head(enigma_MD_adult) 

#dtiFA
enigma_FA_adult <- read.delim("enigma_dtiFA_adult.txt",header=F,stringsAsFactors=F) #txt file without colnames
head(enigma_MD_adult) 

#Compile and save ====
enigma_stradl <- list(eng_subcort = enigma_subcort_adult, 
                  eng_corthick = enigma_corthick_adult_avg, 
                  eng_corSA = enigma_corSA_adult_avg, 
                  eng_dtiMD = enigma_MD_adult, 
                  eng_dtiFA = enigma_FA_adult)

saveRDS(enigma_stradl, "enigma_stradl.rds")
