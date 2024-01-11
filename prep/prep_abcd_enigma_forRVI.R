## PREP ENIGMA DATA FOR RVI ##

library(tidyverse);library(dplyr);library(broom)


#Load abcd imaging
abcd201_img <- readRDS("abcd201_img_forRVI.rds") #list
abcd_subcort <- abcd201_img[["abcd_subcort"]]
abcd_corthick <- abcd201_img[["abcd_corthick"]]
abcd_corSA <- abcd_img[["abcd_corSA"]]
abcd_dtiMD <- abcd_img[["abcd_dtiMD"]]
abcd_dtiFA <- abcd_img[["abcd_dtiFA"]]

#Load rename function 
#rename enigma data to be the same as abcd - is there a more effective way????
source("abcd_rename_enigma_FUN.R")

#Subcortical ====
#https://pubmed.ncbi.nlm.nih.gov/26122586/ 
enigma_subcort_adult <- read.delim("enigma_adult_subcort_cohend.txt",header=F,stringsAsFactors=F) #txt file without colnames
enigma_subcort_adult <- enigma_subcort_adult[2:8,] #remove latven cuz no matching in stradl, remove ICV too cuz not "regional"?

#Cortical ====
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5444023/
enigma_corthick_adult <- read.delim("enigma_adult_corticalthick_cohend.txt", header=F,stringsAsFactors=F) 
enigma_corthick_adult$V2 <- gsub("?^'","-", enigma_corthick_adult$V2)
head(enigma_corthick_adult) 

enigma_corthick_adult <- renameFUN(enigma_corthick_adult) %>% 
  select(annotations, V2) %>% 
  rename(., V1 = annotations) %>% 
  mutate_at(colnames(.)[2], as.numeric)

enigma_corSA_adult <- read.delim("enigma_adult_corSA_cohend.txt", header=F,stringsAsFactors=F) 
enigma_corSA_adult$V2 <- gsub("?^'","-", enigma_corSA_adult$V2)
head(enigma_corSA_adult)

enigma_corSA_adult <- renameFUN(enigma_corSA_adult) %>% 
  select(annotations, V2) %>% 
  rename(., V1 = annotations)%>% 
  mutate_at(colnames(.)[2], as.numeric)

#average values across both hemispheres
avg_hemi_enigma_FUN <- function(dat){
  ls <- gsub("([a-z]{2}\\.)(.+)","\\2", dat$V1)
  names <- ls[!duplicated(ls)] %>% as.vector() #get unique names
  new <- pivot_wider(data,names_from=V1, values_from=V2)
  new_df <- as.data.frame(sapply(names,function(x) (rowSums(new[,grep(paste0("\\b",x,"\\b"),colnames(new))],na.rm=T)/2))) %>% 
    rownames_to_column()
  colnames(new_df) <- c("V1","V2")
  return(new_df)
}

enigma_corSA_adult_avg <- avg_hemi_enigma_FUN(enigma_corSA_adult)
enigma_corthick_adult_avg <- avg_hemi_enigma_FUN(enigma_corthick_adult)

#DTI ====

#MD
enigma_MD_adult <- read.delim("enigma_dtiMD_adult.txt",header=F,stringsAsFactors=F) #txt file without colnames
enigma_MD_adult <- enigma_MD_adult[enigma_MD_adult$V1 %in% colnames(abcd_dtiMD),] #subset columns common to enigma and abcd

#FA
enigma_FA_adult <- read.delim("enigma_dtiFA_adult.txt",header=F,stringsAsFactors=F) #txt file without colnames
enigma_FA_adult <- enigma_FA_adult[enigma_FA_adult$V1 %in% colnames(abcd_dtiFA),] #subset columns common to enigma and abcd

#Compile and save ====
enigma201 <- list(eng_subcort = enigma_subcort_adult, 
                  eng_corthick = enigma_corthick_adult_avg, 
                  eng_corSA = enigma_corSA_adult_avg, 
                  eng_dtiMD = enigma_MD_adult, 
                  eng_dtiFA = enigma_FA_adult)

saveRDS(enigma201, "enigma_abcd201.rds")
