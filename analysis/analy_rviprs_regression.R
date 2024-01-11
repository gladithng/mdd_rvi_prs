## STAT ANALYSIS ABCD ##

#script here shown for abcd baseline, repeat for abcd 2-year and GS (not shown)

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(dplyr);library(broom);library(pbapply);library(nlme)

#Load data ==== 

#rvi
rvi <- readRDS("abcd201_enigma_RVI.rds")

#prs
prs <- read.delim("ABCD3.0_imputed_whiteonlywparents_MDD-PRS_210727.all.score", header=T, sep="") %>% .[,-1] #did not use scoresum for PRS ->values all close to 0
colnames(prs)[1] <- "src_subject_id"
colnames(prs)[2:11] <- gsub("X", "pT_",colnames(prs)[2:11])

#covs
mdd_covs <- readRDS("abcd201_mdd_covs_allsubj.rds")

#combine - merging ensure that only subj with mdd, rvi and prs are included in the analysis
reg_dat <- merge(rvi, prs, by = "src_subject_id") %>% 
  merge(., mdd_covs, by = "src_subject_id")

#Regression ==== 

#define functions - glm for stradl and lme for abcd

glm_FUN <- function(ls, dat){
  
  run_model <- function(ls, dat){ 
    
    depvar <- as.character(ls[1])
    indpvar <- as.character(ls[2])
    cov <- as.character(ls[3])
    mod_type <- as.character(ls[4])
    
    if(mod_type == "logistic"){ #for categorical outcomes like lifetime MDD in gs
      eq <- paste0(depvar,'~',cov,'+','scale(',indpvar,')')
      fit <- glm(as.formula(eq), data = dat, na.action = na.exclude, family="binomial")
      table <- tidy(fit)
      tarv <- nrow(table)
      stats <- table[tarv,-1]
      mod_result <- data.frame(mod_name = paste0(depvar,'~', indpvar), stats)
    }
    
    else{ #for continuous outcomes like totalqids in gs
      eq <- paste0('scale(',depvar,')~', cov,'+ scale(', indpvar, ')')
      fit <- lm(as.formula(eq), data = dat, na.action = na.exclude)
      tb <- tidy(fit) 
      stats <- tail(tb,1)
      mod_result = data.frame(mod_name = paste0(depvar,'~',indpvar), stats[,c(2:5)])
    }
  }
    
    tmp_result<- pbapply(ls.models, 1, run_model, data = dat)
    result_table <- matrix(unlist(tmp_result), ncol = 5, byrow = T)
    result_table <- data.frame(result_table)
    colnames(result_table) <- c('mod_name','beta','sd','tval','pval')
    
    #run fdr correction
    ls.factor <- unique(ls.models$p_batch)
    result_table$p_corrected <- 99999
    for(f in ls.factor){
      loc <- grep(f,ls.models$p_batch) #get the row number of the row that correspond to ls.factor aka p_batch number
      result_table$p_corrected[loc] <- p.adjust(result_table$p.value[loc], method ='fdr')}
    
    return(result_table)
  }
    
mixedlme_FUN <- function(ls, dat){
  
  run_model <- function(ls,dat){
    
    depvar <- as.character(ls[1])
    indpvar <- as.character(ls[2])
    cov <- as.character(ls[3])

    eq <- paste0('scale(',depvar,')~',cov,'+scale(',indpvar,')')
    fit <- lme(as.formula(as.character(eq)), data = dat, na.action = na.exclude,
               random = ~1|site_id, method = "ML",control=lmeControl(opt = "optim"))
    
    #DOUBLE CHECK THIS PART TO MAKE SURE GLMER OUTPUT CAN BE ORGANISED THIS WAY#
    table <- summary(fit)$tTable
    tarv <- nrow(table)
    stats <- table[tarv,c(1,2,4,5)]
    mod_result = data.frame(mod_name = paste0(depvar,'~',indpvar),t(stats))
  }
  
  tmp_result<- pbapply(ls.models, 1, run_model, data = dat)
  result_table <- matrix(unlist(tmp_result), ncol = 5, byrow = T)
  result_table <- data.frame(result_table)
  colnames(result_table) <- c('mod_name','beta','sd','tval','pval')

  #run fdr correction
  ls.factor <- unique(ls.models$p_batch)
  result_table$p_corrected <- 99999
  for(f in ls.factor){
    loc <- grep(f,ls.models$p_batch) #get the row number of the row that correspond to ls.factor aka p_batch number
    result_table$p_corrected[loc] <- p.adjust(result_table$p.value[loc], method ='fdr')}
  
  return(result_table)
}

#specify ls - for abcd
ls.dep <- colnames(reg_dat)[grep("depress_r$", colnames(reg_dat))] #cbcl
ls.factor <- colnames(reg_dat)[grep("^rvi|^pT", colnames(reg_dat))] #rvi, prs
ls.models <- expand.grid(dep = ls.dep, indp = ls.factor, covs = NA, stringsAsFactors = F)
ls.models$covs[grep("^rvi", ls.models$indp)] <- paste0(c('sex','age','parent_edu','fam_income','I(age^2)'),collapse='+') 
ls.models$covs[grep("^pT", ls.models$indp)] <- colnames(reg_dat)[grep("age|sex|edu|income|PC_|plate",colnames(reg_dat))] %>% 
  paste0(.,collapse='+') %>% paste0('+I(age^2)')

ls.models$p_batch <- 99999
ls.models$p_batch[grep('rvi_subcort_adult',ls.models$factor)] <- 1
ls.models$p_batch[grep('rvi_corSA_adult',ls.models$factor)] <- 2
ls.models$p_batch[grep('rvi_corthick_adult',ls.models$factor)] <- 3
ls.models$p_batch[grep('rvi_dtiMD_adult',ls.models$factor)] <- 4
ls.models$p_batch[grep('rvi_dtiFA_adult',ls.models$factor)] <- 5
ls.models$p_batch[grep('rvi_wholebrain',ls.models$factor)] <- 7
ls.models$p_batch[grep('pT_0.001',ls.models$factor)] <- "a"
ls.models$p_batch[grep('pT_0.01',ls.models$factor)] <- "b"
ls.models$p_batch[grep('pT_0.05',ls.models$factor)] <- "c"
ls.models$p_batch[grep('pT_0.1',ls.models$factor)] <- "d"
ls.models$p_batch[grep('pT_0.5',ls.models$factor)] <- "e"
ls.models$p_batch[grep('pT_1',ls.models$factor)] <- "f"

#additionally, specify modtype for GS
# ls.models$mod[grep("lifetime", ls.models$dep)] <- "logistic" #for gs model
# ls.models$mod[grep("qids", ls.models$dep)] <- "cont" #for gs model

#run model 
lme_abcd201 <- mixedlme_FUN(ls.models, reg_dat) #site as a random factor
# glm_gs <- glm_FUN(ls.models, reg_dat)

saveRDS(reg_dat, "abcd201_reg_dat.rds")
saveRDS(lme_abcd201, "abcd201_lme_results.rds")