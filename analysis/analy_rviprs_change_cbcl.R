## Test if change in RVI predicts a change in depressive symptoms in ABCD ##

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(dplyr);library(broom);library(pbapply);library(nlme)

#Load data 
abcd201 <- readRDS("abcd201_reg_dat.rds")
colnames(abcd201)[2:ncol(abcd201)] <- paste0("baseline_",colnames(abcd201)[2:ncol(abcd201)])

abcd3 <- readRDS("abcd3_reg_dat.rds")
colnames(abcd3)[2:ncol(abcd3)] <- paste0("twoyear_",colnames(abcd3)[2:ncol(abcd3)])

comb_dat <- inner_join(abcd201, abcd3) %>% 
  mutate(years_gap = twoyear_age - baseline_age)
comb_dat$years_gap <- as.numeric(comb_dat$years_gap)

#residualise 2year cbcl score 
fit <- lme(twoyear_cbcl_scr_dsm5_depress_r ~ baseline_cbcl_scr_dsm5_depress_r +
             years_gap + baseline_sex + baseline_parent_edu + baseline_fam_income, 
           random = ~1|baseline_site_id, method = "ML",control=lmeControl(opt = "optim"),
           data = comb_dat, na.action = na.exclude) 

resid <- data.frame(src_subject_id = comb_dat$src_subject_id, resid = resid(fit)) 
colnames(resid)[2] <- "cbcl_depress_resid"

#combine predictors and resid dat
#if you regressed site at the start using mixed model, you can just run a general linear model using the residuals
glm_dat <- left_join(resid, comb_dat[, grep("src_subject_id|baseline_rvi|baseline_PC|baseline_pT_|genotype_plate",colnames(comb_dat))])

#regress 2year cbcl residuals on baseline predictors 
ls.dep <- colnames(glm_dat)[grep("cbcl",colnames(glm_dat))] 
ls.factor <- colnames(glm_dat)[grep("rvi|pT",colnames(glm_dat))]
ls.models <- expand.grid(dep = ls.dep, indp = ls.factor, covs = "", stringsAsFactors = F)
ls.models$covs[grep('baseline_pT_',ls.models$indp)] = colnames(glm_dat)[grep("PC_|genotype_plate_no",colnames(glm_dat))] %>% 
  paste0(.,collapse='+')
ls.models$mod <- "cont"

ls.models$p_batch <- 99999
ls.models$p_batch[grep('rvi_subcort_adult',ls.models$factor)] <- 1
ls.models$p_batch[grep('rvi_corSA_adult',ls.models$factor)] <- 2
ls.models$p_batch[grep('rvi_corthick_adult',ls.models$factor)] <- 3
ls.models$p_batch[grep('rvi_dtiMD_adult',ls.models$factor)] <- 4
ls.models$p_batch[grep('rvi_dtiFA_adult',ls.models$factor)] <- 5
ls.models$p_batch[grep('rvi_wholebrain',ls.models$factor)] <- 6
ls.models$p_batch[grep('pT_0.001',ls.models$factor)] <- "a"
ls.models$p_batch[grep('pT_0.01',ls.models$factor)] <- "b"
ls.models$p_batch[grep('pT_0.05',ls.models$factor)] <- "c"
ls.models$p_batch[grep('pT_0.1',ls.models$factor)] <- "d"
ls.models$p_batch[grep('pT_0.5',ls.models$factor)] <- "f"
ls.models$p_batch[grep('pT_1',ls.models$factor)] <- "g"

#run regression
source("glm_FUN.R")
glm_changeDS <- glm_FUN(ls.models, glm_dat)

