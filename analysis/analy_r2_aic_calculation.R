## DETERMINING VARIANCE EXPLAINED AND MODEL FIT ## 

#script here is for abcd baseline, same for abcd 2year and GS (not shown)

rm(list=ls()) 
options(scipen=999)
library(tidyverse);library(dplyr);library(broom);library(performance);
library(nlme);library(pbapply);library(MuMIn)

#Load data 
abcd201 <- readRDS("abcd201_reg_dat.rds")
# abcd3 <- readRDS("abcd3_reg_dat.rds")
# stradl <- readRDS("stradl_reg_dat.rds")

#Standardise first
abcd201_std <-  abcd201 %>% 
  mutate_at(vars(sex),as.factor) %>% 
  mutate(across(where(is.numeric), scale)) %>% 
  as.data.frame()

#Derive R2 and AIC 
#abcd: conditional/marginal r2 for lme 
#stradl: mcfadden pseudo r2 and normal r2
stradl_r2_AIC_FUN <- function(ls, dat){ 
  
  depvar <- ls[1]
  rvi <- ls[2]
  prs <- ls[3]
  covs <- ls[4]
  covs_rvi <- ls[5] #this doesnt include genetic covariates
  
  if(depvar == "lifetime"){
    full_mod <- glm(paste0(depvar,'~',paste0(covs,collapse="+"),'+',prs,'+',rvi),
                    data=dat,na.action=na.exclude,family="binomial")
    prs_mod <- glm(paste0(depvar,'~',paste0(covs,collapse="+"),'+',prs),
                   data=dat,na.action=na.exclude,family="binomial")
    rvi_mod <- glm(paste0(depvar,'~',paste0(covs_rvi,collapse="+"),'+',rvi),
                   data=dat,na.action=na.exclude,family="binomial")
    null_mod <- glm(paste0(depvar,'~',paste0(covs,collapse="+")),
                    data=dat,na.action=na.exclude,family="binomial")
    null_rvi_mod <- glm(paste0(depvar,'~',paste0(covs_rvi,collapse="+")),
                        data=dat,na.action=na.exclude,family="binomial")
    
    #calculate the pseudo r2 for each model using macfadden r2
    #note: mcfadden gives you the proportion of decreased/increased deviance relative to the null model as you add predictors 
    #not exactly comparable to the proportion of variance explained as in linear regression
    
    #adjusted - treat covs as null model 
    #to measure how much model improves from covs-only model after adding predictors
    k_rvi = 1; k_full = 2
    prs_adjr2 <- 1-((logLik(prs_mod)-2)/logLik(null_mod)) #2 cuz 1 prs variable + intercept
    rvi_adjr2 <- 1-((logLik(rvi_mod)-k_rvi-1)/logLik(null_rvi_mod)) #2 cuz 1 rvi variable + intercept
    full_adjr2 <- 1-((logLik(full_mod)-k_full-1)/logLik(null_mod))
    
    #non-adjusted
    prs_r2 <- 1-((logLik(prs_mod))/logLik(null_mod)) 
    rvi_r2 <- 1-((logLik(rvi_mod))/logLik(null_rvi_mod)) 
    full_r2 <- 1-((logLik(full_mod))/logLik(null_mod))
    
    #change in r2 rel to covs only model 
    #covs (all), covs(all)+prs, covs(all)+PRS+RVI
    #basically looking at the sequential addition of predictors 
    covs_dr2 <- logLik(null_mod)
    prs_dr2 <- 1-((logLik(prs_mod))/logLik(null_mod)) 
    full_dr2 <- 1-((logLik(full_mod))/logLik(prs_mod))
  }
  
  else{ #TotalQIDS
    full_mod <- lm(paste0(depvar,'~',paste0(covs,collapse="+"),'+',prs,'+',rvi),
                   data=dat, na.action=na.exclude)
    prs_mod <- lm(paste0(depvar,'~',paste0(covs,collapse="+"),'+',prs),
                  data=dat, na.action=na.exclude)
    rvi_mod <- lm(paste0(depvar,'~',paste0(covs_rvi,collapse="+"),'+',rvi),
                  data=dat, na.action=na.exclude)
    null_mod <- lm(paste0(depvar,'~',paste0(covs,collapse="+")),
                   data=dat, na.action=na.exclude)
    null_rvi_mod <- lm(paste0(depvar,'~',paste0(covs_rvi,collapse="+")),
                       data=dat, na.action=na.exclude)
    
    #adjusted r2 using full-null
    prs_adjr2 <- r2(prs_mod)[[2]]-r2(null_mod)[[2]] 
    rvi_adjr2 <- r2(rvi_mod)[[2]]-r2(null_rvi_mod)[[2]]
    full_adjr2 <- r2(full_mod)[[2]]-r2(null_mod)[[2]]
    
    #non-adjusted r2 using full-null
    prs_r2 <- r2(prs_mod)[[1]]-r2(null_mod)[[1]] 
    rvi_r2 <- r2(rvi_mod)[[1]]-r2(null_rvi_mod)[[1]]
    full_r2 <- r2(full_mod)[[1]]-r2(null_mod)[[1]]
    
    #change in r2 rel to covs only model 
    #covs (all), covs(all)+prs, covs(all)+PRS+RVI
    covs_dr2 <- r2(null_mod)[[1]]
    prs_dr2 <- r2(prs_mod)[[1]]-r2(null_mod)[[1]]
    full_dr2 <- r2(full_mod)[[1]]-r2(prs_mod)[[1]]
  }
  
  #combine results into df
  results_adjr2 <- data.frame(full_model = full_adjr2,
                              rvi_model = rvi_adjr2,
                              prs_model = prs_adjr2,
                              null_model = NA,
                              null_rvi_model = NA, 
                              parameter = "adjR2")
  
  results_r2 <- data.frame(full_model = full_r2,
                           rvi_model = rvi_r2,
                           prs_model = prs_r2,
                           null_model = NA,
                           null_rvi_model = NA, 
                           parameter = "R2")
  
  results_dr2 <- data.frame(full_model = full_dr2,
                            rvi_model = NA,
                            prs_model = prs_dr2,
                            null_model = covs_dr2,
                            null_rvi_model = NA, 
                            parameter = "dR2")
  
  results.AIC <- data.frame(full_model = AIC(full_mod),
                            rvi_model = AIC(rvi_mod),
                            prs_model = AIC(prs_mod),
                            null_model = AIC(null_mod),
                            null_rvi_model = AIC(null_rvi_mod),
                            parameter = "AIC")
  
  results <- rbind(results_adjr2, results_r2, results_dr2,results_AIC)
  results$depvar <- depvar
  results$rvi <- rvi
  results <- results %>% select(depvar,rvi,parameter,everything())
  
  return(results)
}

abcd_r2_AIC_FUN <- function(ls, dat){
  
  depvar <- ls[1]
  rvi <- ls[2]
  prs <- ls[3]
  covs <- ls[4]
  covs_rvi <- ls[5] #this doesnt include genetic covariates
  
  full_mod <- lme(as.formula(paste0(depvar,'~',covs,'+',prs,'+',rvi)),data = dat, 
                  na.action = na.exclude, method = "ML", random = ~1|site_id)
  prs_mod <- lme(as.formula(paste0(depvar,'~',covs,'+',prs)), data = dat, 
                 na.action = na.exclude, method = "ML", random = ~1|site_id)
  rvi_mod <- lme(as.formula(paste0(depvar,'~',covs_rvi,'+',rvi)), data = dat, 
                 na.action = na.exclude, method = "ML", random = ~1|site_id)
  null_mod <- lme(as.formula(paste0(depvar,'~',covs)),data = dat, 
                  na.action = na.exclude, method = "ML", random = ~1|site_id)
  null_rvi_mod <- lme(as.formula(paste0(depvar,'~',covs_rvi)), data = dat, 
                      na.action = na.exclude, method = "ML", random = ~1|site_id)
  
  #use conditional r2 to look at var of random and fixed effects 
  prs_cr2 <- MuMIn::r.squaredGLMM(prs_mod)[[2]]-r.squaredGLMM(null_mod)[[2]] 
  rvi_cr2 <- r.squaredGLMM(rvi_mod)[[2]]-r.squaredGLMM(null_rvi_mod)[[2]]
  full_cr2 <- r.squaredGLMM(full_mod)[[2]]-r.squaredGLMM(null_mod)[[2]]
  
  #use marginal r2 to look at var of fixed effects (preferred)
  prs_mr2 <- MuMIn::r.squaredGLMM(prs_mod)[[1]]-r.squaredGLMM(null_mod)[[1]] 
  rvi_mr2 <- r.squaredGLMM(rvi_mod)[[1]]-r.squaredGLMM(null_rvi_mod)[[1]]
  full_mr2 <- r.squaredGLMM(full_mod)[[1]]-r.squaredGLMM(null_mod)[[1]]
  
  #change in r2 (dr2: delta r2) relative to covs only model 
  #covs (all), covs(all)+prs, covs(all)+PRS+RVI
  covs_dr2 <- r.squaredGLMM(null_mod)[[2]]
  prs_dr2 <- MuMIn::r.squaredGLMM(prs_mod)[[2]]-r.squaredGLMM(null_mod)[[2]] 
  full_dr2 <- r.squaredGLMM(full_mod)[[2]]-r.squaredGLMM(prs_mod)[[2]]
  
  #combine results into df
  results_cr2 <- data.frame(full_model = full_cr2,
                              rvi_model = rvi_cr2,
                              prs_model = prs_cr2,
                              null_model = NA,
                              null_rvi_model = NA, 
                              parameter = "cR2")
  
  results_mr2 <- data.frame(full_model = full_mr2,
                            rvi_model = rvi_mr2,
                            prs_model = prs_mr2,
                            null_model = NA,
                            null_rvi_model = NA, 
                            parameter = "mR2")
  
  results_dr2 <- data.frame(full_model = full_dr2,
                            rvi_model = NA,
                            prs_model = prs_dr2,
                            null_model = covs_dr2,
                            null_rvi_model = NA, 
                            parameter = "dR2")
  
  results_AIC <- data.frame(full_model = AIC(full_mod), #can also use model_performance(full_mod)
                            rvi_model = AIC(rvi_mod),
                            prs_model = AIC(prs_mod),
                            null_model = AIC(null_mod),
                            null_rvi_model = AIC(null_rvi_mod),
                            parameter = "AIC")
  
  results <- rbind(results_cr2, results_mr2, results_dr2,results_AIC)
  results$depvar <- depvar
  results$rvi <- rvi
  results <- results %>% select(depvar,rvi,parameter,everything())
  
  return(results)
} 

ls.rvi <- colnames(abcd201_std)[grep("^rvi", colnames(abcd201_std))]
ls.prs <- "p_0.1" #selected threshold
ls.dep <- "cbcl_scr_dsm5_depress_r"
ls_models <- expand.grid(dep = ls.dep, rvi = ls.rvi,prs = ls.prs,covs = " ", stringsAsFactors = F)
ls_models$covs_rvi <- paste0(c("age","sex","I(age^2)", "parent_edu","fam_income"),collapse='+')
ls_models$covs <- colnames(abcd201_std)[grep("age|sex|edu|income|PC_|plate",colnames(abcd201_std))] %>% 
  paste0(.,collapse='+') %>% paste0('+I(age^2)')

tmp <- pbapply(ls_models, 1, abcd_r2_AIC_FUN, dat = abcd201_std)
abcd201_results <- as.data.frame(do.call(rbind, lapply(tmp, as.data.frame)))


