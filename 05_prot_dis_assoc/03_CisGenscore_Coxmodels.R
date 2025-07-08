#!/usr/bin/env Rscript
#######################################################
#### Cox models with proteins for protein var proj ####
####       KD, for Maik Pietzner   04/04/2024      ####
#######################################################

rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

## correct directory
setwd("/sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction")

## packages needed
require(data.table)
require(doMC)
require(rms)
require(survival)
require(tidyverse)
require(arrow)
require(magrittr)

## import phecode, associated cohort (sex-specific?)
phe    <- args[1]
cohort <- args[2]
# phe    <- "250.2"
# cohort <- "Both"

cat("run Cox-models with", phe, "in", cohort)

#----------------------------#
##-- import relevant data --##
#----------------------------#

## --> define covariates to be imputed <-- ##
covs            <- c("f.21022.0.0", "f.31.0.0")
## covariates (adopt to new covariate data set)
ukb.cov         <- read_parquet("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                            col_select = c("f.eid",covs)) %>% set_colnames(c("f.eid","age","sex")) %>% na.omit()
covs            <- c("age", "sex")

ukb.cov$sex <- as.character(ukb.cov$sex)

## decide whether or not the outcome is sex-specific
if(cohort != "Both"){
  ukb.cov <- ukb.cov %>% filter(sex == cohort)
  ## formula for adjustment
  covs    <- covs[-which(covs == "sex")]
}

ukb.cov <- as.data.table(ukb.cov)

## convert to input format for glmnet
if(length(covs) != 1){
  tmp             <- model.matrix.lm(as.formula(paste0("~ ", paste(covs, collapse = " + "))), ukb.cov, na.action = "na.omit")
  ## add
  ukb.cov         <- data.table(ukb.cov[, "f.eid", with=F], tmp[,-1])
  ## edit names
  names(ukb.cov)  <- gsub("[[:punct:]]| ", "_", names(ukb.cov))
  names(ukb.cov)  <- gsub("f_eid", "f.eid", names(ukb.cov))
  ## rename covs
  covs            <- names(ukb.cov)[-1]
}

## --> import protein label <-- ##
## import
lab.prot                <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/51_protein_variance_updated/02_feature_selection/input/Protein.variance.summary.UKB.Olink.20250612.txt")

## --> import protein levels cis genetic scores <-- ##
ukb.prot <- fread("input/UKB.olink.proteins.cis.genetic_scores.20250612.txt")

ukb.prot <- ukb.prot %>% 
  select(c(IID, starts_with("cis.score")))

## rename to fit to protein names
names(ukb.prot) <- sub("^cis.score.", "", names(ukb.prot))

## get new feature vector (differs for genetic scores)
prots.exp       <- intersect(lab.prot$id, names(ukb.prot))

## --> relevant phecode <-- ##
## update to new set of phecodes
ukb.phe        <- fread("input/UKB.olink.phecodes.fullcohort.20250516.txt", select = c("f.eid", paste0("bin_", phe), paste0("surv_", phe)))
## import baseline date
tmp.dat        <- read_parquet("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                               col_select = c("f.eid", "f.53.0.0"))
names(tmp.dat) <- c("f.eid", "baseline_date")
## add
ukb.phe        <- merge(ukb.phe, tmp.dat)

## --> combine everything <-- ##
## merge
ukb.dat        <- merge(ukb.cov, ukb.phe, by="f.eid")
ukb.dat        <- merge(ukb.dat, ukb.prot, by.x = "f.eid", by.y = "IID")
## delete what is no longer needed
rm(ukb.cov); rm(ukb.prot); rm(ukb.phe); gc(reset=T)

## --> drop prevalent/early cases <-- ##
## drop prevalent cases
ukb.dat        <- ukb.dat[ !is.na(eval(as.name(paste0("bin_", phe)))) ]

## compute follow-up time (compute as years)
ukb.dat[, fol := (eval(as.name(paste0("surv_", phe))) - as.IDate(baseline_date))/365.25]
ukb.dat$fol <- pmin(ukb.dat$fol, 10)
summary(ukb.dat[eval(as.name(paste0("bin_", phe))) == 1]$fol)
## drop cases within first six months
ukb.dat        <- ukb.dat[ !(eval(as.name(paste0("bin_", phe))) == 1 & fol <= .5)]


#-----------------------------------------#
##-- proceed only if at least 200 cases --##
#-----------------------------------------#
if(sum(ukb.dat[, paste0("bin_", phe), with=F]) >= 200){
  
  ## --> rank inverse normal transform proteins <-- ##
  ukb.dat              <- as.data.frame(ukb.dat)
  ukb.dat[, prots.exp] <- apply(ukb.dat[, prots.exp], 2, function(x){
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  })
  ukb.dat              <- as.data.table(ukb.dat)
  
  ## adjustment
  if(cohort == "Both"){
    adj <- paste(covs, collapse = " + ")
  }else{
    adj <- covs
  }
  
  ## intermediate output
  cat("run regression with ", nrow(ukb.dat), " samples in total\n")
  
  ## do in parallel
  registerDoMC(2)
  
  ## run testing
  res.surv <- mclapply(prots.exp, function(x){
    
    print(x)
    
    ## run model: implement cox-prop test and possible step function for time intervals
    ff   <- tryCatch(
      {
        coxph(as.formula(paste0("Surv(fol, bin_", phe,") ~ ", x, " + ", adj)), ukb.dat, ties = "breslow")
      }, error=function(e){
        return(NA)
      }, warning=function(w){
        return(NA)
      })
    
    ## proceed only, if model converged
    if(length(ff) > 1){
      
      ## cox-prop test
      ff.p <- cox.zph(ff)
      ## time-varying effect
      ff.s <- survSplit(as.formula(paste0("Surv(fol, bin_", phe,") ~ ", x, " + ", adj)), ukb.dat, cut=c(2,5,10), episode = "tgroup")
      ## fit the model (need to redefine names of variables)
      ff.s <-  tryCatch(
        {
          coxph(as.formula(paste0("Surv(tstart, fol,","bin_", phe,") ~ ", x, ":strata(tgroup) + ", adj)), ff.s, ties = "breslow")
        }, error=function(e){
          return(NA)
        })
      
      ## proceed only if model converged
      if(length(ff.s) > 1){
        ## compute summary for storage
        ff   <- summary(ff)
        ff.s <- summary(ff.s)
        ## depend on outcome
        tmp <- data.frame(id=x, phecode=phe, beta=ff$coefficients[1,1], se=ff$coefficients[1,3], pval=ff$coefficients[1,5], nevent=ff$nevent, nall=sum(!is.na(ukb.dat[, ..x])),
                          p.resid.prot=ff.p$table[x, 3], p.resid.overall=ff.p$table["GLOBAL", 3])
        ## add estimates for each strata
        jj  <- grep("strata", rownames(ff.s$coefficients), value=T)
        for(j in 1:length(jj)){
          ## estimate
          tmp[, paste0("beta.", j, ".ti")] <- ff.s$coefficients[jj[j],1]
          ## standard error
          tmp[, paste0("se.", j, ".ti")] <- ff.s$coefficients[jj[j],3]
          ## estimate
          tmp[, paste0("pval.", j, ".ti")] <- ff.s$coefficients[jj[j],5]
        }
        ## return
        return(tmp)
        
      }
    }
  }, mc.cores = 2) 
  ## combine into one
  res.surv <- rbindlist(res.surv, fill = T)
  
  ## store the results
  fwrite(res.surv, paste("output/res.surv.cis", phe, "txt", sep="."), sep="\t", row.names = F, na = NA)
  
}else{
  
  cat("not enough cases for", phe)
  
  ## create dummy output to track results later on
  res.surv <- data.table(id=NA, phecode=phe)
  
  ## store the results
  fwrite(res.surv, paste("output/res.surv.cis", phe, "txt", sep="."), sep="\t", row.names = F, na = NA)
}
