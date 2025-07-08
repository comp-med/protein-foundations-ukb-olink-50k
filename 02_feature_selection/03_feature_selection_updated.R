#!/usr/bin/env Rscript

## script to run feature selection for a given trait and a set of explanatory variables
## Maik Pietzner 06/06/2025
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
options(rgl.useNULL = TRUE)
# print(R.Version())

## correct directory
setwd("<path>/02_feature_selection/")

## packages needed
require(data.table)
require(doMC)
# require(caret)
require(glmnet)
require(stabs)
require(readxl)
require(susieR)
require(heplots)

## import Olink variable to be tested
var.olink <- args[1]
## t.b.i. be flexible for sub populations
pop       <- args[2]

cat("Performing variance decomposition for", var.olink, "\n")
cat("--------------------------------------------------\n")

#-----------------------------#
##-- import phenotype data --##
#-----------------------------#

## Olink data (include relevant confounder)
ukb.olink                     <- fread("../01_phenotype_prep/data/UKB.Olink.var.decomp.20250128.txt", select = c("eid", var.olink, "sample_age"))
## omit missing values
ukb.olink                     <- na.omit(ukb.olink)
## imputed data set of phenotypes
ukb.phe                       <- fread("../01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt")
## label for explanatory variables
lab.phe                       <- fread("../01_phenotype_prep/data/UKB.labels.variance.decomp.20280129.txt")
## reduce to entries in the annotation file
ukb.phe                       <- ukb.phe[, c("f.eid", lab.phe$short_name), with=F]

## convert to data frame to ease coding
ukb.phe                       <- as.data.frame(ukb.phe)
## choose most frequent factor as a reference for all other traits
for(j in c(lab.phe[ type %in% c("character", "factor")]$short_name, "centre")){
  ## get frequency
  jj           <- table(ukb.phe[, j])
  ## redefine
  ukb.phe[, j] <- factor(ukb.phe[, j], levels = names(jj[order(jj, decreasing = T)]))
}
## convert back
ukb.phe                       <- as.data.table(ukb.phe)

## create second set with extended labels
lab.ext                       <- lapply(1:nrow(lab.phe), function(x){
  ## check whether factor
  if(lab.phe$type[x] == "character" | lab.phe$short_name[x] == "centre" | lab.phe$type[x] == "factor"){
    ## get the relevant variable
    n <- lab.phe$short_name[x]
    ## get all levels
    l <- levels(as.factor(unlist(ukb.phe[, ..n])))
    print(l)
    ## replace odd characters
    l <- gsub("[[:punct:]]| ", "_", l)
    ## extend labels accordingly
    return(data.table(lab.phe[x,], short_name_new = paste0(n, l[-1]),
                      label_new = paste(lab.phe$label[x], l[-1]), reference = l[1]))
  }else{
    ## extend labels accordingly
    return(data.table(lab.phe[x,], short_name_new = gsub("[[:punct:]]| ", "_", lab.phe$short_name[x]), label_new = lab.phe$label[x], reference = NA))
  }
})
lab.ext                       <- do.call(rbind, lab.ext)

## recode factors to dummy variables
var.phe                       <- names(ukb.phe)[-1]
## convert to model matrix (recode factors)
tmp                           <- model.matrix(as.formula(paste(" ~", paste(var.phe, collapse = " + "))), ukb.phe)
## edit some names
colnames(tmp)                 <- gsub("[[:punct:]]| ", "_", colnames(tmp))
## create
ukb.phe                       <- data.table(ukb.phe[, "f.eid"], tmp[,-1])
## define variable list
var.phe                       <- names(ukb.phe)[-1]
## clear some space
rm(tmp); gc(reset=T)

## merge with protein data
ukb.phe                       <- merge(ukb.phe, ukb.olink, by.x = "f.eid", by.y = "eid")

## add ten random variables
for(j in 1:10){
  ukb.phe[, paste0("rand_", j)] <- rnorm(nrow(ukb.phe), 0, 1)
}

#-----------------------------#
##--  subset to population --##
#-----------------------------#

## currently only accounts for sex and ethnicity
if(pop != "All"){
  if(pop == "Female"){
    ukb.phe <- ukb.phe[ sexMale == 0]
  }else if(pop == "Male"){
    ukb.phe <- ukb.phe[ sexMale == 1]
  }else if(pop == "EUR"){
    ## invert all other ancestries  "AFR", "AMR", "CSA", "EAS", "MID"
    ukb.phe <- ukb.phe[ popAFR == 0 & popCSA == 0 & popMID == 0 & popEAS == 0 & popMID == 0]
  }else{
    ukb.phe <- ukb.phe[ eval(as.name(paste0("pop", pop))) == 1]
  }
}

#-----------------------------#
##--  import genotype data --##
#-----------------------------#

## get all genetic scores computed
pqtl <- dir("<path>/01_phenotype_prep/genetics/")
pqtl <- grep("genetic", pqtl, value = T)

## import if score exists
if(paste0(var.olink, ".genetic") %in% pqtl){
  
  ## import
  snps    <- fread(paste0("<path>/01_phenotype_prep/genetics/", var.olink, ".genetic"))
  ## drop NAs (due to X-chromsome having less observations)
  snps    <- na.omit(snps)
  
  ## add available scores to the data set
  ukb.phe <- merge(ukb.phe, snps[, c("IID", grep("score", names(snps), value=T)), with=F], by.x="f.eid", by.y="IID")
  
  ## add scores to variable list to be used in variable selection
  var.phe <- c(var.phe, grep("score", names(snps), value=T))
  
}

#----------------------------------------------#
##--           split before selection       --##
#----------------------------------------------#

## create 70/30 split of the data for feature selection and percentage estimation
set.seed(42)
id.feat <- sample(ukb.phe$f.eid, round(nrow(ukb.phe)*.7))
id.val  <- ukb.phe$f.eid[ !(ukb.phe$f.eid %in% id.feat)]

## write those to file to possibly use them again later
write.table(id.feat, paste("tmpdir/feature.set", pop, var.olink, "txt", sep="."), sep = "\t", row.names=F)

#----------------------------------------------#
##--           stability selection          --##
#----------------------------------------------#

## run in parallel
registerDoMC(10)

## run selection
set.seed(42)
res.stab <- stabsel(ukb.phe[ f.eid %in% id.feat, c(var.phe, paste0("rand_", 1:10), "sample_age"), with=F], 
                    unlist(ukb.phe[ f.eid %in% id.feat, var.olink, with=F]),
                    fitfun = glmnet.lasso, cutoff = 0.75,
                    PFER = 1, mc.cores=10, B=500)

## get the selected variables
var.sel  <- names(res.stab$selected)

## prune for possible redundant variables (sometimes picked)
dup.col   <- duplicated(as.list(ukb.phe[, ..var.sel, with=F]))
var.sel   <- var.sel[!(dup.col)]

## set to missing if random variables got selected
if(length(grep("rand_", var.sel)) > 0){
  var.sel <- c()
}

## selection matrix
write.table(res.stab$phat, paste("output_updated/variable.matrix", pop, var.olink, "txt", sep="."), sep = "\t", row.names=T)

#----------------------------------------------#
##--       compute explained variance       --##
#----------------------------------------------#

cat("Found", length(var.sel), "informative features \n")
cat("--------------------------------------------------\n")


## --> map original variables <-- ##

## lasso
var.sel <- unique(sapply(var.sel, function(x){
  if(x %in% lab.ext$short_name_new){
    return(lab.ext$short_name[ which(lab.ext$short_name_new == x)])
  }else{
    return(x)
  }
})) 
print(var.sel)

if(length(var.sel) == 0){
  var.sel <- c()
}

cat("Perform variance decomposition \n")
cat("--------------------------------------------------\n")

## --> import original data set (due to splitting of categorical variables ) <-- ##

## imputed data set of phenotypes
ukb.var                       <- fread("../01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt", select=c("f.eid", unique(c(var.sel))))
## create new combined data set
ukb.var                       <- merge(ukb.var, ukb.olink, by.x = "f.eid", by.y = "eid")
## convert to data table
ukb.var                       <- as.data.table(ukb.var)

## add genetic scores if needed
# jj                            <- grep("\\.score", unique(c(var.sel, cred.sel)), value=T)
jj                            <- grep("\\.score", unique(c(var.sel)), value=T)

if(length(jj) > 0){
  ## add to the data
  ukb.var <- merge(ukb.var, ukb.phe[, c("f.eid", jj), with=F])
}

## get variables that possibly need recoding
kk                            <- lab.phe[ (type %in% c("character", "factor") | short_name == "centre") & short_name %in% var.sel]$short_name
## process only if any
if(length(kk) > 0){
  ## convert to data frame to ease coding
  ukb.var                       <- as.data.frame(ukb.var)
  ## choose most frequent factor as a reference for all other traits
  for(j in kk){
    ## get frequency
    jj           <- table(ukb.var[, j])
    ## redefine
    ukb.var[, j] <- factor(ukb.var[, j], levels = names(jj[order(jj, decreasing = T)]))
  }
  ## convert back
  ukb.var                       <- as.data.table(ukb.var)
}



## --> compute variance composition <-- ##

## create 200 bootstraps of the validation data
set.seed(42)
id.boot   <- lapply(1:200, function(x) sample(id.val, ceiling(length(id.val)*.9)))

## --> LASSO <-- ##

## run only if at least one variable has been selected
if(length(var.sel) > 0){
  ## run estimation across subsets
  var.lasso <- mclapply(id.boot, function(x){
    ## run
    tmp <- tryCatch({
      etasq(lm(paste( var.olink, "~", paste(var.sel, collapse = " + ")), ukb.var[ f.eid %in% x ]), partial = T)
    }, 
    error = function(e){
      return(data.frame(Residuals = NA))
    })
    return(as.data.table(t(tmp)))
  }, mc.cores = 10)
  ## combine 
  var.lasso <- rbindlist(var.lasso, fill = T)
  ## replace resdiuals
  var.lasso[, Residuals := 1- apply(var.lasso[, -ncol(var.lasso), with = F], 1, sum, na.rm=T)]
  ## store number of iterations, in case there are massive problems with the bootstrapping
  n.boot    <- nrow(var.lasso)
  ## create aggregated measures for each variable
  var.lasso <- apply(var.lasso, 2, function(x) c(quantile(x, probs=c(.025,.25,.5,.75,.975), na.rm = T), mean(x, na.rm=T), sd(x, na.rm=T)))
  ## add column indicating summary measure
  var.lasso  <- as.data.table(var.lasso)
  var.lasso[, summary.measure := c("p.025", "p.25", "p.50", "p.75", "p.975", "p.mean", "p.sd")]
  ## add variable
  var.lasso[, var := var.olink]
  var.lasso[, n := nrow(ukb.phe)]
  var.lasso[, n.boot := n.boot]
}else{
  ## store results
  var.lasso  <- data.frame(var=var.olink, n=nrow(ukb.phe), Residuals=1)
}

#----------------------------------------------#
##--          write results to file         --##
#----------------------------------------------#

## Explained variance
write.table(var.lasso, paste("output_updated/lasso.explained.variance", pop, var.olink, "txt", sep="."), sep = "\t", row.names=F)
