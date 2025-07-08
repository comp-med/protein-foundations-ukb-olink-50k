#!/usr/bin/env Rscript

## script to perform imputation for phenotypes in UK Biobank
## Maik Pietzner 28/01/2025
rm(list=ls())

## get the arguments from the command line
# args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

## correct directory
setwd("<path>")

## packages needed
require(data.table)
require(doMC)
require(miceRanger)

#----------------------------#
##-- import relevant data --##
#----------------------------#

## import UKB data
ukb.comb <- fread("data//UKB.variables.variance.decomp.for.imputation.20250128.txt")
## define variables relevant for imputation
lab.set  <- fread("data/UKB.labels.variance.decomp.20250128.txt")
## update missing percentage
lab.set[, miss.per := sapply(short_name, function(x) nrow(ukb.comb[is.na(eval(as.name(x)))]))/nrow(ukb.comb)*100]

## convert some date variables
ukb.comb <- as.data.table(ukb.comb)
## take only hour of the day
ukb.comb[, Time_of_blow_measurement := hour(Time_of_blow_measurement) ]
## convert to years from very first visit
ukb.comb[, Date_of_attending_assessment_centre := as.numeric(as.IDate(Date_of_attending_assessment_centre) - min(as.IDate(Date_of_attending_assessment_centre)))/365.25]
## time blood draw
ukb.comb[, Time_blood_draw := hour(Time_blood_draw) + minute(Time_blood_draw)/60]

## get the data type for each feature
lab.set[, type := sapply(short_name, function(x) class(unlist(ukb.comb[, eval(x), with=F])))]
## convert to data frame to ease coding
ukb.comb <- as.data.frame(ukb.comb)

## convert character variables to factors
for(j in lab.set[ type == "character" ]$short_name){
  ukb.comb[, j] <- as.factor(ukb.comb[, j])
}

## convert logical to factor
for(j in lab.set[ type == "logical" ]$short_name){
  ukb.comb[, j] <- as.factor(ukb.comb[, j])
}

## convert integrer to numeric
for(j in lab.set[ type == "integer" ]$short_name){
  ukb.comb[, j] <- as.numeric(ukb.comb[, j])
}

## recode some others as well
ukb.comb$sex     <- as.factor(ukb.comb$sex)
ukb.comb$smoking <- as.factor(ukb.comb$smoking)
ukb.comb$alcohol <- as.factor(ukb.comb$alcohol)

cat("Need to impute", nrow(lab.set[ miss.per > 0]), "variables\n")

#----------------------------#
##--  perform imputation  --##
#----------------------------#

## create five new data sets
ukb.imp <- miceRanger(ukb.comb, m = 5, returnModels = F, verbose = T)
gc(reset=T)

## get the imputed data sets
ukb.imp <- completeData(ukb.imp)
gc(reset=T)

## write to file
for(j in 1:5){
  fwrite(ukb.imp[[j]], paste0("data/UKB.variables.variance.decomp.imputed.dataset.", j,".20280128.txt"), sep="\t", row.names = F, na = NA)
}
