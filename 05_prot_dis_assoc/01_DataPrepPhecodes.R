###################################
####          Data Prep        ####
#### Kamil Demircan 08/11/2024 ####
###################################

rm(list=ls())
setwd("/sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction/input/")
options(stringsAsFactors = F)
# load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(rms)
require(readxl)
require(arrow)
require(tidyverse)
require(gprofiler2)
require(stringr)

####################################
####    for proteomic subset    ####
####################################
ukb.comb                 <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/51_protein_variance_updated/01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt")

## import
ukb.phe        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes_20250120/final/UKB.phecode.first.occurrence.long.format.20250122.txt.gz", sep="\t", header = T)
names(ukb.phe)[1] <- "f.eid"
names(ukb.phe)[3] <- "date"
names(ukb.phe)[4] <- "resource"

## subset to people with protein data
ukb.phe        <- ukb.phe[ f.eid %in% ukb.comb$f.eid ]
## convert to wide format
ukb.phe        <- dcast(ukb.phe, f.eid ~ phecode, value.var = c("date", "resource"))
## keep only what is of immediate interest
ukb.phe        <- ukb.phe[, c("f.eid", grep("date", names(ukb.phe), value=T)), with=F]
## import death dates for censoring
ukb.dod        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/death_certificates/death.txt")
## add date of death (keep only first instance)
ukb.phe        <- merge(ukb.phe, ukb.dod[ ins_index == 0 , c("eid", "date_of_death")], all.x=T, by.x="f.eid", by.y="eid")

## import label
lab.phe <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes_20250120/final/UKB.phecode.label.first.occurrence.20250122.txt", sep="\t", header=T)
## add identifier to match with data set
lab.phe[, id := paste0("date_", phecode)]
## subset to what is in the data
lab.phe        <- lab.phe[ id %in% names(ukb.phe) ]
# ## add sex-specific coding
# lab.tmp        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes/public/phecode_definitions1.2.csv")
# ## create variable to map to UKBB data set
# lab.tmp[, id := paste0("date_", phecode)]
# ## add sex-specific information to the data
# lab.phe        <- merge(lab.phe, lab.tmp[, c("id", "sex")], by="id")
# ## replace missing ones
lab.phe[, sex := ifelse(sex == "", "Both", sex)]

## import baseline date
tmp.dat        <- read_parquet("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                               col_select = c("f.eid", "f.53.0.0"))
names(tmp.dat) <- c("f.eid", "baseline_date")

## convert to data frame to ease coding
ukb.phe        <- as.data.frame(ukb.phe)
## add non-cases and baseline data
ukb.phe        <- merge(ukb.phe, tmp.dat[, c("f.eid", "baseline_date")], all.x=T)

## convert death date
ukb.phe$date_of_death <- as.IDate(ukb.phe$date_of_death, format = "%d/%m/%Y")
ukb.phe$baseline_date <- as.IDate(ukb.phe$baseline_date)
## add sex to count cases (does also bring back participants w/o EHR but proteomic data)
ukb.phe               <- merge(ukb.phe, ukb.comb[, c("f.eid", "sex")], all=T)

## create case and date variable: choose
for(j in lab.phe$id){
  ## define case status (binary)
  ukb.phe[, gsub("date", "bin", j)]  <- ifelse(is.na(ukb.phe[, j]), 0, 
                                               ifelse(ukb.phe[, j] > ukb.phe[, "baseline_date"], 1, NA))
  ## define new date
  ukb.phe[, gsub("date", "surv", j)] <- as.IDate(ifelse(is.na(ukb.phe[, j]) & is.na(ukb.phe[, "date_of_death"]), as.IDate("2020-08-31"), 
                                                        ifelse(is.na(ukb.phe[, j]) & !is.na(ukb.phe[, "date_of_death"]), ukb.phe[, "date_of_death"], ukb.phe[, j])))
}

## compute case numbers
lab.phe[, inc_cases := apply(lab.phe[, c("id", "sex")], 1, function(x){
  if(x[2] == "Both"){
    sum(ukb.phe[, gsub("date", "bin", x[1])] == 1, na.rm=T)
  }else{
    sum(ukb.phe[ukb.phe$sex == x[2], gsub("date", "bin", x[1])] == 1, na.rm=T)
  }
})]

## at least 200 incident cases
## remove congenital anomalies (1), injuries & poisonings (24), no category (10)
lab.phe <-
  lab.phe[ inc_cases >= 200 & !category %in% c("injuries & poisonings", "congenital anomalies") & category != "NULL"]

## keep only what is really needed
ukb.phe <- ukb.phe[, c("f.eid", "baseline_date", paste0("bin_", lab.phe$phecode), paste0("surv_", lab.phe$phecode))]
gc(reset=T)

lab.phe <- lab.phe %>% select(phecode,sex)
lab.phe$sex <- as.factor(lab.phe$sex)
fwrite(lab.phe, "input_phecodes.txt", sep="\t", col.names = F, row.names = F, na = NA, quote = F)

## save
fwrite(ukb.phe, "UKB.olink.phecodes.20250516.txt", sep="\t", row.names = F, na = NA)

####################################
####       for full cohort      ####
####################################
ukb.comb                 <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/51_protein_variance_updated/01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt")

## import
ukb.phe        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes_20250120/final/UKB.phecode.first.occurrence.long.format.20250122.txt.gz", sep="\t", header = T)
names(ukb.phe)[1] <- "f.eid"
names(ukb.phe)[3] <- "date"
names(ukb.phe)[4] <- "resource"

## to exclude
ids_include <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/sample_inclusion/qc_pass_EUR_panukbb_all_unrelated.ids")

## subset to people with protein data
ukb.phe        <- ukb.phe[ f.eid %in% ids_include$V1 ]
## convert to wide format
ukb.phe        <- dcast(ukb.phe, f.eid ~ phecode, value.var = c("date", "resource"))
## keep only what is of immediate interest
ukb.phe        <- ukb.phe[, c("f.eid", grep("date", names(ukb.phe), value=T)), with=F]
## import death dates for censoring
ukb.dod        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/death_certificates/death.txt")
## add date of death (keep only first instance)
ukb.phe        <- merge(ukb.phe, ukb.dod[ ins_index == 0 , c("eid", "date_of_death")], all.x=T, by.x="f.eid", by.y="eid")

## import label
lab.phe <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes_20250120/final/UKB.phecode.label.first.occurrence.20250122.txt", sep="\t", header=T)
## add identifier to match with data set
lab.phe[, id := paste0("date_", phecode)]
## subset to what is in the data
lab.phe        <- lab.phe[ id %in% names(ukb.phe) ]
# ## add sex-specific coding
# lab.tmp        <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes/public/phecode_definitions1.2.csv")
# ## create variable to map to UKBB data set
# lab.tmp[, id := paste0("date_", phecode)]
# ## add sex-specific information to the data
# lab.phe        <- merge(lab.phe, lab.tmp[, c("id", "sex")], by="id")
# ## replace missing ones
lab.phe[, sex := ifelse(sex == "", "Both", sex)]

## import baseline date
tmp.dat        <- read_parquet("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                               col_select = c("f.eid", "f.53.0.0"))
names(tmp.dat) <- c("f.eid", "baseline_date")

## convert to data frame to ease coding
ukb.phe        <- as.data.frame(ukb.phe)
## add non-cases and baseline data
ukb.phe        <- merge(ukb.phe, tmp.dat[, c("f.eid", "baseline_date")], all.x=T)

## convert death date
ukb.phe$date_of_death <- as.IDate(ukb.phe$date_of_death, format = "%d/%m/%Y")
ukb.phe$baseline_date <- as.IDate(ukb.phe$baseline_date)
## add sex to count cases (does also bring back participants w/o EHR but proteomic data)
ukb.phe               <- merge(ukb.phe, ukb.comb[, c("f.eid", "sex")], all=T)

## create case and date variable: choose
for(j in lab.phe$id){
  ## define case status (binary)
  ukb.phe[, gsub("date", "bin", j)]  <- ifelse(is.na(ukb.phe[, j]), 0, 
                                               ifelse(ukb.phe[, j] > ukb.phe[, "baseline_date"], 1, NA))
  ## define new date
  ukb.phe[, gsub("date", "surv", j)] <- as.IDate(ifelse(is.na(ukb.phe[, j]) & is.na(ukb.phe[, "date_of_death"]), as.IDate("2020-08-31"), 
                                                        ifelse(is.na(ukb.phe[, j]) & !is.na(ukb.phe[, "date_of_death"]), ukb.phe[, "date_of_death"], ukb.phe[, j])))
}

## compute case numbers
lab.phe[, inc_cases := apply(lab.phe[, c("id", "sex")], 1, function(x){
  if(x[2] == "Both"){
    sum(ukb.phe[, gsub("date", "bin", x[1])] == 1, na.rm=T)
  }else{
    sum(ukb.phe[ukb.phe$sex == x[2], gsub("date", "bin", x[1])] == 1, na.rm=T)
  }
})]

lab.phe <- fread("input_phecodes.txt")

## keep only what is really needed
ukb.phe <- ukb.phe[, c("f.eid", "baseline_date", paste0("bin_", lab.phe$V1), paste0("surv_", lab.phe$V1))]
gc(reset=T)

## save
fwrite(ukb.phe, "UKB.olink.phecodes.fullcohort.20250516.txt", sep="\t", row.names = F, na = NA)


