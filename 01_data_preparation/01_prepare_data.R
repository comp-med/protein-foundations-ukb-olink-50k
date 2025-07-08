######################################################
#### Prepare comprehensive phenotypic data in UKB ####
#### Maik Pietzner                     28/01/2025 ####
######################################################

rm(list=ls())
setwd("<path>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(arrow)
require(readxl)
require(Hmisc)
require(tidyverse)
require(doMC)
require(igraph)

##################################
#### import basic information ####
##################################

## An example list of columns:
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.20116.0.0", "f.1558.0.0", paste0("f.22009.0.", 1:10), paste0("f.20003.0.", 0:47))
## names to assign
cl.names        <- c("f.eid", "age", "sex", "baseline_date", "centre", "bmi", "smoking", "alcohol", paste0("pc", 1:10), paste0("self_med", 0:47))

## import data from the main release
ukb.dat         <- read_parquet("<path>", col_select = cl.select)
## change names
names(ukb.dat)  <- cl.names

## define factor label for smoking
ukb.dat$smoking <- as.character(ukb.dat$smoking)
ukb.dat$smoking <- ifelse(ukb.dat$smoking == "Prefer not to answer", "Previous", ukb.dat$smoking)
ukb.dat$smoking <- factor(ukb.dat$smoking, levels = c("Never", "Previous", "Current"))

## define factor level for alcohol
ukb.dat$alcohol <- as.character(ukb.dat$alcohol)
ukb.dat$alcohol <- ifelse(ukb.dat$alcohol == "Prefer not to answer", "Three or four times a week", ukb.dat$alcohol)
ukb.dat$alcohol <- factor(ukb.dat$alcohol, levels = c("Never", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week",
                                                      "Daily or almost daily"))
## convert to data table
ukb.dat         <- as.data.table(ukb.dat)

##################################
####  self-report medication  ####
##################################

## import medication mapping
med.map <- data.frame(read_excel("<path>/UK_medication_to_ATC_code_mapping.xlsx", skip=1))
## kepp only what is needed
med.map <- med.map[, c("Category", "Coding.a", "Medication.ATC.code", "Drug.name")]
## expand for ATC codes
med.map <- lapply(1:nrow(med.map), function(x){
  
  ## get all ATC codes separately
  atc <- strsplit(med.map$Medication.ATC.code[x], " \\|")[[1]]
  ## report back
  return(data.frame(category=med.map$Category[x], ukbb.coding=med.map$Coding.a[x], atc_code=atc, drug.name=med.map$Drug.name[x]))
  
})
## combine into one dataframe
med.map <- do.call(rbind, med.map)

## combine with UKB coding to identify whether systemic or local
med.tmp <- read.delim("<path>/coding4.tsv")
## add to medication map
med.map <- merge(med.map, med.tmp, by.x="ukbb.coding", by.y="coding")

#-------------------------#
##--  create data set  --##
#-------------------------#

## generate long format from ukb medical records
tmp      <- melt(ukb.dat[, c("f.eid", paste0("self_med", 0:47))], id.vars = "f.eid")
tmp      <- tmp[ !is.na(value)]
## add ATC codes from mapping table
tmp      <- merge(tmp, unique(med.map[, c("ukbb.coding", "atc_code")]), by.x="value", by.y="ukbb.coding", allow.cartesian = T)
## reduce to one ATC code each
tmp      <- tmp[ order(f.eid, atc_code, value)]
## create indicator
tmp[, ind := 1:.N, by=c("f.eid", "atc_code")]
tmp      <- tmp[ ind == 1]
## convert to wide format
tmp      <- dcast(tmp, f.eid ~ atc_code)
## n = 361,408 participants with at least one mapping ATC code
## create data set
ukb.self <- merge(ukb.dat, tmp, all.x=T, by="f.eid")

## subset medication ma to those
med.map  <- as.data.table(subset(med.map, atc_code %in% names(ukb.self)))
## get all unique ATC codes
atc.self <- data.table(atc_code=unique(med.map$atc_code))
## n = 1014 medications

## recode columns
ukb.self                      <- as.data.frame(ukb.self)
ukb.self[, atc.self$atc_code] <- apply(ukb.self[, atc.self$atc_code], 2, function(x){
  x <- ifelse(!is.na(x), 1, 0)
  return(x)
})
## convert back to data table
ukb.self                      <- as.data.table(ukb.self)
## clean some space
gc(reset=T)

#-------------------------#
##--  count reporting  --##
#-------------------------#

## count how many people report intake of each ATC-code
atc.self[, count.self := sapply(atc.self$atc_code, function(x) nrow(ukb.self[ eval(as.name(x)) == 1]))]

## add column whether prescribed in both sexes
atc.self[, n.men := sapply(atc.self$atc_code, function(x) nrow(ukb.self[ eval(as.name(x)) == 1 & sex == "Male"]))]
atc.self[, n.women := sapply(atc.self$atc_code, function(x) nrow(ukb.self[ eval(as.name(x)) == 1 & sex == "Female"]))]
atc.self[, sex := ifelse(n.men != 0 & n.women != 0, "both", ifelse(n.men != 0 & n.women == 0, "men", "women"))]

##################################
####   augmented phenotypes   ####
##################################

## import main data dictonary
lab.main <- fread("<path>")

#------------------------#
##-- technical variab --##
#------------------------#

# time of blood draw, fasting time, number of draws (field IDs 68, 74, and 3166),
# month of visit (field ID 53), OLink plate ID, OLink well ID; Olink PCs?

## import variables to be selected
tech.var        <- fread("variables_techincal_20230414.txt") 
## create short names
tech.var[, column_name := paste0("f.", id, ".0.0")]
tech.var[, short_name := gsub(" |\\(|\\)|,|-", "_", label)]
## indicate whether already in the data
tech.var[, released := column_name %in% lab.main$id.ukbb]

## import data from the main release
ukb.tech        <- read_parquet("<path>", col_select = c("f.eid", tech.var$column_name))
## change names
names(ukb.tech) <- c("f.eid", tech.var$short_name)
ukb.tech        <- as.data.table(ukb.tech)

## create month as another variable (allow to behave non-linearly)
ukb.tech[, month := month(Date_of_attending_assessment_centre)]
# ## create polynomials
# ukb.tech[, month_2 := month^2]
# ukb.tech[, month_3 := month^3]

## add to the technical variable list
tech.var        <- rbind(tech.var, data.table(id = NA, label = c("Month"), category = "Technical", column_name = NA, 
                                              short_name = c("month"), released = T))
## sample age has already been generated in the Olink data set!

#------------------------#
##-- body composition --##
#------------------------#

## import variables to be selected
body.var          <- fread("variables_body_composition_20230414.txt")
## create short names
body.var[, column_name := paste0("f.", id, ".0.0")]
body.var[, short_name := gsub(" |\\(|\\)|,|-", "_", label)]

## import data from the main release
ukb.body          <- read_parquet("<path>", col_select = c("f.eid", body.var$column_name, "f.31.0.0", "f.21022.0.0"))
## change names
names(ukb.body)   <- c("f.eid", body.var$short_name, "sex", "age")

## create waist-to-hip ratio
ukb.body          <- as.data.table(ukb.body)
ukb.body[, whr := Waist_circumference/Hip_circumference]
## compute values using UKB data
ukb.body[, hip := Hip_circumference]
ukb.body[, waist := Waist_circumference]

## add BMI
ukb.body[, bmi := Weight/(Standing_height/100)^2]

## --> create estimated effects using approximations <-- ##
## (https://www.medrxiv.org/content/10.1101/2020.12.16.20248330v1)

## create intermediate values
ukb.body[, Ht2ImLegs := ((Standing_height/100)^2)/(Impedance_of_leg__left_ + Impedance_of_leg__right_)]
ukb.body[, Ht2ImArms := ((Standing_height/100)^2)/(Impedance_of_arm__left_ + Impedance_of_arm__right_)]

## import model coefficient (sex-specific, ST3 of the paper)
dexa.coefficients <- fread("body_composition_estimation_Powell2020.txt", na.strings = "")
## convert some to be numeric
dexa.coefficients[, Ht2ImArms.m := as.numeric(gsub(",", "", Ht2ImArms.m))]
dexa.coefficients[, Constant.m := as.numeric(gsub(",", "", Constant.m))]
## create shorter variable
dexa.coefficients[, short_name := gsub(" |\\(|\\)|,|-", "_", variable)]
## replace missing values with zeros to ease computation
dexa.coefficients[ is.na(dexa.coefficients)] <- 0

## men
men          <- ukb.body[ sex == "Male"] 
men[, Constant.m := 1]
men[, height := Standing_height]
men          <- cbind(men[, "f.eid"], as.matrix(men[, c("age", "height", "Weight", "waist", "hip", "Impedance_of_whole_body", "Ht2ImLegs", "Ht2ImArms", "Constant.m"), with=F]) %*% 
                        t(as.matrix(dexa.coefficients[, c("age.m", "height.m", "weight.m", "waist.m", "hip.m", "impedance.m", "Ht2ImLegs.m", "Ht2ImArms.m", "Constant.m")])))
## convert
men          <- as.data.table(men)
names(men)   <- c("f.eid", dexa.coefficients$short_name)

## women
women        <- ukb.body[ sex == "Female"] 
women[, Constant.w := 1]
women[, height := Standing_height]
women        <- cbind(women[, "f.eid"], as.matrix(women[, c("age", "height", "Weight", "waist", "hip", "Impedance_of_whole_body", "Ht2ImLegs", "Ht2ImArms", "Constant.w"), with=F]) %*% 
                        t(as.matrix(dexa.coefficients[, c("age.w", "height.w", "weight.w", "waist.w", "hip.w", "impedance.w", "Ht2ImLegs.w", "Ht2ImArms.w", "Constant.w")])))
## convert
women        <- as.data.table(women)
names(women) <- c("f.eid", dexa.coefficients$short_name)

## combine into one (N.B.: negative values are not set to missing)
ukb.body     <- merge(ukb.body[, c("f.eid", "Standing_height", "whr", "hip", "waist", "Weight", "bmi")], rbind(men, women), by="f.eid")
## add to the label file
body.var     <- data.table(id=NA, label=c("Height", "Waist-to-Hip ratio", "Hip circumference", "Waist circumference", "Body mass index", "Weight", dexa.coefficients$variable), 
                           category="Body composition", 
                           column_name=NA, 
                           short_name=c("Standing_height", "whr", "hip", "waist", "bmi", "Weight", dexa.coefficients$short_name))

## calculate missing values
body.var[, miss.per := sapply(short_name, function(x) nrow(ukb.body[is.na(eval(as.name(x)))]))/nrow(ukb.body)*100]

## do some clearning
rm(men); rm(women); gc(reset = T)

#------------------------#
##--   socioeconomic  --##
#------------------------#

## import variables to be selected
soec.var        <- fread("variables_socioeconomic_20230414.txt")
## create short names
soec.var[, column_name := paste0("f.", id, ".0.0")]
soec.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]
## indicate whether already in the data
soec.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
soec.var        <- soec.var[ released == T]

## import data from the main release
ukb.soec        <- read_parquet("<path>", col_select = c("f.eid", soec.var$column_name))
## change names
names(ukb.soec) <- c("f.eid", soec.var$short_name)

## convert to useful data sets
table(ukb.soec$Leisure_social_activities, useNA = "always")
table(ukb.soec$Age_completed_full_time_education, useNA = "always") ## too many NAs
table(ukb.soec$Qualifications, useNA = "always") 
table(ukb.soec$Average_total_household_income, useNA = "always")
table(ukb.soec$Duration_of_moderate_activity, useNA = "always") ## poorly measures with many negative values, and almost 100k missing

#------------------------#
##--    biomarkers    --##
#------------------------#

## import QC'ed biomarker levels
ukb.bio         <- read_parquet("<path>")
## import label
bio.var         <- fread("<path>")
## create short name for cleaned variables
bio.var[, short_name := ifelse(outlier_transformation == "log10", paste0("log_", short_name), short_name)]
bio.var[, short_name := paste0(short_name, "_cleaned")]
## convert to the same format as other variable data sets
bio.var         <- bio.var[, c("UKBB_ID", "id", "short_name", "Biomarker")]
names(bio.var)  <- c("id", "column_name", "short_name", "label")
## add some other variables
bio.var[, released := TRUE]
bio.var[, category := "Biomarker"]

## subset the data
ukb.bio         <- ukb.bio[, c("f.eid", bio.var$short_name)]
## some cleaning
gc(reset = T)

#------------------------#
##--     diseases     --##
#------------------------#

## import phecode assignment
ukb.phe            <- fread("<path>", sep="\t", header = T)
## convert to wide format
ukb.phe            <- dcast(ukb.phe, eid ~ phecode, value.var = "diag.date")
## edit names
names(ukb.phe)[-1] <- paste0("date_", names(ukb.phe)[-1]) 

## import label
lab.phe            <- fread("<path>", sep="\t", header=T)
## add identifier to match with data set
lab.phe[, id := paste0("date_", phecode)]
## replace missing ones
lab.phe[, sex := ifelse(sex == "", "Both", sex)]

## convert to data frame to ease coding
ukb.phe            <- as.data.frame(ukb.phe)
## add non-cases and baseline data
ukb.phe            <- merge(ukb.phe, ukb.dat[, c("f.eid", "baseline_date")], all=T, by.x = "eid", by.y = "f.eid")
## convert date
ukb.phe$baseline_date <- as.IDate(ukb.phe$baseline_date)

## create case and date variable: choose
for(j in lab.phe$id){
  ## define case status (binary)
  ukb.phe[, gsub("date", "bin", j)]  <- ifelse(is.na(ukb.phe[, j]) | ukb.phe[, j] > ukb.phe[, "baseline_date"], 0, 1)
}
## add sex to count cases
ukb.phe            <- merge(ukb.dat[, c("f.eid", "sex")], ukb.phe, by.x = "f.eid", by.y = "eid")

## convert to data frame to ease coding
ukb.phe            <- as.data.frame(ukb.phe)
## compute case numbers
lab.phe[, cases := apply(lab.phe[, .(id, sex)], 1, function(x){
  if(x[2] == "Both"){
    print(gsub("date", "bin", x[1]))
    sum(ukb.phe[, gsub("date", "bin", x[1])] == 1, na.rm=T)
  }else{
    sum(ukb.phe[ukb.phe$sex == x[2], gsub("date", "bin", x[1])] == 1, na.rm=T)
  }
})]

## at least 100 cases overall
lab.phe <- lab.phe[ cases >= 100]
## n = 1198

## keep only what is really needed
ukb.phe <- ukb.phe[, c("f.eid", "baseline_date", paste0("bin_", lab.phe$phecode))]
gc(reset=T)

#------------------------#
##--       diet       --##
#------------------------#

## import variables to be selected
diet.var        <- fread("variables_diet_20230428.txt")
## create short names
diet.var[, column_name := paste0("f.", id, ".0.0")]
diet.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]
## indicate whether already in the data
diet.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
diet.var        <- diet.var[ released == T]

## import data from the main release
ukb.diet        <- read_parquet("<path>", col_select = c("f.eid", diet.var$column_name))
## change names
names(ukb.diet) <- c("f.eid", diet.var$short_name)

## overview
summary(ukb.diet)

#------------------------#
##-- blood cell count --##
#------------------------#

## import variables to be selected
blood.var        <- fread("variables_blood_cell_20230428.txt")
## create short names
blood.var[, column_name := paste0("f.", id, ".0.0")]
blood.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]
## indicate whether already in the data
blood.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
blood.var        <- blood.var[ released == T]

## import data from the main release
ukb.blood        <- read_parquet("<path>", col_select = c("f.eid", blood.var$column_name))
## change names
names(ukb.blood) <- c("f.eid", blood.var$short_name)

## --> do some data cleaning (drop extreme values) <-- ##

## draw simple histograms
pdf("../graphics/Distribution.UKBB.blood.cell.counts.pdf", width = 6.3, height = 6.3)
par(mar=c(1.5,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), mfrow=c(4,4), lwd=.5)
## plot histograms
for(j in blood.var$short_name){
  hist(ukb.blood[ , j], breaks=100, lwd=.2, xlab=j, main = "", col=colorspace::adjust_transparency("#00A4CC", .4))
}
dev.off()

## apply outlier pruning
ukb.blood[, blood.var$short_name] <- apply(ukb.blood[, blood.var$short_name], 2, function(x){
  ## create values of interest
  meds <- median(x, na.rm=T)
  mads <- mad(x, na.rm=T)
  ## test for above or below
  x    <- ifelse(x < (meds - 5*mads ) | x > (meds + 5*mads), NA, x)
  return(x)
})

## draw simple histograms again
pdf("../graphics/Distribution.UKBB.blood.cell.counts.outlier.omitted.pdf", width = 6.3, height = 6.3)
par(mar=c(1.5,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), mfrow=c(4,4), lwd=.5)
## plot histograms
for(j in blood.var$short_name){
  hist(ukb.blood[ , j], breaks=100, lwd=.2, xlab=j, main = "", col=colorspace::adjust_transparency("#00A4CC", .4))
}
dev.off()

## omit two variables from the analysis
blood.var <- blood.var[ !(short_name %in% c("Nucleated_red_blood_cell_count", "Nucleated_red_blood_cell_percentage"))]

#------------------------#
##--     pulmonary    --##
#------------------------#

## import variables to be selected
pulm.var        <- fread("variables_pulmonary_20230503.txt")
## create short names
pulm.var[, column_name := paste0("f.", id, ".0.0")]
pulm.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]
## indicate whether already in the data
pulm.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
pulm.var        <- pulm.var[ released == T]

## import data from the main release
ukb.pulm        <- read_parquet("<path>", col_select = c("f.eid", pulm.var$column_name))
## change names
names(ukb.pulm) <- c("f.eid", pulm.var$short_name)
summary(ukb.pulm)

## some have many missing values

#------------------------#
##--  cardiovascular  --##
#------------------------#

## import variables to be selected
card.var        <- fread("variables_cardio_20230428.txt")
## create short names
card.var[, column_name := paste0("f.", id, ".0.0")]
card.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]
## indicate whether already in the data
card.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
card.var        <- card.var[ released == T]

## import data from the main release
ukb.card        <- read_parquet("<path>", col_select = c("f.eid", card.var$column_name))
## change names
names(ukb.card) <- c("f.eid", card.var$short_name)
summary(ukb.card)

#------------------------#
##--        bone      --##
#------------------------#

## import variables to be selected
bone.var        <- fread("variables_bone_20230503.txt")
## create short names
bone.var[, column_name := paste0("f.", id, ".0.0")]
bone.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]
## indicate whether already in the data
bone.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
bone.var        <- bone.var[ released == T]

## import data from the main release
ukb.bone        <- read_parquet("<path>", col_select = c("f.eid", bone.var$column_name))
## change names
names(ukb.bone) <- c("f.eid", bone.var$short_name)
summary(ukb.bone)

#------------------------#
##--      health      --##
#------------------------#

## import variables to be selected
health.var        <- fread("variables_health_20230428.txt")
## create short names
health.var[, column_name := paste0("f.", id, ".0.0")]
health.var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]
## indicate whether already in the data
health.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
health.var        <- health.var[ released == T]

## import data from the main release
ukb.health        <- read_parquet("<path>", col_select = c("f.eid", health.var$column_name))
## change names
names(ukb.health) <- c("f.eid", health.var$short_name)
summary(ukb.health)

#------------------------#
##--     pollution    --##
#------------------------#

## import variables to be selected
poll.var        <- fread("variables_pollution_20230428.txt")
## create short names
poll.var[, column_name := paste0("f.", id, ".0.0")]
poll.var[, short_name := gsub(" |\\(|\\)|,|-|\\/|;", "_", gsub(", automated reading", "", label))]
## indicate whether already in the data
poll.var[, released := column_name %in% lab.main$id.ukbb]
## drop those
poll.var        <- poll.var[ released == T]

## import data from the main release
ukb.poll        <- read_parquet("<path>", col_select = c("f.eid", poll.var$column_name))
## change names
names(ukb.poll) <- c("f.eid", poll.var$short_name)
summary(ukb.poll)
## subsequent years are highly correlated, do no further processing

##################################
####    combined data set     ####
##################################

#------------------------------#
##--     combine labels     --##
#------------------------------#

## create combined data set of labels
lab.set    <- rbind(bio.var, blood.var, body.var, bone.var, card.var, diet.var, health.var, poll.var, pulm.var, soec.var, tech.var, fill=T)

## given to at least 50 patients
tmp        <- atc.self[ count.self >= 50 ] ## n = 704
## add label column, add first entry for drug name
tmp[, label := sapply(atc_code, function(x) paste(unique(med.map$drug.name[ which(med.map$atc_code == x)]), collapse = "|"))]
tmp[, category := "Drugs"]
## keep to be aligned with lab values
tmp        <- tmp[, c("atc_code", "label", "category")]
## rename
names(tmp) <- c("short_name", "label", "category")
## add 
lab.set    <- rbind(lab.set, tmp, fill=T)

## add prevalent diseases
tmp        <- lab.phe[, c("id", "phenotype")]
tmp$cat    <- "Diseases"
names(tmp) <- c("short_name", "label", "category")
## add
lab.set    <- rbind(lab.set, tmp, fill=T)

## replace names
lab.set[, short_name := gsub("^date_", "bin_", short_name)]

#------------------------------#
##--    combine data set    --##
#------------------------------#

## combine the relevant phenotypic data
ukb.comb   <- list(ukb.bio, ukb.blood, ukb.body, ukb.bone, ukb.card, ukb.dat[, c("f.eid", "age", "sex", "centre", "smoking", "alcohol", paste0("pc", 1:10)), with=F],
                   ukb.diet, ukb.health, ukb.poll, ukb.pulm, ukb.soec, ukb.tech) %>% reduce(full_join, by="f.eid")
## add medications
ukb.comb   <- merge(ukb.comb, ukb.self[, c("f.eid", atc.self[ count.self >= 50 ]$atc_code ), with=F], by="f.eid", all.x=T)
## add diseases
ukb.comb   <- merge(ukb.comb, ukb.phe, by="f.eid", all.x=T)
ukb.comb   <- as.data.table(ukb.comb)
## delete what is no longer needed
rm(list=c("ukb.bio", "ukb.blood", "ukb.body", "ukb.bone", "ukb.card", "ukb.diet", "ukb.health", "ukb.poll", "ukb.pulm", "ukb.soec", "ukb.tech", "ukb.phe", "ukb.self"))
gc(reset = T)

## import ancestry information (loop over all)
ukb.anc    <- fread("<path>")
## add to the data (this will reduce the size quite a lot)
ukb.comb   <- merge(ukb.comb, ukb.anc[, c("EID.44448", "pop")], by.x="f.eid", by.y="EID.44448")

## add certain variables to lab.set
lab.set  <- rbind(data.table(id=NA, column_name=NA, short_name=c("age", "sex", "centre", "smoking", "alcohol", "pop"), 
                             label=c("Age", "Sex", "Centre", "Smoking", "Alcohol consumption", "Ancestry"), released=T,
                             category=c(rep("Basic demographics", 5), "Ancestry"), miss.per=NA), lab.set)

## compute missing values
lab.set[, miss.per := sapply(short_name, function(x) nrow(ukb.comb[is.na(eval(as.name(x)))]))/nrow(ukb.comb)*100]
## drop variables with >50% missing values
lab.set    <- lab.set[ miss.per <= 50 ]
## n = 2054

## compute missing values by person
eid.miss <- apply(ukb.comb[, lab.set$short_name, with=F], 1, function(x) sum(is.na(x))/nrow(lab.set))
eid.miss <- data.frame(eid=ukb.comb$f.eid, miss.per = eid.miss)
## drop participants with >50% missing values
eid.miss <- subset(eid.miss, miss.per <= .5)
## one participant with all missing values --> drop
ukb.comb <- ukb.comb[ f.eid %in% eid.miss$eid ]

## write to file for imputation
fwrite(ukb.comb[, c("f.eid", lab.set$short_name), with=F], "UKB.variables.variance.decomp.20250128.txt", sep = "\t", row.names = F, na = NA)

## add type of variable
lab.set[, type := sapply(short_name, function(x) class(unlist(ukb.comb[, eval(x), with=F])))]

## corresponding labels
fwrite(lab.set, "UKB.labels.variance.decomp.20250128.txt", sep = "\t", row.names = F, na = NA)

###################################
####     import Olink data     ####
###################################

## import protein data
ukb.prot      <- fread("<path>")
## import label
lab.prot      <- fread("<path>")

#------------------------#
##--    do some QC    --##
#------------------------#

## count missing values by participant (start here, since one entire batch of proteins is missing for a subset of individuals)
eid.miss.prot           <- apply(ukb.prot[, lab.prot$id, with=F], 1, function(x) sum(is.na(x))/nrow(lab.prot))
eid.miss.prot           <- data.frame(eid=ukb.prot$eid, miss.per = eid.miss.prot)
## drop participants with >50% missing values
ukb.prot                <- ukb.prot[ eid %in% subset(eid.miss.prot, miss.per <= .5)$eid]

## same for proteins
lab.prot[, miss.per := sapply(id, function(x) nrow(ukb.prot[is.na(eval(as.name(x)))]))/nrow(ukb.prot)*100]
## drop proteins with >20% missing values
lab.prot                <- lab.prot[ miss.per <= 20 ]
## n = 2919 proteins

## apply normalization
ukb.prot                <- as.data.frame(ukb.prot)
ukb.prot[, lab.prot$id] <- apply(ukb.prot[, lab.prot$id], 2, function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
})

#----------------------------#
##--        export        --##
#----------------------------#

## export data to run association testing
fwrite(ukb.prot, "UKB.Olink.var.decomp.20250128.txt", sep="\t", row.names=F, na = NA)
## export list of olink variables to be used as outcome
write.table(lab.prot$id, "olink.variable.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(lab.prot, "Olink.proteins.variance.decompostion.20250128.txt", sep="\t", row.names=F)

#----------------------------#
##-- phenotypic data imp  --##
#----------------------------#

## subset to people with protein data
ukb.comb                <- ukb.comb[ f.eid %in% ukb.prot$eid ]
rm(ukb.dat); gc(reset=T)

## store to perform imputation of missing values
fwrite(ukb.comb, "UKB.variables.variance.decomp.for.imputation.20250128.txt", sep="\t", row.names=F, na = NA)

##################################
####    import imputed data   ####
##################################

## import
ukb.phe  <- fread("../../01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt")

#--------------------------------------#
##-- among purely numeric variables --##
#--------------------------------------#

## compute simple correlation
cont.cor                                   <- lab.set[ type %in% c("numeric", "integer")]$short_name
## exclude test centre from the list
cont.cor                                   <- cont.cor[-which(cont.cor == "centre")]
## compute simple correlation
cont.cor                                   <- Rfast::cora(ukb.phe[, cont.cor, with=F])
## drop redundancies
cont.cor[ lower.tri(cont.cor, diag = T)]   <- NA
## convert into table
cont.cor                                   <- reshape2::melt(cont.cor, na.rm=T)
cont.cor                                   <- as.data.table(cont.cor)
## convert to character
cont.cor[, Var1 := as.character(Var1)]
cont.cor[, Var2 := as.character(Var2)]

#--------------------------------------#
##--     among factor variables     --##
#--------------------------------------#

## use advice from: https://rpubs.com/hoanganhngo610/558925 (contigency coefficient)
# DescTools::ContCoef(ukb.phe$Spread_type, ukb.phe$centre)
require(DescTools)

## compute at scale
fact.cor                                   <- c(lab.set[ !(type %in% c("numeric", "integer"))]$short_name, "centre")
## get all possible combinations (scale to 0 to 1)
fact.cor                                   <- PairApply(ukb.phe[, fact.cor, with=F], ContCoef, correct = TRUE)
## drop redundancies
fact.cor[ lower.tri(fact.cor, diag = T)]   <- NA
## convert into table
fact.cor                                   <- reshape2::melt(fact.cor, na.rm=T)
fact.cor                                   <- as.data.table(fact.cor)
## convert to character
fact.cor[, Var1 := as.character(Var1)]
fact.cor[, Var2 := as.character(Var2)]

## don't seem to be much issues; only drop 'Lamb_mutton_intake'

#--------------------------------------#
##--        mixed data types        --##
#--------------------------------------#

## compute explained variance by pair (to allow categorical variables to contribute as well);
mix.cor                                    <- expand.grid(Var1 = unique(c(fact.cor$Var1, fact.cor$Var2)),
                                                          Var2 = unique(c(cont.cor$Var1, cont.cor$Var2)), 
                                                          stringsAsFactors = F)
## loop through
registerDoMC(12)
mix.cor                                    <- mclapply(1:nrow(mix.cor), function(x){
  
  ## get pair to best tested
  t1 <- mix.cor$Var1[x]
  t2 <- mix.cor$Var2[x]
  
  ## report
  cat("testing", t1, "and", t2, "\n")
  cat("\n-----------------------------\n")
  
  ## run linear regression model 
  ff <- summary(lm(paste(t2, t1, sep="~"), ukb.phe))
  
  ## return results
  return(data.table(Var1 = t1, Var2 = t2, value=sqrt(ff$adj.r.squared)))
  
  
}, mc.cores = 12)
## combine
mix.cor                                    <- rbindlist(mix.cor)
## only 'sex' and testosterone are possibly 'problematic' 

#--------------------------------------#
##--      prune cont. variables     --##
#--------------------------------------#

## map to network (chosen, as it separates lean mass traits from other anthro measures that would be otherwise lost)
cont.net <- graph_from_data_frame(cont.cor[ value >= .85])
## get all separate components
cont.net <- components(cont.net)$membership
## convert to data frame
cont.net <- data.table(short_name=names(cont.net), id.group=cont.net)
## add missingness
cont.net <- merge(cont.net, lab.set)
## order accordingly
cont.net <- cont.net[ order(id.group, miss.per)]
## add indicator
cont.net[, ind := 1:.N, by="id.group"]

## do some manual tweaking to keep 'biological interpretability'
cont.net$ind[ cont.net$short_name == "ldl_chol_cleaned"]      <- 1
cont.net$ind[ cont.net$short_name == "chol_cleaned"]          <- 2
cont.net$ind[ cont.net$short_name == "bmi"]                   <- 1
cont.net$ind[ cont.net$short_name == "waist"]                 <- 2
cont.net$ind[ cont.net$short_name == "Total_Lean_Mass__g_"]   <- 1
cont.net$ind[ cont.net$short_name == "Android_Lean_Mass__g_"] <- 2

#--------------------------------------#
##--   final feature set for pred.  --##
#--------------------------------------#

## generate new feature set
lab.final <- lab.set[ !(short_name %in% c(cont.net[ ind != 1]$short_name, "Lamb_mutton_intake"))]
## n = 1,876 feature

## corresponding labels
fwrite(lab.final, "UKB.labels.variance.decomp.20280129.txt", sep = "\t", row.names = F, na = NA)

#--------------------------------------#
##--  prepare extended feature set  --##
#--------------------------------------#

## convert to data frame to ease coding
ukb.phe                       <- as.data.frame(ukb.phe)
## choose most frequent factor as a reference for all other traits
for(j in c(lab.set[ type %in% c("character", "factor")]$short_name, "centre")){
  ## get frequency
  jj           <- table(ukb.phe[, j])
  ## redefine
  ukb.phe[, j] <- factor(ukb.phe[, j], levels = names(jj[order(jj, decreasing = T)]))
}
## convert back
ukb.phe                       <- as.data.table(ukb.phe)

## create second set with extended labels
lab.ext                       <- lapply(1:nrow(lab.set), function(x){
  ## check whether factor
  if(lab.set$type[x] == "character" | lab.set$short_name[x] == "centre" | lab.set$type[x] == "factor"){
    ## get the relevant variable
    n <- lab.set$short_name[x]
    ## get all levels
    l <- levels(as.factor(unlist(ukb.phe[, ..n])))
    print(l)
    ## replace odd characters
    l <- gsub("[[:punct:]]| ", "_", l)
    ## extend labels accordingly
    return(data.table(lab.set[x,], short_name_new = paste0(n, l[-1]),
                      label_new = paste(lab.set$label[x], l[-1]), reference = l[1]))
  }else{
    ## extend labels accordingly
    return(data.table(lab.set[x,], short_name_new = gsub("[[:punct:]]| ", "_", lab.set$short_name[x]), label_new = lab.set$label[x], reference = NA))
  }
})
lab.ext                       <- do.call(rbind, lab.ext)

## export
fwrite(lab.ext[ short_name %in% lab.final$short_name], "UKB.labels.variance.decomp.factor.extended.20250129.txt", sep = "\t", row.names = F, na = NA)

##################################
####     ST1 - demographics   ####
##################################

## import function to do so
source("../scripts/create_tab_1.R")

#------------------------------#
##-- compute summary by sex --##
#------------------------------#

## all
st1 <- lapply(c("All", "Female", "Male", "EUR", "AFR", "CSA"), function(x){
  ## create table in relevant subgroup
  if(x == "All"){
    tmp <- create.table1(ukb.phe, lab.final$short_name)
  }else if(x %in% c("Female", "Male")){
    tmp <- create.table1(ukb.phe[ sex == x], lab.final$short_name)
  }else{
    tmp <- create.table1(ukb.phe[ pop == x], lab.final$short_name)
  }
  ## return results
  return(data.table(strata=x, tmp))
})
## combine and convert
st1 <- rbindlist(st1)
## convert 
st1 <- dcast(st1, var ~ strata, value.var = c("N", "statistic"))
## add label
lab.final[, srt := 1:nrow(lab.final)]
st1 <- merge(lab.final, st1, by.x = "short_name", by.y = "var")

## write to file
write.table(st1, "Supplementary.Table.1.UKB.summary.20250129.txt", sep="\t", row.names = F)

##################################
####  prepare genetic scores  ####
##################################

#--------------------------#
##--    import stats    --##
#--------------------------#

## import Supplemental Table with all cis/trans pQTLs (from: https://www.nature.com/articles/s41586-023-06592-6)
pqtl                     <- as.data.table(readxl::read_excel("<path>/41586_2023_6592_MOESM3_ESM.xlsx", 11, skip = 4))
## subset to protein of interest
pqtl[, protein_id := sapply(`UKBPPP ProteinID`, function(x) tolower(strsplit(x, ":")[[1]][1]))]

## compute MarkerName to ease merging
pqtl[, MarkerName := sapply(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, function(x){
  ## split
  x <- strsplit(x, ":")[[1]]
  ## return MarkerName
  return(paste0("chr", as.numeric(ifelse(x[1] == "X", 23, x[1])), ":", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_")))
})]
## extract effect allele
pqtl[, effect_allele := sapply(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, function(x) strsplit(x, ":")[[1]][4])]

## import UKB file to map possible missing rsIDs
ukb.snps                 <- fread("<path>")
## subset to IDs needed
ukb.snps                 <- ukb.snps[ MarkerName %in% pqtl$MarkerName]
gc(reset=T)
## replace rsIDs accordingly
pqtl                     <- merge(pqtl, unique(ukb.snps[, .(MarkerName, rsid)]), by="MarkerName", all.x=T)
## replace rsIDs if needed
pqtl[, rsID := ifelse(rsID == "-", rsid, rsID)]
## drop SNPs still not found
pqtl                     <- pqtl[ rsID != "-"]
## subset to proteins included
pqtl                     <- pqtl[ protein_id %in% lab.prot$id]

#--------------------------#
##--     obtain SNPs    --##
#--------------------------#

## write to file by Chromosome: chunk in batches!
for(j in 1:23){
  print(j)
  ## get the list of snps
  snps <- unique(pqtl[ CHROM == j]$rsID)
  ## chunk
  snps <- split(snps, ceiling(seq_along(snps)/100))
  print(length(snps))
  ## write to file
  for(k in 1:length(snps)){
    write.table(snps[[k]], paste("snps", ifelse(j == 23, "X", j), k,"txt", sep = "."), col.names = F, row.names = F, quote = F)
  }
}

## obtain input list
tmp <- dir()
tmp <- grep("snps", tmp, value = T)
tmp <- data.table(file=tmp)
tmp[, chr := sapply(file, function(x) strsplit(x, "\\.")[[1]][2])]
tmp[, chunk := sapply(file, function(x) strsplit(x, "\\.")[[1]][3])]
## reorder by chunk
tmp <- tmp[order(chunk)]
## write to file
write.table(tmp, "files.obtain.snps.txt", sep = "\t", row.names = F, col.names = F, quote = F)

## --> gather SNP info <-- ##

## import snp info
snp.info <- grep("info", dir("../genetics/"), value = T)
snp.info <- lapply(snp.info, function(x){
  ## import
  tmp <- fread(paste0("../genetics/", x))
  ## return along wiht identifier
  tmp[, block := x]
  return(tmp)
})
## combine
snp.info      <- rbindlist(snp.info)
## create MarkerName 
snp.info[, MarkerName := paste0("chr", ifelse(chromosome == "X", 23, chromosome), ":", position, "_", pmin(alleleA, alleleB), "_", pmax(alleleA, alleleB))]
## preserve order to be able to map SNPs to position in dosage files
snp.info[, srt := 1:.N, by = "block"]

## check whether SNPs have not been found
pqtl[ !(rsID %in% snp.info$rsid)]
## all found!!

## --> compute scores <-- ##

### what is missing
jj            <- grep("genetic", dir("../genetics/"), value = T) 
jj            <- gsub("\\.genetic", "", jj)
jj            <- unique(pqtl[!(protein_id %in% jj)]$protein_id)

## loop through and write to file
prot.genetics <- lapply(jj, function(x){
  
  print(x)
  
  ## get relevant SNPs
  tmp         <- pqtl[ protein_id == x]
  ## merge with tmp data (allele B is effect allele in the requested data)
  tmp         <- merge(tmp, snp.info, by="MarkerName")
  ## recode effect according to the data set
  tmp[, weight_pgs := ifelse(alleleB == effect_allele, BETA, -BETA)]
  
  ## take care of multiallelic variants
  ii          <- table(tmp$MarkerName)
  tmp         <- tmp[ MarkerName %in% names(ii[ ii == 1])]
  
  if(nrow(tmp) > 0){
    
    ## sort
    tmp         <- tmp[ order(block, srt)]
    
    ## import the SNP dosages
    snp.dosage  <- mclapply(unique(tmp$block), function(k){
      ## import
      foo        <- fread(paste0("../genetics/", gsub("info", "dosage.transpose", k)), select = c(1, tmp[ block == k]$srt+1))
      ## assign names
      names(foo) <- c("IID", tmp[ block == k]$MarkerName)
      return(foo)
    }, mc.cores = 5)
    ## combine
    snp.dosage  <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), snp.dosage)
    
    ## compute cis score
    if("cis" %in% tmp$`cis/trans`){
      ## get the variables needed
      ii                        <- tmp[ `cis/trans` == "cis"]$MarkerName
      snp.dosage[, "cis.score"] <- as.matrix(snp.dosage[, ii, with=F]) %*% tmp[ `cis/trans` == "cis"]$weight
    }
    
    ## compute trans score
    if("trans" %in% tmp$`cis/trans`){
      ## get the variables needed
      ii                          <- tmp[ `cis/trans` == "trans"]$MarkerName
      snp.dosage[, "trans.score"] <- as.matrix(snp.dosage[, ii, with=F]) %*% tmp[ `cis/trans` == "trans"]$weight
    }
    
    ## rename
    names(snp.dosage) <- gsub("score", paste0("score.", x), names(snp.dosage))
    gc(reset=T)
    ## write to file
    fwrite(snp.dosage[, c("IID", grep("score", names(snp.dosage), value=T)), with=F], paste0("../genetics//", x, ".genetic"), sep = "\t", row.names = F, na = NA)
  }
})

## --> sanity checking <-- ##

## get proteins with scores
jj            <- grep("genetic", dir("../genetics/"), value = T) 

## loop
prot.sanity   <- mclapply(jj, function(x){
  
  ## import genetic
  prot.genetic <- fread(paste0("../genetics/", x))
  
  ## get the protein of interest
  prot         <- gsub("\\.genetic", "", x)
  
  ## merge with protein data
  prot.genetic <- merge(ukb.prot[, c("eid", prot)], prot.genetic, by.x = "eid", by.y = "IID")
  
  ## perform association analysis
  res          <- lapply(names(prot.genetic)[-c(1:2)], function(k){
    
    ## run simple linear model
    ff <- summary(lm(paste(prot, k, sep = " ~ "), prot.genetic))$coefficients[2,, drop=F]
    ## return results
    return(data.table(var=k, ff))
    
  })
  
  return(rbindlist(res))

}, mc.cores = 10)
## combine
prot.sanity   <- rbindlist(prot.sanity)

## for some cis-pQTLs, there is a discrepancy of association...

############################################################################
####                       REVISION - MSB 24/06/2025                    ####
############################################################################

#-----------------------------------#
##-- downsample EUR to match AFR --##
#-----------------------------------#

## create columns to ease matching
ukb.phe[, popAFR := ifelse(pop == "AFR", 1, 0)]
table(ukb.phe$popAFR)

## create subset of relevant populations
tmp          <- ukb.phe[ pop %in% c("EUR", "AFR"), .(f.eid, pop, popAFR, age, sex)]

## match EUR people to AFR
set.seed(42)
controls.AFR <- MatchIt::matchit(as.formula(as.factor(popAFR) ~ age), 
                                 data = tmp,
                                 method = "nearest", exact = ~ sex,
                                 distance = "euclidean", ratio = 4, replace = F)

## get the selected samples
controls.AFR <- apply(controls.AFR$match.matrix, 2, as.numeric)
controls.AFR <- as.data.table(reshape2::melt(controls.AFR))
## assign EID: Var2 contains the batch of matched controls.AFR
controls.AFR[, f.eid := sapply(value, function(x) tmp$f.eid[x])]

## rename
names(controls.AFR) <- gsub("Var2", "controls.AFR", names(controls.AFR))

## add the information to the ukb data
ukb.phe      <- merge(ukb.phe, controls.AFR[, .(f.eid, controls.AFR)], all.x = T)

## cross check
table(ukb.phe$pop, ukb.phe$controls.AFR)

#-----------------------------------#
##-- downsample EUR to match CSA --##
#-----------------------------------#

## create columns to ease matching
ukb.phe[, popCSA := ifelse(pop == "CSA", 1, 0)]
table(ukb.phe$popCSA)

## create subset of relevant populations
tmp          <- ukb.phe[ pop %in% c("EUR", "CSA"), .(f.eid, pop, popCSA, age, sex)]

## match EUR people to CSA
set.seed(42)
controls.CSA <- MatchIt::matchit(as.formula(as.factor(popCSA) ~ age), 
                                 data = tmp,
                                 method = "nearest", exact = ~ sex,
                                 distance = "euclidean", ratio = 4, replace = F)

## get the selected samples
controls.CSA <- apply(controls.CSA$match.matrix, 2, as.numeric)
controls.CSA <- as.data.table(reshape2::melt(controls.CSA))
## assign EID: Var2 contains the batch of matched controls.CSA
controls.CSA[, f.eid := sapply(value, function(x) tmp$f.eid[x])]

## rename
names(controls.CSA) <- gsub("Var2", "controls.CSA", names(controls.CSA))

## add the information to the ukb data
ukb.phe      <- merge(ukb.phe, controls.CSA[, .(f.eid, controls.CSA)], all.x = T)
## cross check
table(ukb.phe$pop, ukb.phe$controls.CSA)

#-----------------------------------#
##-- store file to run analysis  --##
#-----------------------------------#

## store to perform imputation of missing values
fwrite(ukb.phe, "UKB.variables.variance.decomp.imputed.for.revision.20250624.txt", sep="\t", row.names=F, na = NA)


