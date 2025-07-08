######################################################
#### perform additional regression analysis       ####
#### Maik Pietzner                     04/06/2025 ####
######################################################

rm(list=ls())
setwd("<path>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(arrow)
require(doMC)
require(stringr)

#############################################
####      import relevant data sets      ####
#############################################

## Olink data (include relevant confounder)
ukb.olink                     <- fread("../../01_phenotype_prep/data/UKB.Olink.var.decomp.20250128.txt")
## imputed data set of phenotypes
ukb.phe                       <- fread("../../01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt")
## label for explanatory variables
lab.phe                       <- fread("../../02_feature_selection/input/UKB.labels.variance.decomp.20250604.txt")
## reduce to entries in the annotation file
ukb.phe                       <- ukb.phe[, c("f.eid", intersect(lab.phe$short_name, names(ukb.phe))), with=F]

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
## merge with protein data
ukb.phe                       <- merge(ukb.phe, ukb.olink, by.x = "f.eid", by.y = "eid")

#----------------------------#
##--    Label data sets   --##
#----------------------------#

## factor extended labels
lab.ext                       <- fread("../../01_phenotype_prep/data/UKB.labels.variance.decomp.factor.extended.20250129.txt")
## import protein labels
prot.order                    <- fread("../../02_feature_selection/input/Protein.variance.summary.UKB.Olink.20250612.txt")
## import variance explained
res.var.sex                   <- fread("../../02_feature_selection/input/Results.variance.decomposition.UKB.Olink.20250612.txt")

################################################
####         import genetic scores          ####
################################################

## get proteins that are missing
jj                       <- dir("../../01_phenotype_prep/genetics/")
jj                       <- grep(".genetic", jj, value=T)

## import again
registerDoMC(5)
prot.genetics            <- mclapply(jj, function(x){
  ##  import
  tmp <- fread(paste0("../../01_phenotype_prep/genetics/", x))
  ## subset to IIDs in proteomic data
  if(nrow(tmp) > 0) return(tmp[ IID %in% ukb.phe$f.eid])
}, mc.cores=5)
## combine
prot.genetics            <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), prot.genetics)

## write to file for later use
fwrite(prot.genetics, "UKB.olink.proteins.genetic_scores.20250612.txt", sep = "\t", row.names = F, na = NA)

## add to the overall data
ukb.phe                  <- merge(ukb.phe, prot.genetics, by.x = "f.eid", by.y = "IID")

## delete what is no longer needed
rm(ukb.olink); rm(prot.genetics); gc(reset=T)

#############################################
####     perform standard regression     ####
#############################################

#--------------------------------------#
##-- joint model with all variables --##
#--------------------------------------#

## run in parallel
registerDoMC(5)

## run regression models across all selected variabels
res.joint <- expand.grid(protein = prot.order$id, population = c("All", "Male", "Female"), stringsAsFactors = F)
## loop through
res.joint <- mclapply(1:nrow(res.joint), function(x){
  
  print(res.joint[x,])
  
  ## get the relevant variables
  vars <- res.var.sex[ var == res.joint$protein[x] & population == res.joint$population[x] & variable != "Residuals"]$variable
  
  ## drop if residuals
  if(length(vars) > 0){
    
    ## account for genetic scores
    vars <- gsub("score$", paste0("score.", res.joint$protein[x]), vars)
    ## create formula
    ff   <- paste0(res.joint$protein[x], "~", paste(vars, collapse = " + "))
    
    ## run regression model
    if(res.joint$population[x] == "All"){
      ff <- summary(lm(ff, ukb.phe))$coefficients
    }else{
      ff <- summary(lm(ff, ukb.phe[ sex == res.joint$population[x]]))$coefficients
    }
    ## return results
    return(data.table(res.joint[x,], short_name_new=row.names(ff), ff))
  }
  
}, mc.cores = 5)
## combine
res.joint <- rbindlist(res.joint, fill = T)
## drop intercept
res.joint <- res.joint[ short_name_new != "(Intercept)"]
## edit names
res.joint[, short_name_new := gsub("[[:punct:]]| ", "_", short_name_new)]

## write to file
tmp       <- merge(res.var.sex[ variable != "Residuals"], unique(lab.ext[, .(short_name, short_name_new, label_new, reference)]), 
                   by.x = "variable", by.y = "short_name", allow.cartesian = T, all.x = T)
## accout for genetic scores
tmp[, short_name_new := ifelse(is.na(short_name_new), variable, short_name_new)]

## make score names comparable
res.joint[, short_name_edited := sapply(short_name_new, function(x){
  ## only for genetic scores
  if(substr(x, 1, 4) == "cis_"){
    return("cis.score")
  }else if(substr(x, 1, 6) == "trans_"){
    return("trans.score")
  }else{
    return(x)
  }
})]

## add joint results
tmp       <- merge(res.joint, tmp, by.x = c("protein", "population", "short_name_edited"), by.y = c("var", "population", "short_name_edited"), all = T)

## write to file
write.table(tmp, "Results.joint.models.explained.variance.UKB.sex.20250620.txt", sep = "\t", row.names = F)

#############################################
####       test for sex-interaction      ####
#############################################

## identify characteristics that may differ by sex
res.var.sex.split                           <- merge(res.var.sex[ population == "Female"], res.var.sex[ population == "Male"], by=c("var", "variable"), suffixes = c(".female", ".male"), all = T)
## replace NA with 0 to ease downstream computation
res.var.sex.split[is.na(res.var.sex.split)] <- 0
## indicate whether there are possible differences (non-overlapping confidence intervals)
res.var.sex.split[, sex.diff := ifelse(p.025.female > p.975.male | p.025.male > p.975.female, 1, 0)]
## compute difference
res.var.sex.split[, p.50.diff := p.50.female - p.50.male]

#-------------------------------------------# 
##--  test for a sex-interaction effect  --##
#-------------------------------------------# 

## create set to start with
res.sex.interaction <- res.var.sex.split[ variable != "Residuals" & sex.diff == 1, .(var, variable)]

## run in parallel
registerDoMC(10)

## loop through each examples
res.sex.interaction <- mclapply(1:nrow(res.sex.interaction), function(x){
  
  ## replace if cis of trans score
  if(res.sex.interaction$variable[x] %in% c("cis.score", "trans.score")){
    expo <- paste(res.sex.interaction$variable[x], res.sex.interaction$var[x], sep=".")
  }else{
    expo <- res.sex.interaction$variable[x]
  }
  
  ## create formula
  ff   <- paste0(res.sex.interaction$var[x], "~", expo, "*sex")
  print(ff)
  ## run regression model
  ff   <- lm(ff, ukb.phe)
  ## get estimates, but be careful to flag instances w/o interaction term due to sex-specificity
  ff.s <- summary(ff)$coefficients
  ## add empty rows if needed
  if(nrow(ff.s) < length(ff$coefficients)){
    ff.s            <- rbind(ff.s, matrix(data=NA, length(ff$coefficients)-nrow(ff.s), ncol(ff.s)))
    ## set new rownames
    row.names(ff.s) <- names(ff$coefficients)
  }
  rm(ff); gc(reset=T)
  ## return results
  res.int   <- data.table(res.sex.interaction[x,], short_name_new=row.names(ff.s), ff.s)
  ## subset to interaction
  res.int   <- res.int[ grep(":", short_name_new)]
  res.int[, short_name_new := gsub(":sexMale", "", short_name_new, fixed = T)]
  
  ## run analysis in each sex
  ff        <- paste0(res.sex.interaction$var[x], "~", expo)
  ## women
  ff.w      <- summary(lm(ff, ukb.phe[ sex == "Female"]))$coefficients
  ## men
  ff.m      <- summary(lm(ff, ukb.phe[ sex == "Male"]))$coefficients
  
  ## add results
  res.women <- data.table(short_name_new=row.names(ff.w), ff.w)
  res.men   <- data.table(short_name_new=row.names(ff.m), ff.m)
  ## combine (keep track of variable present in only one sex!)
  res.comb  <- merge(res.women, res.men, by="short_name_new", suffixes = c(".female", ".male"), all=T)
  ## merge with interaction results
  res.comb  <- merge(res.int, res.comb, by="short_name_new")
  ## return
  return(res.comb)
  
}, mc.cores = 10)
## combine
res.sex.interaction        <- rbindlist(res.sex.interaction, fill = T)
## edit names
res.sex.interaction[, short_name_new := gsub("[[:punct:]]| ", "_", short_name_new)]
save.image()
## edit column names
names(res.sex.interaction) <- c("short_name_new", "var", "variable", "beta.inter", "se.inter", "tval.inter", "pval.inter", "beta.female", "se.female", "tval.female", "pval.female",
                                "beta.male", "se.male", "tval.male", "pval.male")

## write results to file
write.table(res.sex.interaction[ pval.inter < .05/nrow(res.sex.interaction) ], "Results.sex.interaction.UKB.olink.feature.selection.20250617.txt", sep="\t", row.names = F)

#############################################
####       test for pop-interaction      ####
#############################################

## import results
res.var.pop                                 <- fread("../../02_feature_selection/input/Results.variance.decomposition.UKB.Olink.ancestry.20250616.txt")
## add rank for each variable
tmp                                         <- res.var.pop[ variable != "Residuals"]
tmp                                         <- tmp[ order(population, var, -p.50)]
## create indicator
tmp[, variable.rank := 1:.N, by=c("population", "var")]
## add to main table
res.var.pop                                 <- merge(res.var.pop, tmp[, .(population, var, variable, variable.rank)], by=c("var", "variable", "population"),
                                                     all.x=T)
## identify characteristics that may differ by sex
res.var.pop.split                           <- dcast(res.var.pop, var + variable ~ population, sep=".", value.var = c("n", "n.boot", "variable.rank", grep("p\\.", names(res.var.pop), value=T)))
## replace NA with 0 to ease downstream computation
res.var.pop.split[is.na(res.var.pop.split)] <- 0
## keep only variables within sample size
res.var.pop.split[, variable.rank.non.EUR := pmax(variable.rank.AFR, variable.rank.CSA, na.rm=T)]
summary(res.var.pop.split)
## maximum rank 13

## indicate whether there are possible differences (non-overlapping confidence intervals)
res.var.pop.split[, pop.diff := ifelse(p.025.EUR > p.975.AFR | p.025.AFR > p.975.EUR | p.025.EUR > p.975.CSA | p.025.CSA > p.975.EUR, 1, 0)]
## compute difference
res.var.pop.split[, p.50.diff.EUR.AFR := p.50.EUR - p.50.AFR]
res.var.pop.split[, p.50.diff.EUR.CSA := p.50.EUR - p.50.CSA]

#-------------------------------------------# 
##--  test for a pop-interaction effect  --##
#-------------------------------------------# 

## create set to start with
res.pop.interaction <- res.var.pop.split[ variable != "Residuals" & pop.diff == 1 & (variable.rank.non.EUR != 0 | variable.rank.EUR < variable.rank.non.EUR), .(var, variable)]

## run in parallel
registerDoMC(10)

## loop through each examples
res.pop.interaction <- mclapply(1:nrow(res.pop.interaction), function(x){
  
  ## replace if cis of trans score
  if(res.pop.interaction$variable[x] %in% c("cis.score", "trans.score")){
    expo <- paste(res.pop.interaction$variable[x], res.pop.interaction$var[x], sep=".")
  }else{
    expo <- res.pop.interaction$variable[x]
  }
  
  ## create formula
  ff   <- paste0(res.pop.interaction$var[x], "~", expo, "*pop")
  print(ff)
  ## run regression model (restrict to populations of interest)
  ff   <- lm(ff, ukb.phe[ pop %in% c("EUR", "AFR", "CSA")])
  ## get estimates, but be careful to flag instances w/o interaction term due to pop-specificity
  ff.s <- summary(ff)$coefficients
  ## add empty rows if needed
  if(nrow(ff.s) < length(ff$coefficients)){
    ff.s            <- rbind(ff.s, matrix(data=NA, length(ff$coefficients)-nrow(ff.s), ncol(ff.s)))
    ## set new rownames
    row.names(ff.s) <- names(ff$coefficients)
  }
  rm(ff); gc(reset=T)
  ## return results
  res.int   <- data.table(res.pop.interaction[x,], short_name_new=row.names(ff.s), ff.s)
  ## subset to interaction
  res.int   <- res.int[ grep(":pop", short_name_new)]
  ## split by population
  res.int[, pop := sapply(short_name_new, function(x) str_extract(x, "(?<=pop)(AFR|CSA)"))]
  res.int[, short_name_new := gsub(":popCSA|:popAFR", "", short_name_new, fixed = F)]
  res.int   <- reshape(res.int, idvar = c("var", "variable", "short_name_new"), timevar = "pop", direction = "wide")
  
  ## run analysis in each pop
  ff        <- paste0(res.pop.interaction$var[x], "~", expo)
  ## EUR
  ff.e      <- summary(lm(ff, ukb.phe[ pop == "EUR"]))$coefficients
  ## AFR
  ff.a      <- summary(lm(ff, ukb.phe[ pop == "AFR"]))$coefficients
  ## men
  ff.c      <- summary(lm(ff, ukb.phe[ pop == "CSA"]))$coefficients
  
  ## add results
  res.eur   <- data.table(short_name_new=row.names(ff.e), ff.e)
  res.afr   <- data.table(short_name_new=row.names(ff.a), ff.a)
  res.csa   <- data.table(short_name_new=row.names(ff.c), ff.c)
  ## combine (keep track of variable present in only one pop!)
  res.comb  <- merge(res.csa, res.afr, by="short_name_new", suffixes = c(".CSA", ".AFR"), all=T)
  res.comb  <- merge(res.comb, res.eur, all=T)
  ## merge with interaction results
  res.comb  <- merge(res.int, res.comb, by="short_name_new", suffixes = c(".inter", ".strata"))
  ## return
  return(res.comb)
  
}, mc.cores = 10)
## combine
res.pop.interaction        <- rbindlist(res.pop.interaction, fill = T)
## edit names
res.pop.interaction[, short_name_new := gsub("[[:punct:]]| ", "_", short_name_new)]
save.image()
## edit column names
names(res.pop.interaction) <- c("short_name_new", "var", "variable", 
                                "beta.inter.AFR", "se.inter.AFR", "tval.inter.AFR", "pval.inter.AFR",
                                "beta.inter.CSA", "se.inter.CSA", "tval.inter.CSA", "pval.inter.CSA",
                                "beta.CSA", "se.CSA", "tval.CSA", "pval.CSA",
                                "beta.AFR", "se.AFR", "tval.AFR", "pval.AFR",
                                "beta.EUR", "se.EUR", "tval.EUR", "pval.EUR")

## write results to file
write.table(res.pop.interaction[  pval.inter.AFR < .05/nrow(res.pop.interaction) | pval.inter.CSA < .05/nrow(res.pop.interaction) ], 
            "Results.ancestry.interaction.UKB.olink.feature.selection.20260617.txt", sep="\t", row.names = F)
write.table(res.pop.interaction, "Results.ancestry.interaction.UKB.olink.feature.selection.all.20260617.txt", sep="\t", row.names = F)


#############################################
####      test of different lead SNPs    ####
#############################################

## import list of proteins with evidence for differential cis-effects
cis.prot             <- fread("../../02_feature_selection/input/Proteins.cis.pQTLs.differential.20250617.txt")
## add protein position (build38 from UKB-PPP consortium)
ppp.lab              <- as.data.table(readxl::read_excel("../../02_feature_selection/input/41586_2023_6592_MOESM3_ESM (1).xlsx", 4, skip = 2))
## need to manually tweak some
tmp                  <- ppp.lab[ `Assay Target` %in% toupper(cis.prot$var)]
## FUT3_FUT5
ii                   <- which(tmp$`Assay Target` == "FUT3_FUT5")
tmp$`Gene CHROM`[ii] <- 19
tmp$`Gene start`[ii] <- 5842888
tmp$`Gene end`[ii]   <- 5870540
## FUT3_FUT5
ii                   <- which(tmp$`Assay Target` == "MICB_MICA")
tmp$`Gene CHROM`[ii] <- 6
tmp$`Gene start`[ii] <- 31399784
tmp$`Gene end`[ii]   <- 31511124

## generate file to query variants
write.table(tmp, "cis.pQTL.query.ancestries.txt", sep="\t", row.names = F, quote=F)

#--------------------------------#
##-- import top SNPs for each --##
#--------------------------------#

## get strongest cis-pQTL for each ancestry
cis.ancestry <- lapply(ppp.lab[ `Assay Target` %in% toupper(cis.prot$var)]$`Assay Target`, function(x){
  print(x)
  ## loop over the different ancestries
  res <- lapply(c("EUR", "AFR", "CSA"), function(c){
    ## import
    tmp <- fread(paste(c, x, "cis.region.txt", sep = "."))
    ## get top
    tmp <- tmp[ which.max(tmp$LOG10P),]
    ## return
    return(data.table(tmp, ancestry = c, Assay = x))
  })
  ## return
  return(rbindlist(res))
}) 
## combine
cis.ancestry <- rbindlist(cis.ancestry, fill = T)
## compute minor allele frequency
cis.ancestry[, maf := ifelse(A1FREQ > .5, 1-A1FREQ, A1FREQ)]

## --> map to UKB release <-- ##

## create new MarkerName
cis.ancestry[, pos_hg19 := as.numeric(gsub("^[^:]+:([^:]+).*", "\\1", ID))]
cis.ancestry[, MarkerName := paste(CHROM, pos_hg19, pmin(ALLELE0, ALLELE1), pmax(ALLELE0, ALLELE1), sep="_")]

## import variant statistics -- needs to be recomputed
tmp          <- mclapply(1:nrow(cis.ancestry), function(x){
  
  ## get required variant info (QC'ed)
  snp.stat <- fread(paste0("/sc-resources/ukb/data/bulk/genetic/imputed/variant_qc/", cis.ancestry$ancestry[x],"/ukb_imp_", cis.ancestry$ancestry[x], 
                           "_chr", ifelse(cis.ancestry$CHROM[x] == 23, "X", cis.ancestry$CHROM[x]),"_snpstat.out"), 
                    select = c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB"))
  ## create MarkerName to match
  snp.stat[, MarkerName := paste(ifelse(chromosome == "X", 23, chromosome), position, pmin(alleleA, alleleB), pmax(alleleA, alleleB), sep="_")]
  
  ## add to the variant needed
  return(data.table(cis.ancestry[x,], snp.stat[ MarkerName == cis.ancestry$MarkerName[x]]))
  
  
}, mc.cores = 10)
## combine
tmp          <- rbindlist(tmp)
## remap
cis.ancestry <- tmp

## delete one MarkerName due to poor coding
cis.ancestry[, MarkerName := NULL]

## write to file to check for LD and explained variance
write.table(cis.ancestry, "Ancestry.differential.cis.pQTL.top.cis.pQTL.lookup.20250617.txt", sep = "\t", row.names = F)

#--------------------------------#
##--      import results      --##
#--------------------------------#

## import extended cis-ancestry results
jj <- dir("../output/")
## subset to what of interest
jj <- grep("Ancestry.ld.overlap", jj, value = T)

## import
tmp          <- rbindlist(mclapply(jj, function(x){ fread(paste0("../output/", x))}, mc.cores = 10))
## reassign
cis.ancestry <- tmp

## import interaction results
jj <- dir("../output/")
## subset to what of interest
jj <- grep("Ancestry.interaction", jj, value = T)
## all found

## import
cis.inter    <- rbindlist(mclapply(jj, function(x){ fread(paste0("../output/", x))}, mc.cores = 10))

## import explained variance across ancestries
jj <- dir("../output/")
## subset to what of interest
jj <- grep("Ancestry.explained.variance", jj, value = T)
## all found

## import
cis.expl     <- rbindlist(mclapply(jj, function(x){ 
  ## import
  tmp <- fread(paste0("../output/", x))
  ## add protein target
  tmp[, var := gsub("Ancestry.explained.variance.cis.pQTL\\.|\\.txt", "", x)]
  ## return
  return(tmp)
}, mc.cores = 10))

## write to file to consolidate with feature selection results
write.table(cis.ancestry, "Results.cis.pQTL.ancestries.explained.var.top.ancestral.20250618.txt", sep="\t", row.names = F)
write.table(cis.inter, "Results.cis.pQTL.ancestries.interaction.all.ancestral.20250618.txt", sep="\t", row.names = F)
write.table(cis.expl, "Results.cis.pQTL.ancestries.explained.var.all.ancestral.20250618.txt", sep="\t", row.names = F)
