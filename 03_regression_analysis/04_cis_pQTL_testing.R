#!/usr/bin/env Rscript

## script to extract SNPs to test for differential effects of pQTLs across ancestries
## Maik Pietzner 17/06/2025
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

## correct directory
setwd("<path>/03_regression_analysis/")

## packages needed
require(data.table)
require(doMC)

## import Olink variable to be tested
var.olink <- tolower(args[1])

# ## for testing
# var.olink <- 'il13ra1'

cat("Performing cross-ancestry variance decomposition for", var.olink, "\n")
cat("-------------------------------------------------------------------\n")


#----------------------------------#
##--     import protein data    --##
#----------------------------------#

## import protein data
ukb.olink   <- fread("../01_phenotype_prep/data/UKB.Olink.var.decomp.20250128.txt", select = c("eid", var.olink))
## omit missing values
ukb.olink   <- na.omit(ukb.olink)

## create file for SNP calling
write.table(ukb.olink$eid, paste("tmpdir/ukb.proteomics", var.olink, "sample", sep = "."), row.names = F, col.names = F, quote = F)

#----------------------------------#
##--        generic stats       --##
#----------------------------------#

## import lookup
lookup     <- fread("input/Ancestry.differential.cis.pQTL.top.cis.pQTL.lookup.20250617.txt")
## subset to what is of interest
lookup     <- lookup[ Assay == toupper(var.olink)]

## import cis-regions
cis.region <- lapply(c("EUR", "AFR", "CSA"), function(x){
  ## import
  tmp <- fread(paste0("input/", x, ".", toupper(var.olink), ".cis.region.txt"))
  ## add ancestry
  return(data.table(tmp, ancestry = x))
  })
## combine - N.B. GENPOS is build38!
cis.region <- rbindlist(cis.region)

## convert
cis.region <- dcast(cis.region, CHROM + GENPOS + ID + ALLELE0 + ALLELE1 ~ ancestry, value.var = c("A1FREQ", "BETA", "SE", "LOG10P"), sep = ".")

## create pos in hg19
cis.region[, pos_hg19 := as.numeric(gsub("^[^:]+:([^:]+).*", "\\1", ID))]
cis.region[, MarkerName := paste(CHROM, pos_hg19, pmin(ALLELE0, ALLELE1), pmax(ALLELE0, ALLELE1), sep="_")]

## import UKB rsID mapping
ukb.rsid   <- lapply(c("EUR", "AFR", "CSA"), function(x){
  
  ## get required variant info (QC'ed)
  snp.stat <- fread(paste0("<path>", x,"/ukb_imp_", x, 
                           "_chr", ifelse(lookup$CHROM[1] == 23, "X", lookup$CHROM[1]), "_snpstat.out"), 
                    select = c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB"))
  ## create MarkerName to match
  snp.stat[, MarkerName := paste(ifelse(chromosome == "X", 23, chromosome), position, pmin(alleleA, alleleB), pmax(alleleA, alleleB), sep="_")]
  
  ## add to the variant needed
  return(data.table(snp.stat[ MarkerName %in% cis.region$MarkerName]))
  
})
## combine
ukb.rsid   <- unique(rbindlist(ukb.rsid)) 

## add to cis region
cis.region <- merge(cis.region, ukb.rsid, by = "MarkerName")

## export list of SNPs to be queried from UKB genotype files
write.table(unique(cis.region$rsid), paste("tmpdir/snps", var.olink, "txt", sep="."), col.names = F, row.names = F, quote = F)

## run scripts to extract from autosomes
system(paste("./scripts/obtain_snps.sh", var.olink, ifelse(lookup$CHROM[1] == 23, "X", lookup$CHROM[1])))

## import SNPs
snps        <- fread(paste0("tmpdir/snp.", var.olink, ".dosage.transpose"))
## import SNP information
snps.info   <- fread(paste0("tmpdir/snp.", var.olink, ".info"))
## create MarkerName 
snps.info[, MarkerName := paste(ifelse(chromosome == "X", 23, chromosome), position, pmin(alleleA, alleleB), pmax(alleleA, alleleB), sep = "_")]
snps.info[, snpid := paste0("chr", snps.info$MarkerName)]
## assign names to SNP file
names(snps) <- c("eid", snps.info$snpid)
## remove duplications
jj          <- table(snps.info$snpid)
jj          <- names( jj[jj>1])
## omit from the data
snps        <- snps[, !(names(snps) %in% jj), with=F]
## also from snp info
snps.info   <- snps.info[ !(snpid %in% jj)]

#----------------------------------#
##-- import relevant phenotypes --##
#----------------------------------#

## import imputed ukb data
ukb.phe     <- fread("../01_phenotype_prep/data/UKB.variables.variance.decomp.imputed.dataset.1.20280128.txt", select = c("f.eid", "pop", "sex", "age"))

## clear some space
gc(reset=T)

#----------------------------------#
##--    interaction testing     --##
#----------------------------------#

## bring the relevant data together
ukb.phe     <- merge(ukb.phe, ukb.olink, by.x = "f.eid", by.y = "eid")
## SNPs selected across ancestries
ukb.phe     <- merge(ukb.phe, snps[, c("eid", paste0("chr", unique(lookup$MarkerName))), with=F], by.x = "f.eid", by.y = "eid")
# ## define EUR as reference group
ukb.phe[, pop := factor(pop, levels = c("EUR", "AFR", "CSA"))]

## --> evidence for interaction effects <-- ##

## run through each SNP
registerDoMC(3)

## loop through each examples
res.pop <- data.table(MarkerName = paste0("chr", unique(lookup$MarkerName)), var = var.olink)
res.pop <- lapply(1:nrow(res.pop), function(x){
  
  ## replace if cis of trans score
  expo <- res.pop$MarkerName[x]
  
  ## create formula
  ff   <- paste0(res.pop$var[x], " ~ ", expo, "*pop")
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
  res.int   <- data.table(res.pop[x,], short_name=row.names(ff.s), ff.s)
  ## subset to interaction
  res.int   <- res.int[ grep(":pop", short_name)]
  ## split by population
  res.int[, pop := sapply(short_name, function(x) stringr::str_extract(x, "(?<=pop)(AFR|CSA)"))]
  res.int[, short_name := gsub(":popCSA|:popAFR", "", short_name, fixed = F)]
  res.int   <- reshape(res.int, idvar = c("var", "short_name"), timevar = "pop", direction = "wide")
  
  ## run analysis in each pop
  ff        <- paste0(res.pop$var[x], "~", expo)
  ## EUR
  ff.e      <- summary(lm(ff, ukb.phe[ pop == "EUR"]))$coefficients
  ## AFR
  ff.a      <- summary(lm(ff, ukb.phe[ pop == "AFR"]))$coefficients
  ## men
  ff.c      <- summary(lm(ff, ukb.phe[ pop == "CSA"]))$coefficients
  
  ## add results
  res.eur   <- data.table(MarkerName=row.names(ff.e), ff.e)
  res.afr   <- data.table(MarkerName=row.names(ff.a), ff.a)
  res.csa   <- data.table(MarkerName=row.names(ff.c), ff.c)
  ## combine (keep track of variable present in only one pop!)
  res.comb  <- merge(res.csa, res.afr, by="MarkerName", suffixes = c(".CSA", ".AFR"), all=T)
  res.comb  <- merge(res.comb, res.eur, all=T)
  ## merge with interaction results
  res.comb  <- merge(res.int, res.comb, by.x = "short_name", by.y = "MarkerName", suffixes = c(".inter", ".strata"))
  ## return
  return(res.comb)
  
})
## combine
res.pop        <- rbindlist(res.pop, fill = T)
## edit column names
names(res.pop) <- c("MarkerName", "var",  
                    "MarkerName.AFR", "beta.inter.AFR", "se.inter.AFR", "tval.inter.AFR", "pval.inter.AFR",
                    "MarkerName.CSA", "beta.inter.CSA", "se.inter.CSA", "tval.inter.CSA", "pval.inter.CSA",
                    "beta.CSA", "se.CSA", "tval.CSA", "pval.CSA",
                    "beta.AFR", "se.AFR", "tval.AFR", "pval.AFR",
                    "beta.EUR", "se.EUR", "tval.EUR", "pval.EUR")

## write to file
write.table(res.pop, paste("output/Ancestry.interaction.cis.pQTL", var.olink, "txt", sep = "."), row.names = F, sep = "\t")

#----------------------------------#
##--   compute overlap LD-ovl.  --##
#----------------------------------#

## add population variable to SNP data
snps   <- merge(snps, ukb.phe[, .(f.eid, pop)], by.x = "eid", by.y = "f.eid")

## compute R2 across all SNPs by population!
res.r2 <- lapply(c("EUR", "AFR", "CSA"), function(x){
  
  ## compute squared correlation
  tmp <- cor(snps[ pop == x, snps.info$snpid, with=F], snps[ pop == x, paste0("chr", unique(lookup$MarkerName)), with=F])^2
  ## return
  return(data.table(ancestry=x, MarkerName = snps.info$MarkerName, tmp))
  
})
## combine into one large file
res.r2 <- rbindlist(res.r2)

## --> define whether ancestry SNPs have LD-backbones in common <-- ##

## loop through all ancestries, whether there is a subset of SNPs shared across haplotypes
lookup[, EUR.r2 := apply(lookup[, .(ancestry, MarkerName)], 1, function(x){
  ## get LD block for EUR top signal
  jj <- res.r2[ ancestry == "EUR" & eval(as.name(paste0("chr", lookup[ancestry == "EUR"]$MarkerName))) > .1]$MarkerName
  ## get LD block for top SNP for pop
  ii <- res.r2[ ancestry == x[1] & eval(as.name(paste0("chr", x[2]))) > .1]$MarkerName
  ## report the size of overlap
  return(length(intersect(jj,ii)))
})]

## loop through all ancestries, whether there is a subset of SNPs shared across haplotypes
lookup[, CSA.r2 := apply(lookup[, .(ancestry, MarkerName)], 1, function(x){
  ## get LD block for EUR top signal
  jj <- res.r2[ ancestry == "CSA" & eval(as.name(paste0("chr", lookup[ancestry == "EUR"]$MarkerName))) > .1]$MarkerName
  ## get LD block for top SNP for pop
  ii <- res.r2[ ancestry == x[1] & eval(as.name(paste0("chr", x[2]))) > .1]$MarkerName
  ## report the size of overlap
  return(length(intersect(jj,ii)))
})]

## loop through all ancestries, whether there is a subset of SNPs shared across haplotypes
lookup[, AFR.r2 := apply(lookup[, .(ancestry, MarkerName)], 1, function(x){
  ## get LD block for EUR top signal
  jj <- res.r2[ ancestry == "AFR" & eval(as.name(paste0("chr", lookup[ancestry == "EUR"]$MarkerName))) > .1]$MarkerName
  ## get LD block for top SNP for pop
  ii <- res.r2[ ancestry == x[1] & eval(as.name(paste0("chr", x[2]))) > .1]$MarkerName
  ## report the size of overlap
  return(length(intersect(jj,ii)))
})]

#----------------------------------#
##--   compute explained var.   --##
#----------------------------------#

## loop through all unique SNPs
res.var <- lapply(unique(lookup$MarkerName), function(x){
  ## compute explained variance for each ancestry
  r2 <- lapply(c("EUR", "CSA", "AFR"), function(c){
    summary(lm(paste0( var.olink, " ~ ", paste0("chr", x)), ukb.phe[ pop == c,]))$adj.r.squared
  })
  ## return the results
  return(data.table(MarkerName=x, ancestry=c("EUR", "CSA", "AFR"), r2=unlist(r2)))
})
## combine
res.var <- rbindlist(res.var) 

## add to the look-up data
lookup  <- merge(lookup, res.var)

## write to file
write.table(lookup, paste("output/Ancestry.ld.overlap.cis.pQTL", var.olink, "txt", sep = "."), row.names = F, sep = "\t")
write.table(res.var, paste("output/Ancestry.explained.variance.cis.pQTL", var.olink, "txt", sep = "."), row.names = F, sep = "\t")
