######################################################
#### Protein - Drug - Disease - pQTL network      ####
#### Maik Pietzner                     17/06/2025 ####
######################################################

## set-up
rm(list=ls())
setwd("<path>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(colorspace)
require(igraph)
require(readxl)
library(arrow)

###########################################
####    import the relevant results    ####
###########################################

#--------------------------------#
##--    variance explained    --##
#--------------------------------#

## import explained variance across all, and sexes
res.var.sex <- fread("../../02_feature_selection/input/Results.variance.decomposition.UKB.Olink.20250612.txt")
## import explained varaince across three different ethnicities
res.var.pop <- fread("../../02_feature_selection/input/Results.variance.decomposition.UKB.Olink.ancestry.20250616.txt")
## import protein label (this also contains colouring code - 'cluster.cl' and position in the umap plot - 'umap.1' & 'umap.2')
prot.lab    <- fread("../../02_feature_selection/input/Protein.variance.summary.UKB.Olink.20250612.txt")
## import label for phenotypes (contains label, category and colour code)
phe.lab     <- fread("../../02_feature_selection/input/UKB.labels.variance.decomp.20250612.txt")

#--------------------------------#
##--         pQTL data        --##
#--------------------------------#

## import Supplemental Table with all cis/trans pQTLs (from: https://www.nature.com/articles/s41586-023-06592-6)
pqtl     <- as.data.table(read_excel("../../02_feature_selection/input/41586_2023_6592_MOESM3_ESM (1).xlsx", 11, skip = 4))
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
pqtl[, other_allele := sapply(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, function(x) strsplit(x, ":")[[1]][3])]
## import UKB file to map possible missing rsIDs
ukb.snps <- fread("<path>")
## subset to IDs needed
ukb.snps <- ukb.snps[ MarkerName %in% pqtl$MarkerName]
gc(reset=T)
## replace rsIDs accordingly
pqtl[, rsID := sapply(MarkerName, function(x){
  ## test whether included in the data
  if(x %in% ukb.snps$MarkerName){
    ukb.snps$rsid[which(ukb.snps$MarkerName == x)]
  }else{
    return("")
  }
})]
## drop SNPs still not found
pqtl        <- pqtl[ rsID != ""]
## convert to single column
pqtl[, rsID := sapply(rsID, function(x) x[1])]
## create simplified table to map alleles
pqtl[, alleles := paste(pmin(other_allele, effect_allele), pmax(other_allele, effect_allele), sep="_")]

## subset to pQTLs that map to proteins associated with drugs or disease
pqtl        <- pqtl[ protein_id %in% unique(res.var.sex[ variable %in% phe.lab[ category %in% c("Diseases", "Drugs")]$short_name]$var)]

## export to generate a list of proxy variables
for(j in 1:23){
  ## write all rsIDs into file
  write.table(unique(pqtl[ CHROM == j]$rsID), paste("snp.list", j, "txt", sep="."), col.names = F, row.names = F, quote = F)
}

#----------------------------#
##--  import proxy list   --##
#----------------------------#

## import: N.B.: SNP pairs are only listed once!
r2.proxies <- rbindlist(lapply(1:23, function(x) fread(paste0("ld.proxies.", x,".vcor"))))
## edit naming
names(r2.proxies) <- gsub("#", "", names(r2.proxies))

#-------------------------------------#
##--      lift over to build 38    --##
#-------------------------------------#

## import file
lift.over  <- fread("<path>", select = c("chr_hg19", "pos_hg19", "pos_hg38"))
lift.over  <- unique(lift.over[, .(chr_hg19, pos_hg19, pos_hg38)])
## align chromosome naming
lift.over[, chr_hg19 := gsub("chr", "", chr_hg19)]

## map to proxy data
r2.proxies <- merge(r2.proxies, lift.over, by.x = c("CHROM_A", "POS_A"), by.y = c("chr_hg19", "pos_hg19"), all.x=T) 
r2.proxies <- merge(r2.proxies, lift.over, by.x = c("CHROM_B", "POS_B"), by.y = c("chr_hg19", "pos_hg19"), all.x=T, suffixes = c("_A", "_B")) 

## drop
rm(lift.over); gc(reset = T)

#-----------------------------#
##--   create r2 groups    --##
#-----------------------------#

## get all unique variants from regional results; based on r2>0.6
pqtl.var                 <- unique(pqtl[, .(MarkerName, rsID, CHROM, `GENPOS (hg38)`, effect_allele, other_allele)])
## n = 9587 unique variants

## position in build 38
pqtl.var[, pos_hg19 := sub(".*:(\\d+)_.*", "\\1", MarkerName)]

## add regional sentinels
tmp        <- pqtl.var[, .(CHROM, `GENPOS (hg38)`, rsID, pos_hg19)]
tmp        <- cbind(tmp, tmp)
tmp[, UNPHASED_R2 := 1]
names(tmp) <- c("CHROM_A", "pos_hg38_A", "ID_A", "POS_A", "CHROM_B", "POS_B", "ID_B", "pos_hg38_B", "UNPHASED_R2")
## combine with proxies
r2.proxies <- rbind(r2.proxies, tmp)
## exchange chromosome coding
r2.proxies[, CHROM_A := as.numeric(ifelse(CHROM_A == "X", 23, CHROM_A))]
r2.proxies[, CHROM_B := as.numeric(ifelse(CHROM_B == "X", 23, CHROM_B))]

## convert to graph (use max to allow for edge weights)
ld.sub     <- graph_from_data_frame(r2.proxies[ UNPHASED_R2 >= .6, .(ID_A, ID_B, UNPHASED_R2)])
## get all separate components
ld.sub     <- components(ld.sub)$membership
## convert to data frame
ld.sub     <- data.table(rsID=names(ld.sub), R2.group=ld.sub)

## add to results
pqtl       <- merge(pqtl, ld.sub, by = "rsID", all.x=T)

###########################################
####        overlap GWAS catalog       ####
###########################################

#-------------------------------------#
##--       import GWAS catalog     --##
#-------------------------------------#

## import latest release (22/05/2025)
gwas.catalogue           <- fread("<path>/input/alternative")

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- gwas.catalogue[ !is.na(`OR or BETA`) & is.finite(`OR or BETA`) ]
## generate risk allele and drop everything w/o this information
gwas.catalogue[, riskA := sapply(`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])] 
gwas.catalogue[, riskA := trimws(riskA, which = "b")] 
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- gwas.catalogue[ PVALUE_MLOG > 7.3 ]
## N = 420,745 entries

## create another entry to possible merge on (careful, genome build 38 mapping)
gwas.catalogue[, snp.id := paste0(ifelse(CHR_ID == "X", 23, CHR_ID), ":", CHR_POS)]

#-------------------------------------#
##-- map variants to the GWAS cat. --##
#-------------------------------------#

## do in parallel
registerDoMC(10)

## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
res.gwas                 <- mclapply(1:nrow(pqtl.var), function(x){
  
  print(x)
  
  ## get all possible proxies (careful; does not include the SNP itself)
  snp               <- r2.proxies[ (ID_A == pqtl.var$rsID[x] | ID_B == pqtl.var$rsID[x]) & UNPHASED_R2 >= .8]
  
  ## convert back to easy format
  snp[, lead.rsID := pqtl.var$rsID[x]]
  snp[, proxy.rsID := ifelse(ID_A == pqtl.var$rsID[x], ID_B, ID_A)]
  snp[, pos.proxy := ifelse(ID_A == pqtl.var$rsID[x], pos_hg38_B, pos_hg38_A)]
  ## subset to minimum needed
  snp              <- snp[, .(lead.rsID, proxy.rsID, CHROM_A, pos.proxy, UNPHASED_R2)]
  
  ## create snp id to optimize merging multiple mappings
  snp[, snp.id := paste0(CHROM_A, ":", pos.proxy)]
  
  ## create two versions of mapping
  snp.rsid          <- merge(snp, gwas.catalogue, by.x="proxy.rsID", by.y="SNPS")
  snp.pos           <- merge(snp, gwas.catalogue, by = "snp.id")
  ## edit
  snp.rsid$snp.id   <- snp.rsid$snp.id.x
  snp.rsid$snp.id.x <- snp.rsid$snp.id.y <- NULL
  snp.pos$SNPS      <- snp.pos$proxy.rsID
  snp.rsid$SNPS     <- snp.rsid$proxy.rsID
  ## combine
  snp               <- unique(rbind(snp.rsid, snp.pos))
  
  ## prepare return
  if(nrow(snp) > 0){
    ## sort 
    snp              <- snp[order(TRAIT, lead.rsID, -UNPHASED_R2)]
    ## create indicator
    snp[, ind := 1:.N, by=c("TRAIT", "lead.rsID")]
    ## keep only one finding per trait
    snp              <- snp[ind == 1]
    
    ## report summary back
    snp              <- data.table(rsid.gwas=paste(sort(unique(snp$proxy.rsID)), collapse = "||"),
                                   trait_reported=paste(sort(unique(snp$TRAIT)), collapse = "||"),
                                   mapped_trait=paste(sort(unique(snp$MAPPED_TRAIT)), collapse = "||"),
                                   mapped_trait_efo=paste(sort(unique(snp$MAPPED_TRAIT_URI)), collapse = "||"),
                                   study_id=paste(sort(unique(snp$`STUDY ACCESSION`)), collapse = "||"),
                                   source_gwas=paste(sort(unique(snp$PUBMEDID)), collapse = "||"),
                                   num_reported=nrow(snp))
    
  }else{
    
    ## report summary back
    snp              <- data.table(rsid.gwas="",
                                   trait_reported="",
                                   mapped_trait="",
                                   mapped_trait_efo="",
                                   study_id="",
                                   source_gwas="",
                                   num_reported=0)
  }
  
  ## return data set
  return(data.table(pqtl.var[x,], snp))
}, mc.cores=10)
## combine 
res.gwas                 <- rbindlist(res.gwas)

## add r2 groups
res.gwas                 <- merge(res.gwas, ld.sub, by = "rsID", all.x=T)

###########################################
####       map ATC codes to targets    ####
###########################################

#---------------------------#
##--   map ATC to Chembl --##
#---------------------------#

## import Chembl to ATC map
atc.chembl <- fread("<path>/atc_to_chembl_id_34_map.tsv")

## add to phe.lab (N.B.: some ATC codes are higher level --> match all of those)
phe.lab[, chembl.id := sapply(short_name, function(x){
  ## if short code
  if(nchar(x) < 7){
    ## get all possible codes
    return(paste(unique(atc.chembl[ substr(level5, 1, nchar(x)) == x]$chembl_id), collapse = "|"))
  }else{
    return(paste(unique(atc.chembl[ level5 == x]$chembl_id), collapse = "|"))
  }}
)]

## how many are missing
nrow(phe.lab[ category == "Drugs" & chembl.id == ""])
## n = 77 missing --> fill in manually

## write to file
write.table(phe.lab[ category == "Drugs" & chembl.id == ""], "missing.chemblids.atc.codes.txt", sep="\t", row.names=F)
## import manually annotated ones based on search at Chembl, OpenTargets, and Drug Bank
tmp.chembl <- fread("missing.chemblids.atc.codes.annotated.20250619.txt")
## replace
phe.lab[, chembl.id := apply(phe.lab[, .(short_name, chembl.id)], 1, function(x) ifelse(x[1] %in% tmp.chembl$short_name, tmp.chembl[ short_name == x[1]]$chembl.id, x[2]))]
## missing ones are not mappable

#----------------------------#
##-- map Chembl to target --##
#----------------------------#

## path to results
ot.drugs <- dir("<path>")
ot.drugs <- grep("part", ot.drugs, value=T)
ot.drugs <- lapply(ot.drugs, function(x) read_parquet(paste0("<path>", x)))
ot.drugs <- lapply(ot.drugs, as.data.table)
ot.drugs <- rbindlist(ot.drugs)
## transform some columns to make life easier
ot.drugs[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "|"))]
ot.drugs[, crossReferences := sapply(crossReferences, function(x) paste(x, collapse = "|"))]
ot.drugs[, tradeNames := sapply(tradeNames, function(x) paste(x, collapse = "|"))]
ot.drugs[, childChemblIds := sapply(childChemblIds, function(x) paste(x, collapse = "|"))]
ot.drugs[, linkedDiseases.rows := sapply(linkedDiseases.rows, function(x) paste(x, collapse = "|"))]
ot.drugs[, linkedTargets.rows := sapply(linkedTargets.rows, function(x) paste(x, collapse = "|"))]

## add gene targets
phe.lab[, drug.target.ensembl := sapply(chembl.id, function(x){
  ## do only if any
  if(x != ""){
    ## separate out
    x  <- strsplit(x, "\\|")[[1]]
    ## find all matching entries
    ii <- lapply(unique(ot.drugs[ id %in% x]$linkedTargets.row), function(c) strsplit(c, "\\|")[[1]])
    return(paste(unique(unlist(ii)), collapse = "|"))
  }else{
    return("")
  }
})]

## get all protein coding genes: from https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/ build 38
tmp.genes <- rtracklayer::readGFF("<path>/Homo_sapiens.GRCh38.112.gtf.gz")
tmp.genes <- as.data.table(tmp.genes)
## unique protein coding genes
tmp.genes <- unique(tmp.genes[ gene_biotype == "protein_coding", .(seqid, start, end, gene_id, gene_name)])
tmp.genes[, chromosome := as.numeric(ifelse(seqid == "X", 23, as.character(seqid)))]
tmp.genes <- tmp.genes[ !is.na(chromosome)]
tmp.genes <- tmp.genes[, .(chrom=paste(unique(chromosome), collapse = "|"),
                           start=paste(start, collapse = "|"),
                           end=paste(end, collapse = "|")), 
                       by=c("gene_id", "gene_name")]

## add gene names
phe.lab[, drug.target.gene := sapply(drug.target.ensembl, function(x){
  ## do only if any
  if(x != ""){
    ## separate out
    x  <- strsplit(x, "\\|")[[1]]
    ## find all matching entries
    return(paste(unique(tmp.genes[ gene_id %in% x]$gene_name), collapse = "|"))
  }else{
    return("")
  }
})]

###########################################
####       map disease to targets      ####
###########################################

#------------------------------#
##-- phecode to disease ids --##
#------------------------------#

## import all entries and generate one large table
ot.diseases <- dir("../ot_download//diseases/")
ot.diseases <- grep("part", ot.diseases, value=T)
ot.diseases <- lapply(ot.diseases, function(x) read_parquet(paste0("../ot_download//diseases/", x)))
ot.diseases <- lapply(ot.diseases, as.data.table)
ot.diseases <- rbindlist(ot.diseases)
## transform some columns to make life easier
ot.diseases[, dbXRefs := sapply(dbXRefs, function(x) paste(x, collapse = "|"))]
ot.diseases[, directLocationIds := sapply(directLocationIds, function(x) paste(x, collapse = "|"))]
ot.diseases[, obsoleteTerms := sapply(obsoleteTerms, function(x) paste(x, collapse = "|"))]
ot.diseases[, parents := sapply(parents, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasBroadSynonym := sapply(synonyms.hasBroadSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasExactSynonym := sapply(synonyms.hasExactSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasNarrowSynonym := sapply(synonyms.hasNarrowSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasRelatedSynonym := sapply(synonyms.hasRelatedSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, descendants := sapply(descendants, function(x) paste(x, collapse = "|"))]
ot.diseases[, children := sapply(children, function(x) paste(x, collapse = "|"))]
ot.diseases[, therapeuticAreas := sapply(therapeuticAreas, function(x) paste(x, collapse = "|"))]
ot.diseases[, indirectLocationIds := sapply(indirectLocationIds, function(x) paste(x, collapse = "|"))]
ot.diseases[, ancestors := sapply(ancestors, function(x) paste(x, collapse = "|"))]
## n = 28,327 diseases

## harmonize codings
ot.diseases[, dbXRefs := gsub("\\.", "", dbXRefs)]

#-----------------------------#
##-- map OT ID to phecodes --##
#-----------------------------#

## import phecode mapping
lab.phecodes <- fread("<path>")

## create an unpacked set of ot.diseases (inflate)
tmp.ot       <- rbindlist(mclapply(1:nrow(ot.diseases), function(x) return(data.table(id=ot.diseases$id[x], db.ref = strsplit(ot.diseases$dbXRefs[x], "\\|")[[1]]))))

## add to phecode label
phe.lab[, ot.disease.mapping := sapply(short_name, function(x){
  
  ## get relevant code
  jj      <- which(lab.phecodes$id == gsub("bin_", "date_", x))
  
  ## do only if any mapping
  if(length(jj) > 0){
    ## get all relevant codes
    icd10   <- sapply(strsplit(lab.phecodes$icd10[jj], "\\|")[[1]], function(c) paste0("ICD10:", c))
    snomed1 <- sapply(strsplit(lab.phecodes$snomed[jj], "\\|")[[1]], function(c) paste0("SNOMEDCT:", c))
    snomed2 <- sapply(strsplit(lab.phecodes$snomed[jj], "\\|")[[1]], function(c) paste0("SCTID:", c))
    
    ## return relevant mappings
    return(paste(unique(tmp.ot[ db.ref %in% c(icd10, snomed1, snomed2)]$id), collapse = "|"))
  }else{
    return("")
  }
  
})]

## import phecode to EFO assignment
phecode.efo <- fread("<path>")

###########################################
####        create simple network      ####
###########################################

#--------------------------#
##--  drug --> protein  --##
#--------------------------#

## create simple edge data set
edges.drug.protein    <- res.var.sex[ variable %in% phe.lab[ category == "Drugs"]$short_name & p.025 > 0, .(p.50.max = max(p.50)), by=c("variable", "var")]
## add attributes (e.g., whether already a target)
edges.drug.protein[, drug.target := apply(edges.drug.protein, 1, function(x){
  ## get all targets of the drug
  tmp <- phe.lab[ short_name == x[1]]$drug.target.ensembl
  ## search in the protein label
  if(tmp != ""){
    nrow(prot.lab[ ensembl_gene_id %in% grep(tmp, prot.lab$ensembl_gene_id, value = T) & id == x[2]])
  }else{
    ## return blank
    return("")
  }
})]
## how many unique drugs
length(unique(edges.drug.protein$variable)) ## n = 194
length(unique(edges.drug.protein$var)) ## n = 1,494

#--------------------------#
##--  snp --> protein   --##
#--------------------------#

## create simple set, but use R2 groups as anchor; include only proteins related to drugs or diseases (n = 1718)
edges.snp.protein     <- unique(res.var.sex[ variable %in% phe.lab[ category %in% c("Drugs", "Diseases")]$short_name & p.025 > 0, .(p.50.max = max(p.50)), by=c("variable", "var")]$var)
edges.snp.protein     <- pqtl[ protein_id %in% edges.snp.protein, .(R2.group, protein_id, `Bioinfomatic annotated gene`, `cis/trans`, `log10(p)`)]
## still n = 18,355 edges
length(unique(edges.snp.protein$R2.group)) ## n = 5,330
length(unique(edges.snp.protein$protein_id)) ## n = 1,731

#--------------------------#
##--disease --> protein --##
#--------------------------#

## tmp data to ease mapping on whether proteins are drug targets for associated diseases
tmp.drugs            <- rbindlist(lapply(1:nrow(ot.drugs), function(x){
  print(x)
  return(expand.grid(id=ot.drugs$id[x], 
                     linkedTargets=strsplit(ot.drugs$linkedTargets.rows[x], "\\|")[[1]],
                     linkedDiseases=strsplit(ot.drugs$linkedDiseases.rows[x], "\\|")[[1]],
                     stringsAsFactors = F))
}))

## create set of edges based on explained variance
edges.disease.protein <- res.var.sex[ variable %in% phe.lab[ category == "Diseases"]$short_name & p.025 > 0, .(p.50.max = max(p.50)), by=c("variable", "var")]
## add attributes (e.g., whether already a target)
edges.disease.protein[, drug.target := apply(edges.disease.protein, 1, function(x){
  ## get all disease mappings 
  tmp    <- phe.lab[ short_name == x[1]]$ot.disease.mapping
  ## get the ensembl gene assignment for the protein
  esembl <- strsplit(prot.lab[ id == x[2]]$ensembl_gene_id, "\\|")[[1]]
  ## map across the OT data
  if(tmp != ""){
    ## get all disease mappings 
    tmp    <- strsplit(tmp, "\\|")[[1]]
    ## return
    return(nrow(tmp.drugs[ linkedTargets %in% esembl & linkedDiseases %in% tmp]))
  }else{
    ## return blank
    return("")
  }
})]

#--------------------------#
##--   snp --> disease  --##
#--------------------------#

## map EFO terms more broadly to phecodes
icd10.efo             <- read.delim("<path>/UK_Biobank_master_file.tsv", sep= "\t", fill = T)
icd10.efo             <- as.data.table(icd10.efo)

## add to label
phe.lab[, efo.mapping := sapply(short_name, function(x){
  
  ## test if looking at a phecode
  if(gsub("bin_", "date_", x) %in% lab.phecodes$id){
    ## get the relevant icd-10 codes
    icd10 <- strsplit(lab.phecodes[ which(id == gsub("bin_", "date_", x))]$icd10, "\\|")
    ## shrink to three characters
    icd10 <- sapply(icd10, function(c) substr(c, 1, 3))
    print(icd10)
    ## return all matching EFO terms
    return(paste(icd10.efo[ICD10_CODE.SELF_REPORTED_TRAIT_FIELD_CODE %in% icd10]$MAPPED_TERM_URI, collapse = "|"))
  }else{
    return("")
  }
  
})]


## create a list of SNP - outcome associations based on GWAS catalog look-up
edges.snp.disease     <- res.gwas[ num_reported > 0]
## summarize by R2 group
edges.snp.disease     <- rbindlist(lapply(unique(edges.snp.disease$R2.group), function(x){
  
  ## get all relevant findings
  traits <- unlist(lapply(edges.snp.disease[R2.group == x]$mapped_trait_efo, function(c){
    strsplit(c, "\\|\\||, ")[[1]]
  }))
  
  ## keep only those consistently linked to all members of the r2 group
  traits <- table(traits)
  traits <- names(traits[ traits >= nrow(edges.snp.disease[R2.group == x])])
  
  ## return info
  return(data.table(R2.group = x, mapped_trait_efo = traits))
  
})) 
## contain only entries mapping to phecodes
edges.snp.disease[, efo.term := gsub(".*/(.*)", "\\1", mapped_trait_efo)]
## generate all terms mapping
tmp                   <- unique(c(unlist(lapply(phe.lab$ot.disease.mapping, function(x) strsplit(x, "\\|")[[1]])),
                                  unlist(lapply(phe.lab$efo.mapping, function(x) strsplit(x, "\\|")[[1]]))))
edges.snp.disease     <- edges.snp.disease[ efo.term %in% tmp ]
## n = 2,884
## add back disease names ('phecodes'); create dummy data frame first to ease mapping
tmp                   <- phe.lab[ category == "Diseases" ]
tmp                   <- rbindlist(lapply(1:nrow(tmp), function(x) data.table(short_name=tmp$short_name[x], 
                                                                              efo.all=unique(c(strsplit(tmp$ot.disease.mapping[x], "\\|")[[1]],
                                                                                               strsplit(tmp$efo.mapping[x], "\\|")[[1]])))))
## now add phecodes back in
edges.snp.disease     <- rbindlist(lapply(1:nrow(edges.snp.disease), function(x){
  ## get the relevant EFO code
  return(data.table(edges.snp.disease[x,], short_name = unique(tmp[ efo.all %in% edges.snp.disease$efo.term[x]]$short_name)))
}))
## expand for multiple mappings
edges.snp.disease     <- unique(edges.snp.disease)

## restrict to diseases included in the var explained data
edges.snp.disease     <- edges.snp.disease[ short_name %in% res.var.sex[ variable %in% phe.lab[ category == "Diseases"]$short_name & p.025 > 0]$variable ]

## prun many to many mappings
edges.snp.disease[, ind := 1:.N, by = c("R2.group", "short_name")]
## drop 
edges.snp.disease     <- edges.snp.disease[ ind == 1]
edges.snp.disease[, ind := NULL]

#--------------------------#
##--   snp --> effector --##
#--------------------------#

## create simple edges of SNP to effector gene
edges.snp.gene        <- unique(pqtl[ , .(R2.group, `Bioinfomatic annotated gene`)]) 
## consolidate effector gene by locus
edges.snp.gene        <- rbindlist(lapply(unique(edges.snp.gene$R2.group), function(x){
  
  ## get the relevant genes
  tmp <- edges.snp.gene[ R2.group == x]$`Bioinfomatic annotated gene`
  ## drop missing assignments
  tmp <- tmp[!(tmp == "-")]
  ## table
  if(length(tmp) > 0){
    ## table results
    tmp <- table(tmp)
    ## return results
    return(data.table(R2.group=x, effecto.gene=names(tmp[which.max(tmp)])))
  }else{
    ## return blank
    return(data.table(R2.group=x, effecto.gene="-"))
  }
  
  
}))

#--------------------------#
##--   drug --> target  --##
#--------------------------#

## create edges based on expanded label mapping
edges.drug.gene       <- phe.lab[ drug.target.gene != "",  .(short_name, drug.target.gene)]
## inflate
edges.drug.gene       <- rbindlist(lapply(1:nrow(edges.drug.gene), function(x) return(data.table(short_name = edges.drug.gene$short_name[x],
                                                                                                 target.geme = strsplit(edges.drug.gene$drug.target.gene[x], "\\|")[[1]]))))
#--------------------------#
##-- bring it together  --##
#--------------------------#

head(edges.drug.protein, 1)
head(edges.disease.protein, 1)
head(edges.snp.protein, 1)
head(edges.snp.disease, 1)

## --> tweak each data set <-- ##

## DRUG --> PROTEIN
edges.drug.protein[, anno.3 := ""]
names(edges.drug.protein)    <- c("node1", "node2", "anno.1", "anno.2", "anno.3")
edges.drug.protein[, type := "drug-protein"]

## DISEASE --> PROTEIN
edges.disease.protein[, anno.3 := ""]
names(edges.disease.protein) <- c("node1", "node2", "anno.1", "anno.2", "anno.3")
edges.disease.protein[, type := "disease-protein"]

## SNP --> PROTEIN (make sure to keep variable type of annotation columns consistent; anno.1 == 'numeric')
names(edges.snp.protein)     <- c("node1", "node2", "anno.2", "anno.3", "anno.1")
edges.snp.protein[, type := "snp-protein"]

## SNP --> DISEASE
edges.snp.disease[, anno.1 := -1]
names(edges.snp.disease)     <- c("node1", "anno.2", "anno.3", "node2", "anno.1")
edges.snp.disease[, type := "snp-disease"]

## SNP --> GENE
edges.snp.gene[, anno.1 := -1]
edges.snp.gene[, anno.2 := ""]
edges.snp.gene[, anno.3 := ""]
names(edges.snp.gene)        <- c("node1", "node2", "anno.1", "anno.2", "anno.3")
edges.snp.gene[, type := "snp-gene"]

## DRUG --> GENE
edges.drug.gene[, anno.1 := -1]
edges.drug.gene[, anno.2 := ""]
edges.drug.gene[, anno.3 := ""]
names(edges.drug.gene)        <- c("node1", "node2", "anno.1", "anno.2", "anno.3")
edges.drug.gene[, type := "drug-gene"]

## --> combine <-- ##

## combine all into one large edge file
edges.network                <- rbind(edges.drug.protein, 
                                      edges.disease.protein, 
                                      edges.snp.disease, 
                                      edges.snp.protein,
                                      edges.snp.gene,
                                      edges.drug.gene)

## prune some redundant edges
edges.network                <- edges.network[ order(node1, node2, -anno.1, anno.3)]
edges.network[, ind := 1:.N, by=c("node1", "node2")]
edges.network                <- edges.network[ ind == 1]

## --> node label data <-- ##

## create node labels
nodes.network                <- data.table(node = unique(c(edges.drug.protein$node1, edges.drug.protein$node2,
                                                           edges.disease.protein$node1, edges.disease.protein$node2,
                                                           edges.snp.disease$node1, edges.snp.disease$node2,
                                                           edges.snp.protein$node1, edges.snp.protein$node2,
                                                           edges.snp.gene$node1, edges.snp.gene$node2,
                                                           edges.drug.gene$node1, edges.drug.gene$node2)))
## define the type of node
nodes.network[, type := sapply(node, function(x){
  if(x %in% phe.lab$short_name){
    return(phe.lab[ short_name == x]$category)
  }else if(x %in% prot.lab$id){
    return("Protein")
  }else if(x %in% pqtl$R2.group){
    return("SNP")
  }else{
    return("Gene")
  }
})]

## add a label
nodes.network[, label := sapply(node, function(x){
  if(x %in% phe.lab$short_name){
    return(phe.lab[ short_name == x]$label)
  }else if(x %in% prot.lab$id){
    return(paste0("P::", prot.lab[ id == x]$Assay))
  }else if(x %in% pqtl$R2.group){
    ## get all possible genes
    tmp <- table(pqtl[ R2.group == x]$rsID)
    return(names(tmp[which.max(tmp)]))
  }else{
    return(paste0("G::", x))
  }
})]

#--------------------------#
##-- prune the network  --##
#--------------------------#

## get edge type
nodes.network[, edge.set := sapply(node, function(x) paste(sort(unique(edges.network[ node1 == x | node2 == x]$type)), collapse = "|"))]

###########################################
####          create network           ####
###########################################

## create the network from a data frame
network.graph <- graph_from_data_frame(edges.network, directed = F, vertices = nodes.network)

## drop nodes with node degree of one recursively
n.old <- 0; n.new <- 1
while(n.old != n.new){
  ## assign
  n.old         <- vcount(network.graph)
  ## drop nodes
  network.graph <- subgraph(network.graph, which(degree(network.graph) > 1))
  ## report
  print(ecount(network.graph))
  print(vcount(network.graph))
  ## assign
  n.new         <- vcount(network.graph)
}
ecount(network.graph) ## 30,068
vcount(network.graph) ## 6,591

#----------------------#
##--  circle length --##
#----------------------#

## write the network to file
save(network.graph, file = "Protein.SNP.Drug.Disease.network.20250619.RData")
## write edges and nodes to file as well
write.table(as_data_frame(network.graph, what = "edges"), "Protein.SNP.Drug.Disease.network.edges.20250619.txt", sep="\t", row.names=F, quote = F)
## add colour code
tmp <- as_data_frame(network.graph, what = "vertices")
## create colour vector
col.network        <- c("#A6CEE3", "#1F78B4", "#E31A1C", "#FDBF6F", "#FF7F00")
names(col.network) <- c("Diseases", "Drugs", "Protein", "SNP", "Gene")
## add 
tmp$colour         <- col.network[ tmp$type]
write.table(tmp, "Protein.SNP.Drug.Disease.network.nodes.20250619.txt", sep="\t", row.names=F, quote = F)

#----------------------------------------#
##-- Protein - Disease - SNP triplets --##
#----------------------------------------#

## Protein --> Disease --> SNP
prot.snp.disease <- merge(edges.disease.protein, edges.snp.disease, 
                          by.x="node1", by.y="node2", suffixes = c(".protein", ".snp"),
                          allow.cartesian = T)
## SNP --> Protein 
prot.snp.disease <- merge(prot.snp.disease, edges.snp.protein, 
                          by.x=c("node2", "node1.snp"), 
                          by.y=c("node2", "node1"), 
                          suffixes = c(".protein", ".disease"),
                          allow.cartesian = T)
## add labels for diseases
prot.snp.disease <- merge(prot.snp.disease, phe.lab, by.x = "node1", by.y = "short_name")

## collapse and perform enrichment analysis
prot.snp.disease.trans <- prot.snp.disease[ anno.3 == "trans", .(common.proteins = paste(unique(node2), collapse = "|"),
                                                                 number.shared = length(unique(node2))), 
                                            by = c("node1", "node1.snp")]
## add label
prot.snp.disease.trans  <- merge(prot.snp.disease.trans, nodes.network, by.x = "node1", by.y = "node")
prot.snp.disease.trans[, node1.snp := as.character(node1.snp)]
prot.snp.disease.trans  <- merge(prot.snp.disease.trans, nodes.network, by.x = "node1.snp", by.y = "node", suffixes = c(".disease", ".snp"))

## test for enrichment
registerDoMC(10)
prot.snp.disease.trans  <- mclapply(1:nrow(prot.snp.disease.trans), function(x){
  
  ## get disease
  disease <- prot.snp.disease.trans$node1[x]
  ## get SNP
  snp     <- prot.snp.disease.trans$node1.snp[x]
  ## proteins common to both
  prot <- strsplit(prot.snp.disease.trans$common.proteins[x], "\\|")[[1]]
  
  ## --> compute enrichment statistics <-- ##
  
  ## define protein back ground for population
  prot.back <- prot.lab$id
  
  ## protein related to disease and SNP
  d1        <- length(prot)
  ## protein not related to disease but SNP
  d2        <- sum( !(unique(pqtl[ R2.group == snp]$protein_id) %in% prot))
  ## protein related to disease but not SNP
  d3        <- sum( !(unique(res.var.sex[ variable == disease ]$var) %in% prot))
  ## protein not related to disease and not SNP
  d4        <- length(prot.back[ !(prot.back %in% unique(pqtl[ R2.group == snp]$protein_id) | prot.back %in% unique(res.var.sex[ variable == disease]$var))]) 
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.table(prot.snp.disease.trans[x,], or=enr$estimate, pval=enr$p.value, 
                    d1=d1, d2=d2, d3=d3, d4=d4))
  
}, mc.cores = 10)
## combine again
prot.snp.disease.trans  <- rbindlist(prot.snp.disease.trans)

## add GWAS associations more broadly to each SNP
prot.snp.disease.trans[, snp.gwas := sapply(node1.snp, function(x){
  ## get all possible GWAS consequences
  traits <- unlist(lapply(res.gwas[R2.group == x]$mapped_trait, function(c){
    strsplit(c, "\\|\\||, ")[[1]]
  }))
  
  ## keep only those consistently linked to all members of the r2 group
  traits <- table(traits)
  traits <- names(traits[ traits >= nrow(res.gwas[R2.group == x])])
  
  ## return info
  return(paste(traits, collapse = "|"))
  
  
})]
## add pQTL gene annotation
prot.snp.disease.trans[, pqtl.gene := sapply(node1.snp, function(x){
  ## get all possible GWAS consequences
  genes <- sort(table(pqtl[ R2.group == x]$`Bioinfomatic annotated gene`), decreasing = T)
  ## return info
  return(paste(names(genes), collapse = "|"))
})]

#----------------------------------------#
##--  Protein - Drug - SNP triplets   --##
#----------------------------------------#

## drug target engagement markers: requires drug target to gene mapping as well; SNP -> gene
prot.drug.gene.snp <- merge(edges.drug.protein, edges.drug.gene, by.x = "node1", by.y = "node1", suffixes = c(".protein", ".gene"),
                            allow.cartesian = T)
prot.drug.gene.snp <- merge(prot.drug.gene.snp, edges.snp.gene, by.x = "node2.gene", by.y = "node2", suffixes = c(".gene", ".snp"))
## final mapping back to protein
prot.drug.gene.snp <- merge(prot.drug.gene.snp, edges.snp.protein, 
                            by.x = c("node1.snp", "node2.protein"), 
                            by.y = c("node1", "node2"), 
                            suffixes = c(".snp", ".protein.out"))
## write to file
write.table(prot.drug.gene.snp, "Drug.Target.Engagement.Examples.KG.20250619.txt", sep="\t", row.names=F)

#----------------------------------------#
##--     Drug - SNP - similarities    --##
#----------------------------------------#

## compute drug/SNP overlap for all possible pairs
drug.prot.snp      <- expand.grid(drug=nodes.network[ type == "Drugs"]$node, 
                                  snps=nodes.network[ type == "SNP"]$node) 

# ## Drug --> Protein --> SNP
drug.prot.snp      <- merge(edges.drug.protein, edges.snp.protein,
                            by.x="node2", by.y="node2", suffixes = c(".drug", ".snp"),
                            allow.cartesian = T)
## aggregate by drug - SNP combination
drug.prot.snp      <- drug.prot.snp[, .(common.proteins = paste(node2, collapse = "|"),
                                        number.shared = length(node2)), 
                                    by = c("node1.drug", "node1.snp")]
## add label
drug.prot.snp      <- merge(drug.prot.snp, nodes.network, by.x = "node1.drug", by.y = "node")
drug.prot.snp[, node1.snp := as.character(node1.snp)]
drug.prot.snp      <- merge(drug.prot.snp, nodes.network, by.x = "node1.snp", by.y = "node", suffixes = c(".drug", ".snp"))

## test for enrichment
registerDoMC(10)
drug.prot.snp      <- mclapply(1:nrow(drug.prot.snp), function(x){
  
  ## get drug
  drug <- drug.prot.snp$node1.drug[x]
  ## get SNP
  snp  <- drug.prot.snp$node1.snp[x]
  ## proteins common to both
  prot <- strsplit(drug.prot.snp$common.proteins[x], "\\|")[[1]]
  
  ## --> compute enrichment statistics <-- ##
  
  ## define protein back ground for population
  prot.back <- prot.lab$id
  
  ## protein related to drug and SNP
  d1        <- length(prot)
  ## protein not related to drug but SNP
  d2        <- sum( !(unique(pqtl[ R2.group == snp]$protein_id) %in% prot))
  ## protein realted to drug but not SNP
  d3        <- sum( !(unique(res.var.sex[ variable == drug]$var) %in% prot))
  ## protein not realted to drug and not SNP
  d4        <- length(prot.back[ !(prot.back %in% unique(pqtl[ R2.group == snp]$protein_id) | prot.back %in% unique(res.var.sex[ variable == drug]$var))]) 
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.table(drug.prot.snp[x,], or=enr$estimate, pval=enr$p.value, 
                    d1=d1, d2=d2, d3=d3, d4=d4))
  
}, mc.cores = 10)
## combine again
drug.prot.snp      <- rbindlist(drug.prot.snp)

## add GWAS associations more broadly to each SNP
drug.prot.snp[, snp.gwas := sapply(node1.snp, function(x){
  ## get all possible GWAS consequences
  traits <- unlist(lapply(res.gwas[R2.group == x]$mapped_trait, function(c){
    strsplit(c, "\\|\\||, ")[[1]]
  }))
  
  ## keep only those consistently linked to all members of the r2 group
  traits <- table(traits)
  traits <- names(traits[ traits >= nrow(res.gwas[R2.group == x])])
  
  ## return info
  return(paste(traits, collapse = "|"))
  
  
})]
## add pQTL gene annotation
drug.prot.snp[, pqtl.gene := sapply(node1.snp, function(x){
  ## get all possible GWAS consequences
  genes <- sort(table(pqtl[ R2.group == x]$`Bioinfomatic annotated gene`), decreasing = T)
  ## return info
  return(paste(names(genes), collapse = "|"))
})]
save.image()

## how often is each medication linked to a SNP
sort(table(drug.prot.snp[ number.shared > 1 & pval < .05/nrow(drug.prot.snp)]$label.drug))

## write to file
write.table(drug.prot.snp, "Results.Drug.Protein.SNP.triangles.enrichment.20250619.txt", sep="\t", row.names = F)
write.table(prot.snp.disease.trans, "Results.Disease.Protein.SNP.triangles.enrichment.20250619.txt", sep="\t", row.names = F)
