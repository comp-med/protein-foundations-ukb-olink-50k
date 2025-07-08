######################################################
#### Variance Decomposition UKB proteins          ####
#### Maik Pietzner                     29/01/2025 ####
######################################################

rm(list=ls())
setwd("<path>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(rms)
require(umap)
require(basicPlotteR)
require(colorspace)
require(gprofiler2)
require(readxl)
require(pscl)
require(circlize)
require(igraph)
require(vioplot)

##############################################
####          import label data           ####
##############################################

## import labels for phenotypes
lab.phe  <- fread("../../01_phenotype_prep/data/UKB.labels.variance.decomp.20280129.txt")
## label proteins 
lab.prot <- fread("../../01_phenotype_prep/data/Olink.proteins.variance.decompostion.20250128.txt")  

## create file to run analysis by protein and strata
tmp.sex  <- expand.grid(olink.id=lab.prot$id, population=c("All", "Female", "Male"), stringsAsFactors = F)
tmp.pop  <- expand.grid(olink.id=lab.prot$id, population=c("EUR", "AFR", "CSA"), stringsAsFactors = F) 
## needs to be split in multiple
write.table(tmp.sex, "variance.decomp.sex.proteins", sep="\t", row.names = F, col.names = F, quote = F)
write.table(tmp.pop, "variance.decomp.pop.proteins", sep="\t", row.names = F, col.names = F, quote = F)

##############################################
####            results by sex            ####
##############################################

## get the output
jj <- dir("../output_updated//")
## get only the variance components
jj <- grep("variance", jj, value=T)
## only variables related to All, Female, Male
jj <- grep("All|Female|Male", jj, value=T)

## --> import <-- ##

## import extended mapping to correct for naming difference when only one variable (a factor) was selected to be important
lab.ext     <- fread("../../01_phenotype_prep/data/UKB.labels.variance.decomp.factor.extended.20250129.txt")

## run import in parallel
registerDoMC(10)

## import
res.var.sex.lasso <- mclapply(grep("lasso", jj, value=T), function(x){
  ## import 
  tmp        <- fread(paste0("../output_updated//", x))
  ## add population
  tmp[, population := strsplit(x, "\\.")[[1]][4]]
  ## edit names
  names(tmp) <- sapply(names(tmp), function(x) ifelse(x %in% lab.ext$short_name_new, lab.ext$short_name[ which(lab.ext$short_name_new == x)], x))
  ## account for genetic scores
  names(tmp) <- gsub(paste0("\\.", tmp$var[1]), "", names(tmp))
  return(tmp)
}, mc.cores = 10)
## combine
res.var.sex.lasso <- rbindlist(res.var.sex.lasso, fill = T)
## look for possibly poorly run estimates
summary(res.var.sex.lasso$n.boot)

## replace missing summary measure with p.50!
res.var.sex.lasso[, summary.measure := ifelse(is.na(summary.measure), "p.50", summary.measure)]

## add cis/trans score to the set of labels
lab.phe     <- plyr::rbind.fill(lab.phe, data.table(short_name = c("sample_age", "cis.score", "trans.score"), 
                                                    label = c("Sample age", "cis-pQTL score", "trans-pQTL score"), released=T,
                                                    category = c("Technical", rep("Genetic", 2))))

## look at the distribution of well explained assays
hist(1-res.var.sex.lasso[ summary.measure == "p.50" & population == "All"]$Residuals)
## substantial proportion not well explained!

## transform in more efficient shape
res.var.sex.long <- melt.data.table(res.var.sex.lasso, id.vars = c("var", "summary.measure", "n.boot", "population", "n"))
## reshape again by summary measure
res.var.sex.long <- dcast(res.var.sex.long, var + variable + population + n + n.boot  ~ summary.measure, sep = ".", value.var = "value")
## drop missing values
res.var.sex.long <- res.var.sex.long[ !is.na(p.50) ]
## avoid coding as a factor
res.var.sex.long[, variable := as.character(variable)]

#---------------------------------#
##-- look at specific examples --##
#---------------------------------#

## introduce new colour gradient
cat.col    <- data.table(category=sort(unique(lab.phe$category)), 
                         cl=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#D37295", "#A0CBE8", "#F1CE63", "#D4A6C8", "#8CD17D", "#86BCB6"))

## add to label
lab.phe$cl <- NULL
lab.phe    <- merge(lab.phe, cat.col)
lab.phe    <- as.data.table(lab.phe)

## simple summary

## create simple pie chart
cat.col[, number.feature := sapply(category, function(x) nrow(lab.phe[ category == x]))]

## store results
fwrite(res.var.sex.long, "Results.UKB.Olink.variance.explained.all.men.women.20250612.txt", sep="\t", row.names = F, na = NA)
fwrite(lab.phe, "Labels.UKB.Olink.variance.explained.all.men.women.20250612.txt", sep="\t", row.names = F, na = NA)

## draw barchart
pdf("../graphics/Barchart.UKB.features.20241017.pdf", width = 3.15, height = 3.15)
## create plot with axis break
par(mar=c(1.5,6,.5,.1), lwd=.5, mgp=c(.6,0,0), tck=-.01, bty="n", xaxs="i", yaxs="i", cex.axis=.5, cex.lab=.5)
layout(matrix(1:2, 1, 2), widths = c(.7,.3))

## --> first part <-- ##

## empty plot
plot(c(0, 30), c(-.5, nrow(cat.col)+.5), type="n", xlab="", ylab="",
     yaxt="n", xlim=c(0,30), ylim=rev(c(-.5, nrow(cat.col)+.5)))
## add each category
rect(0, 1:nrow(cat.col)-.4, cat.col$number.feature, 1:nrow(cat.col)+.4, lwd=.5, col=cat.col$cl)

## add label
text(0, 1:nrow(cat.col), pos=2, xpd=NA, labels=stringr::str_to_sentence(cat.col$category), cex=.5)

## --> second part <-- ##
par(mar=c(1.5,.1,.5,.5))

## empty plot
plot(c(35, max(cat.col$number.feature)), c(-.5, nrow(cat.col)+.5), type="n", xlab="Number of features", ylab="",
     yaxt="n", xlim=c(35, max(cat.col$number.feature)), ylim=rev(c(-.5, nrow(cat.col)+.5)))
## add each category
rect(0, 1:nrow(cat.col)-.4, cat.col$number.feature, 1:nrow(cat.col)+.4, lwd=.5, col=cat.col$cl)
## close device
dev.off()

#---------------------------------#
##--        overview           --##
#---------------------------------#

## run in parallel
registerDoMC(10)

## add overall characteristics
tmp     <- mclapply(1:nrow(lab.phe), function(x){
  
  ## get the characteristic of interest
  ii <- lab.phe$short_name[x]
  
  print(ii)
  
  ## do only for variables actually selected
  if(ii %in% res.var.sex.long$variable){
    ## do by population
    res <- lapply(c("All", "Female", "Male"), function(k){
      ## count
      ic <- nrow(res.var.sex.long[ population == k & variable == ii, .(p.50)])
      ## summary
      is <- quantile(res.var.sex.long[ population == k & variable == ii, .(p.50)], probs = c(.25,.5,.75), na.rm = T)
      ## range
      ir <- range(res.var.sex.long[ population == k & variable == ii, .(p.50)], na.rm=T)
      ## mean
      ie <- mean(unlist(res.var.sex.long[ population == k & variable == ii, .(p.50)]), na.rm=T)
      ## top protein
      tp <- res.var.sex.long[population == k & variable == ii]$var[which.max(res.var.sex.long[population == k & variable == ii]$p.50)]
      ## store everything
      return(data.table(lab.phe[x,], i.count=ic, i.25=is[1], i.50=is[2], i.mean=ie, i.75=is[3], 
                        i.min=ir[1], i.max=ir[2], top.protein=tp, population=k))
    })
    ## combine
    res <- rbindlist(res, fill = T)
    
  }else{
    return(lab.phe[x,])
  }
  
}, mc.cores=10)
## combine and add missing values for variables not selected
tmp     <- rbindlist(tmp, fill = T)
## drop missing values
tmp     <- tmp[ !is.na(i.25)]
## spread by population
tmp     <- dcast(tmp, short_name + label + category + id + released + miss.per + type + cl ~ population, sep=".", value.var = c("top.protein", grep("i\\.", names(tmp), value=T)))
## define as new reduced label set
lab.red <- tmp
lab.red <- as.data.table(lab.red)

## overall contribution stratified by technical factors, non-modifiable, and modifiable
lab.red[, type.factor := ifelse(category == "Technical" | short_name == "centre", "technical", 
                                ifelse(category %in% c("Genetic", "Ancestry") | short_name == "sex", "non-modifiable", "modifiable"))]
table(lab.red$type.factor, useNA = "always")

#---------------------------------#
##--           Ordering        --##
#---------------------------------#

## add category sorting to align with colour legend
cat.col[, cat.srt := 1:nrow(cat.col)]
lab.red <- merge(lab.red, cat.col, by=c("category", "cl"))
## sort accordingly
lab.red <- lab.red[order(cat.srt, -i.50.All)]

## --> explained variance <-- ## 

## establish concise ordering of proteins
prot.order       <- merge(lab.prot, res.var.sex.long[ population == "Female" & variable == "Residuals", .(var, p.50)], by.x="id", by.y="var", all.x=T) 
prot.order       <- merge(prot.order, res.var.sex.long[ population == "Male" & variable == "Residuals", .(var, p.50)], by.x="id", by.y="var", all.x=T, suffixes = c(".female", ".male")) 
prot.order       <- merge(prot.order, res.var.sex.long[ population == "All" & variable == "Residuals", .(var, p.50)], by.x="id", by.y="var", all.x=T) 
## replace missing values accordingly
prot.order[, p.50.female := ifelse(is.na(p.50.female), 1, p.50.female)]
prot.order[, p.50.male := ifelse(is.na(p.50.male), 1, p.50.male)]
prot.order[, p.50 := ifelse(is.na(p.50), 1, p.50)]
## assign an ordering
prot.order       <- prot.order[ order(p.50) ]
prot.order[, prt.srt := 1:nrow(prot.order)]
## add to explained variance
res.var.sex.long <- merge(res.var.sex.long, prot.order[, c("id", "prt.srt")], by.x="var", by.y="id")
## add ordering by proteins
tmp              <- res.var.sex.long[ variable != "Residuals"]
tmp              <- tmp[ order(population, var, -p.50)]
## create indicator
tmp[, variable.rank := 1:.N, by=c("population", "var")]
## add to main table
res.var.sex.long <- merge(res.var.sex.long, tmp[, .(population, var, variable, variable.rank)], by=c("var", "variable", "population"), all.x=T)

#---------------------------------#
##--            UMAP           --##
#---------------------------------#

## create tmp data
tmp.umap                   <- as.matrix(res.var.sex.lasso[ population == "All" & summary.measure == "p.50", lab.red$short_name, with=F])
tmp.umap[ is.na(tmp.umap)] <- 0
rownames(tmp.umap)         <-  res.var.sex.lasso[ population == "All" & summary.measure == "p.50"]$var
## drop completely zero columns
jj                         <- apply(tmp.umap, 1, function(x) sum(x != 0))
ii                         <- apply(tmp.umap, 2, function(x) sum(x != 0))
## subset
tmp.umap                   <- tmp.umap[names(jj[ jj > 0]),  names(ii[ ii > 0])]
## compute umap mapping
set.seed(42)
tmp.umap                   <- umap(tmp.umap)$layout
tmp.umap                   <- data.frame(id=rownames(tmp.umap), tmp.umap)
## edit names
names(tmp.umap)            <- c("id", "umap.1", "umap.2")
## add protein label
tmp.umap                   <- merge(tmp.umap, lab.prot)
## add explained variance 
tmp.umap                   <- as.data.table(tmp.umap)
tmp.umap                   <- merge(tmp.umap, res.var.sex.long[ population == "All" & variable == "Residuals", .(var, p.50)], 
                                    by.x = "id", by.y = "var")
## add top explained variable
tmp.umap                   <- merge(tmp.umap, res.var.sex.long[ population == "All" & variable.rank == 1, .(var, variable)], 
                                    by.x = "id", by.y = "var", all.x=T)
## add some information on labels
tmp.umap                   <- merge(tmp.umap, lab.red[, c("short_name", "label", "category", "cl")],
                                    by.x = "variable", by.y = "short_name", all.x=T)
## add remaining variables
tmp.umap                  <- merge(tmp.umap, res.var.sex.lasso[population == "All" & summary.measure == "p.50"],
                                   by.x=c("id", "p.50"), by.y=c("var", "Residuals"))
tmp.umap[is.na(tmp.umap)] <- 0
tmp.umap                  <- as.data.frame(tmp.umap)

#####################################################################
####                     cluster protein targets                 ####
#####################################################################

## add umap projection to check clustering
prot.order                    <- merge(prot.order, as.data.table(tmp.umap[, c("id", "umap.1", "umap.2")]), all.x=T)

#------------------------#
##--  create cluster  --##
#------------------------#

## --> try k-means clustering <-- ##

## robust sparse k-means: 
kmeans_object           <- as.data.frame(tmp.umap[, lab.red$short_name])
rownames(kmeans_object) <- tmp.umap$id
colnames(kmeans_object) <- lab.red$short_name
## drop completely zero columns
jj                      <- apply(kmeans_object, 1, function(x) sum(x != 0))
ii                      <- apply(kmeans_object, 2, function(x) sum(x != 0))
## subset
kmeans_object           <- kmeans_object[names(jj[ jj > 0]),  names(ii[ ii > 0])]

## drop rare explanatory variables
ii                      <- apply(kmeans_object, 2, function(x) sum(x != 0))
## subset
kmeans_object           <- kmeans_object[,  names(ii[ ii >= 10])]


## --> run through different k and plot original umap <-- ##

## do assessment across multiple rounds using random seeds
set.seed(42)
# determine overall SS
wss <- (nrow(kmeans_object)-1)*sum(apply(kmeans_object,2,var))

for (i in 2:25){
  ## print progress
  print(i)
  ## allow for multiple random starting points to achieve stable solution
  fit.KM                                   <- kmeans(kmeans_object, centers=i, nstart = 100)
  ## get within cluster sums of squares
  wss[i]                                   <- sum(fit.KM$withinss)
  ## store assignment
  kmeans_object[, paste("km", i, sep=".")] <- fit.KM$cluster
} 

## add UMAP candidates
kmeans_object$id        <- row.names(kmeans_object)
kmeans_object           <- merge(kmeans_object, prot.order[, .(id, umap.1, umap.2)], all.x=T)

## plot for using original umap
pdf("../graphics/Kmeans.variance.explained.simple.20250603.pdf", width = 6.3, height = 6.3)
par(mar=c(1.5,1.5,.5,.5), mgp=c(.4,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="r", yaxs="r", mfrow=c(3,3))

## loop through
for(j in 2:25){
  ## do the actual plot
  plot(umap.2 ~ umap.1, kmeans_object, pch=20, col=kmeans_object[, paste("km", j, sep=".")]+1, 
       xlab="UMAP 1", ylab="UMAP 2", xaxt="n", yaxt="n",
       cex=.4)
  ## add axis
  axis(1, lwd=.5); axis(2, lwd=.5)
  ## add legend
  legend("topleft", lty=0, pch=20, legend=1:j, col = (1:j)+1, cex=.5, bty="n")
}

dev.off()

## --> go for 8 as the cut-off (mimimizes WSS) + visual inspection <-- ##

## add to protein order
kmeans_object <- as.data.table(kmeans_object)
prot.order    <- merge(prot.order, kmeans_object[, .(id, km.8)], all.x=T)

## --> create cluster profiles <-- ##

## create variables to ease coding
for(j in 1:8){
  kmeans_object[, paste0("cluster.", j)] <- ifelse(kmeans_object$km.8 == j, 1, 0)
}

## compute cluster marker based on simple testing
registerDoMC(10)
cluster.marker <- mclapply(1:max(prot.order$km.8, na.rm = T), function(x){
  
  ## run zero inflation model on pseudo counts for each variable
  res <- lapply(lab.red[ i.count.All >= 10]$short_name, function(k){
    
    ## create formula
    zero.test <- paste0("round(", k, "*1e6) ~ cluster.", x, " | cluster.", x)
    ## run model
    zero.test <- summary(zeroinfl(as.formula(zero.test), dist = "negbin", data = kmeans_object))
    ## return results
    return(data.table(short_name=k, km.8 = x, 
                      beta.negbin = zero.test$coefficients$count[2,1], se.negbin = zero.test$coefficients$count[2,2], 
                      zval.negbin = zero.test$coefficients$count[2,3], pval.negbin = zero.test$coefficients$count[2,4],
                      beta.zero = zero.test$coefficients$zero[2,1], se.zero = zero.test$coefficients$zero[2,2], 
                      zval.zero = zero.test$coefficients$zero[2,3], pval.zero = zero.test$coefficients$zero[2,4]))
    
  })
  ## combine
  res <- rbindlist(res, fill = T)
  return(res)
  
}, mc.cores=10)
## combine results
cluster.marker <- rbindlist(cluster.marker)
## add label
cluster.marker <- merge(cluster.marker , lab.red[, .(short_name, label, category, cl)])

## get marker feature
View(cluster.marker[ (pval.negbin < .05/nrow(cluster.marker) | pval.zero < .05/nrow(cluster.marker)) & km.8 == 8]) 
View(cluster.enrich[ or > 1 & pval < .05/nrow(cluster.enrich) & km.8 == 8])

#----------------------------------------#
##--         pathway enrichment       --##
#----------------------------------------#

## import gene assignment from the Nature flagship paper
olink.lab           <- as.data.table(read_excel("41586_2023_6592_MOESM3_ESM (1).xlsx", 4, skip = 2))
## rename
names(olink.lab)    <- c("ukb.ppp.id", "olink.id", "assay", "panel", "hgnc_symbol", "uniprot", "chr", "gene_start", "gene_end", "dilution", "below.LOD", "CV", "fail.per", "icc")
## create identifier to map to selected proteins
olink.lab[, prot.id := gsub("-", "_", tolower(assay))]

## create background list of genes
prot.back           <- na.omit(unique(unlist(lapply(olink.lab$hgnc_symbol, function(x){
  x <- strsplit(x, ";")[[1]]
  return(x)
}))))
## n = 2938 protein coding genes


## perform pathway enrichment for each factor
enrich.pathway.cluster <- lapply(1:8, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- prot.order[ km.8 == x ]$id
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## only if any
  if(length(prot.diff) > 0){
    
    ## run g:Profiler
    enrich    <- gost(query = prot.diff,
                      organism = "hsapiens", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "fdr", 
                      domain_scope = "annotated", custom_bg = prot.back, 
                      numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
    ## keep only what is needed
    enrich    <- as.data.table(enrich$result)
    
    ## if anything was found
    if(nrow(enrich) > 0){
      ## add fold change
      enrich[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]
      ## return results
      return(data.table(km.8 = x, enrich))
    }
  }
})
## combine
enrich.pathway.cluster <- rbindlist(enrich.pathway.cluster, fill=T)
## nothing screaming, but some matching results

#----------------------------------------#
##--          tissue enrichment       --##
#----------------------------------------#

## import HPA assignment from Burulça
hpa.mapping        <- unique(fread("<path>"))

## understand possible tissue patterns
table(hpa.mapping$`RNA tissue specific nTPM`)
## get all possibly named tissues
tmp                <- lapply(1:nrow(hpa.mapping), function(x){
  print(x)
  ## get all possible tissues
  jj <- strsplit(hpa.mapping$`RNA tissue specific nTPM`[x], ";")[[1]]
  ## proceed only if any evidence
  if(length(jj) > 0){
    ## get the relevant data
    jj <- do.call(rbind, lapply(jj, function(k) strsplit(k, ": ")[[1]]))
    ## convert to data set
    jj <- data.table(tissue = jj[, 1], ntpm = as.numeric(jj[, 2]), ind = 1)
    jj <- reshape(jj, timevar = "tissue", idvar = "ind", direction = "wide", sep=".")
    ## return
    return(data.table(Gene=hpa.mapping$Gene[x], jj))
  }
})
## combine
tmp                <- rbindlist(tmp, fill = T)
## edit names
names(tmp)         <- gsub("ntpm\\.", "", names(tmp))
names(tmp)         <- gsub(" ", "_", names(tmp))
tmp$ind            <- NULL
## get tissues
tissues            <- names(tmp)[-1]
## add some info
tmp                <- merge(tmp, hpa.mapping[, c("Gene", "Ensembl", "RNA tissue specificity")], all=T)
hpa.tissue         <- tmp

## do in parallel
registerDoMC(10)

## test whether proteins in a cluster are specific to a given tissue
enrich.tissue.cluster <- mclapply(1:8, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- prot.order[ km.8 == x ]$id
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## test for the enrichment across different tissues
  enr       <- lapply(tissues, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.tissue[ Gene %in% prot.diff & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(tissue=k, km.8=x, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  enr       <- rbindlist(enr, fill=T)
  return(enr)
}, mc.cores=10)
## combine results
enrich.tissue.cluster <- rbindlist(enrich.tissue.cluster, fill = T)

#----------------------------------------#
##--       cell-type enrichment       --##
#----------------------------------------#

## understand possible tissue patterns
table(hpa.mapping$`RNA single cell type specificity`)
## get all possibly named tissues
tmp                   <- lapply(1:nrow(hpa.mapping), function(x){
  print(x)
  ## get all possible tissues
  jj <- strsplit(hpa.mapping$`RNA single cell type specific nTPM`[x], ";")[[1]]
  ## proceed only if any evidence
  if(length(jj) > 0){
    ## get the relevant data
    jj <- do.call(rbind, lapply(jj, function(k) strsplit(k, ": ")[[1]]))
    ## convert to data set
    jj <- data.table(tissue = jj[, 1], ntpm = as.numeric(jj[, 2]), ind = 1)
    jj <- reshape(jj, timevar = "tissue", idvar = "ind", direction = "wide", sep=".")
    ## return
    return(data.table(Gene=hpa.mapping$Gene[x], jj))
  }
})
## combine
tmp                   <- rbindlist(tmp, fill = T)
## edit names
names(tmp)            <- gsub("ntpm\\.", "", names(tmp))
names(tmp)            <- gsub(" ", "_", names(tmp))
tmp$ind               <- NULL
## get tissues
cell.type             <- names(tmp)[-1]
## add some info
tmp                   <- merge(tmp, hpa.mapping[, c("Gene", "Ensembl", "RNA single cell type specificity")], all=T)
hpa.cell              <- tmp

## do in parallel
registerDoMC(10)

## test whether proteins in clusters are enriched among certain cell types
enrich.cell.cluster <- mclapply(1:8, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- prot.order[ km.8 == x ]$id
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## test for the enrichment across different tissues
  enr       <- lapply(cell.type, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.cell[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.cell[ !(Gene %in% prot.diff) & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.cell[ Gene %in% prot.diff & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.cell[ !(Gene %in% prot.diff) & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(cell.type=k, km.8=x, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.cell[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  enr       <- rbindlist(enr, fill=T)
  return(enr)
}, mc.cores=10)
## combine results
enrich.cell.cluster <- rbindlist(enrich.cell.cluster, fill = T)

#----------------------------------#
##--     cluster assignment     --##
#----------------------------------#

## import cluster assignment
cluster.assignment <- fread("Kmeans_cluster_label_20250612.csv")
## order and create colour vector
cluster.assignment <- cluster.assignment[ order(order)]
## add colour code
cluster.assignment[, cl.cluster := c("#A05486", "#606c38", "#9AA677", "#D8E4B6", "#007AA3", "#de0a26",  "#f4a261", "#1477d2")]
## add to protein data
prot.order         <- merge(prot.order, cluster.assignment, by.x="km.8", by.y="cluster", all.x=T)

## add number of proteins by cluster
cluster.assignment <- merge(cluster.assignment, prot.order[, .(count=.N), by="km.8"], by.x="cluster", by.y="km.8")

## write to file
write.table(cluster.marker, "Protein.Cluster.Marker.Phenotypes.20250612.txt", sep="\t", row.names=F)

#####################################################################
####         pathway/tissue enrichment analysis by factor        ####
#####################################################################

#----------------------------------------#
##--         pathway enrichment       --##
#----------------------------------------#

## perform pathway enrichment for each factor
res.enr.pathway.all <- lapply(lab.red$short_name, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- res.var.sex.long[population == "All" & variable == x ]$var
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## only if any
  if(length(prot.diff) > 0){
    
    ## run g:Profiler
    enrich    <- gost(query = prot.diff,
                      organism = "hsapiens", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "fdr", 
                      domain_scope = "annotated", custom_bg = prot.back, 
                      numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
    ## keep only what is needed
    enrich    <- as.data.table(enrich$result)
    
    ## if anything was found
    if(nrow(enrich) > 0){
      ## add fold change
      enrich[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]
      ## return results
      return(data.table(short_name = x, enrich))
    }
  }
})
## combine
res.enr.pathway.all <- rbindlist(res.enr.pathway.all, fill=T)
## add label
res.enr.pathway.all <- merge(res.enr.pathway.all, lab.red[, .(short_name, label)])

#----------------------------------------#
##--          tissue enrichment       --##
#----------------------------------------#

## do in parallel
registerDoMC(10)

## test whether proteins explained by certain factors are enriched for 
## expression in certain tissues
res.enr.tissue.all <- mclapply(lab.red$short_name, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- res.var.sex.long[population == "All" & variable == x ]$var
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## test for the enrichment across different tissues
  enr       <- lapply(tissues, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.tissue[ Gene %in% prot.diff & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(tissue=k, short_name=x, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  enr       <- rbindlist(enr, fill=T)
  return(enr)
}, mc.cores=10)
## combine results
res.enr.tissue.all <- rbindlist(res.enr.tissue.all, fill = T)
## add some label
res.enr.tissue.all <- merge(res.enr.tissue.all, lab.red[, .(short_name, label, category, cl, cat.srt, type.factor, i.count.All)])

## prune for terms selected at least 5 time
res.enr.tissue.all <- res.enr.tissue.all[ i.count.All > 4]

#----------------------------------------#
##--       cell-type enrichment       --##
#----------------------------------------#

## do in parallel
registerDoMC(10)

## test whether proteins explained by certain factors are enriched for 
## expression in certain tissues
res.enr.cell.type.all <- mclapply(lab.red$short_name, function(x){
  
  print(x)
  
  ## get the proteins of interest
  prot.diff <- res.var.sex.long[population == "All" & variable == x ]$var
  prot.diff <- unique(unlist(lapply(olink.lab[ prot.id %in% prot.diff]$hgnc_symbol, function(x) strsplit(x, ";")[[1]])))
  
  ## test for the enrichment across different tissues
  enr       <- lapply(cell.type, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.cell[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.cell[ !(Gene %in% prot.diff) & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.cell[ Gene %in% prot.diff & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.cell[ !(Gene %in% prot.diff) & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(cell.type=k, short_name=x, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.cell[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  enr       <- rbindlist(enr, fill=T)
  return(enr)
}, mc.cores=10)
## combine results
res.enr.cell.type.all <- rbindlist(res.enr.cell.type.all, fill = T)
## add some label
res.enr.cell.type.all <- merge(res.enr.cell.type.all, lab.red[, .(short_name, label, category, cl, cat.srt, type.factor, i.count.All)])
## prune for terms selected at least 5 time
res.enr.cell.type.all <- res.enr.cell.type.all[ i.count.All > 4]

## write enrichment results to file
write.table(res.enr.cell.type.all, "Results.cell.type.enrichment.explained.variance.20250612.txt", sep="\t", row.names=F)
write.table(res.enr.tissue.all, "Results.tissue.enrichment.explained.variance.20250612.txt", sep="\t", row.names=F)

#----------------------------------------#
##--   map to entire Olink platform   --##
#----------------------------------------#

## convolut table
tmp.tissue <- melt(hpa.tissue, id.vars = c("Gene", "Ensembl", "RNA tissue specificity"), na.rm = T)
tmp.tissue <- tmp.tissue[, .(tissues = paste(variable, collapse = "|")) , by="Gene"]
## add to protein labels
lab.prot   <- merge(lab.prot, tmp.tissue, by.x="Assay", by.y="Gene", all.x=T)

## same for cell types
tmp.cells  <- melt(hpa.cell, id.vars = c("Gene", "Ensembl", "RNA single cell type specificity"), na.rm = T)
tmp.cells  <- tmp.cells[, .(cells = paste(variable, collapse = "|")) , by="Gene"]
## add to protein labels
lab.prot   <- merge(lab.prot, tmp.cells, by.x="Assay", by.y="Gene", all.x=T)

## need tissue assignment for all genes in HPA
hpa.all    <- fread("<path>")
## process some columns
hpa.all[, `RNA single cell type specific nTPM` := gsub("[{}\"]", "", `RNA single cell type specific nTPM`)]
hpa.all[, `RNA tissue specific nTPM` := gsub("[{}\"]", "", `RNA tissue specific nTPM`)]

## --> cell-type data set <-- ##

## understand possible cell patterns
table(hpa.all$`RNA single cell type specificity`)
## get all possibly named tissues
tmp                   <- lapply(1:nrow(hpa.all), function(x){
  print(x)
  ## get all possible tissues
  if(hpa.all$`RNA single cell type specific nTPM`[x] != "null"){
    jj <- strsplit(hpa.all$`RNA single cell type specific nTPM`[x], ", ")[[1]]
    ## proceed only if any evidence
    if(length(jj) > 0){
      ## get the relevant data
      jj <- do.call(rbind, lapply(jj, function(k) strsplit(k, ": ")[[1]]))
      ## convert to data set
      jj <- data.table(tissue = jj[, 1], ntpm = as.numeric(jj[, 2]), ind = 1)
      jj <- reshape(jj, timevar = "tissue", idvar = "ind", direction = "wide", sep=".")
      ## return
      return(data.table(Gene=hpa.all$Gene[x], jj))
    }
  }
})
## combine
tmp                   <- rbindlist(tmp, fill = T)
## edit names
names(tmp)            <- gsub("ntpm\\.", "", names(tmp))
names(tmp)            <- gsub(" ", "_", names(tmp))
tmp$ind               <- NULL
## add some info
tmp                   <- merge(tmp, hpa.all[, c("Gene", "Ensembl", "RNA single cell type specificity")], all=T)
hpa.cell.all          <- tmp

## 'specific' expression
hpa.cell.all[, cell.specific := ifelse(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched"), 1, 0)]

## add column whether included in Olink
hpa.cell.all[, olink.included := ifelse(Gene %in% hpa.cell$Gene, 1, 0)]

## get all possibly named tissues
tmp                   <- lapply(1:nrow(hpa.all), function(x){
  print(x)
  ## get all possible tissues
  if(hpa.all$`RNA tissue specific nTPM`[x] != "null"){
    jj <- strsplit(hpa.all$`RNA tissue specific nTPM`[x], ", ")[[1]]
    ## proceed only if any evidence
    if(length(jj) > 0){
      ## get the relevant data
      jj <- do.call(rbind, lapply(jj, function(k) strsplit(k, ": ")[[1]]))
      ## convert to data set
      jj <- data.table(tissue = jj[, 1], ntpm = as.numeric(jj[, 2]), ind = 1)
      jj <- reshape(jj, timevar = "tissue", idvar = "ind", direction = "wide", sep=".")
      ## return
      return(data.table(Gene=hpa.all$Gene[x], jj))
    }
  }
})
## combine
tmp                   <- rbindlist(tmp, fill = T)
## edit names
names(tmp)            <- gsub("ntpm\\.", "", names(tmp))
names(tmp)            <- gsub(" ", "_", names(tmp))
tmp$ind               <- NULL
## add some info
tmp                   <- merge(tmp, hpa.all[, c("Gene", "Ensembl", "RNA tissue specificity")], all=T)
hpa.tissue.all        <- tmp

## 'specific' expression
hpa.tissue.all[, tissue.specific := ifelse(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched"), 1, 0)]
table(hpa.tissue.all$tissue.specific)

## add whether captured by Olink
hpa.tissue.all[, olink.included := ifelse(Gene %in% hpa.tissue$Gene, 1, 0)]
table(hpa.tissue.all$olink.included)

#####################################################################
####            some numbers/figures for the manuscript          ####
#####################################################################

#---------------------------------------#
##--          general numbers        --##
#---------------------------------------#

## how many protein targets correlated with at least one in *All* UKB
nrow(prot.order[ p.50 != 1])
## average contribution
mean(1-prot.order[ p.50 != 1]$p.50)
range(1-prot.order[ p.50 != 1]$p.50)*100
## how many factors
length(unique(res.var.sex.long[population == "All" & variable != "Residuals"]$variable))

## how many factors per protein target
prot.order <- merge(prot.order, dcast(res.var.sex.long[ variable != "Residuals", .(count = length(variable)), by = c("population", "var")], var ~ population, value.var = "count"), 
                    all.x = T, by.x="id", by.y="var")

## summary of factors contributing
summary(prot.order$All)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1.00   12.00   20.00   18.75   26.00   37.00      66 

# ## contribution of frequently selected factors
# View(lab.red[ i.count.All > (nrow(prot.order)/3)])
# 
# ## 'specific' factors
# View(lab.red[ i.count.All < 80 & i.50.All > .01])

## for how many of the protein targets >1% explained variance, did only one factor contribute most
prot.order <- merge(prot.order, dcast(res.var.sex.long[ variable != "Residuals", .(max = max(p.50)), by = c("population", "var")], var ~ population, value.var = "max"), 
                    all.x = T, by.x="id", by.y="var", suffixes = c(".count", ".explained.top"))
# ## count
# nrow(prot.order[ (1-p.50) > .05]) ## n = 2320
# nrow(prot.order[ (1-p.50) > .05 & (2*All.explained.top >= (1-p.50))]) ## n = 742

## add strongest contributing factor
tmp        <- res.var.sex.long[ variable != "Residuals" ]
tmp        <- tmp[order(population, var, -p.50)]
tmp[, ind := 1:.N, by=c("population", "var")]
tmp        <- tmp[ ind == 1]
## add
prot.order <- merge(prot.order, dcast(tmp, var ~ population, value.var = "variable"), all.x = T, by.x="id", by.y="var")
## look into it
sort(table(prot.order[ (1-p.50) > .05 & (2*All.explained.top >= (1-p.50))]$All))
## compute sum: 1 + 1 + 1 + 1 + 1...

## how many proteins are most explained by general axis
nrow(prot.order[ 2*All.explained.top >= (1-p.50) & All %in% lab.red[ i.count.All > (nrow(prot.order)/3)]$short_name])
## n = 808

## --> compute contribution by each category <-- ##

## add type of variable
res.var.sex.long <- merge(res.var.sex.long, lab.red[, .(short_name, type.factor)], by.x = "variable", by.y = "short_name", all.x=T)
prot.order       <- merge(prot.order, dcast(res.var.sex.long[ population == "All" & variable != "Residuals", .(p50.category = sum(p.50)), by=c("var", "type.factor")], var ~ type.factor), 
                          by.x = "id", by.y = "var", all.x=T)

## compare overall results
wilcox.test(prot.order$modifiable, prot.order$`non-modifiable`)$p.value
## 7.495093e-47

## how often does modifiable outweight non-modifiable
sum(ifelse(is.na(prot.order$modifiable), 0, prot.order$modifiable) > 2*ifelse(is.na(prot.order$`non-modifiable`), 0, prot.order$`non-modifiable`))
## n = 1632

## general summary
summary(prot.order[, .(modifiable, `non-modifiable`, technical)])

## how many protein targets best explained by technical variation
nrow(prot.order[ (1-p.50)/2 < technical])
## n = 15

## contribution cluster
summary(res.var.sex.long[ variable == "centre" & var %in% prot.order[ km.8 == 2]$id & population == "All"]$p.50)
summary(res.var.sex.long[ variable == "Platelet_crit" & var %in% prot.order[ km.8 == 2]$id & population == "All"]$p.50)

#---------------------------------------#
##--         Secretion status        --##
#---------------------------------------#

## add status on whether or not protein is secreted
prot.order <- merge(prot.order, hpa.mapping[, .(Gene, `Biological process`, `Molecular function`, `RNA tissue specificity`, `RNA single cell type specificity`, `Subcellular location`,
                                                `Subcellular main location`, `Secretome location`)], by.x = "Assay", by.y = "Gene", all.x=T)
## create color coding for secretome
# sec.col    <- data.table(`Secretome location` = rev(names(sort(table(prot.order$`Secretome location`)))), sec.cl = c("grey70", "#E63946", "#1D3557", "#F4A261", "#2A9D8F", "#9B5DE5", "#F15BB5", "#00B4D8", "#FFB703", "#06D6A0", "#118AB2"))
sec.col    <- data.table(`Secretome location` = rev(names(sort(table(prot.order$`Secretome location`)))), sec.cl = RColorBrewer::brewer.pal(11, "Paired"))
sec.col    <- data.table(`Secretome location` = rev(names(sort(table(prot.order$`Secretome location`)))), sec.cl = c("#A6CEE3", "#E31A1C", "#FDBF6F", "#B2DF8A", "#33A02C", "#1F78B4",  "#FB9A99", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"))
## add to prot.order
prot.order <- merge(prot.order, sec.col, by = "Secretome location", all.x=T)
## add counts
sec.col[, number.proteins := sapply(`Secretome location`, function(x) nrow(prot.order[ `Secretome location` == x]))]
sec.col$`Secretome location`[1] <- "Not secreted"

## create pie chart
pdf("../graphics/Olink.proteins.secretome.20250604.pdf", width = 3.15, height = 3.15)
par(mar=rep(.5,4), lwd=.5)
pie(sec.col$number.proteins, col = lighten(sec.col$sec.cl, .1), labels = sec.col$`Secretome location`, init.angle = 317, clockwise = T)
dev.off()

#---------------------------------------#
##--      some results to files      --##
#---------------------------------------#

## write explained variance to file
fwrite(res.var.sex.long, "Results.variance.decomposition.UKB.Olink.20250612.txt", sep="\t", row.names = F, na=NA)
## protein summary
fwrite(prot.order, "Protein.variance.summary.UKB.Olink.20250612.txt", sep="\t", row.names = F, na=NA)

## introduce new colour gradient
colors     <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#D37295", "#A0CBE8", "#F1CE63", "#D4A6C8", "#8CD17D", "#86BCB6")
cat.col$cl <- colors
## change accordingly
lab.red[, cl := NULL]
lab.red    <- merge(lab.red, cat.col)
## change accordingly
lab.phe[, cl := NULL]
lab.phe    <- merge(lab.phe, cat.col)

## label set
fwrite(lab.phe, "UKB.labels.variance.decomp.20250612.txt", sep="\t", row.names = F, na=NA)

#---------------------------------------#
##--            Figure 1             --##
#---------------------------------------#

## store as png
pdf("../graphics/Figure_1.pdf", width = 6.3, height = 6.3)
## graphical parameters
par(mar=c(1.5,1.5,.5,.5), mgp=c(.6,0,0), tck=-.01, cex.axis=.5, cex.lab=.5, lwd=.5, xaxs="i", yaxs="i", bty="n")
layout(matrix(c(1,3,3,2,3,3,4,4,5), 3, 3, byrow = T), heights = c(1/6,1/2,1/3))

## -->  Secretome status  <-- ##

par(mar=rep(.5,4), lwd=.5)
pie(sec.col$number.proteins, col = lighten(sec.col$sec.cl, .1), labels = sec.col$`Secretome location`, init.angle = 317, clockwise = T, cex=.5,
    radius = .9)

## --> Number of features <-- ##

par(mar=c(1.5,1.5,.5,.5))

## empty plot
plot(c(0, sqrt(max(cat.col$number.feature))), c(-.5, nrow(cat.col)+.5), type="n", xlab="Number of variables", ylab="",
     yaxt="n",  ylim=rev(c(-.5, nrow(cat.col)+.5)), xaxt="n")
## add axis
axis(1, lwd=.5, at=sqrt(c(1,5,10,50,100,500,1000)), labels = c(1,5,10,50,100,500,1000))
## order by the number of features
cat.col <- cat.col[order(-number.feature)]
## add each category
rect(0, 1:nrow(cat.col)-.4, sqrt(cat.col$number.feature), 1:nrow(cat.col)+.4, lwd=.5, col=cat.col$cl)
## add label
text(0, 1:nrow(cat.col), pos=4, xpd=NA, labels=stringr::str_to_sentence(cat.col$category), cex=.5)

## -->       UMAP       <-- ##

## adopt parameters
par(mar=c(2.5,2.5,.5,.5), mgp=c(.6,0,0), tck=-.01, cex.axis=.6, cex.lab=.6, lwd=.5, xaxs="r", yaxs="r", bty="o")

## do the actual plot
plot(umap.2 ~ umap.1, prot.order, pch=20, col=prot.order$cl.cluster, 
     xlab="UMAP 1", ylab="UMAP 2", xaxt="n", yaxt="n",
     cex=.5)
## add axis
axis(1, lwd=.5); axis(2, lwd=.5)

## get plotting coordinates
pm <- par("usr")

## add legend
cluster.assignment <- cluster.assignment[order(order)]
legend(pm[1]+(pm[2]-pm[1])*.7, pm[3]+(pm[4]-pm[3])*.7, lty=0, pch=22, pt.bg=cluster.assignment$cl.cluster, pt.lwd=.3, pt.cex=.9,
       cex=.6, legend=cluster.assignment$label, lwd=.1, bty="n")

## --> Overall variance <-- ##

## adjust plotting parameters
par(mar=c(1,1.5,.5,.5), mgp=c(.6,0,0), bty="o", xaxs="i", yaxs="i")

## empty plot
plot(c(.5,nrow(prot.order)+.5), c(0,.85), type="n", xlab="Ordered protein targets", ylab="Explained variance [%]", xaxt="n", yaxt="n")
## add y-axis
axis(2, lwd=.5, at=c(0,.2,.4,.6,.8), labels = c(0,20,40,60,80))
## add type
pm <- par("usr")
## add label
mtext("Ordered protein targets", 1, cex=.5)

## draw each protein
for(j in prot.order$prt.srt){
  ## get the relevant data
  tmp <- res.var.sex.long[ var ==  prot.order$id[which(prot.order$prt.srt == j)] & population == "All" & variable != "Residuals"]
  if(nrow(tmp) > 0){
    ## reorder to fit label set
    tmp <- tmp[order(match(tmp$variable, lab.red$short_name))]
    ## draw rectanlges
    rect(j-.5, c(0, cumsum(tmp$p.50)[-nrow(tmp)]), j+.5, cumsum(tmp$p.50), lwd=.1, border = NA,
         col=sapply(tmp$variable, function(x) lab.red$cl[which(lab.red$short_name == x)]))
  }
}

## add legend
cat.col <- cat.col[order(cat.col$cat.srt)]
legend("topright", pt.cex=.8, pch=22, pt.lwd = .3, lty=0, bty="n",
       legend = cat.col$category, pt.bg = cat.col$cl, cex=.6, ncol=2)


## --> Violin plots <-- ##

## adopt parameters
par(mar=c(1.5,1.5,.5,.5), xaxs="r", yaxs="r")


## violin plot (make sure to use the very same bandwith to estimate density)
vioplot(prot.order[, .(modifiable, `non-modifiable`, technical)]*100, ylab="Explained variance [%]", xlab="",
        names=c("modifiable", "non-modifiable", "technical"), col=adjust_transparency("grey60", .5), lwd=.4, h=1.755,
        side="both")

## close device
dev.off()

#---------------------------------------#
##--               SF2               --##
#---------------------------------------#

pdf("../graphics/UMAP.explained.variance.selected.cluster.20250604.pdf", width = 6.3, height = 6.3)
par(mar=c(1.5,1.5,.5,.5), mgp=c(.4,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="r", yaxs="r")
layout(matrix(c(1,1,2,3,1,1,2,4:12), 4, 4, byrow = T))

#----------------------------#
##--     overall plot     --##
#----------------------------#

## do the actual plot
plot(umap.2 ~ umap.1, prot.order, pch=20, col=prot.order$cl.cluster, 
     xlab="UMAP 1", ylab="UMAP 2", xaxt="n", yaxt="n",
     cex=.4)
## add axis
axis(1, lwd=.5); axis(2, lwd=.5)
## get plotting coordinates
pm <- par("usr")

## add legend
cluster.assignment <- cluster.assignment[order(order)]
legend(pm[1]+(pm[2]-pm[1])*.7, pm[3]+(pm[4]-pm[3])*.7, lty=0, pch=22, pt.bg=cluster.assignment$cl.cluster, pt.lwd=.3, pt.cex=.9,
       cex=.6, legend=cluster.assignment$label, lwd=.1, bty="n")



#----------------------------#
##-- Proteins by cluster  --##
#----------------------------#

## adapt parameters
par(xaxs="i", yaxs="i")
## empty plot
plot(c(0,930), c(.5,nrow(cluster.assignment)+.5), type="n", xlab="Number of proteins", ylab="", xaxt="n", yaxt="n",
     ylim=rev(c(.5,nrow(cluster.assignment)+.5)))
## add axis
axis(1, lwd=.5)
## add rectangle
rect(0, cluster.assignment$order-.4, cluster.assignment$count, cluster.assignment$order+.4, border="grey30", lwd=.2, 
     col=adjust_transparency(cluster.assignment$cl.cluster, .5))
## add label
text(0, cluster.assignment$order, cex=.5, labels=cluster.assignment$label, col="grey10", pos=4, offset=.1)


#----------------------------#
##--     single plots     --##
#----------------------------#


## adapt parameters
par(xaxs="r", yaxs="r")

## go thorugh all variables
for(j in c("cis.score", "trans.score", "sex", "pop", "log_cysc_cleaned", "bin_585.3", "log_alt_cleaned", "log_crp_cleaned","Platelet_crit", "centre")){
  
  ## define the maximum to colour for
  ii <- round(max(tmp.umap[,j], na.rm=T)*10000+1)
  ii <- colorRampPalette(c("grey90", "#00A4CC"))(ii)
  
  ## do the actual plot
  plot(umap.2 ~ umap.1, tmp.umap, pch=20, 
       col=ii[(tmp.umap[,j])*10000+1], 
       xlab="UMAP 1", ylab="UMAP 2", xaxt="n", yaxt="n",
       cex=.4, lwd=.1)
  ## overlay again
  jj <- which(tmp.umap[,j] > 0)
  points(tmp.umap[jj, "umap.1"], tmp.umap[jj, "umap.2"], cex=.4, pch=20, lwd=.1, col=ii[(tmp.umap[jj,j])*10000+1])
  ## overlay again
  jj <- which(tmp.umap[,j] > .01)
  points(tmp.umap[jj, "umap.1"], tmp.umap[jj, "umap.2"], cex=.4, pch=21, lwd=.1, bg=ii[(tmp.umap[jj,j])*10000+1], col="grey30")
  ## add axis
  axis(1, lwd=.5); axis(2, lwd=.5)
  ## add title
  mtext(lab.red$label[which(lab.red$short_name == j)], cex=.3)
  
  ## add legend
  legend("bottomleft", bty="n", cex=.4, lty=0, pch=21, pt.bg="grey80", pt.lwd=.1, pt.cex=.7,
         legend = "Explained variance > 1%")
  
  ## colour gradient
  pm <- par("usr")
  ## length
  l <- seq(pm[1]+(pm[2]-pm[1])*.05, pm[1]+(pm[2]-pm[1])*.35, length.out = length(ii))
  ## rectangle for the colours
  rect(l-(l[2]-l[1])/2, pm[4]-(pm[4]-pm[3])*.05, l+(l[2]-l[1])/2, pm[4]-(pm[4]-pm[3])*.1, border=NA, col=ii)
  ## box
  rect(l[1]-(l[2]-l[1])/2, pm[4]-(pm[4]-pm[3])*.05, l[length(l)]+(l[2]-l[1])/2, pm[4]-(pm[4]-pm[3])*.1, border="black", col=NA, lwd=.3)
  ## add header
  text(pm[1]+(pm[2]-pm[1])*.1, pm[4]-(pm[4]-pm[3])*.03, cex=.4, labels = "Explained variance [%]", pos=4,
       offset = .2)
  ## simple axis
  text(l[round(c(1, c(.2, .4, .6, .8, 1)*length(ii)))], pm[4]-(pm[4]-pm[3])*.11, 
       labels=sprintf("%.2f", c(0, .2, .4, .6, .8, 1)*length(ii)/100), pos=1, cex=.3, offset = .1)
  
}

dev.off()

#---------------------------------------#
##--         Figure 2 - cells        --##
#---------------------------------------#

## import HPA annotation files
hpa.cell.anno         <- fread("HPA.cell.type.annotation.20241017.txt")
hpa.tissue.anno       <- fread("HPA.tissue.annotation.20241017.txt")

## add colour and ordering
res.enr.cell.type.all <- merge(res.enr.cell.type.all, hpa.cell.anno[, .(cell.type, category, colour, label, order)], 
                               by="cell.type", suffixes = c(".feature", ".hpa"), all.x = T)
res.enr.tissue.all    <- merge(res.enr.tissue.all, hpa.tissue.anno[, .(tissue, group, colour, label, order)], 
                               by="tissue", suffixes = c(".feature", ".hpa"), all.x = T)

## loop through each cell type
for(j in unique(res.enr.cell.type.all[ pval < .05/nrow(res.enr.cell.type.all) & or > 1]$label.hpa)){
  
  ## open device
  pdf(paste0("../graphics/tissue.cell.type.enrichment.20250612.",gsub(" ", "_", j),".pdf"), width = 6.3, height = 6.3)
  ## define graphical parameters
  par(mfrow=c(1,1), mar=rep(0,4))
  
  ## clear any previous plottings
  circos.clear()
  
  ## create input
  tmp               <- res.enr.cell.type.all[ pval < .05/nrow(res.enr.cell.type.all) & or > 1, .(label.hpa, label.feature, cl, order, cat.srt)]
  ## dummy
  tmp[, value := 1]
  ## new order
  tmp                <- tmp[, .(label.hpa, label.feature, value, order, cl, cat.srt)]
  # tmp                <- tmp[ order(order, cat.srt)]
  ## define vector to customize the plot
  chord.order        <- rev(c(hpa.cell.anno[ label %in% tmp$label.hpa]$label, lab.phe[ label %in% tmp$label.feature]$label))
  ## colour
  chord.color        <- rev(c(rep("grey80", length(unique(tmp$label.hpa))), lab.phe$cl[lab.phe$label %in% tmp$label.feature])) 
  ## assign names
  names(chord.color) <- chord.order
  
  ## define gap parameters
  circos.par(gap.after = c(rep(.5, length(unique(tmp[[2]]))-1), 3,
                           rep(2, length(unique(tmp[[1]]))-1), 3),
             start.degree = 160)
  
  # 
  ## create diagram
  chordDiagram(tmp,
               ## define order based on grouping of labels
               order = chord.order,
               ## define colour vector
               grid.col = chord.color,
               annotationTrack = "grid",
               annotationTrackHeight = .01,
               ## larger layer for labels
               preAllocateTracks = list(track.height = .35),
               ## colour links according to characteristics, only colour cell-type of interest
               col = ifelse(tmp$label.hpa == j, tmp$cl, adjust_transparency(tmp$cl, .2)))
  
  ## change orientation of labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .3)
    # circos.axis(h = "top", labels.cex = .3, major.tick.length = .01, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  
  dev.off()
}

#---------------------------------------#
##--               SF3               --##
#---------------------------------------#

pdf("../graphics/tissue.tissue.enrichment.20250612.pdf", width = 6.3, height = 6.3)

## define graphical parameters
par(mfrow=c(1,1), mar=rep(0,4))

## clear any previous plottings
circos.clear()

## create input
tmp               <- res.enr.tissue.all[ pval < .05/nrow(res.enr.tissue.all) & or > 1, .(label.hpa, label.feature, cl, order, cat.srt)]
## dummy
tmp[, value := 1]
## new order
tmp                <- tmp[, .(label.hpa, label.feature, value, order, cl, cat.srt)]
# tmp                <- tmp[ order(order, cat.srt)]
## define vector to customize the plot
chord.order        <- rev(c(hpa.tissue.anno[ label %in% tmp$label.hpa]$label, lab.phe[ label %in% tmp$label.feature]$label))
## colour
chord.color        <- rev(c(rep("grey80", length(unique(tmp$label.hpa))), lab.phe$cl[lab.phe$label %in% tmp$label.feature])) 
## assign names
names(chord.color) <- chord.order

## define gap parameters
circos.par(gap.after = c(rep(.5, length(unique(tmp[[2]]))-1), 3,
                         rep(2, length(unique(tmp[[1]]))-1), 3),
           start.degree = 160)

# 
## create diagram
chordDiagram(tmp,
             ## define order based on grouping of labels
             order = chord.order,
             ## define colour vector
             grid.col = chord.color,
             annotationTrack = "grid",
             annotationTrackHeight = .01,
             ## larger layer for labels
             preAllocateTracks = list(track.height = .35),
             ## colour links according to characteristics
             col = tmp$cl)

## change orientation of labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .4)
  # circos.axis(h = "top", labels.cex = .3, major.tick.length = .01, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

##################################################
####  meta-regression protein characteristics ####
##################################################

#------------------------------#
##--      query UniProt     --##
#------------------------------#

## edit some odd namings
tmp.uni                  <- unique(unlist(lapply(lab.prot$UniProt, function(x) strsplit(x, "_")[[1]])))

## query uniprot instead: https://github.com/VoisinneG/queryup
require(queryup)
uni.df                   <- query_uniprot(list(accession_id = tmp.uni, organism_id = "9606"), 
                                          columns = c("id", "accession", "gene_names", "xref_proteomes", "fragment", "organelle", "length", "mass",
                                                      "ft_non_cons", "ft_non_std", "ft_non_ter", "sequence", "ft_act_site", "ft_binding",
                                                      "cc_activity_regulation", "cc_catalytic_activity", "cc_cofactor", "ft_dna_bind",
                                                      "kinetics", "ph_dependence", "redox_potential", "temp_dependence",
                                                      "cc_subunit", "cc_pharmaceutical", "ft_intramem", "cc_subcellular_location",
                                                      "ft_topo_dom", "ft_transmem", "ft_chain", "ft_crosslnk", "ft_disulfid", "ft_carbohyd", "ft_init_met",
                                                      "cc_tissue_specificity", "ft_lipid", "ft_mod_res", "ft_peptide", "cc_ptm", "ft_strand",
                                                      "ft_helix", "ft_turn", "ft_coiled", "ft_region", "protein_families", "ft_zn_fing"),
                                          show_progress = FALSE)
## recode to simplify
uni.df[, -c(1:4,7:8, 12)] <- apply(uni.df[, -c(1:4,7:8, 12)], 2, function(x) ifelse(x == "", 0, 1))

## convert to data table
uni.df                    <- as.data.table(uni.df)

## add Assay ID
uni.df[, id := sapply(Entry, function(x){
  ## get matching ids
  return(prot.order[ grep(x, UniProt)]$id[1])
})]
## three proteins are missing

## edit names
names(uni.df)             <- sapply(names(uni.df), function(x){
  if(x %in% return_fields$label){
    return(return_fields$field[ which(return_fields$label == x)])
  }else{
    return(x)
  }
})
## delete one
uni.df$`Subcellular location [CC]` <- NULL

#------------------------------#
##--   additional HPA data  --##
#------------------------------#

## create Protein class data
tmp        <- lapply(1:nrow(hpa.mapping), function(x) return(data.table(Gene=hpa.mapping$Gene[x], protein.class = strsplit(hpa.mapping$`Protein class`[x], ", ")[[1]])))
tmp        <- rbindlist(tmp)
tmp[, ind := 1]
## convert to wide
tmp        <- dcast(tmp, Gene ~ protein.class, fill = 0, value.var = "ind")
## change names
names(tmp) <- gsub(" |-", "_", names(tmp))

#------------------------------#
##--     create data set    --##
#------------------------------#

## merge with protein info
prot.rf        <- merge(prot.order, unique(uni.df), by="id") ## duplication!
## add additional HPA data; this will drop ambigious protein assignments
prot.rf        <- merge(prot.rf, tmp, by.x = "Assay", by.y = "Gene")
## edit names
names(prot.rf) <- gsub(" |-", "_", names(prot.rf))

## define vector of input features
rf.feature     <- c("xref_proteomes", "fragment", "organelle", "length", "mass",
                    "ft_non_cons", "ft_non_std", "ft_non_ter", "ft_act_site", "ft_binding",
                    "cc_activity_regulation", "cc_catalytic_activity", "cc_cofactor", "ft_dna_bind",
                    "kinetics", "ph_dependence", "redox_potential", "temp_dependence",
                    "cc_pharmaceutical", "ft_intramem", 
                    "ft_topo_dom", "ft_transmem", "ft_chain", "ft_crosslnk", "ft_disulfid", "ft_carbohyd", "ft_init_met",
                    "cc_tissue_specificity", "ft_lipid", "ft_mod_res", "ft_peptide", "cc_ptm", "ft_strand",
                    "ft_helix", "ft_turn", "ft_coiled", "ft_region", "protein_families", "ft_zn_fing",
                    names(tmp)[-1], "Secretome_location", "Panel", "chromosome_name", "mv", "md",
                    "std", "skew", "curt", "miss.per", "RNA_tissue_specificity", "RNA_single_cell_type_specificity")

## recode character as factor
rf.feature     <- data.table(var = rf.feature, type = sapply(rf.feature, function(x) class(unlist(prot.rf[, ..x]))))

## recode
prot.rf        <- as.data.frame(prot.rf)

## recode factor based on most common level
for(j in rf.feature[ type == "character"]$var){
  ## get frequency
  jj           <- table(prot.rf[, j])
  ## redefine
  prot.rf[, j] <- factor(prot.rf[, j], levels = names(jj[order(jj, decreasing = T)]))
}
## convert back
prot.rf        <- as.data.table(prot.rf)

#------------------------------#
##--       run Boruta       --##
#------------------------------#

## required package
require(Boruta)

## run feature selection
set.seed(42)
prot.boruta.all <- Boruta(prot.rf[, rf.feature$var, with=F], prot.rf[, p.50], maxRuns = 500)
plot(prot.boruta.all)
## assign last tentative ones
prot.boruta.all <- TentativeRoughFix(prot.boruta.all)

## add median value to protein features
rf.feature[, boruta.median := apply(prot.boruta.all$ImpHistory[, rf.feature$var], 2, function(x) median(x[is.finite(x)]))]
rf.feature[, boruta.selected := prot.boruta.all$finalDecision[ rf.feature$var]]

## add label
rf.feature[, label := sapply(var, function(x){
  if(x %in% return_fields$field){
    return(return_fields$label[ which(return_fields$field == x)])
  }else{
    return("")
  }
})] 

## replace infinite values
prot.boruta.all$ImpHistory[!is.finite(prot.boruta.all$ImpHistory)] <- NA

## add colour code fro groupings
col.rf          <- data.table(category = unique(rf.feature$category), cl = c("#ffaf24", "#c2364c", "#575fa1"))
rf.feature      <- merge(rf.feature, col.rf)

## order by strength
rf.feature      <- rf.feature[order(boruta.median)]

#------------------------------#
##--  Supplementary Figure  --##
#------------------------------#

## write to PDF
pdf("../graphics/Boruta.variance.explained.20250612.pdf", width = 6.3, height = 6.3)
## graphical parameters
par(mar=c(1.5,10,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), lwd=.5, yaxs="i")
## layout
layout(matrix(c(1,2,3,1,4,4,1,5,6,1,7,8), 4, 3, byrow = T), widths = c(.5,.25,.25))

#-------------------------------------#
##--       Variable importance     --##
#-------------------------------------#

## boxplots
boxplot(prot.boruta.all$ImpHistory[, rf.feature$var], cex=.3, medlwd=.7, whisklty=1, staplelty=0,
        col=ifelse(rf.feature$boruta.selected == "Confirmed", rf.feature$cl, adjust_transparency(rf.feature$cl, .2)), 
        horizontal=T, yaxt="n", xaxt="n", xlab = "Variable importance", xlim=c(.5,nrow(rf.feature)+.5))
## add x-axis
axis(1, lwd=.5)

## separate non-informative features by line
abline(h=42.5, lwd=.5, lty=2)

## add label
pm <- par("usr")
text(pm[1], 1:nrow(rf.feature), pos=2, cex=.6, labels=rf.feature$label, xpd=NA)

## add legend
legend("bottomright", bty="n", cex=.7, lty=0, pch=22,
       pt.lwd=.2, pt.bg = c(col.rf$cl, "white"), 
       legend = c(col.rf$category, "Not confirmed"))

#-------------------------------------#
##--            Skewness           --##
#-------------------------------------#

## adjsut plotting parameters
par(mar=c(1.5,1.5,.5,.5), yaxs="r")

## simple plot
plot((1-p.50) ~ skew, prot.rf, pch=20, cex=.3, col=rgb(0,0,0,.2), 
     xlab="Skewness", ylab="Explained variance [%]", 
     xaxt="n", yaxt="n")
## add axis
axis(1, lwd=.5); axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)

#-------------------------------------#
##--            Missing            --##
#-------------------------------------#

## simple plot
boxplot((1-p.50) ~ cut(miss.per, c(seq(0,10,2), 20)), prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", 
        xlab="Missing values [%]", ylab="Explained variance [%]")
## add axis
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)

#-------------------------------------#
##--      Secretome location       --##
#-------------------------------------#

## adjsut plotting parameters
par(mar=c(8,1.5,.5,.5))

## simple plot
boxplot((1-p.50) ~ Secretome_location, prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", xaxt="n",
        xlab="", ylab="Explained variance [%]")
## add axis
axis(1, lwd=.5, at=1:nlevels(prot.rf$Secretome_location), labels = NA)
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)
## add levels
pm <- par("usr")
text(1:nlevels(prot.rf$Secretome_location), pm[3]-(pm[4]-pm[3])*.05, pos=2, offset = 0, srt=60,
     xpd=NA, labels = sapply(levels(prot.rf$Secretome_location), function(x) ifelse(x == "", "Not secreted", x)),
     cex=.6)

#-------------------------------------#
##--             Panel             --##
#-------------------------------------#

## adjsut plotting parameters
par(mar=c(5,1.5,.5,.5))


## simple plot
boxplot((1-p.50) ~ Panel, prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", xaxt="n",
        xlab="", ylab="Explained variance [%]")
## add axis
axis(1, lwd=.5, at=1:nlevels(prot.rf$Panel), labels = NA)
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)
## add levels
pm <- par("usr")
text(1:nlevels(prot.rf$Panel), pm[3]-(pm[4]-pm[3])*.05, pos=2, offset = 0, srt=60,
     xpd=NA, labels = sapply(levels(prot.rf$Panel), function(x) ifelse(x == "", "Not secreted", x)),
     cex=.6)

#-------------------------------------#
##--         Glycosylation         --##
#-------------------------------------#

## simple plot
boxplot((1-p.50) ~ ft_carbohyd, prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", xaxt="n",
        xlab="Glycosylation", ylab="Explained variance [%]")
## add axis
axis(1, lwd=.5, at=1:2, labels = NA)
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)
## add levels
pm <- par("usr")
text(1:2, pm[3]-(pm[4]-pm[3])*.05, pos=2, offset = 0, srt=60,
     xpd=NA, labels = c("no", "yes"),
     cex=.6)


#-------------------------------------#
##--     RNA_tissue_specificity    --##
#-------------------------------------#

## simple plot
boxplot((1-p.50) ~ RNA_tissue_specificity, prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", xaxt="n",
        xlab="", ylab="Explained variance [%]")
## add axis
axis(1, lwd=.5, at=1:nlevels(prot.rf$RNA_tissue_specificity), labels = NA)
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)
## add levels
pm <- par("usr")
text(1:nlevels(prot.rf$RNA_tissue_specificity), pm[3]-(pm[4]-pm[3])*.05, pos=2, offset = 0, srt=60,
     xpd=NA, labels = levels(prot.rf$RNA_tissue_specificity),
     cex=.6)


#-------------------------------------#
##--         Di-sulfite bond       --##
#-------------------------------------#

## simple plot
boxplot((1-p.50) ~ ft_disulfid, prot.rf,
        cex=.3, medlwd=.7, whisklty=1, staplelty=0, 
        yaxt="n", xaxt="n",
        xlab="Disulfide bond", ylab="Explained variance [%]")
## add axis
axis(1, lwd=.5, at=1:2, labels = NA)
axis(2, lwd=.5, at=seq(0,.8,.2), labels = seq(0,.8,.2)*100)
## add levels
pm <- par("usr")
text(1:2, pm[3]-(pm[4]-pm[3])*.05, pos=2, offset = 0, srt=60,
     xpd=NA, labels = c("no", "yes"),
     cex=.6)

## close device
dev.off()


#------------------------------#
##--    numbers manuscript  --##
#------------------------------#

## variance for proteins secreted
aggregate((1-p.50) ~ Secretome_location, prot.rf, summary)

## missing proteins
aggregate((1-p.50) ~ ifelse(miss.per >= 5, 0, 1), prot.rf, summary)
nrow(prot.rf[ miss.per <= 5])

#####################################################################
####                    Results ethnic groups                    ####
#####################################################################

## get the output
jj <- dir("../output_updated///")
## get only the variance components
jj <- grep("variance", jj, value=T)
## restrict to population results
jj <- grep("EUR|AFR|CSA", jj, value=T)

## investigate what has been missing
tmp.pop          <- as.data.table(tmp.pop)
tmp.pop[, missing.lasso := !(paste("lasso.explained.variance", population, olink.id, "txt", sep=".") %in% jj)]

## --> import <-- ##

## run import in parallel
registerDoMC(10)

## import
res.var.pop      <- mclapply(jj, function(x){
  ## import 
  tmp <- fread(paste0("../output_updated/", x))
  ## add population
  tmp[, population := strsplit(x, "\\.")[[1]][4]]
  ## edit names
  names(tmp) <- sapply(names(tmp), function(x) ifelse(x %in% lab.ext$short_name_new, lab.ext$short_name[ which(lab.ext$short_name_new == x)], x))
  ## account for genetic scores
  names(tmp) <- gsub(paste0("\\.", tmp$var[1]), "", names(tmp))
  return(tmp)
}, mc.cores = 10)
## combine
res.var.pop      <- rbindlist(res.var.pop, fill = T)
## look for possibly poorly run estimates
summary(res.var.pop$n.boot)

## replace missing summary measure with p.50!
res.var.pop[, summary.measure := ifelse(is.na(summary.measure), "p.50", summary.measure)]

## transform in more efficient shape
res.var.pop.long <- melt.data.table(res.var.pop, id.vars = c("var", "summary.measure", "n.boot", "population", "n"))
## reshape again by summary measure
res.var.pop.long <- dcast(res.var.pop.long, var + variable + population + n + n.boot  ~ summary.measure, sep = ".", value.var = "value")
## drop missing values
res.var.pop.long <- res.var.pop.long[ !is.na(p.50) ]
## avoid coding as a factor
res.var.pop.long[, variable := as.character(variable)]

#####################################################################
####                     Compare ancestries                      ####
#####################################################################

## create summary by ancestry
prot.ancestry        <- dcast(res.var.pop.long[ variable == "Residuals", .(var, p.50, population)], var ~ population, sep=".", value.var = "p.50")
names(prot.ancestry) <- c("var", "p.50.AFR", "p.50.CSA", "p.50.EUR")
prot.ancestry        <- merge(prot.ancestry, dcast(res.var.pop.long[ variable != "Residuals", .(count = length(variable), max = max(p.50)), by = c("population", "var")], 
                                                   var ~ population, value.var = c("count", "max"), sep="."), all.x=T)
## add strongest contributing factor
tmp                   <- res.var.pop.long[ variable != "Residuals" ]
tmp                   <- tmp[order(population, var, -p.50)]
tmp[, ind := 1:.N, by=c("population", "var")]
tmp                   <- tmp[ ind == 1]
## add
prot.ancestry         <- merge(prot.ancestry, dcast(tmp, var ~ population, value.var = "variable"), all.x = T)

## create summary using equal numbers of variables - AFR
prot.ancestry[, p.50.EUR.num.AFR := apply(prot.ancestry[, .(var, count.AFR)], 1, function(x){
  
  if(!is.na(x[2])){
    ## get the explained variance components
    tmp <- res.var.pop.long[ population == "EUR" & var == x[1] & variable != "Residuals"]
    ## order
    return(1-sum(sort(tmp$p.50, decreasing = T)[1:as.numeric(x[2])], na.rm = T))
  }else{
    return(NA)
  }
  
})]

## create summary using equal numbers of variables - CSA
prot.ancestry[, p.50.EUR.num.CSA := apply(prot.ancestry[, .(var, count.CSA)], 1, function(x){
  
  if(!is.na(x[2])){
    ## get the explained variance components
    tmp <- res.var.pop.long[ population == "EUR" & var == x[1] & variable != "Residuals"]
    ## order
    return(1-sum(sort(tmp$p.50, decreasing = T)[1:as.numeric(x[2])], na.rm = T))
  }else{
    return(NA)
  }
  
})]

## --> compare against European - AFR <-- ## (LAIR1 - https://pubmed.ncbi.nlm.nih.gov/38106031/)
summary(prot.ancestry$p.50.EUR - prot.ancestry$p.50.AFR) 
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.778127 -0.126401 -0.050685 -0.069232 -0.003232  0.576327 
wilcox.test(prot.ancestry$p.50.EUR, prot.ancestry$p.50.AFR)$p.value
# 1.491814e-91

## same number of variables
summary(prot.ancestry$p.50.EUR.num.AFR - prot.ancestry$p.50.AFR) 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.7049 -0.0658 -0.0138 -0.0191  0.0315  0.6868     556 
wilcox.test(prot.ancestry$p.50.EUR.num.AFR, prot.ancestry$p.50.AFR)$p.value
# 1.009336e-66

## --> compare against European - CSA <-- ##
summary(prot.ancestry$p.50.EUR - prot.ancestry$p.50.CSA)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.564842 -0.043555 -0.003299 -0.002837  0.039016  0.494035  
wilcox.test(prot.ancestry$p.50.EUR, prot.ancestry$p.50.CSA)$p.value
# 0.01304485

## same number of variables
summary(prot.ancestry$p.50.EUR.num.CSA - prot.ancestry$p.50.CSA)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.3104 -0.0214  0.0124  0.0125  0.0516  0.3081     458 
wilcox.test(prot.ancestry$p.50.EUR.num.CSA, prot.ancestry$p.50.CSA)$p.value
# 1.587841e-10

#---------------------------------------#
##--         strong differences      --##
#---------------------------------------#

## --> flag extreme outliers for each comparison <-- ##

## EUR vs AFR
prot.ancestry[, outlier.EUR.AFR := ifelse(abs(p.50.EUR.num.AFR - p.50.AFR) > .25, 1, 0)]

## EUR vs CSA
prot.ancestry[, outlier.EUR.CSA := ifelse(abs(p.50.EUR.num.CSA - p.50.CSA) > .25, 1, 0)]

#---------------------------------------#
##--         import genetics         --##
#---------------------------------------#

## import Supplemental Table with all cis/trans pQTLs (from: https://www.nature.com/articles/s41586-023-06592-6)
pqtl     <- as.data.table(read_excel("../input/41586_2023_6592_MOESM3_ESM (1).xlsx", 11, skip = 4))
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

## add protein cluster
prot.ancestry <- merge(prot.ancestry, prot.order[, .(id, category, order, cl.cluster)], by.x="var", by.y="id")

## export for interaction testing
fwrite(res.var.pop.long, "Results.variance.decomposition.UKB.Olink.ancestry.20250616.txt", sep="\t", row.names = F, na=NA)
## protein summary
fwrite(prot.ancestry, "Protein.variance.summary.UKB.Olink.ancestry.20250616.txt", sep="\t", row.names = F, na=NA)

## comprehensive file
tmp <- rbindlist(list(res.var.sex.long, res.var.pop.long), fill = T)
tmp <- merge(tmp, lab.phe[, .(short_name, label)], by.x = "variable", by.y = "short_name")
tmp[, Assay := toupper(var)]
fwrite(tmp[, .(variable, label, Assay, population, n, p.025, p.50, p.975)], "Results.variance.decomposition.UKB.Olink.all.strata.20250617.txt", sep="\t", row.names = F, na=NA)

## --> import interaction results <-- ##
res.pop.interaction <- fread("../../03_regression_analysis/input/Results.ancestry.interaction.UKB.olink.feature.selection.all.20260617.txt")
## process variable names

## add label
res.pop.interaction <- merge(res.pop.interaction, lab.phe, by.x="variable", by.y="short_name")

## write a list of combinations to pull from UKB-PPP stats
res.pop.interaction[ variable == "cis.score" & (pval.inter.AFR < .05/nrow(res.pop.interaction) | pval.inter.CSA < .05/nrow(res.pop.interaction))]
write.table(unique(res.pop.interaction[ variable == "cis.score" & (pval.inter.AFR < .05/nrow(res.pop.interaction) | pval.inter.CSA < .05/nrow(res.pop.interaction)), .(var)]),
            "Proteins.cis.pQTLs.differential.20250617.txt", sep="\t", row.names = F)

## add variance explained estimates
res.pop.interaction <- merge(res.pop.interaction, res.var.pop.long[population == "AFR", .(var, variable, p.50)], by = c("var", "variable"), all.x=T)
res.pop.interaction <- merge(res.pop.interaction, res.var.pop.long[population == "CSA", .(var, variable, p.50)], by = c("var", "variable"), all.x=T, suffixes = c(".AFR", ".CSA"))
res.pop.interaction <- merge(res.pop.interaction, res.var.pop.long[population == "EUR", .(var, variable, p.50)], by = c("var", "variable"), all.x=T)
## how many extreme examples (replace by zero); p.50 is the estimates from EUR
res.pop.interaction[, p.50.AFR := ifelse(is.na(p.50.AFR), 0, p.50.AFR)]
res.pop.interaction[, p.50.CSA := ifelse(is.na(p.50.CSA), 0, p.50.CSA)]
res.pop.interaction[, p.50 := ifelse(is.na(p.50), 0, p.50)]


#-------------------------------#
##--   pull in allele freq   --##
#-------------------------------#

## import allele freq in EUR
registerDoMC(10)
allele.eur          <- mclapply(1:23, function(x){
  ## import
  tmp <- fread(paste0("<path>/ukb_imp_EUR_chr", ifelse(x == 23, "X", x),"_snpstat.out"), 
               select = c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB", "alleleA_frequency", "info"))
  ## return only variants needed
  return(tmp[ rsid %in% pqtl$rsID])
}, mc.cores=10)
## combine
allele.eur          <- rbindlist(allele.eur)

## import allele freq in AFR
registerDoMC(10)
allele.afr          <- mclapply(1:23, function(x){
  ## import
  tmp <- fread(paste0("<path>/ukb_imp_AFR_chr", ifelse(x == 23, "X", x),"_snpstat.out"), 
               select = c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB", "alleleA_frequency", "info"))
  ## return only variants needed
  return(tmp[ rsid %in% pqtl$rsID])
}, mc.cores=10)
## combine
allele.afr          <- rbindlist(allele.afr)

## import allele freq in EUR
registerDoMC(10)
allele.csa          <- mclapply(1:23, function(x){
  ## import
  tmp <- fread(paste0("<path>/ukb_imp_CSA_chr", ifelse(x == 23, "X", x),"_snpstat.out"), 
               select = c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB", "alleleA_frequency", "info"))
  ## return only variants needed
  return(tmp[ rsid %in% pqtl$rsID])
}, mc.cores=10)
## combine
allele.csa          <- rbindlist(allele.csa)

## combine into one large file
allele.pqtl         <- merge(allele.afr, allele.csa, by=c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB"), suffixes = c(".AFR", ".CSA"))
allele.pqtl         <- merge(allele.pqtl, allele.eur, by=c("alternate_ids", "rsid", "chromosome", "position", "alleleA", "alleleB"))
## contains multi-allelic variants...
allele.pqtl <- allele.pqtl[order(rsid, alternate_ids), ]
allele.pqtl[, ind := 1:.N, by="rsid"]
allele.pqtl[ ind > 1]
## create allele mapper
allele.pqtl[, alleles := paste(pmin(alleleA, alleleB), pmax(alleleA, alleleB), sep="_")]

## add to the pQTL data
pqtl                <- merge(pqtl, unique(allele.pqtl), by.x = c("rsID", "alleles"), by.y = c("rsid", "alleles"), all.x=T)
## delete three duplication
tail(sort(table(paste(pqtl$rsID, pqtl$`UKBPPP ProteinID`))))
pqtl[ rsID == "rs200683903" & `UKBPPP ProteinID` == "GPHA2:Q96T91:OID31437:v1"]
pqtl                <- pqtl[-which(pqtl$rsID == "rs200683903" & pqtl$`UKBPPP ProteinID` == "GPHA2:Q96T91:OID31437:v1" & is.na(pqtl$`Distance to gene`))]
pqtl[ rsID == "rs533589698" & `UKBPPP ProteinID` == "SCGB3A1:Q96QR1:OID30577:v1"]
pqtl                <- pqtl[-which(pqtl$rsID == "rs533589698" & pqtl$`UKBPPP ProteinID` == "SCGB3A1:Q96QR1:OID30577:v1" & pqtl$alleleA_frequency.AFR > .9)]
pqtl[ rsID == "rs375092822" & `UKBPPP ProteinID` == "TCL1A:P56279:OID20987:v1"]
pqtl                <- pqtl[-which(pqtl$rsID == "rs375092822" & pqtl$`UKBPPP ProteinID` == "TCL1A:P56279:OID20987:v1" & pqtl$alleleA_frequency > .9)]

## add to interaction results
pqtl[, variable := paste(`cis/trans`, "score", sep=".")]
tmp                 <- pqtl[ variable == "cis.score", .(protein_id, variable, alleleA_frequency, alleleA_frequency.AFR, alleleA_frequency.CSA)]
## make unique by averaging
tmp                 <- tmp[, .(alleleA_frequency=mean(alleleA_frequency), 
                               alleleA_frequency.AFR=mean(alleleA_frequency.AFR),
                               alleleA_frequency.CSA=mean(alleleA_frequency.CSA)), 
                           by=c("protein_id", "variable")]
res.pop.interaction <- merge(res.pop.interaction, tmp, by.x = c("var", "variable"), by.y = c("protein_id", "variable"), all.x = T)
## create differences
res.pop.interaction[, p.50.difference.AFR.EUR := abs(p.50.AFR - p.50)]
res.pop.interaction[, p.50.difference.CSA.EUR := abs(p.50.CSA - p.50)]
res.pop.interaction[, allele.difference.AFR.EUR := abs(alleleA_frequency - alleleA_frequency.AFR)]
res.pop.interaction[, allele.difference.CSA.EUR := abs(alleleA_frequency - alleleA_frequency.CSA)]
## add overall explained variance to normalize effects
res.pop.interaction <- merge(res.pop.interaction, dcast(res.var.pop.long[ variable == "Residuals"], var  ~ population, sep=".", value.var = c("p.50")),
                             by="var")

## add protein cluster colours
res.pop.interaction <- merge(res.pop.interaction, prot.ancestry[, .(var, cl.cluster)])

#-----------------------------------------------#
##--  import results testing for diff. vars  --##
#-----------------------------------------------#

## import results testing for differential causal variants
cis.ancestry <- fread("../../03_regression_analysis/input/Results.cis.pQTL.ancestries.explained.var.top.ancestral.20250618.txt")
cis.inter    <- fread("../../03_regression_analysis/input/Results.cis.pQTL.ancestries.interaction.all.ancestral.20250618.txt")
cis.expl     <- fread("../../03_regression_analysis/input/Results.cis.pQTL.ancestries.explained.var.all.ancestral.20250618.txt")

## look into some extreme changes; add data first
cis.ancestry[, var := tolower(Assay) ]
cis.ancestry <- merge(cis.ancestry, res.var.pop.long[variable == "cis.score", .(var, population, p.50)],
                      by.x = c("var", "ancestry"), by.y = c("var", "population"), all.x=T)
## look for extreme differences
cis.ancestry[ abs(p.50 - r2) > .4]
cis.ancestry[ var == "arsa"]

#-------------------------------#
##--     summary figure      --##
#-------------------------------#

## open device
pdf("../graphics/comparison.feature.selection.ancestry.20250618.pdf", width = 6.3, height = 6.3*(2/3))
## plotting parameters
par(mar=c(1.5,1.5,.5,.5), cex=.5, mgp=c(.6,0,0), tck=-.01, lwd=.5, cex.lab=.5, cex.axis=.5, xaxs="r", yaxs="r")
layout(matrix(c(1,2,3,4,5,5), 2, 3))

## --> EUR vs AFR <-- ##

## empty plot
plot(c(0,.9), c(0,.9), xaxt="n", yaxt="n", xlab="Explained variance [%] - White European*", ylab="Explained variance [%] - British African",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
axis(2, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
## add diagonale
abline(a=0, b=1, lwd=.5)
## add points
points(1-prot.ancestry$p.50.EUR.num.AFR, 1-prot.ancestry$p.50.AFR, cex=.5, col=prot.ancestry$cl.cluster, pch=20)
## annotate strong outliers
jj  <- which(prot.ancestry$outlier.EUR.AFR == 1)
## annotate
text(1-prot.ancestry$p.50.EUR.num.AFR[jj], 1-prot.ancestry$p.50.AFR[jj], cex=.4, labels = toupper(prot.ancestry$var[jj]), xpd=NA)

## --> EUR vs CSA <-- ##

## empty plot
plot(c(0,.9), c(0,.9), xaxt="n", yaxt="n", xlab="Explained variance [%] - White European*", ylab="Explained variance [%] - British Central South Asian",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
axis(2, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
## add diagonale
abline(a=0, b=1, lwd=.5)
## add points
points(1-prot.ancestry$p.50.EUR.num.CSA, 1-prot.ancestry$p.50.CSA, cex=.5, col=prot.ancestry$cl.cluster, pch=20)
## annotate strong outliers
jj  <- which(prot.ancestry$outlier.EUR.CSA == 1)
## annotate
text(1-prot.ancestry$p.50.EUR.num.CSA[jj], 1-prot.ancestry$p.50.CSA[jj], cex=.4, labels = toupper(prot.ancestry$var[jj]), xpd=NA)
## add colour legend
legend("topleft", lty=0, pch=22, pt.bg=cluster.assignment$cl.cluster, pt.lwd=.3, pt.cex=.9,
       cex=.5, legend=cluster.assignment$label, lwd=.1, bty="n")

## --> differences in Allele frequencies <-- ##

## empty plot: EUR vs AFR
plot(c(0,.9), c(0,4), xaxt="n", yaxt="n", xlab="Absolute difference in allele frequency [%] - EUR/AFR", ylab="Fold change in explained variance cis-pQTL - EUR/AFR",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75,1), labels = c(0,25,50,75,100))
axis(2, lwd=.5)
## add points
tmp <-  res.pop.interaction[ p.50 > 0 & p.50.AFR > 0 & variable == "cis.score" & !is.na(allele.difference.AFR.EUR)]
points(tmp$allele.difference.AFR.EUR, tmp$p.50.difference.AFR.EUR/(1-tmp$EUR), cex=.5, pch=20, col=tmp$cl.cluster)
## annotate some
tmp[, plt := (allele.difference.AFR.EUR - (tmp$p.50.difference.AFR.EUR/(1-tmp$EUR)))^2]
tmp <- tmp[ order(-abs(plt))]
tmp <- tmp[1:10]
text(tmp$allele.difference.AFR.EUR, tmp$p.50.difference.AFR.EUR/(1-tmp$EUR), cex=.4, labels = toupper(tmp$var), xpd=NA)

## empty plot: EUR vs CSA
plot(c(0,.9), c(0,4), xaxt="n", yaxt="n", xlab="Absolute difference in allele frequency [%] - EUR/CSA", ylab="Fold change in explained variance cis-pQTL - EUR/CSA",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75,1), labels = c(0,25,50,75,100))
axis(2, lwd=.5)
## add points
tmp <-  res.pop.interaction[ p.50 > 0 & p.50.CSA > 0 & variable == "cis.score" & !is.na(allele.difference.CSA.EUR)]
points(tmp$allele.difference.CSA.EUR, tmp$p.50.difference.CSA.EUR/(1-tmp$EUR), cex=.5, pch=20, col=tmp$cl.cluster)
## annotate some
tmp[, plt := (allele.difference.CSA.EUR - (tmp$p.50.difference.CSA.EUR/(1-tmp$EUR)))^2]
tmp <- tmp[ order(-abs(plt))]
tmp <- tmp[1:10]
text(tmp$allele.difference.CSA.EUR, tmp$p.50.difference.CSA.EUR/(1-tmp$EUR), cex=.4, labels = toupper(tmp$var), xpd=NA)


## --> Interaction testing <-- ##

## adopt plotting range
par(mar=c(1.5,7,1.5,.5), xaxs="i", yaxs="i")

## reduce to effects with differential ancestral variants
jj  <- unique(cis.ancestry[ EUR.r2 == 0 | CSA.r2 == 0 | AFR.r2 == 0]$Assay)
tmp <- res.pop.interaction[ var %in% tolower(jj) & variable == "cis.score"]
## add ethnic-specific estimates
foo <- dcast(cis.ancestry, Assay ~ ancestry, value.var = "r2")
foo[, var := tolower(Assay)]
tmp <- merge(tmp, foo, suffixes = c(".trans", ".specific"))

## order
tmp <- tmp[order(-abs(pmax(abs(tval.inter.AFR), abs(tval.inter.CSA))))]
# tmp <- tmp[1:30]

## empty plotting frame
plot(c(0,.7), c(.5,nrow(tmp)+.5), xlab="Explained variance [%]", ylab="", xaxt="n", yaxt="n", type="n", ylim=rev(c(.5,nrow(tmp)+.5)))
## store parameters
pm <- par("usr")
## display adjusted axis
axis(1, lwd=.5, at=seq(0,.7,.1), labels=seq(0,70,10))

## background best ethnic variant each
rect(0, 1:nrow(tmp)-.45,  tmp$EUR.specific, 1:nrow(tmp)-.15,  col=adjust_transparency("#F9A12EFF", .2), lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)-.15,  tmp$AFR.specific, 1:nrow(tmp)+.15,  col=adjust_transparency("#FC766AFF", .2), lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)+.15,  tmp$CSA.specific, 1:nrow(tmp)+.4,  col=adjust_transparency("#9B4A97FF", .2), lwd=.2, border="grey30")

## add rectangles - each ancestry
rect(0, 1:nrow(tmp)-.45,  tmp$p.50, 1:nrow(tmp)-.15,  col="#F9A12EFF", lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)-.15,  tmp$p.50.AFR, 1:nrow(tmp)+.15,  col="#FC766AFF", lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)+.15,  tmp$p.50.CSA, 1:nrow(tmp)+.4,  col="#9B4A97FF", lwd=.2, border="grey30")
## add names
text(pm[1], 1:nrow(tmp), cex=.5, labels=paste0(tmp$label, " - ", toupper(tmp$var)), xpd=NA, col="grey30", pos=2)

## add legend
legend("bottomright", bty="n", lty=0, pt.cex=1,
       pt.bg = c("#F9A12EFF", "#FC766AFF", "#9B4A97FF", "grey30", "grey90"), pch=22,
       legend = c("EUR", "AFR", "CSA", "trans-ethnic", "ethnic-specific"), cex=.5)

## close device
dev.off()

#####################################################################
####                     Compare the sexes                       ####
#####################################################################

#---------------------------------------#
##--         strong differences      --##
#---------------------------------------#

## --> flag extreme outliers for each comparison <-- ##

## female vs all
tmp <- resid(lm(p.50.female ~ p.50, prot.order, na.action = "na.pass"))
jj  <- which(tmp > mean(tmp) + 4*sd(tmp) | tmp < mean(tmp) - 4*sd(tmp))
prot.order[, outlier.female.all := ifelse(1:nrow(prot.order) %in% jj, 1, 0)]

## male vs all
tmp <- resid(lm(p.50.male ~ p.50, prot.order, na.action = "na.pass"))
jj  <- which(tmp > mean(tmp) + 4*sd(tmp) | tmp < mean(tmp) - 4*sd(tmp))
prot.order[, outlier.male.all := ifelse(1:nrow(prot.order) %in% jj, 1, 0)]

## female vs male
tmp <- resid(lm(p.50.male ~ p.50.female, prot.order, na.action = "na.pass"))
jj  <- which(tmp > mean(tmp) + 4*sd(tmp) | tmp < mean(tmp) - 4*sd(tmp))
prot.order[, outlier.female.male := ifelse(1:nrow(prot.order) %in% jj, 1, 0)]

#---------------------------------------#
##--          Jaccard index          --##
#---------------------------------------#

## compute jaccard index for overlap of selected features
prot.order[, jaccard.sex := sapply(id, function(x){
  ## get variable selection among females
  female       <- res.var.sex.long[variable != "Residuals" & population == "Female" & var == x]$variable
  male         <- res.var.sex.long[variable != "Residuals" & population == "Male" & var == x]$variable
  ## compute index
  intersection <- length(intersect(female, male))
  union        <- length(female) + length(male) - intersection
  return(intersection/union)
})]

## how many have a Jaccard index below .5
nrow(prot.order[ jaccard.sex < .5])/nrow(prot.order)
summary(prot.order$jaccard.sex)


#---------------------------------------#
##--       interaction results       --##
#---------------------------------------#

## import sex-interaction
res.sex.interaction <- fread("../../03_regression_analysis/input/Results.sex.interaction.UKB.olink.feature.selection.20250617.txt")
## add variance explained estimates
res.sex.interaction <- merge(res.sex.interaction, res.var.sex.long[population == "Female", ], by = c("var", "variable"), all.x=T)
res.sex.interaction <- merge(res.sex.interaction, res.var.sex.long[population == "Male", ], by = c("var", "variable"), all.x=T, suffixes = c(".female", ".male"))
## how many extreme examples (replace by zero)
res.sex.interaction[, p.50.female := ifelse(is.na(p.50.female), 0, p.50.female)]
res.sex.interaction[, p.50.male := ifelse(is.na(p.50.male), 0, p.50.male)]

#---------------------------------------#
##--                SF4              --##
#---------------------------------------#

## open device
pdf("../graphics/comparison.feature.selection.sex.all.20250619.pdf", width = 6.3, height = 6.3/2)
## plotting parameters
par(mar=c(1.5,1.5,.5,.5), cex=.5, mgp=c(.6,0,0), tck=-.01, lwd=.5, mfrow=c(1,2), cex.lab=.5, cex.axis=.5)

## --> female vs all <-- ##

## empty plot
plot(c(0,.9), c(0,.9), xaxt="n", yaxt="n", xlab="Explained variance [%] - Both sexes", ylab="Explained variance [%] - Female",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
axis(2, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
## add diagonale
abline(a=0, b=1, lwd=.5)
## add points
points(1-prot.order$p.50, 1-prot.order$p.50.female, cex=.5, col=prot.order$cl.cluster, pch=20)
## add colour legend
legend("topleft", lty=0, pch=22, pt.bg=cluster.assignment$cl.cluster, pt.lwd=.3, pt.cex=.9,
       cex=.5, legend=cluster.assignment$label, lwd=.1, bty="n")
## annotate strong outliers
tmp <- resid(lm(p.50.female ~ p.50, prot.order, na.action = "na.pass"))
jj  <- which(tmp > mean(tmp) + 4*sd(tmp) | tmp < mean(tmp) - 4*sd(tmp))
## annotate
text(1-prot.order$p.50[jj], 1-prot.order$p.50.female[jj], cex=.4, labels = prot.order$Assay[jj], xpd=NA)

## --> male vs all <-- ##

## empty plot
plot(c(0,.9), c(0,.9), xaxt="n", yaxt="n", xlab="Explained variance [%] - Both sexes", ylab="Explained variance [%] - Male",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
axis(2, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
## add diagonale
abline(a=0, b=1, lwd=.5)
## add points
points(1-prot.order$p.50, 1-prot.order$p.50.male, cex=.5, col=prot.order$cl.cluster, pch=20)
## annotate strong outliers
tmp <- resid(lm(p.50.male ~ p.50, prot.order, na.action = "na.pass"))
jj  <- which(tmp > mean(tmp) + 4*sd(tmp) | tmp < mean(tmp) - 4*sd(tmp))
## annotate
text(1-prot.order$p.50[jj], 1-prot.order$p.50.male[jj], cex=.4, labels = prot.order$Assay[jj], xpd=NA)

## close device
dev.off()
## open device
pdf("../graphics/comparison.feature.selection.sex.20250619.pdf", width = 6.3, height = 6.3/2)
## plotting parameters
par(mar=c(1.5,1.5,.5,.5), cex=.5, mgp=c(.6,0,0), tck=-.01, lwd=.5, cex.lab=.5, cex.axis=.5)
## define layout
layout(matrix(c(1,1,3,4,1,2,3,4),2,4, byrow = T), widths = c(.35,.15,.25,.25))

## --> female vs male <-- ##

## empty plot
plot(c(0,.8), c(0,.8), xaxt="n", yaxt="n", xlab="Explained variance [%] - Female", ylab="Explained variance [%] - Male",
     type="n")
## add axis
axis(1, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
axis(2, lwd=.5, at=c(0,.25,.5,.75), labels = c(0,25,50,75))
## add diagonale
abline(a=0, b=1, lwd=.5)
## add points
points(1-prot.order$p.50.female, 1-prot.order$p.50.male, cex=.5, col=prot.order$cl.cluster, pch=20)
## annotate strong outliers
jj  <- which(prot.order$outlier.female.male == 1)
## annotate
text(1-prot.order$p.50.female[jj], 1-prot.order$p.50.male[jj], cex=.5, labels = prot.order$Assay[jj], xpd=NA)
## add colour legend
legend("topleft", lty=0, pch=22, pt.bg=cluster.assignment$cl.cluster, pt.lwd=.3, pt.cex=.9,
       cex=.7, legend=cluster.assignment$label, lwd=.1, bty="n")


## --> Jaccard index <-- ##

## plot
vioplot::vioplot(prot.order$jaccard.sex, lwd=.5, ylab="Jaccard index feature overlap: Female - Male", xlab="", col="grey80",
                 xaxt="n", yaxt="n")
axis(2, lwd=.5)

## --> Sex-specific selection (top XXX) <-- ##

## adopt plotting range
par(mar=c(1.5,7,1,.5), xaxs="i", yaxs="i")

## create data frame
tmp <- res.sex.interaction[ pval.inter < .05/nrow(res.sex.interaction) & abs(p.50.female - p.50.male) > .07]
## add colour coding
tmp <- merge(tmp, cat.col, by="category")
## order
# tmp <- tmp[order(-abs(p.50.female - p.50.male))]
tmp <- tmp[order(-abs(tval.inter))]

## empty plotting frame
plot(c(0,.7), c(.5,nrow(tmp)+.5), xlab="Explained variance [%]", ylab="", xaxt="n", yaxt="n", type="n", ylim=rev(c(.5,nrow(tmp)+.5)))
## store parameters
pm <- par("usr")
## display adjusted axis
axis(1, lwd=.5, at=c(0,.1,.2,.3)+.1, labels=c(0,10,20,30))
## add rectangles - Female
rect(0, 1:nrow(tmp)-.4,  tmp$p.50.female, 1:nrow(tmp)+.0,  col=tmp$cl, lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)-.0,  tmp$p.50.male, 1:nrow(tmp)+.4,  col=tmp$cl, lwd=.2, border="grey30", density = 40)
## add names
text(pm[1], 1:nrow(tmp), cex=.5, labels=paste0(tmp$label, " - ", toupper(tmp$var)), xpd=NA, col="grey30", pos=2)
# ## add legend
# legend("topright", legend=c("Female", "Male"), cex=.5, lty=0, pch=22, bty="n",
#        pt.bg = c("grey", "white"))
## lines dividing proteins
abline(h=seq(.5,nrow(tmp)+.5,1), lwd=.3)

## add a title
mtext("Most sex-differential", cex=.5)


## --> Oxytocin example <-- ##

## create data frame
tmp <- dcast(res.var.sex.long[ var == "oxt"], var + variable ~ population, sep=".", value.var = c("n", "n.boot", "variable.rank", "p.50"))
## add label information
tmp <- merge(tmp, lab.red[, .(short_name, label, cl, cat.srt)], by.x="variable", by.y="short_name")
## order
tmp[is.na(tmp)] <- 0
# tmp <- tmp[order(cat.srt, label)]
tmp <- tmp[ order(-p.50.Female)]

## empty plotting frame
plot(sqrt(c(0,.32)), c(.5,nrow(tmp)+.5), xlab="Explained variance [%]", ylab="", xaxt="n", yaxt="n", type="n", ylim=rev(c(.5,nrow(tmp)+.5)))
## store parameters
pm <- par("usr")
## display adjusted axis
axis(1, lwd=.5, at=sqrt(c(0,.05,.1,.2,.3)), labels=c(0,5,10,20,30))
## add rectangles - Female
rect(0, 1:nrow(tmp)-.4,  sqrt(tmp$p.50.Female), 1:nrow(tmp)+.0,  col=tmp$cl, lwd=.2, border="grey30")
rect(0, 1:nrow(tmp)-.0,  sqrt(tmp$p.50.Male), 1:nrow(tmp)+.4,  col=tmp$cl, lwd=.2, border="grey30", density = 40)
## add names
text(pm[1], 1:nrow(tmp), cex=.5, labels=tmp$label, xpd=NA, col="grey30", pos=2)
## add legend
legend("topright", legend=c("Female", "Male"), cex=.5, lty=0, pch=22, bty="n",
       pt.bg = c("grey", "white"), pt.cex=1)
## lines dividing proteins
abline(h=seq(.5,nrow(tmp)+.5,1), lwd=.3)


## add a title
mtext("Plasma oxytocin", cex=.5)


## close device
dev.off()

## protein summary
fwrite(prot.order, "Protein.variance.summary.UKB.Olink.20250612.txt", sep="\t", row.names = F, na=NA)

##############################################
####            survival models           ####
##############################################

## import results
registerDoMC(20)
res.surv <- dir("<path>")
res.surv <- grep("res.surv", res.surv, value = T)
res.surv <- mclapply(res.surv, function(x){
  ## read in the file
  tmp <- fread(paste0("<path>", x))
  ## add file name
  tmp[, file := x]
  ## return results
  return(tmp)
}, mc.cores = 20)
## combine again
res.surv <- rbindlist(res.surv, fill = T)
## delete things that not have run
res.surv <- res.surv[ !is.na(beta)]
## add type of analysis
res.surv[, type.analysis := gsub(".*\\.((?:adj|cis|trans|top\\d+)(?:\\.top\\d+)?)\\..*", "\\1", file, perl = TRUE)]
res.surv[, type.analysis := gsub(".*\\.txt.*", "min", type.analysis)]
table(res.surv$type.analysis)

#-----------------------------------#
##--          summary            --##
#-----------------------------------#

## subset to what is of interest
tmp.phe  <- lab.phe[ short_name %in% unique(paste0("bin_", res.surv$phecode))] 

## prepare lab for plotting
tmp.phe[, phecode := as.numeric(gsub("bin_", "", short_name))]
tmp.phe  <- tmp.phe[ order(phecode) ]
## add sorting
tmp.phe[, srt := 1:nrow(tmp.phe)]

## get proteins by genomic position (make unique first)
lab.prot[, chr_numeric := sapply(chromosome_name, function(x){
  if(x %in% c("X", "X|X")){
    return(23)
  }else{
    return(as.numeric(strsplit(x, "\\|")[[1]][1]))
  }
})]
## one protein missing: KIR2DL2, needs to be placed manually
lab.prot[, chr_numeric := ifelse(is.na(chr_numeric), 19, chr_numeric)]
lab.prot[, start_unique := sapply(start_position, function(x) as.numeric(strsplit(x, "\\|")[[1]][1]))]
## edit two entries manually
lab.prot$start_unique[ which(lab.prot$Assay == "KIR2DL2")] <- 106543
lab.prot$start_unique[ which(lab.prot$Assay == "PRSS2")]   <- 142760398
## now sort based on unqiue position in the genome
lab.prot  <- lab.prot[ order(chr_numeric, start_unique)]
lab.prot[, p.srt := 1:nrow(lab.prot)]

#-----------------------------------#
##--           Figure            --##
#-----------------------------------#

## add colours
phe.col        <- fread("<path>", sep="\t", header=T)
## add identifier to match with data set
phe.col[, id := paste0("bin_", phecode)]
phe.col        <- phe.col[ id %in% tmp.phe$short_name]
phe.col        <- phe.col[ order(phecode)]
phe.col[, srt := 1:nrow(phe.col)]

## add phecode category to results --> reduces results to less redundant ones as well
res.surv       <- merge(res.surv, phe.col[, .(phecode, category)])

## --> association count for each model <-- ##

## proteins: N.B. This drops all proteins with no association with any of the phecodes for any of the variables: 4.039895e-08
assoc.prot     <- res.surv[pval < .05/(nrow(lab.prot)*nrow(tmp.phe)), .(n.assoc = length(phecode), n.categories = length(unique(category))), by = c("id", "type.analysis")]
## add to the label
assoc.prot     <- dcast(assoc.prot, id ~ type.analysis, value.var = c("n.assoc", "n.categories"), sep=".")
lab.prot       <- merge(lab.prot, assoc.prot, by = "id", all.x=T)

## diseases
assoc.phe      <- res.surv[, .(n.assoc = sum(pval < .05/(nrow(lab.prot)*nrow(tmp.phe)))), by = c("phecode", "type.analysis")]
tmp.phe        <- merge(tmp.phe, dcast(assoc.phe, phecode ~ type.analysis, value.var = "n.assoc"), by = "phecode")

## --> prepare plotting <-- ##

## prepare plotting boundaries
chr.dat        <- do.call(data.frame, aggregate(p.srt ~ chr_numeric, lab.prot, function(x) c(min(x), mean(x), max(x))))
## edit names
names(chr.dat) <- c("chr_numeric", "min", "mid", "max")
## same for diseases
phe.dat        <- do.call(data.frame, aggregate(srt ~ category, phe.col, function(x) c(min(x), mean(x), max(x))))
## edit names
names(phe.dat) <- c("category", "min", "mid", "max")
## add colours
phe.dat        <- merge(phe.dat, unique(phe.col[, c("category", "cl")]))
phe.dat        <- phe.dat[ order(phe.dat$min),]

## start plotting device
pdf("../graphics/Summary.Olink.phecode.incident.diseases.20250701.pdf", width = 6.3, height = 3.15)
## plotting parameters
par(mar=c(.1,1.5,1.5,.1), mgp=c(.6,0,0), cex.axis=.4, cex.lab=.4, tck=-.01, lwd=.5, xaxs="i", yaxs="i", bty="n")
## layout
layout(matrix(1:4, 2, 2), heights = c(.2,.8), widths = c(.9, .1))

#---------------------------------------#
##-- Fraction of proteins persistent --##
#---------------------------------------#

## empty plot
plot(c(.5, nrow(lab.prot)+.5), c(0, 105), type="n", xlab="", ylab="Fraction persistent [%]",
     xaxt="n", yaxt="n")
## add axis
axis(2, lwd=.5)
## rectanlges to divide chromosomes
pm <- par("usr")

## add regional estimates
tmp <- lab.prot[ n.assoc.min > 4]
arrows(tmp$p.srt, 0, tmp$p.srt, (tmp$n.assoc.adj/tmp$n.assoc.min)*100,  length = 0, lwd=.4, col="grey20")

                                                ## add box
rect(pm[1], pm[3], pm[2], pm[4], lwd=.5, xpd=NA)

## add least affected
tmp <- lab.prot[ n.assoc.min*.5 <= n.assoc.adj & n.assoc.min > 4]
## add labels
text(tmp$p.srt, pm[4]+(pm[4]-pm[3])*.05, xpd=NA, labels = tmp$Assay, cex=.4, pos=4, srt=90, offset = 0)

#----------------------#
##--     2D plot    --##
#----------------------#

## adapt plotting parameters
par(mar=c(3.5,1.5,.1,.1))

## empty plot
plot(c(.5, nrow(lab.prot)+.5), c(.5,nrow(tmp.phe)+.5), xaxt="n", yaxt="n",
     xlab="Chromosomal position", ylab="", type="n", ylim=rev(c(.5,nrow(tmp.phe)+.5)))
## add axis
axis(1, lwd=.5, at=chr.dat$mid, labels=1:23)
## get plotting coordinates
pm <- par("usr")
## add phenotype categories
rect(pm[1], phe.dat$min, pm[2], phe.dat$max, col=colorspace::lighten(phe.dat$cl, .6), border=phe.dat$cl, lwd=.1)
abline(v=chr.dat$min, lwd=.3, lty=2, col="white")

## add all points
tmp <- res.surv[ pval < .05/(nrow(lab.prot)*nrow(tmp.phe)) & type.analysis == "min"]
## add plotting position
tmp <- merge(tmp, lab.prot[, .(id, p.srt)], by = "id")
tmp <- merge(tmp, tmp.phe[, .(phecode, srt)], by = "phecode")
## add
points(tmp$p.srt, tmp$srt, cex=.2, lwd=.2, col="grey90", pch=20)

## highlight protein associations with persistent significance
tmp <- res.surv[ pval < .05/(nrow(lab.prot)*nrow(tmp.phe)) & type.analysis == "adj"]
## add plotting position
tmp <- merge(tmp, lab.prot[, .(id, p.srt)], by = "id")
tmp <- merge(tmp, tmp.phe[, .(phecode, srt)], by = "phecode")
## add
points(tmp$p.srt, tmp$srt, cex=.3, lwd=.2, col="black", pch=20)

## add box
rect(pm[1], pm[3], pm[2], pm[4], lwd=.5, xpd=NA)

## add legend
legend(pm[1]-(pm[2]-pm[1])*.02, pm[3]-(pm[4]-pm[3])*.1, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=phe.dat$cl, legend = stringr::str_to_sentence(phe.dat$category),
       cex=.45, ncol=12, pt.cex=1)

#----------------------#
##--   empty plot   --##
#----------------------#

## empty plot
plot(1,1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
## add legend
legend(.1, .9, bty="n", lty=0, pch=22, pt.bg = c("grey90", "black"),
       pt.lwd = .4, legend = c("Age & sex", 
                               "Fully adjusted"), 
       cex=.5, title = paste0("p<", sprintf("%.1e", .05/(nrow(lab.prot)*nrow(tmp.phe)))),
       xpd=NA)

#----------------------#
##--   extend. adj  --##
#----------------------#

## adapt plotting parameters
par(mar=c(3.5,.1,.1,.5))

## at least some associations
tmp <- tmp.phe[ min > 4]

## empty plot
plot(c(0, max(tmp$adj/tmp$min, na.rm=T)*100+2), c(.5, nrow(tmp.phe)+.5), type="n", ylab="", xlab="Fraction persistent [%]",
     yaxt="n", xaxt="n", ylim=rev(c(.5,nrow(tmp.phe)+.5)))
## add axis
axis(1, lwd=.5)
pm <- par("usr")
## add phenotype categories
rect(pm[1], phe.dat$min, pm[2], phe.dat$max, col=colorspace::lighten(phe.dat$cl, .6), border=phe.dat$cl, lwd=.1)

## add regional estimates
arrows(0, tmp$srt, (tmp$adj/tmp$min)*100, tmp$srt, length = 0, lwd=.4, col="grey20")
## add box
rect(pm[1], pm[3], pm[2], pm[4], lwd=.5, xpd=NA)

## close device
dev.off()

## start plotting device
pdf("../graphics/Comparison.Effects.Summary.Olink.phecode.incident.diseases.20250701.pdf", width = 6.3, height = 2.1)
## plotting parameters
par(mar=c(1.5,1.5,.5,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i", mfrow=c(1,3))

#-----------------------------#
##--       min vs all      --##
#-----------------------------#

## create data set
tmp <- dcast(res.surv[ type.analysis %in% c("min", "adj") ], id + phecode ~ type.analysis, value.var = c("beta", "se", "pval"),
             sep = ".")
## subset to what is of interest
tmp <- tmp[ pval.min < .05/(nrow(lab.prot)*nrow(tmp.phe))]
## simple plot
plot(beta.adj ~ beta.min, tmp, cex=.3, pch=20, 
     col=ifelse(tmp$pval.adj < .05/(nrow(lab.prot)*nrow(tmp.phe)), "red2", adjust_transparency("grey80", .3)),
     xlim=c(-.9,1.5), ylim=c(-.9,1.5), type="n",
     xlab=expression(Hazard~ratio[minimal]), ylab=expression(Hazard~ratio[extended]),
     xaxt="n", yaxt="n")
## axis - adjust to represent hazard rations
axis(1, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
axis(2, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
## add polygon for ±20% change
polygon(c(1.5, 1.5, -.9*.8,-.9*1.2), c(1.5/1.2, 1.5/0.8, -.9, -.9),col="grey90", border=NA)
## add line of identify
abline(a=0, b=1, lwd=.5); abline(h=0, lwd=.5, lty=2); abline(v=0, lwd=.5, lty=2)

## add all dots
points(tmp$beta.min, tmp$beta.adj, cex=.3, pch=20, col=adjust_transparency("grey50", .5))

## highlight persisting
tmp <- tmp[ pval.adj < .05/(nrow(lab.prot)*nrow(tmp.phe))]
points(tmp$beta.min, tmp$beta.adj, pch=20,
       cex=ifelse(abs(tmp$beta.adj/tmp$beta.min) < .8  , .3, .4),  
       col=ifelse(abs(tmp$beta.adj/tmp$beta.min) > .8  & sign(tmp$beta.min) == sign(tmp$beta.adj), 
                  "red2", adjust_transparency("orange", .5)))

## add legend
legend("topleft", cex=.5, lty=0, bty="n", pch=20,
       col=c("grey80", "orange", "red2"),
       legend=c("Minimal only", "Both", "Both & <20% effect attenuation"),
       title = paste0("Significance (p<", sprintf("%.1e", .05/(nrow(lab.prot)*nrow(tmp.phe))), "):"))

## annotate some staggering examples: 
tmp <- merge(tmp, phe.col[, .(phecode, phenotype)], by = "phecode")
jj  <- data.table(id = c("klk3", "sftpd", "pdcd1", "col9a1", "itgav", "mmp12", "proc"),
                  phecode = c(185.00, 502.00, 202.00, 740.11, 332.00, 442.11, 571.51))
## loop to assign labels
for(j in 1:nrow(jj)){
  text(tmp$beta.min[ which(tmp$id == jj$id[j] & tmp$phecode == jj$phecode[j])],
       tmp$beta.adj[ which(tmp$id == jj$id[j] & tmp$phecode == jj$phecode[j])],
       cex=.4,
       labels = paste0(toupper(jj$id[j]), " - ", tmp$phenotype[ which(tmp$phecode == jj$phecode[j] & tmp$id == jj$id[j])]),
       xpd=NA)
}

#-----------------------------#
##--       min vs cis      --##
#-----------------------------#

## p-value threshold for cis-scores
p.adj <- .05/(nrow(tmp.phe) * length(unique(res.surv[ type.analysis == "cis"]$id)))

## create data set
tmp <- dcast(res.surv[ type.analysis %in% c("min", "cis") ], id + phecode ~ type.analysis, value.var = c("beta", "se", "pval"),
             sep = ".")
## subset to what is of interest
tmp <- tmp[ pval.min < p.adj | pval.cis < p.adj]
## simple plot
plot(beta.cis ~ beta.min, tmp, cex=.3, pch=20, 
     col=ifelse(tmp$pval.cis < p.adj, "red2", adjust_transparency("grey80", .3)),
     xlim=c(-.9,1.5), ylim=c(-.9,1.5), type="n",
     xlab=expression(Hazard~ratio[minimal]), ylab=expression(Hazard~ratio["cis-pQTL-score"]),
     xaxt="n", yaxt="n")
## axis - cisust to represent hazard rations
axis(1, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
axis(2, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
## add polygon for ±20% change
polygon(c(1.5, 1.5, -.9*.8,-.9*1.2), c(1.5/1.2, 1.5/0.8, -.9, -.9),col="grey90", border=NA)
## add line of identify
abline(a=0, b=1, lwd=.5); abline(h=0, lwd=.5, lty=2); abline(v=0, lwd=.5, lty=2)

## add all dots
points(tmp$beta.min, tmp$beta.cis, cex=.3, pch=20, col=adjust_transparency("grey50", .5))

## highlight persisting
tmp <- tmp[ pval.cis < p.adj &  pval.min < p.adj]
points(tmp$beta.min, tmp$beta.cis, pch=20,
       cex=ifelse(abs(tmp$beta.cis/tmp$beta.min) < .8  , .3, .4),  
       col=ifelse(abs(tmp$beta.cis/tmp$beta.min) > .8 & abs(tmp$beta.cis/tmp$beta.min) < 1.2 & sign(tmp$beta.cis) == sign(tmp$beta.min),
                  "red2", adjust_transparency("orange", .5)))

## add legend
legend("topleft", cex=.5, lty=0, bty="n", pch=20,
       col=c("grey80", "orange", "red2"),
       legend=c("Minimal only", "Both", "Both & <20% effect attenuation"),
       title = paste0("Significance (p<", sprintf("%.1e", p.adj), "):"))
                                                
#-----------------------------#
##--       min vs trans    --##
#-----------------------------#

## p-value threshold for trans-scores
p.adj <- .05/(nrow(tmp.phe) * length(unique(res.surv[ type.analysis == "trans"]$id)))

## create data set
tmp <- dcast(res.surv[ type.analysis %in% c("min", "trans") ], id + phecode ~ type.analysis, value.var = c("beta", "se", "pval"),
             sep = ".")
## subset to what is of interest
tmp <- tmp[ pval.min < p.adj | pval.trans < p.adj]
## simple plot
plot(beta.trans ~ beta.min, tmp, cex=.3, pch=20, 
     col=ifelse(tmp$pval.trans < p.adj, "red2", adjust_transparency("grey80", .3)),
     xlim=c(-.9,1.5), ylim=c(-.9,1.5), type="n",
     xlab=expression(Hazard~ratio[minimal]), ylab=expression(Hazard~ratio["trans-pQTL-score"]),
     xaxt="n", yaxt="n")
## axis - transust to represent hazard rations
axis(1, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
axis(2, lwd=.5, at=log(c(.5,1,.75,1.5,2,3)), labels = c(.5,1,.75,1.5,2,3))
## add polygon for ±20% change
polygon(c(1.5, 1.5, -.9*.8,-.9*1.2), c(1.5/1.2, 1.5/0.8, -.9, -.9),col="grey90", border=NA)
## add line of identify
abline(a=0, b=1, lwd=.5); abline(h=0, lwd=.5, lty=2); abline(v=0, lwd=.5, lty=2)

## add all dots
points(tmp$beta.min, tmp$beta.trans, cex=.3, pch=20, col=adjust_transparency("grey50", .5))

## highlight persisting
tmp <- tmp[ pval.trans < p.adj &  pval.min < p.adj]
points(tmp$beta.min, tmp$beta.trans, pch=20,
       cex=ifelse(abs(tmp$beta.trans/tmp$beta.min) < .8  , .3, .4),  
       col=ifelse(abs(tmp$beta.trans/tmp$beta.min) > .8 & abs(tmp$beta.trans/tmp$beta.min) < 1.2 & sign(tmp$beta.trans) == sign(tmp$beta.min),
                  "red2", adjust_transparency("orange", .5)))

## add legend
legend("topleft", cex=.5, lty=0, bty="n", pch=20,
       col=c("grey80", "orange", "red2"),
       legend=c("Minimal only", "Both", "Both & <20% effect attenuation"),
       title = paste0("Significance (p<", sprintf("%.1e", p.adj), "):"))

## close device
dev.off()

#---------------------------------------#
##--           final table           --##
#---------------------------------------#

## get all the protein disease associations significantly associated in 
## at least one strata
tmp <- unique(res.surv[ (pval <  .05/(nrow(tmp.phe) * nrow(lab.prot)) & !(type.analysis %in% c("cis", "trans"))) | 
                          (pval <  .05/nrow(res.surv[ type.analysis == "cis"]) & type.analysis == "cis") | 
                          (pval <  .05/nrow(res.surv[ type.analysis == "trans"]) & type.analysis == "trans"), .(id, phecode)])
tmp <- merge(tmp, res.surv, allow.cartesian = T, by = c("id", "phecode"))

## create new HR colum
tmp[, hr := paste0(sprintf("%.2f", exp(beta)), " (", sprintf("%.2f", exp(beta - 1.96*se)), ";", sprintf("%.2f", exp(beta + 1.96*se)), ")")]
## expand what is needed
tmp <- merge(tmp, phe.col[, .(phecode, phenotype)], by = "phecode")
tmp <- dcast(tmp, id + phenotype + category ~ type.analysis, value.var = c("hr", "pval", "nevent"))
## write to file
fwrite(tmp, "Results.UKB.proteins.incident.disease.models.20250701.txt", sep = "\t", row.names = F, na = "")

##############################################
####          Enrichment analysis         ####
##############################################

## import function to do so
source("../scripts/enrich_participant_characteristics.R")

#-------------------------------#
##--       fasting study     --##
#-------------------------------#

## import results from fasting study
res.fasting    <- fread("<path>")
## fdr.aov gives the multiple testing correction!

## create time resolved results (do not run in parallele, since this is done within the function)
enr.fasting    <- lapply(c(1:7, 10), function(x){
  print(unique(tolower(res.fasting[ eval(as.name(paste0("pval.exposure", x))) < .05/nrow(res.fasting)]$Assay)))
  ## run the enrichment with the respective set of proteins
  tmp <- enrich.chrarac.prot(res.var.sex.long, lab.phe, unique(tolower(res.fasting[ eval(as.name(paste0("pval.exposure", x))) < .05/nrow(res.fasting)]$Assay)))
  ## report findings, but restrict to overall population
  return(data.table(day=x, tmp[ population == "All"]))
})
## combine
enr.fasting <  - rbindlist(enr.fasting)

#-------------------------------#
##--      ovarian cancer     --##
#-------------------------------#

## import results for plasma proteomics
ovarian.cancer <- as.data.table(read_excel("OvarianCancer_Paper.xlsx", 7))
## keep only what is needed
ovarian.cancer <- ovarian.cancer[, .(X, log2.foldchange., P_value, P_value_adjust)]
## create gene names that matches with Olink
ovarian.cancer[, protein_id := gsub("^.*_", "", X)]

## perform enrichment
enr.ovarian    <- enrich.chrarac.prot(res.var.sex.long[ population == "Female"], lab.phe, 
                                      prot.vec = tolower(ovarian.cancer[ P_value_adjust < .05]$protein_id), 
                                      prot.back = tolower(ovarian.cancer$protein_id))
## smoking as potential confounder;

## SPINT1 explained variance --> strongest predictor
pdf("../graphics/Example.SPINT1.pdf", width = 3.15, height = 3.15)
plot.expl.var(res.var.sex.long[ population == "Female"], "spint1",lab.phe)
dev.off()

#-------------------------------#
##--   CAD prediction model  --##
#-------------------------------#

## import proteins selected
cad.proteins <- as.data.table(read_excel("Protein_risk_model_deCODE_CAD.xlsx", 1))
## SomaLogic background
cad.back     <- as.data.table(read_excel("Protein_risk_model_deCODE_CAD.xlsx", 2))

## perform enrichment
enr.cad      <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, 
                                    prot.vec = tolower(cad.proteins$`Gene Name`), 
                                    prot.back = tolower(cad.back$Gene))
# https://jamanetwork.com/journals/jama/fullarticle/2808522

#-------------------------------#
##--       organ clocks      --##
#-------------------------------#

## import proteins selected for different organ systems
organ.clocks <- as.data.table(read_excel("organ_clock_proteins.xlsx"))
## create protein ID matching to Olink proteins
organ.clocks[, prot.id := tolower(gsub("\\..*", "", Somamer))]

## import SomaScan background (overall)
soma.back    <- as.data.table(read_excel("organ_clock_proteins.xlsx", 2))
soma.back[, prot.id := tolower(gsub("\\..*", "", Somamer))]
## subset to proteins actually used
soma.back    <- soma.back[ Stable_assay == T]

## import sorting by organ
soma.organ   <- as.data.table(read_excel("organ_clock_proteins.xlsx", 3))
## simplify
soma.organ   <- melt(soma.organ, measure.vars = names(soma.organ), na.rm = T)
## add protein identifier
soma.organ[, prot.id := tolower(gsub("\\..*", "", value))]
## careful: Some inconsistencies in the tables

## --> test organ background <-- ##

## perform enrichment by organ
enr.organ <- lapply(unique(organ.clocks$organ), function(x){
  ## run the enrichment
  tmp <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, 
                             prot.vec = soma.organ[ variable == x]$prot.id,
                             prot.back = soma.back$prot.id) 
  ## return results
  return(data.table(organ=x, tmp))
})
## combine results
enr.organ <- rbindlist(enr.organ)

## --> test w/o organ-specific background <-- ##

## perform enrichment by organ
enr.organ.wo <- lapply(unique(organ.clocks$organ), function(x){
  ## run the enrichment
  tmp <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, 
                             prot.vec = organ.clocks[ organ == x]$prot.id,
                             prot.back = soma.back$prot.id) 
  ## return results
  return(data.table(organ=x, tmp))
})
## combine results
enr.organ.wo <- rbindlist(enr.organ.wo)

## perform enrichment by organ
enr.organ.wo.v2 <- lapply(unique(organ.clocks$organ), function(x){
  ## run the enrichment
  tmp <- enrich.chrarac.prot_v2(res.var.sex.long[ population == "All"], lab.phe, 
                                prot.vec = organ.clocks[ organ == x]$prot.id,
                                prot.back = soma.back$prot.id) 
  ## return results
  return(data.table(organ=x, tmp))
})
## combine results
enr.organ.wo.v2 <- rbindlist(enr.organ.wo.v2, fill = T)

## --> test w/ organ-specific background <-- ##

## perform enrichment by organ
enr.organ.w  <- lapply(unique(organ.clocks$organ), function(x){
  ## run the enrichment
  print(x)
  tmp <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, 
                             prot.vec = organ.clocks[ organ == x]$prot.id,
                             prot.back = soma.organ[ variable == x]$prot.id) 
  ## return results
  return(data.table(organ=x, tmp))
})
## combine results
enr.organ.w  <- rbindlist(enr.organ.w, fill = T)

#---------------------------------#
##--         Figure(s)         --##
#---------------------------------#

## start device
pdf("../graphics/Enrichment.Fasting.Ovarian.cancer.CAD.20250619.pdf", width = 3.15, height = 2.5)
## graphical parameters
par(mar=c(6,1.5,3,.5), tck=-.01, cex.axis=.5, cex.lab=.5, mgp=c(.6,0,0), bty="l", lwd=.5, yaxs="i", xaxs="i", mfrow=c(1,1))

## sort and reduce to characteristics enriched at least once
tmp          <- unique(c(enr.ovarian[ pval < .05/nrow(enr.ovarian) & or > 1]$label, 
                         enr.cad[ pval < .05/nrow(enr.cad) & or > 1]$label,
                         enr.fasting[ day == 1 & or > 1 & pval < .05/nrow(enr.fasting)]$label))
tmp          <- rbind(data.table(outcome="OC", enr.ovarian[ label %in% tmp]),
                      data.table(outcome="CAD", enr.cad[ label %in% tmp]),
                      data.table(outcome="Fasting", enr.fasting[ day == 1 & label %in% tmp]), fill = T)
tmp          <- tmp[ order(cat.srt, variable),]
tmp[, srt := 1:length(label), by="outcome"]

## define colour gradient for p-values
p.col        <- colorRampPalette(c("white", "#CB625F"))(90)

## vector to ease plotting
p.vec        <- 1:3
names(p.vec) <- c("OC", "CAD", "Fasting")

## empty plot
plot(c(.5,max(tmp$srt)+.5),  c(.5,3.5), ylab="", xlab="", type="n", xaxt="n", yaxt="n")

## devide fasting and refeeding
pm           <- par("usr")

## add dots
points(tmp$srt, p.vec[as.character(tmp$outcome)],  
       cex=ifelse(is.finite(tmp$or), log10(tmp$or), log10(100))*2,
       # cex=ifelse(tmp$pval < .05/nrow(enr.fasting), 
       #                              ifelse(tmp$or > 5, 2, tmp$or/4), .3), 
       pch=ifelse(tmp$pval < .05/min(c(nrow(enr.cad), nrow(enr.ovarian))), 21, NA),
       bg=p.col[ceiling(-log10(tmp$pval)*10)], xpd=NA)

## add names
text(1:max(tmp$srt), 0, cex=.4, labels=tmp[ outcome == "OC" ]$label, pos=2, xpd=NA, srt=60, offset = 0)

## add axis
axis(2, lwd=.5, at=p.vec, labels = names(p.vec))

## p-value gradient
l <- seq(pm[1]+(pm[2]-pm[1])*.15, pm[1]+(pm[2]-pm[1])*.45, length.out = length(p.col))
## rectangle for the colours
rect(l-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.2, l+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.3, border=NA, col=p.col, xpd=NA)
## box
rect(l[1]-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.2, l[length(l)]+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.3, border="black", col=NA, lwd=.3, xpd=NA)
## add header
text(pm[1]+(pm[2]-pm[1])*.05, pm[4]+(pm[4]-pm[3])*.45, cex=.4, labels = expression(-log[10]("p-value")), pos=4,
     offset = .2, xpd=NA)
## simple axis
text(l[round(c(1, c(.2, .4, .6, .8, 1)*length(p.col)))], pm[4]+(pm[4]-pm[3])*.19, 
     labels=sprintf("%.1f", c(0, .2, .4, .6, .8, 1)*length(p.col)/10), pos=1, cex=.3, offset = .1, xpd=NA)

## odds ratio
legend(pm[1]+(pm[2]-pm[1])*.5, pm[4]+(pm[4]-pm[3])*.8, bty="n", cex=.5, xpd=NA, ncol=4, 
       lty=0, pt.cex=log10(c(1, 5, 10, 50))*2, pch=21,
       legend = c(1, 5, 10, 50), title = "Odds ratio",
       x.intersp = 1.5)

## close device
dev.off()

## start device
pdf("../graphics/Enrichment.Organ.clocks.20250619.pdf", width = 3.15, height = 5)
## graphical parameters
par(mar=c(3,7,3,.5), tck=-.01, cex.axis=.5, cex.lab=.5, mgp=c(.6,0,0), bty="l", lwd=.5, yaxs="i", xaxs="i")

## sort and reduce to characteristics enriched at least once
tmp          <- enr.organ.wo[ label %in% enr.organ.wo[ pval < .05/nrow(enr.organ.wo) & or > 1]$label & organ %in% enr.organ.wo[ pval < .05/nrow(enr.organ.wo) & or > 1]$organ]
tmp          <- tmp[ order(cat.srt, variable),]
tmp[, srt := 1:length(label), by="organ"]

## define colour gradient for p-values
p.col        <- colorRampPalette(c("white", "#CB625F"))(130)

## vector to ease plotting
p.vec        <- 1:length(unique(tmp$organ))
names(p.vec) <- unique(tmp$organ)

## empty plot
plot(c(.5,length(p.vec)+.5), c(.5,max(tmp$srt)+.5), ylab="", xlab="", type="n", xaxt="n", yaxt="n")

## devide fasting and refeeding
pm           <- par("usr")

## add dots
points(p.vec[as.character(tmp$organ)], tmp$srt,   
       cex=ifelse(is.finite(tmp$or), log10(tmp$or), log10(100))*1.5,
       # cex=ifelse(tmp$pval < .05/nrow(enr.fasting), 
       #                              ifelse(tmp$or > 5, 2, tmp$or/4), .3), 
       pch=ifelse(tmp$pval < .05/.05/nrow(enr.organ.wo), 21, NA),
       bg=p.col[ceiling(-log10(tmp$pval)*10)], xpd=NA)

## add names
text(.5, 1:max(tmp$srt), cex=.4, labels=tmp[ organ == "Liver" ]$label, pos=2, xpd=NA, offset = .1)

## add axis
axis(1, lwd=.5, at=p.vec, labels = names(p.vec), las=2)

## p-value gradient
l <- seq(pm[1]-(pm[2]-pm[1])*.4, pm[1]+(pm[2]-pm[1])*.2, length.out = length(p.col))
## rectangle for the colours
rect(l-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.05, l+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.07, border=NA, col=p.col, xpd=NA)
## box
rect(l[1]-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.05, l[length(l)]+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.07, border="black", col=NA, lwd=.3, xpd=NA)
## add header
text(pm[1]-(pm[2]-pm[1])*.5, pm[4]+(pm[4]-pm[3])*.08, cex=.4, labels = expression(-log[10]("p-value")), pos=4,
     offset = .2, xpd=NA)
## simple axis
text(l[round(c(1, c(.2, .4, .6, .8, 1)*length(p.col)))], pm[4]+(pm[4]-pm[3])*.04, 
     labels=sprintf("%.1f", c(0, .2, .4, .6, .8, 1)*length(p.col)/10), pos=1, cex=.3, offset = .1, xpd=NA)

## odds ratio
legend(pm[1]+(pm[2]-pm[1])*.3, pm[4]+(pm[4]-pm[3])*.1, bty="n", cex=.5, xpd=NA, ncol=4, 
       lty=0, pt.cex=log10(c(1, 5, 10, 20))*1.5, pch=21,
       legend = c(1, 5, 10, 20), title = "Odds ratio",
       x.intersp = 1.5)

## close device
dev.off()

#-------------------------------#
##--    HPA atlas results    --##
#-------------------------------#

## import results
hpa.disease     <- fread("blood_pea_disease_de.tsv")

## define protein background and align with format of UKB proteins
prot.hpa        <- tolower(unique(hpa.disease$Gene))
## n = 1162

## count how many sig findings
prot.hpa.sig.1  <- hpa.disease[, .(num.sig = sum(`p-value adjusted` < .05)), 
                               by=c("Disease", "Control")]
## dcast
prot.hpa.sig.1  <- dcast(prot.hpa.sig.1, Disease ~ Control, value.var = "num.sig")

## flag unspecific findings
prot.hpa.sig.2  <- hpa.disease[, .(num.sig = sum(`p-value adjusted` < .05)), 
                               by=c("Gene", "Control")]
## dcast
prot.hpa.sig.2  <- dcast(prot.hpa.sig.2, Gene ~ Control, value.var = "num.sig")

## look into pleiotropic proteins first
enr.hpa.prot    <- lapply(seq(5,50,5), function(x){
  ## apply for different thresholds
  tmp <- enrich.chrarac.prot(res.var.sex.long, lab.phe, 
                             prot.vec = tolower(prot.hpa.sig.2[ Healthy >= x]$Gene), 
                             prot.back = prot.hpa)
  ## return results
  return(data.table(num.prot=x, tmp[ population == "All"]))
})
## combine again
enr.hpa.prot     <- rbindlist(enr.hpa.prot)

## --> no restriction of proteins <-- ##

## create time resolved results (do not run in parallel, since this is done within the function)
enr.hpa.disease <- unique(hpa.disease[, .(Disease, Control, Class)])
enr.hpa.disease <- lapply(1:nrow(enr.hpa.disease), function(x){
  ## progress
  cat("Testing differential proteins for", enr.hpa.disease$Disease[x], "against", enr.hpa.disease$Control[x],"\n")
  ## get set of proteins
  prots <- tolower(unique(hpa.disease[ Disease == enr.hpa.disease$Disease[x] & Control == enr.hpa.disease$Control[x] & `p-value adjusted` < .05 & abs(logFC) > .5]$Gene))
  ## which proteins
  print(prots)
  ## run the enrichment with the respective set of proteins
  tmp <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, prot.vec = prots, prot.back = prot.hpa)
  cat("\n-----------------------------------\n")
  ## report findings, but restrict to overall population
  return(data.table(enr.hpa.disease[x,], tmp[ population == "All"]))
})
## combine
enr.hpa.disease <- rbindlist(enr.hpa.disease)

#-------------------------------#
##--    compare vs healthy   --##
#-------------------------------#

## sort and reduce to characteristics enriched at least once
tmp          <- enr.hpa.disease[ label %in% enr.hpa.disease[ pval < .05/(nrow(enr.hpa.disease)/3) & or > 1 & Control == "Healthy"]$label &
                                   Control == "Healthy" & Disease %in% enr.hpa.disease[ pval < .05/(nrow(enr.hpa.disease)/3) & or > 1 & Control == "Healthy"]$Disease]
tmp          <- tmp[ order(Class, cat.srt, variable),]
tmp[, srt := 1:length(label), by="Disease"]

## define colour gradient for p-values
p.col        <- colorRampPalette(c("white", "#CB625F"))(32)

## vector to ease plotting
p.vec        <- 1:length(unique(tmp$Disease))
names(p.vec) <- unique(tmp$Disease)

## start device
pdf("../graphics/Enrichment.HPA.Disease.Healthy.20250619.pdf", width = 6.3, height = 6)
## graphical parameters
par(mar=c(7,7,2,.5), tck=-.01, cex.axis=.5, cex.lab=.5, mgp=c(.6,0,0), bty="l", lwd=.5, yaxs="i", xaxs="i")

## empty plot
plot(c(.5,max(tmp$srt)+.5),  c(.5,length(p.vec)+.5), ylab="", xlab="", type="n", xaxt="n", yaxt="n")

## devide fasting and refeeding
pm           <- par("usr")
rect(pm[1], 1:length(p.vec)-.5, pm[2], 1:length(p.vec)+.5, col=c("grey90", "white"), border = NA)
## separate clumns as well
abline(v=1:max(tmp$srt)+.5, lwd=.5, lty=2, col="white")

## add dots
points(tmp$srt, p.vec[as.character(tmp$Disease)],  
       cex=ifelse(is.finite(tmp$or), log10(tmp$or), log10(100)),
       pch=ifelse(tmp$pval < .05/nrow(enr.fasting), 21, NA),
       bg=p.col[ceiling(-log10(tmp$pval))], xpd=NA)

## add names
text(1:max(tmp$srt), 0, cex=.4, labels=tmp[ Disease == tmp$Disease[1] ]$label, pos=2, xpd=NA, srt=60, offset = 0)

## add axis
text(pm[1], 1:length(p.vec), cex=.4, labels = names(p.vec), pos=2, xpd=NA)

## p-value gradient
l <- seq(pm[1]+(pm[2]-pm[1])*.15, pm[1]+(pm[2]-pm[1])*.45, length.out = length(p.col))
## rectangle for the colours
rect(l-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.05, l+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.08, border=NA, col=p.col, xpd=NA)
## box
rect(l[1]-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.05, l[length(l)]+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.08, border="black", col=NA, lwd=.3, xpd=NA)
## add header
text(pm[1]+(pm[2]-pm[1])*.02, pm[4]+(pm[4]-pm[3])*.07, cex=.4, labels = expression(-log[10]("p-value")), pos=4,
     offset = .2, xpd=NA)
## simple axis
text(l[round(c(1, c(.2, .4, .6, .8, 1)*length(p.col)))], pm[4]+(pm[4]-pm[3])*.04, 
     labels=sprintf("%.1f", c(0, .2, .4, .6, .8, 1)*length(p.col)), pos=1, cex=.3, offset = .1, xpd=NA)

## odds ratio
legend(pm[1]+(pm[2]-pm[1])*.5, pm[4]+(pm[4]-pm[3])*.1, bty="n", cex=.5, xpd=NA, ncol=5, 
       lty=0, pt.cex=log10(c(1, 5, 10, 20, 50)), pch=21,
       legend = c(1, 5, 10, 20, 50), title = "Odds ratio",
       x.intersp = 1.5)

## close device
dev.off()

##########################################################################################################
##########################################################################################################
####                                    REVISION - MSB 23/06/2025                                     ####
##########################################################################################################
##########################################################################################################


#################################################
####          comparison Fenland study       ####
#################################################

#---------------------------------#
##-- import relevant data sets --##
#---------------------------------#

## import results
fenland.soma          <- as.data.table(readxl::read_excel("Carrasco-Zanini_2024_NatMetab.xlsx", skip = 2, sheet = "Supp. Table 4"))
## split to have unique UniProt identifier
fenland.soma          <- rbindlist(lapply(1:nrow(fenland.soma), function(x) data.table(fenland.soma[x, ], 
                                                                                       uniprot.wide = strsplit(fenland.soma$UniProt[x], "\\|")[[1]])))
## n = 5048

## import correlation coefficients SomaScan v4 and Olink Explore, from EPIC or deCODE
cor.soma.olink.decode <- as.data.table(readxl::read_excel("SomaLogic_Olink_comparison_Eldjarn_deCODE_2023.xlsx"))

## make SeqId comparable
cor.soma.olink.decode[, SeqId := gsub("SeqId\\.", "", seqid)]

#---------------------------------#
##--   merge with UKB results  --##
#---------------------------------#

## add possible matching SomaScan entry
prot.comparison        <- prot.order
prot.comparison[, SeqId := sapply(prot.order$UniProt, function(x){
  ## split if needed
  x <- strsplit(x, "_")[[1]]
  ## search for matching entries in Fenland SomaScan data
  return(paste(sort(unique(fenland.soma[uniprot.wide %in% x]$SeqId)), collapse = "|"))
  
})]

## split by SomaScan ID
prot.comparison        <- rbindlist(lapply(1:nrow(prot.comparison), function(x) data.table(prot.comparison[x, ], 
                                                                                           SeqId.wide = strsplit(prot.comparison$SeqId[x], "\\|")[[1]])))
## n = 1,996 (drops non-mapping proteins)

## add explained variance
prot.comparison        <- merge(prot.comparison, unique(fenland.soma[, -c("uniprot.wide"), with=F]), 
                                by.x = "SeqId.wide", by.y = "SeqId", suffixes = c(".olink", ".soma"))

## overall correlation
cor.test(prot.comparison$`Total variance explained`, 1-prot.comparison$p.50)$p.value

#---------------------------------#
##--  factors poor correlation --##
#---------------------------------#

## --> cross-platform correlation coefficients <-- ##

## add correlation coefficient
prot.comparison        <- merge(prot.comparison, 
                                unique(cor.soma.olink.decode[!is.na(olink_smpnorm_corr), .(SeqId, olink_smpnorm_corr)]), 
                                all.x = T, by.x = "SeqId.wide", by.y = "SeqId")
## take care of redundancies (take highest correlation coefficient)
prot.comparison        <- prot.comparison[ order(Assay, SeqId.wide, -olink_smpnorm_corr)] 
prot.comparison[, ind := 1:.N, by = c("SeqId.wide", "Assay")]
prot.comparison        <- prot.comparison[ ind == 1]

## --> UKB categories  <-- ##

## add category: Residuals do not have a category!
res.var.sex.long       <- merge(res.var.sex.long, lab.phe[, .(short_name, category)], by.x = "variable", by.y = "short_name", all.x = T) 

## add variance explained by UKB categories: fill missing values with zeros; drops Olink proteins w/o explained variance
prot.comparison        <- merge(prot.comparison, 
                                dcast(res.var.sex.long[ population == "All" & variable != "Residuals", .(explained.category = sum(p.50)),
                                                  by = c("var", "category")], var ~ category, sep = ".", fill = 0),
                                by.x = "id", by.y = "var")

## compute difference in explained variance
prot.comparison[, diff.expl.var.olink.soma := ((1-p.50)*100)-`Total variance explained`]

## edit names to remove white spaces (needs to be considered when remapping later)
names(prot.comparison)  <- gsub("[^A-Za-z0-9._]", "_", names(prot.comparison))

## create fileds for fenland
prot.comparison[, Genetic.fenland := apply(prot.comparison[, .(cis_pQTL_score, trans_pQTL_score)], 1, sum, na.rm=T)]
prot.comparison[, Basic_demographics.fenland := apply(prot.comparison[, .(Age, Sex, Alcohol_intake, current_smoker, ever_smoker)], 1, sum, na.rm=T)]
prot.comparison[, Biomarker.fenland := apply(prot.comparison[, .(eGFR, HDL, LDL, Triglycerides, Total_cholesterol, Vitamin_C, HbA1c, Fasting_glucose, Insulin, CRP, ALP, ALT, FT4, TSH, FT3, Bilirubin)], 1, sum, na.rm=T)]
prot.comparison[, Body_composition.fenland := apply(prot.comparison[, .(Total_fat_mass, Total_lean_mass, Arm_fat_mass, Legs_fat_mass, BMI, WHR, Peripheral_fat, Subcutaneous_fat, Visceral_fat)], 1, sum, na.rm=T)]
prot.comparison[, Bone.fenland := apply(prot.comparison[, .(Total_bone_mass)], 1, sum, na.rm=T)]
prot.comparison[, Cardiovascular.fenland := apply(prot.comparison[, .(DBP, SBP)], 1, sum, na.rm=T)]
prot.comparison[, Diet.fenland := apply(prot.comparison[, .(MDS, DASH_diet_score)], 1, sum, na.rm=T)]
prot.comparison[, Technical.fenland := apply(prot.comparison[, .(Appointment_Date, Appointment_Date__spring_, Appointment_Date__summer_, Appointment_Date__winter_, Proteomic_PC_1__0.005__dilution_proteins_,
                                                                 Proteomic_PC_1__20__dilution_proteins_, Proteomic_PC_1__0.5__dilution_proteins_)], 1, sum, na.rm=T)]
prot.comparison[, Genetic.scores.fenland := apply(prot.comparison[, .(FG_GRS, FI_GRS, BMI_GRS, WHR_GRS, Hip_specific_GRS, Waist_specific_GRS, T2D_GRS,
                                                                      eGFR_GRS, DBP_GRS, SBP_GRS, CAD_GRS)], 1, sum, na.rm=T)]

## write to file
write.table(prot.comparison, "Explained.variance.comparison.UKB.olink.Fenland.soma.20250703.txt", sep = "\t", row.names = F)


## test whether differences in explained variance might be explained by specific fields in UKB
res.diff.var.olink.soma <- lapply(c(gsub("[^A-Za-z0-9._-]", "_", unique(lab.phe$category)), "Genetic.fenland", "Basic_demographics.fenland",
                                         "Biomarker.fenland", "Body_composition.fenland", "Bone.fenland", "Cardiovascular.fenland",
                                         "Diet.fenland", "Technical.fenland"), function(x){
  ## compute regression model, add factor to make estimates comparable
  if(length(grep("fenland", x) > 0)){
    ff <- summary(lm(paste0("diff.expl.var.olink.soma ~ I(",  x, ") + olink_smpnorm_corr"), prot.comparison))$coefficients
  }else{
    ff <- summary(lm(paste0("diff.expl.var.olink.soma ~ I(",  x, "*100) + olink_smpnorm_corr"), prot.comparison))$coefficients
  }
  ## return results
  return(data.table(expl.factor = x, ff[2,,drop=F]))
})
## combine
res.diff.var.olink.soma <- rbindlist(res.diff.var.olink.soma)
## write to file
write.table(res.diff.var.olink.soma, "Results.variance.comparison.UKB.olink.Fenland.soma.20250623.txt", sep = "\t", row.names = F)

#---------------------------------#
##-- enrichment prediction work--##
#---------------------------------#

## import relevant table (proteins selected)
ukb.prediction         <- data.table(read_excel("UKB_olink_prediction.xlsx", 8, skip=3))
## keep only results relevant to entire platform
ukb.prediction         <- ukb.prediction[, 1:12]
## edit names
names(ukb.prediction ) <- c("Disease", "Specialty", "c.index", "cil", "ciu", "delta.ci", "delta.ci.cil", "delta.ci.ciu", "n.incident", "n.proteins", "sig.model", "proteins")

## run enrichment
enr.prediction         <- lapply(1:nrow(ukb.prediction), function(x){
  ## get the proteins needed for enrichment
  prots <- tolower(strsplit(ukb.prediction$proteins[x], " , ")[[1]])
  ## get what is needed
  prots <- gsub("_.*$", "", prots)
  print(prots)
  ## perform enrichment
  prots <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe,
                               prot.vec = prots)
  ## return
  return(data.table(ukb.prediction[x, .(Disease, Specialty, proteins)], prots))
  
})
## combine
enr.prediction         <- rbindlist(enr.prediction)

## add model estimates
enr.prediction         <- merge(enr.prediction, ukb.prediction)

## export results
write.table(enr.prediction[ pval < .05/nrow(enr.prediction[ Disease %in% ukb.prediction[ sig.model == "yes"]$Disease]) & Disease %in% ukb.prediction[ sig.model == "yes"]$Disease], "Results.enrichment.sparse.protein.signatures.20250703.txt", sep = "\t", row.names = F)

## --> how many instances, in which we observe enrichment for significant models <-- ##

## compute only for models with improvement: p-value = 3.346496e-06
View(enr.prediction[ or > 1 & pval < .05/nrow(enr.prediction[ Disease %in% ukb.prediction[ sig.model == "yes"]$Disease]) & Disease %in% ukb.prediction[ sig.model == "yes"]$Disease])
View(enr.prediction[ or < 1 & pval < .05/nrow(enr.prediction[ Disease %in% ukb.prediction[ sig.model == "yes"]$Disease]) & Disease %in% ukb.prediction[ sig.model == "yes"]$Disease])

## how many disease
length(unique(enr.prediction[ or > 1 & pval < .05/nrow(enr.prediction[ Disease %in% ukb.prediction[ sig.model == "yes"]$Disease]) & Disease %in% ukb.prediction[ sig.model == "yes"]$Disease]$Disease))

## --> how many instances, in which we observe enrichment for significant models <-- ##

## look at pleiotropic proteins
pleio.prediction       <- table(unlist(lapply(ukb.prediction$proteins, function(x){
  ## get the proteins needed for enrichment
  prots <- tolower(strsplit(x, " , ")[[1]])
  ## get what is needed
  prots <- gsub("_.*$", "", prots)
  return(prots)
})))
## proteins selected 5 times or more
pleio.prediction       <- names(pleio.prediction[ pleio.prediction >= 5])
pleio.prediction       <- enrich.chrarac.prot(res.var.sex.long[ population == "All"], lab.phe, prot.vec = pleio.prediction )

#---------------------------------#
##--      subsampling EUR      --##
#---------------------------------#

## create file for jub submission
write.table(expand.grid(olink.id=lab.prot$id,
                        population = c("AFR", "CSA"), 
                        fold = 1:4), 
            "EUR.subsample.ancestry.txt", sep = "\t", row.names = F, col.names = F, quote = F)

## --> import results <-- ##

## import
jj               <- dir("../output_updated/")
## pick out the subsampling sets
jj               <- grep("\\.\\d\\.", jj, value = T)
jj               <- grep("lasso", jj, value = T)
## one is missing

## find what is missing
tmp.sampling     <- expand.grid(olink.id=lab.prot$id,
                                population = c("AFR", "CSA"), 
                                fold = 1:4, stringsAsFactors = F)
tmp.sampling     <- as.data.table(tmp.sampling)
tmp.sampling[, missing := paste("lasso.explained.variance", population, fold, olink.id, "txt", sep = ".") %in% jj]
tmp.sampling[ missing == F]
## rerun manually: slmap - AFR - 4

## --> import <-- ##

## run import in parallel
registerDoMC(10)

## import
res.var.sub      <- mclapply(jj, function(x){
  ## import 
  tmp <- fread(paste0("../output_updated/", x))
  ## add population
  tmp[, population := strsplit(x, "\\.")[[1]][4]]
  ## subsample
  tmp[, fold := strsplit(x, "\\.")[[1]][5]]
  ## edit names
  names(tmp) <- sapply(names(tmp), function(x) ifelse(x %in% lab.ext$short_name_new, lab.ext$short_name[ which(lab.ext$short_name_new == x)], x))
  ## account for genetic scores
  names(tmp) <- gsub(paste0("\\.", tmp$var[1]), "", names(tmp))
  return(tmp)
}, mc.cores = 10)
## combine
res.var.sub      <- rbindlist(res.var.sub, fill = T)
## look for possibly poorly run estimates
summary(res.var.sub$n.boot)

## replace missing summary measure with p.50!
res.var.sub[, summary.measure := ifelse(is.na(summary.measure), "p.50", summary.measure)]

## transform in more efficient shape
res.var.sub.long <- melt.data.table(res.var.sub, id.vars = c("var", "summary.measure", "n.boot", "population", "fold", "n"))
## reshape again by summary measure
res.var.sub.long <- dcast(res.var.sub.long, var + variable + population + n + n.boot + fold  ~ summary.measure, sep = ".", value.var = "value")
## drop missing values
res.var.sub.long <- res.var.sub.long[ !is.na(p.50) ]
## avoid coding as a factor
res.var.sub.long[, variable := as.character(variable)]

## --> compute summary to compare with ancestral results <-- ##

## create summary by ancestry: median across the four folds
tmp              <- dcast(res.var.sub.long[ variable == "Residuals", .(p.50 = median(p.50, na.rm = T)), by = c("var", "population")], var ~ population, sep=".", value.var = "p.50")
## average by fold
names(tmp)       <- c("var", "p.50.EUR.sub.AFR", "p.50.EUR.sub.CSA")
prot.ancestry    <- merge(prot.ancestry, tmp, all.x = T)
