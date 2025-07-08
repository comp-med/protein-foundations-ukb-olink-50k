#################################################################
## function to test for enrichment of participant characterisitcs
## among a set of proteins given 

enrich.chrarac.prot <- function(res.var, lab, prot.vec, min.per=0, prot.back=NA){
  
  ## 'res.var'   -- results from variance
  ## 'lab'       -- label for phenotypes
  ## 'prot.vec'  -- vector of protein targets to test for enrichment for
  ##                N.B. needs to map to protein names in 'res.var' 
  ## 'min.per'   -- require a minimum of explained variance by certain factors to be considered
  ##                as a fraction between 0 and 1 (1 = 100%)
  ## 'prot.back' -- optional lsit of proteins to be used as background
  
  ## drop associations not passing certain explained variance threshold
  res.var <- res.var[ p.50 > min.per & p.025 > 0]
  
  ## subset to background if required
  if(!is.na(prot.back[1])){
    res.var <- res.var[ var %in% prot.back]
    cat(length(unique(res.var$var)), "proteins remained after applying background filltering\n")
  }
  
  ## subset protein input to those in the data set
  prot.vec <- prot.vec[ prot.vec %in% res.var$var] 
  
  ## report how many different
  cat(length(prot.vec), "differentially expressed proteins remained after applying filltering\n")
  
  ## define variables to test for across different populations
  vars     <- res.var[ variable != "Residuals", .(variable.freq=length(var)), by=c("population", "variable")]
  
  # print(vars)
  
  ## reduce to at least five variables
  vars     <- vars[ variable.freq >= 5]
  
  ## proceed only if any
  if(nrow(vars) > 0){
    ## perform enrichment by variable and population
    registerDoMC(30)
    ## run in parallel
    res      <- mclapply(1:nrow(vars), function(x){
      
      ## variable of interest
      var.id    <- vars$variable[x]
      ## population of interest
      pop.id    <- vars$population[x]
      
      ## define protein back ground for population
      prot.back <- unique(res.var[ variable != "Residuals" & population == pop.id ]$var)
      
      ## protein of interest and selected for the variable of interest
      d1        <- nrow(res.var[ variable == var.id & population == pop.id & var %in% prot.vec])
      ## protein not of interest and selected for the variable of interest
      d2        <- nrow(res.var[ variable == var.id & population == pop.id & !(var %in% prot.vec)])
      ## protein of interest but not selected
      d3        <- length(prot.vec) - d1
      ## protein not of interest and not selected
      d4        <- length(prot.back[ !(prot.back %in% prot.vec | prot.back %in% res.var[ variable == var.id & population == pop.id]$var)]) 
      
      ## test for enrichment
      enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
      
      ## return information needed
      return(data.table(population=pop.id, variable=var.id, or=enr$estimate, pval=enr$p.value, 
                        intersection=paste(res.var[ variable == var.id & population == pop.id & var %in% prot.vec]$var, collapse = "|"), 
                        d1=d1, d2=d2, d3=d3, d4=d4))
      
    }, mc.cores = 30)
    ## return list
    res         <- rbindlist(res, fill = T)
    ## add label
    res         <- merge(res, lab, by.x = "variable", by.y = "short_name")
    ## return results
    return(res)
  }else{
    cat("Not enough proteins to test, please check input files\n")
  }
  


}