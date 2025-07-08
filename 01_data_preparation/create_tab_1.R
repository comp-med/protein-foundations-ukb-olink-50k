#################################
## function to compute simple
## table 1 characteristics

create.table1 <- function(dat, feat){
  
  ## 'dat'  -- data set containing all needed data
  ## 'feat' -- which feature to display
  
  ## create object to store results
  tab1 <- data.frame(var=feat)
  
  ## convert to data frame
  dat  <- as.data.frame(dat)
  
  ## now compute either mean and sd or %-fraction of cases
  tab1[,2:3] <- t(sapply(tab1$var, function(x){
    print(x)
    if(length(unique(dat[,x]))>5 & !(is.factor(dat[,x]))){
      ## compute mean and sd
      mw  <- mean(dat[,x], na.rm=T)
      sa  <- sd(dat[,x], na.rm=T)
      n   <- sum(!is.na(dat[,x]))
      ## now create a column
      dig <- ifelse(mw <= 5, 2, ifelse(mw <= 200, 1, 0))
      cl  <- paste(round(mw,dig), " (", round(sa,dig), ")", sep="")
      return(c(n, cl))
    }else{
      n   <- sum(!is.na(dat[,x]))
      # ne  <- sum(dat[,x] == 1, na.rm=T)
      # return(c(n, round(ne/n*100,1)))
      ne <- table(dat[,x])
      ne <- sapply(ne, function(x) paste0(x, " (", sprintf("%.1f", x/n*100), "%)"))
      print(ne)
      ne <- paste(paste(names(ne), ne), collapse = "|")
      # print(paste(names(x), x, collapse = " "))
      return(c(n, ne))
    }
  }))
  ## adapt names
  names(tab1)[-1] <- c("N", "statistic") 
  return(tab1)
}