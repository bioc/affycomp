read.spikein <- function(filename){
#######################################################
###prep spike in exprSet
#######################################################
  s <- read.csv(filename,check.names=FALSE,row.names=1)
  samplenames <- colnames(s)
  ##remove the .cel if its there
  samplenames <- sub("\\.gz$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.Z$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.cel$","",samplenames,ignore.case=TRUE)
  colnames(s) <- samplenames
  ##read phenodata
  data(spikein.phenodata)
  ##putit in order
  s <- s[,rownames(pData(spikein.phenodata))]
  s <- new("exprSet",exprs=as.matrix(s),phenoData=spikein.phenodata)
  s <- exprset.log(s) ##take log
  return(s)
}

read.dilution <- function(filename){
#######################################################
###prep dilution exprSet
#######################################################
  d <- read.csv(filename,check.names=FALSE,row.names=1)
  
  samplenames <- colnames(d)
  ##remove the .cel if its there
  samplenames <- sub("\\.gz$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.Z$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.cel$","",samplenames,ignore.case=TRUE)
  colnames(d) <- samplenames
  ##read phenodata
  data(dilution.phenodata)
  ##putit in order
  d <- d[,rownames(pData(dilution.phenodata))]
  d <- new("exprSet",exprs=as.matrix(d),phenoData=dilution.phenodata)
  d <- exprset.log(d) ##take log
  return(d)
}


exprset.log <- function(exprset){
    e <- exprs(exprset)
    e <- log2(e)
    o <- abs(e)==Inf | is.na(e)
    e[o] <- min(e[!o])
    exprset@exprs <- e
    return(exprset)
  }
