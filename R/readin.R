read.spikein <- function(filename,cdfName=c("hgu95a","hgu133a"),
                         remove.xhyb=TRUE){
#######################################################
###prep spike in ExpressionSet
#######################################################
  cdfName <- match.arg(cdfName)
  s <- read.csv(filename,check.names=FALSE,row.names=1)
  samplenames <- colnames(s)
  ##remove the .cel if its there
  samplenames <- sub("\\.gz$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.Z$","",samplenames,ignore.case=TRUE)
  samplenames <- sub("\\.cel$","",samplenames,ignore.case=TRUE)
  colnames(s) <- samplenames
  ##read phenodata
  if(cdfName=="hgu95a"){
    data(spikein.phenodata)
    pd <- spikein.phenodata
  }
  if(cdfName=="hgu133a"){
    data(hgu133a.spikein.phenodata)
    pd <- hgu133a.spikein.phenodata  
  }
  ##putit in order
  s <- s[,rownames(pData(pd))]
  s <- new("ExpressionSet",exprs=as.matrix(s),phenoData=pd)
  s <- exprset.log(s) ##take log
  if(remove.xhyb & cdfName=="hgu133a") s <- remove.hgu133a.xhyb(s)
  return(s)
}

read.newspikein <- function(filename)
  read.spikein(filename,cdfName="hgu133a")
                                                  
read.dilution <- function(filename){
#######################################################
###prep dilution ExpressionSet
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
  d <- new("ExpressionSet",exprs=as.matrix(d),phenoData=dilution.phenodata)
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
