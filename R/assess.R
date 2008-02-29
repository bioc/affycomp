affycomp <- function(d,s,method.name=NULL,verbose=TRUE,return.it=TRUE){
  if(is.null(method.name)) method.name <- "New expression measure"
  l <- assessAll(d,s,method.name=method.name,verbose=verbose)
  data(mas5.assessment)
  affycompPlot(mas5.assessment,l)
  tmp <- affycompTable(mas5.assessment,l)
  print(format(tmp,digits=2))
  if(return.it) return(l)
  else return(NULL)
}

assessSpikeIn <- function(s,method.name=NULL,verbose=TRUE){
  if(verbose) cat("Performing 6 assessments that will take a few minutes")
  tmp1 <- assessFC(s,method.name=method.name)
  if(verbose) cat("...")
  tmp2 <- assessFC2(s,method.name=method.name)
  if(verbose) cat(".")
  tmp3 <- assessMA(s,method.name=method.name)
  if(verbose) cat(".")
  tmp4 <- assessSignal(s,method.name=method.name)
  if(verbose) cat(".\n")
  return(list(MA=tmp3,Signal=tmp4,FC=tmp1,FC2=tmp2,what="SpikeIn",
              method.name=method.name))
}

assessAll <- function(d,s,method.name=NULL,verbose=TRUE){
  if(verbose) cat("Performing 9 assessments that will take a few minutes\nWe start with 3 on dilution data")
  tmp1 <- assessDilution(d,method.name=method.name)
  cat("...\n")
  tmp2 <- assessSpikeIn(s,verbose=verbose,method.name=method.name)
  tmp <- c(Dilution=list(tmp1),tmp2)
  tmp["what"] <- "All"
  return(tmp)
}

assessDilution <- function(exprset,method.name=NULL){
  require(splines,quietly = TRUE)

  e <- exprs(exprset)
  pdata <- pData(exprset)
  ##GET R^2: depends on the right order
  o <- c()
  for(j in 1:12) o <- rbind(o,expand.grid( ( (j-1)*5+1 ):( j*5),( (j-1)*5+1 ):( j*5)))
  o <- o[o[,1]< o[,2],] ##o gives us all replicate pairs
  o <- as.matrix(o)
  rownames(o) <- NULL
  R2 <- apply(o,1,function(x) (cor((e[,x]))[1,2])^2)

  ##average per concentration group
  tmp <- cbind(c(1.25,0,2.5,0,5,0,7.5,0,10,0,20,0),
               c(0,1.25,0,2.5,0,5,0,7.5,0,10,0,20))
  m <- apply(tmp,1,function(x){
    o <- which(pdata[,1]==x[1] & pdata[,2]==x[2])
    rowMeans(e[,o])
})
  ##sd per concentration group
  s <- apply(tmp,1,function(x){
    o <- which(pdata[,1]==x[1] & pdata[,2]==x[2])
    apply(e[,o],1,sd)
  })

  ##we will "normalize" using the spike-in
  spikedin <- colnames(pdata)[-c(1:2,ncol(pdata))]
  for(j in 1:2){ ##1 is CNS and 2 liver
    Index <- tmp[,j]==0 ##first for liver
    k <- colMeans(m[spikedin, Index]) ## these should be the same
    k <- k-mean(k)##so let's make them the same
    m[,Index] <- sweep(m[,Index],2,k)
  }

  ##regression coef for each gene for liver and cns. get slope

  o <- which(tmp[,1]!=0)
  x <- log2(tmp[o,1])
  x <- x-mean(x)
  x2 <- sum(x^2)
  beta1 <- apply(m[,o],1,function(y) sum((y-mean(y))*x)/x2)##fast regression
  o <- which(tmp[,2]!=0)
  x <- log2(tmp[o,2])
  x <- x-mean(x)
  x2 <- sum(x^2)
  beta2 <- apply(m[,o],1,function(y) sum((y-mean(y))*x)/x2)

  ##consistency
  fc1.25 <- m[,1]-m[,2]
  fc20 <- m[,11]-m[,12]

  ##things to keep: curves, etc..
  ##sd vs mean plot
  x <- as.vector(m)
  y <- as.vector(s)
  smooth1 <- lm(y~ns(x,7))
  x1 <- sort(x)[seq(1,length(x),length=100)]
  y1 <- smooth1$fitted[order(x)][seq(1,length(x),length=100)]
  ##took out this plot. will keep actual points instead.
  ##beta vs mean plot
  x <- rowMeans(m)
  x <- c(x,x)
  y <- c(beta1,beta2)
  smooth1 <- lm(y~ns(x,7))
  x2 <- sort(x)[seq(1,length(x),length=100)]
  y2 <- smooth1$fitted[order(x)][seq(1,length(x),length=100)]

  list(R2=R2,sdplotx=x1,sdploty=y1,slopeplotx=x,slopeploty=y,
       mediansd=median(s),medianbeta=median( c(beta1,beta2)),
       fc1.25=fc1.25,fc20=fc20,
       consitency=cor(fc1.25,fc20)^2,
       two.fold.discrepancy=sum(abs(fc1.25-fc20)>1),
       three.fold.discrepancy=sum(abs(fc1.25-fc20)>log2(3)),
       slopesmoothx=x2,slopesmoothy=y2,what="Dilution",method.name=method.name)
}

assessSignal <- function(exprset,method.name=NULL){
  e <- exprs(exprset)
  pdata <- pData(exprset)
  genenames <- colnames(pdata)
  y <- as.vector(t(e[match(genenames,geneNames(exprset)),]))
  names(y) <- rep(colnames(pdata),nrow(pdata))
  x <- log2(as.vector(as.matrix(pdata)))
  names(x) <- names(y)
  tmp <- tapply(y,x,mean)
  fit1 <- lm(y~x,subset=x>-Inf)
  list(slope=fit1$coef[2],R2=summary(fit1)$r.squared,
       plotx=x,ploty=y,linex=names(tmp),liney=as.numeric(tmp),what="Signal",method.name=method.name)
}

assessFC2 <- function(exprset,method.name=NULL){
  e <- exprs(exprset)
  pdata <- pData(exprset)

  ##this is bad, but works
  if(ncol(e)==59)
    WHICHSPIKEIN <- "HGU95A"
  else{
    if(ncol(e)==42){
      WHICHSPIKEIN <- "HGU133A"
    }
    else
      stop("Not the right number of columns in expression matrix\n")
  }



  if(WHICHSPIKEIN=="HGU95A"){
    e <- e[,-grep("2353d99hpp_av08",colnames(e))] ##take d_08 because it has no pair
    o1 <- c(2 ,4 ,6 ,8 ,10,12,17,18,19,20)
    o1 <- c(o1,o1+20,c(42,44,46,48,50,55,56,57)) ##took out 58 and 54.. 54 is bad
    o2 <- c(1,3,5,7, 9,11,13,14,15,16)
    o2 <- c(o2,o2+20,c(41,43,45,47,49,51,52,53))##took out 54.see above
  }
  else{
    o1 <- seq(2,42,2)
    o2 <- seq(1,41,2)
  }

  m <- e[,o1]-e[,o2]
  a <- (e[,o1] + e[,o2])/2

  fc2 <- apply(m,2,function(x){
    tp <- sum(abs(x[match(colnames(pdata),names(x))]) >= 1) ##true positives
    fp <- sum(abs(x[-match(colnames(pdata),names(x))]) >= 1) ##false positives
    c(fp=fp,tp=tp)
  })

  rocs <- apply(m,2,function(x){
    x <- sort(-abs(x))
    y <- rep(0,length(x))
    y[match(colnames(pdata),names(x))] <- 1
    return(cumsum(y)) ##this is number of true positive
  })
  tp <- rowMeans(rocs)
  fp <- seq(along=tp)-tp ##total calls minus true positives

  N <- ncol(pdata)
  return(list(fc2=t(fc2),m=m,a=a,fp=fp,tp=tp,
              area=c(a10=mean(tp[fp<10]/N),
                a15=mean(tp[fp<15]/N),
                a25=mean(tp[fp<25]/N),
                a100=mean(tp[fp<100]/N)),what="FC2",method.name=method.name))
}


assessFC <- function(exprset,method.name=NULL){
  e <- exprs(exprset)
  pdata <- pData(exprset)

  ##this is not good but works..
  if(ncol(e)==59)
    WHICHSPIKEIN <- "HGU95A"
  else{
    if(ncol(e)==42)
      WHICHSPIKEIN <- "HGU133A"
    else
      stop("Not the right number of columns in expression matrix\n")
  }


  genenames <- colnames(pdata)
  N <- length(genenames)
  M <- nrow(e) - N
  intended <- array(0,dim=c(570,N,2)) ##570 is a maximum.. later reduced
  observed <- matrix(0,570,N)
  fc2 <- matrix(0,570,2)
  probs <- c(0,25/M,100/M,.25,.75,1-100/M,1-25/M,1)
  quantiles <- matrix(0,570,length(probs))
  pdata <- as.matrix(pdata)
  spikeindex <- match(genenames,rownames(e))
  roc <- vector(mode="numeric",length=nrow(e))
  Count <- 0

  if(WHICHSPIKEIN=="HGU95A") J <- 20 else J <- 14
  for(i in 1:(J-1)){
    for(j in (i+1):J){
      i1 <- pdata[i,]
      i2 <- pdata[j,]
      if(!all(i1-i2==0)){ ##if not reps lets do it
        Count <- Count + 1
        intended[Count,,1] <- i1
        intended[Count,,2] <- i2
        m <- e[,j]-e[,i]
        quantiles[Count,] <- quantile(m[-spikeindex],prob=probs)
        fc2[Count,1] <- sum(abs(m[-spikeindex]) >= 1)
        fc2[Count,2] <- sum(abs(m[spikeindex]) >= 1)
        observed[Count,] <- m[genenames]
        m <- sort(-abs(m))
        y <- rep(0,length(m))
        y[match(genenames,names(m))] <- 1
        roc <- roc + cumsum(y) ##this is number of true positive
      }
    }
  }


  for(i in (J+1):(2*J-1)){
    for(j in (i+1):(2*J)){
      i1 <- pdata[i,]
      i2 <- pdata[j,]
      if(!all(i1-i2==0)){ ##if not reps lets do it
        Count <- Count + 1
        intended[Count,,1] <- i1
        intended[Count,,2] <- i2
        m <- e[,j]-e[,i]
        quantiles[Count,] <- quantile(m[-spikeindex],prob=probs)
        fc2[Count,1] <- sum(abs(m[-spikeindex]) >= 1)
        fc2[Count,2] <- sum(abs(m[spikeindex]) >= 1)
        observed[Count,] <- m[genenames]
        m <- sort(-abs(m))
        y <- rep(0,length(m))
        y[match(genenames,names(m))] <- 1
        roc <- roc + cumsum(y) ##this is number of true positive
      }
    }
  }

  if(WHICHSPIKEIN=="HGU95A"){
    J1 <- 41; J2 <- 58
  }
  else{
    J1 <- 29; J2 <- 42
  }

  for(i in J1:(J2-1)){
    for(j in (i+1):J2){
      i1 <- pdata[i,]
      i2 <- pdata[j,]
      if(!all(i1-i2==0) & i!=54 & j!=54){ ##if not reps lets do it
        Count <- Count + 1
        intended[Count,,1] <- i1
        intended[Count,,2] <- i2
        m <- e[,j]-e[,i]
        quantiles[Count,] <- quantile(m[-spikeindex],prob=probs)
        fc2[Count,1] <- sum(abs(m[-spikeindex]) >= 1)
        fc2[Count,2] <- sum(abs(m[spikeindex]) >= 1)
        observed[Count,] <- m[genenames]
        m <- sort(-abs(m))
        y <- rep(0,length(m))
        y[match(genenames,names(m))] <- 1
        roc <- roc + cumsum(y) ##this is number of true positive
      }
    }
  }

  intended <- intended[1:Count,,]
  observed <- observed[1:Count,]
  fc2 <- fc2[1:Count,]
  quantiles <- quantiles[1:Count,]
  quantiles <- colMeans(quantiles)
  tp <- roc/Count
  fp <- seq(along=tp)-tp ##total calls minus true positives

  colnames(observed) <- genenames
  dimnames(intended) <- list(NULL,genenames,NULL)
  names(quantiles) <- c("lowest","lowest25","lowest100","25",
                     "75","highest100","highest25","highest")


  intended.log.ratios <- log2(intended[,,2]/intended[,,1])
  x <- as.vector(intended.log.ratios)
  y <- as.vector(observed)
  Index <- as.vector(intended[,,2])<=2 & as.vector(intended[,,1])<=2 &
  as.vector(intended[,,2])>0 & as.vector(intended[,,1])>0 ##small signal

  N <- ncol(pdata)
  list(signal=intended,
       intended.log.ratios=intended.log.ratios,
       observed.log.ratios=observed,
       quantiles=quantiles,
       fc2=fc2,
       tp=tp,fp=fp,
       area=c(a10=mean(tp[fp<10]/N),
         a15=mean(tp[fp<15]/N),
         a25=mean(tp[fp<25]/N),
         a100=mean(tp[fp<100]/N)),
       slope=lm(y~x,subset=abs(x)<Inf)$coef[2],
       low.signal.slope= lm(y~x,subset=Index)$coef[2],
       index.low.signal=Index,what="FC",method.name=method.name)
}


assessMA<- function(exprset,method.name=NULL){
  e <- exprs(exprset)
  pdata <- pData(exprset)

  if(ncol(e)==59)
    WHICHSPIKEIN <- "HGU95A"
  else{
    if(ncol(e)==42)
      WHICHSPIKEIN <- "HGU133A"
    else
      stop("Not the right number of columns in expression matrix\n")
  }

  N <- nrow(e)
  genenames <- colnames(pdata)
  spikeindex <- match(genenames,rownames(e))
  M <- length(genenames)

  thearray <- 1
  if(WHICHSPIKEIN=="HGU95A"){
    thearray <- 1
    otherarrays <-  c(2:12,13,17)
  }
  else{
    thearray <- 1
    otherarrays <- 2:14#c(1:4,6:14)
  }
  L <- length(otherarrays)
  m <- matrix(0,N,L)
  a <- matrix(0,N,L)
  intended <- matrix(0,M,L)

  for(i in 1:L){
    m[,i] <- e[,otherarrays[i]] - e[,thearray]
    a[,i] <- (e[,otherarrays[i]]+e[,thearray])/2
    intended[,i] <- log2(as.numeric(pdata[otherarrays[i],])/as.numeric(pdata[thearray,]))
  }
  rownames(m) <- rownames(e)
  rownames(a) <- rownames(e)
  rownames(intended) <- genenames
  colnames(m) <- colnames(e)[otherarrays]
  colnames(a) <- colnames(e)[otherarrays]
  colnames(intended) <- colnames(e)[otherarrays]

  list(m=m,a=a,intended=intended,index=spikeindex,what="MA",method.name=method.name)
}

assessSD <- function(exprset,method.name=NULL,logx=FALSE){
  e <- exprs(exprset)
  se <- assayDataElement(exprset,"se.exprs")
  pdata <- pData(exprset)

  tmp <- cbind(c(1.25,0,2.5,0,5,0,7.5,0,10,0,20,0),
               c(0,1.25,0,2.5,0,5,0,7.5,0,10,0,20))

  m <- apply(tmp,1,function(x){
    o <- which(pdata[,1]==x[1] & pdata[,2]==x[2])
    rowMeans(e[,o])
  })
  ##sd per concentration group
  observed <- apply(tmp,1,function(x){
    o <- which(pdata[,1]==x[1] & pdata[,2]==x[2])
    apply(e[,o],1,sd)
  })

  nominal <- apply(tmp,1,function(x){
    o <- which(pdata[,1]==x[1] & pdata[,2]==x[2])
    sqrt(rowMeans(se[,o]^2))
})

  tmp1 <- log2(as.vector(nominal))
  tmp2 <- log2(as.vector(observed))
  y <- tmp1-tmp2
  x <- as.vector(m)
  if(logx) x <- log2(x)

  list(average.log.expression=x,log.ratio=y,what="SD",
       corr=cor(tmp1,tmp2),method.name=method.name)
}









