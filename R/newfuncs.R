######
##ASSESSMENTS
#######
remove.hgu133a.xhyb <- function(s,bp=c("200","150","100")){
  bp <- match.arg(bp)
  bp <- paste("bp",bp,sep="")
  data(hgu133a.spikein.xhyb)
  out <- match(hgu133a.spikein.xhyb[[bp]],geneNames(s))

  if(any(is.na(out))) warning("No probesets removed. Make sure you are using the hgu133a spikein exprSet or that you didnt already remove the xbybs.")
  else s <- s[-out,]
  return(s)
}

assessSpikeIn2 <- function(s,method.name=NULL,verbose=TRUE){
  if(ncol(exprs(s))==59){
    s <- s[,c(1:13,17,21:33,37)]
    cat("Using only a subset of the spike in data\n")
  }
  if(verbose) cat("Performing 3 assessments that will take a few seconds")
  tmp1 <- assessSpikeInSD(s,method.name=method.name)
  if(verbose) cat(".")
  tmp2 <- assessLS(s,method.name=method.name)
  if(verbose) cat(".")
  tmp3 <- assessMA2(s,method.name=method.name)
  if(verbose) cat(".\n")
  return(list(MA2=tmp3,SpikeInSD=tmp1,LS=tmp2,what="SpikeIn2",
              method.name=method.name))
}

assessSpikeInSD <- function(exprset,method.name=NULL,span=1/3){
  require("splines",quietly = TRUE)
  genenames <- colnames(pData(exprset))
  spikein <-match(genenames,geneNames(exprset))
  y <- esApply(exprset[-spikein,],1,sd)
  x <- esApply(exprset[-spikein,],1,mean)
  smooth1 <- loess(y~x,span=span,family="gaussian",degree=1)
  x2 <- sort(x)[seq(1,length(x),length=100)]
  y2 <- smooth1$fitted[order(x)][seq(1,length(x),length=100)]

  list(x=x,y=y,xsmooth=x2,ysmooth=y2,loess=smooth1,method.name=method.name,what="SpikeInSD")
}

assessLS <-  function(exprset,method.name=NULL){
  e <- exprs(exprset)
  pdata <- pData(exprset)
  genenames <- colnames(pdata)
  y <- as.vector(t(e[match(genenames,geneNames(exprset)),]))
  names(y) <- rep(colnames(pdata),rep(nrow(pdata),ncol(pdata)))
  x <- log2(as.vector(as.matrix(pdata)))
  names(x) <- names(y)
  fit1 <- lm(y~x,subset=x>-Inf)
  tmp <- tapply(y,x,mean) ##this is the line we plot

  nc <- unique(x)
  localslope <- rep(0,length(nc)-1)
  localr2 <- rep(0,length(nc)-1)
  for(i in 1:(length(nc)-1)){
    Index1 <- which(x==nc[i])
    Index2 <- which(x==nc[i+1])
    localslope[i] <- mean(y[Index2])-mean(y[Index1])
    localr2[i]<-cor(c(rep(nc[i+1],length(Index2)),
                      rep(nc[i],length(Index1))),
                    c(y[Index2],y[Index1]))^2
  }

  tmplist <- list()
  tmplist[[1]]   <- which(x <= 0)
  tmplist[[2]]  <- which(x >= 1 & x <= 5)
  tmplist[[3]]   <- which(x>=6)

  threeslopes <- matrix(0,length(tmplist),2) #2 columns: slope and R^2
  for(i in 1:length(tmplist)){
    xx <- x[tmplist[[i]]]
    yy <- y[tmplist[[i]]]
    tmpfit <-lm(yy~xx,subset=xx>-Inf)
    threeslopes[i,] <- c(tmpfit$coef[2],summary(tmpfit)$r.squared)
  }

  list(slope=fit1$coef[2],R2=summary(fit1)$r.squared,
       plotx=x,ploty=y,linex=names(tmp),liney=as.numeric(tmp),
       localslopes=localslope,localr2s=localr2,method.name=method.name,
       low.slope=threeslopes[1,1],
       med.slope=threeslopes[2,1],
       high.slope=threeslopes[3,1],
       low.R2=threeslopes[1,2],
       med.R2=threeslopes[2,2],
       high.R2=threeslopes[3,2],
       what="LS")
}

assessMA2 <- function (exprset, method.name = NULL) 
{
    mat <- exprs(exprset)
    if (ncol(mat) == 28) {
        WHICHSPIKEIN <- "HGU95A"
        NCOMP <- 91 * 2
        I <- c(1:13, 15:27)
    }
    else {
        if (ncol(mat) == 42) {
            WHICHSPIKEIN <- "HGU133A"
            NCOMP <- 91 * 3
            I <- c(1:13, 15:27, 29:41)
        }
        else stop("Not the right number of columns in expression matrix\n")
    }
    pdata <- pData(exprset)
    genenames <- colnames(pdata)
    spikein <- match(genenames, geneNames(exprset))
    quants <- matrix(0, nrow(mat) - length(spikein), NCOMP)
    m <- matrix(0, nrow(mat), NCOMP)
    a <- matrix(0, nrow(mat), NCOMP)
    intended <- matrix(0, length(spikein), NCOMP)
    denom <- matrix(0, length(spikein), NCOMP)
    num <- matrix(0, length(spikein), NCOMP)
    count <- 0
    for (i in I) {
        for (j in (i + 1):(14 * ceiling(i/14))) {
            count <- count + 1
            fc <- mat[, j] - mat[, i]
            quants[, count] <- sort(fc[-spikein])
            a[, count] <- (mat[, i] + mat[, j])/2
            m[, count] <- fc
            num[, count] <- as.numeric(pdata[j, ])
            denom[, count] <- as.numeric(pdata[i, ])
            intended[, count] <- log2(as.numeric(pdata[j, ])/as.numeric(pdata[i, ]))
        }
    }
    set.seed(12345)
    flip = which(sample(c(FALSE, TRUE), NCOMP, replace = T))
    m[, flip] = -m[, flip]
    intended[, flip] = -intended[, flip]
    N2D = num[, flip]
    num[, flip] = denom[, flip]
    denom[, flip] = N2D
    Index <- spikein
    lowIndex <- which(as.vector(num) <= 2 & as.vector(denom) <= 
        2)
    medIndex <- which(as.vector(num) >= 4 & as.vector(denom) >= 
        4 & as.vector(num) <= 32 & as.vector(denom) <= 32 & abs(log2(as.vector(num)/as.vector(denom))) <= 
        2)
    highIndex <- which(as.vector(num) >= 64 & as.vector(denom) >= 
        64 & abs(log2(as.vector(num)/as.vector(denom))) <= 2)
    spikes <- as.vector(m[Index, ])
    nulls <- apply(quants, 1, median)
    tp.low <- sort(abs(spikes[lowIndex]), decreasing = TRUE)
    fp.low <- sapply(tp.low, function(k) sum(abs(nulls) >= k))
    tp.low <- seq(along = tp.low)/length(tp.low)
    tp.med <- sort(abs(spikes[medIndex]), decreasing = TRUE)
    fp.med <- sapply(tp.med, function(k) sum(abs(nulls) >= k))
    tp.med <- seq(along = tp.med)/length(tp.med)
    tp.high <- sort(abs(spikes[highIndex]), decreasing = TRUE)
    fp.high <- sapply(tp.high, function(k) sum(abs(nulls) >= 
        k))
    tp.high <- seq(along = tp.high)/length(tp.high)
    rownames(intended) <- genenames
    SUB <- sample(1:length(m[-Index, ]), length(nulls))
    m <- c(m[Index, ], m[-Index, ][SUB])
    a <- c(a[Index, ], a[-Index, ][SUB])
    names(m)<- c(rep(rownames(mat)[Index], NCOMP), rep(NA,length(m)-length(Index)*NCOMP))
    names(a)<-names(m)
    return(list(qs = nulls, m = m, a = a, spikein = spikein, 
        intended = intended, num = num, denom = denom, fp.low = fp.low, 
        tp.low = tp.low, fp.med = fp.med, tp.med = tp.med, fp.high = fp.high, 
        tp.high = tp.high, method.name = method.name, what = "MA2"))
}


######
## FIGURES
#######
affycomp.compfig4c <- function(l,method.names=as.character(1:length(l)),
                               add.legend=TRUE,rotate=TRUE, main="Figure 4c"){
  xs <- c()
  ys <- c()

  for(i in seq(along=l)){
    if(rotate){
      x <- as.numeric(l[[i]]$linex)[-c(1,2)]
      y <- l[[i]]$localslopes[-1]-1 ##take out the ones related to -Inf
      ylab="Bias"
    }
    else{
      x <- as.numeric(l[[i]]$linex)
      y <- l[[i]]$liney
      In<- which(abs(x)<Inf)
      x<- x[In];y<-y[In]
     ylab<-"Observed log expression"
    }
    xs <- rbind(xs,x)
    ys <- rbind(ys,y)
  }
  if(length(l)==1){ xs <- matrix(xs,nrow=1);ys <- matrix(ys,nrow=1); }
  plot(0,0,type="n",ylim=range(ys),xlim=range(xs),ylab=ylab,
       xlab="Log nominal concentration",main=main,las=1)

  for(i in seq(along=l))
    lines(xs[i,],ys[i,],col=i+1,lty=i,lwd=3)

  if(add.legend){
    if(rotate){
      legend(max(xs)*.7,max(ys)*.95,method.names,lwd=2,col=seq(along=l)+1,lty=seq(along=l))
      abline(h=0)
    }
    else
      legend(min(xs),max(ys)*.95,method.names,lwd=2,col=seq(along=l)+1,lty=seq(along=l))
  }
}
affycomp.figure4c <- function(l,rotate=TRUE, main="Figure 4c")
  affycomp.compfig4c(list(l),method.names=l$method.name,add.legend=FALSE,
                     main=main,rotate=rotate)

affycomp.compfig2b <- function(l,method.names=as.character(1:length(l)),
                              add.legend=TRUE,main="Figure 2b"){
  xs <- c()
  ys <- c()
  for(i in seq(along=l)){
    x <- l[[i]]$xsmooth
    y <- l[[i]]$ysmooth
    xs <- rbind(xs,x)
    ys <- rbind(ys,y)
  }
  plot(0,0,type="n",ylim=range(ys),xlim=range(xs),xlab="Average log expression",ylab="Log expression SD", main=main,las=1)

  for(i in seq(along=l))
    lines(xs[i,],ys[i,],col=i+1,lty=i,lwd=3)

  if(add.legend)
    legend(max(xs)*.7,max(ys)*.95,method.names,lwd=2,col=seq(along=l)+1,lty=seq(along=l))
}
affycomp.figure2b <- function(l,main="Figure 2b")
  affycomp.compfig2b(list(l),method.name=l$method.name,add.legend=FALSE,
                     main=main)

affycomp.compfig5cde <- function(l,method.names=as.character(1:length(l)),
                                 add.legend=TRUE,main="Figure 5c",maxfp=100,
type=c("low","med","high")){
  type <- match.arg(type)
  tp.type <- paste("tp.",type,sep="")
  fp.type <- paste("fp.",type,sep="")
  N <- length(l)
  FP <- vector(mode="list",length=N)
  TP <- vector(mode="list",length=N)
  MIN <- 10
  MAX <- 0
  for(i in 1:N){
    x <- l[[i]][[fp.type]]
    y <- l[[i]][[tp.type]]
    Index <- x<=maxfp
    if(any(Index)){
      MIN <- min(MIN,min(y[Index]))
      MAX <- max(MAX,max(y[Index]))
      TP[[i]] <- c(0,y[Index])
      FP[[i]] <- c(0,x[Index])
    }
    else{ ##there are no true positives!
      TP[[i]] <- c(0,0)
      FP[[i]] <- c(0,maxfp)
      MIN <- 0
    }
  }
  XLIM <- c(0,maxfp)
  YLIM <- c(MIN,MAX)
  plot(1,1,xlab="False Positives",ylab="True Positives",type="n",xlim=XLIM,ylim=YLIM,main=main,las=1)
  for(i in 1:N)
    lines(FP[[i]],TP[[i]],col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.6*(XLIM[2]-XLIM[1]),
           YLIM[1]+.3*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=1)

}

affycomp.compfig5c <- function(l,method.names=as.character(1:length(l)),
                              add.legend=TRUE,main="Figure 5c",maxfp=100)
  affycomp.compfig5cde(l,method.names=method.names,add.legend=add.legend,
                       main=main,maxfp=100,type="low")
affycomp.figure5c <- function(l,main="Figure 5c",maxfp=100)
  affycomp.compfig5c(list(l),method.names=l$method.name,main=main,maxfp=maxfp,
                     add.legend=FALSE)

affycomp.compfig5d <- function(l,method.names=as.character(1:length(l)),
                               add.legend=TRUE,main="Figure 5d",maxfp=100)
  affycomp.compfig5cde(l,method.names=method.names,add.legend=add.legend,
                       main=main,maxfp=100,type="med")
affycomp.figure5d <- function(l,main="Figure 5d",maxfp=100)
  affycomp.compfig5d(list(l),method.names=l$method.name,main=main,maxfp=maxfp,
                     add.legend=FALSE)

affycomp.compfig5e <- function(l,method.names=as.character(1:length(l)),
                              add.legend=TRUE,main="Figure 5e",maxfp=100)
  affycomp.compfig5cde(l,method.names=method.names,add.legend=add.legend,
                        main=main,maxfp=100,type="high")
affycomp.figure5e <- function(l,main="Figure 5e",maxfp=100)
  affycomp.compfig5e(list(l),method.names=l$method.name,main=main,maxfp=maxfp,
                     add.legend=FALSE)


###MvA plot
affycomp.figure1b <- function(l,main="Figure 1b",xlim=NULL,ylim=NULL,cex=.85,
                              all=FALSE){
  x <- l$a
  y <- l$m
  Index <- which(!is.na(names(x)))
  xx <- as.vector(x[Index])
  yy <- as.vector(y[Index])
  x <- as.vector(x[-Index])
  y <- as.vector(y[-Index])

  if(all){
    fc <- as.character(l$intended)
    colors <- abs(l$intended)
    colors[colors>11] <- 11
    Colors <- rev((c(topo.colors(8)[1:6],rev(heat.colors(8))[4:8])))

    if(is.null(ylim))
      ylim <- range(c(y,yy),na.rm=TRUE,finite=TRUE)
    if(is.null(xlim))
      xlim <- range(c(x,xx),na.rm=TRUE,finite=TRUE)

    plot(x,y,pch=".",xlim=xlim,ylim=ylim,main=main,xlab="A",ylab="M",
       las=1)
    o <- abs(y)>1
    points(x[o],y[o],pch=".",col="red")
    o1 <- fc=="Inf"
    o2 <- fc=="-Inf"
    text(xx[!o1 & !o2],yy[!o1&!o2],fc[!o1 & !o2],col=Colors[colors[!o1 & !o2]])
    text(xx[o1],yy[o1],expression(infinity),col="black")
    text(xx[o2],yy[o2],expression(-infinity),col="black")

  }
  else{
    num <- as.vector(l$num); denom <- as.vector(l$denom)
    tmplist <- list()
    tmplist[[1]] <- which(as.vector(num) <= 2 & as.vector(denom) <= 2)
    tmplist[[2]] <- which(as.vector(num) >= 4 & as.vector(denom) >= 4 &
                          as.vector(num) <= 32 & as.vector(denom) <= 32 &
                         abs(log2(as.vector(num)/as.vector(denom)))<=2)
    tmplist[[3]] <- which(as.vector(num) >= 64 & as.vector(denom) >= 64&
                         abs(log2(as.vector(num)/as.vector(denom)))<=2)


    if(is.null(ylim))
      ylim <- range(c(y,yy[unlist(tmplist)]),na.rm=TRUE,finite=TRUE)
    if(is.null(xlim))
      xlim <- range(c(x,xx[unlist(tmplist)]),na.rm=TRUE,finite=TRUE)


    plot(x,y,pch=".",xlim=xlim,ylim=ylim,main=main,xlab="A",ylab="M",
         las=1)

    colors <- c("blue","darkgreen","red")
    FC <- as.character(l$intended)
    for(i in 1:length(tmplist)){
      fc <- FC[ tmplist[[i]] ]
      xxx <- xx[ tmplist[[i]] ]
      yyy <- yy[ tmplist[[i]] ]
      o1 <- fc=="Inf"
      o2 <- fc=="-Inf"
      text(xxx[!o1 & !o2],yyy[!o1 & !o2],fc[!o1 & !o2],col=colors[i],cex=cex)
      text(xxx[o1],yyy[o1],expression(infinity),col=colors[i],cex=cex)
      text(xxx[o2],yyy[o2],expression(-infinity),col=colors[i],cex=cex)
    }
    o <- abs(y)>1
    points((x)[o],(y)[o],pch=".",col="red")
  }

}


#########
## TABLES
#########

tableOverallSNR <- function(...,assessment.list=NULL,method.names=NULL,
                            ngenes=12626){
  if(is.null(assessment.list)) l<-list(...) else l <- assessment.list
  N <- length(l)

  if(is.null(method.names)){
    method.names <- vector(mode="character",len=N)
    for(i in 1:N){
      tmp <- l[[i]]$method.name
      if(is.null(tmp)) method.names[i] <- i else method.names[i] <- tmp
    }
  }

  results <- matrix(0,length(l),6)
  colnames(results) <- c("slope","R2","medianSD","Null FC IQR","Null FC 99.9%","Rank")
  rownames(results) <- method.names
  for(i in 1:N){
    results[i,1] <- l[[i]]$LS$slope
    results[i,2] <- l[[i]]$LS$R2
    results[i,3] <- median(l[[i]]$SpikeInSD$y)
    results[i,4] <-  IQR(l[[i]]$MA$qs)
    results[i,5] <- quantile(abs(l[[i]]$MA$qs),prob=c(0.999))
    qs <- l[[i]]$MA2$qs
    tmp <- abs(qs-results[i,1])
    tmp <- which(tmp==min(tmp))/length(qs)
    results[i,6] <- round(ngenes-tmp*ngenes)
  }
  return(results)
}

tableLS <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,length(l[[1]]$localslopes)+6,N)
  colnames(results) <- method.names
  tmp <-  as.character((2^as.numeric(l[[1]]$linex)))
  rownames(results) <- c("low.slope","med.slope","high.slope",
                         "low.R2","med.R2","high.R2",
                         paste(tmp[-1],tmp[-length(tmp)],sep=":"))

  for(i in 1:N){
    results[1:6,i] <- unlist(l[[i]][10:15])
    results[-(1:6),i] <- l[[i]]$localslopes
  }
  return(results)
}

tableSpikeInSD <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,4,N)
  colnames(results) <- method.names
  rownames(results) <- c("25% SD",
                         "Median SD",
                         "75% SD",
                         "99% SD")
  for(i in 1:N)
    results[,i] <- quantile(l[[i]]$y,prob=c(0.25,0.5,0.75,0.99))
  return(results)
}

tableMA2 <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,7,N)
  colnames(results) <- method.names
  rownames(results) <- c("null log-fc IQR",
                         "null log-fc 99%",
                         "null log-fc 99.9%",
                         "low AUC",
                         "med AUC",
                         "high AUC",
                         "weighted avg AUC")

  myauc <- function(mat){
    y <- tapply(c(0,mat[,2]),c(0,mat[,1]),max)
    x <- as.numeric(names(y))
    y <- y[x<100]
    x <- x[x<100]
    return(mean(approx(x,y,0:100,method="constant",rule=2)$y))
  }
  for(i in 1:N){
    results[1,i] <- IQR(l[[i]]$qs)
    results[2:3,i] <- quantile(abs(l[[i]]$qs),prob=c(0.99,0.999))
    results[4,i] <- myauc(cbind(l[[i]]$fp.low,l[[i]]$tp.low))
    results[5,i] <- myauc(cbind(l[[i]]$fp.med,l[[i]]$tp.med))
    results[6,i] <- myauc(cbind(l[[i]]$fp.high,l[[i]]$tp.high))
    results[7,i] <- sum(c(0.75,0.24,0.01)*results[4:6,i])
  }
  return(results)
}
tableRanks <- function(...,assessment.list=NULL,method.names=NULL,
                            ngenes=12626,rank=TRUE){
  if(is.null(assessment.list)) l<-list(...) else l <- assessment.list
  N <- length(l)

  if(is.null(method.names)){
    method.names <- vector(mode="character",len=N)
    for(i in 1:N){
      tmp <- l[[i]]$method.name
      if(is.null(tmp)) method.names[i] <- i else method.names[i] <- tmp
    }
  }

  results <- matrix(0,length(l[[1]]$LS$localslopes),N)
  colnames(results) <- method.names
  tmp <-  as.character((2^as.numeric(l[[1]]$LS$linex)))
  rownames(results) <- paste(tmp[-1],tmp[-length(tmp)],sep=":")

  for(i in 1:N)
    results[,i] <- l[[i]]$LS$localslopes

  if(!rank) return(results)
  else{
    for(i in 1:N){
       qs <- l[[i]]$MA2$qs
       results[,i] <- sapply(results[,i],function(x){
         tmp <- abs(qs-x)
         tmp <- which(tmp==min(tmp))/length(qs)
         round(ngenes-tmp*ngenes)+1
       })
     }
    return(results)
  }
}














