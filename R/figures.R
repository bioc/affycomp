affycompPlot <- function(...,assessment.list=NULL,method.names=NULL,
                         figure1.xlim=c(-4,15),figure1.ylim=c(-10,12),
                         figure1b.xlim=c(-2,14),figure1b.ylim=c(-6,5),
                         figure6a.xlim=c(-12,12),figure6a.ylim=c(-12,12),
                         figure6b.xlim=c(-3,3),figure6b.ylim=c(-6,6)){
  if(is.null(assessment.list)) l<-list(...) else l <- assessment.list
  N <- length(l)

  if(is.null(method.names)){
    method.names <- vector(mode="character",len=N)
    for(i in 1:N){
      tmp <- l[[i]]$method.name
      if(is.null(tmp)) method.names[i] <- i else method.names[i] <- tmp
    }
  }  

  for(i in 1:N){
    tmp <- l[[i]]
    if(tmp$what=="All") tmp <- tmp[1:5]
    else{
      if(tmp$what=="SpikeIn")
        tmp <- tmp[1:4]
      else{
        if(tmp$what=="SpikeIn2")
          tmp <- tmp[1:3]
        else
          tmp <- list(tmp)
      }
    }
    l[[i]] <- tmp
  }
  
  if(N==1) affycomp.figures(l[[1]])
  else{
    if(is.null(method.names)) method.names <- 1:N
    affycomp.compfigs(l,method.names=method.names,
                      figure1.xlim=figure1.xlim,figure1.ylim=figure1.ylim,
                      figure1b.xlim=figure1b.xlim,figure1b.ylim=figure1b.ylim,
                      figure6a.xlim=figure6a.xlim,figure6a.ylim=figure6a.ylim,
                      figure6b.xlim=figure6b.xlim,figure6b.ylim=figure6b.ylim)
  }
}

affycomp.figures <- function(l){
  n <- length(l)
  
  scn <- prod(par("mfrow"))
  ask <- dev.interactive()
  which.plot <- 0

  ##next block is just so that figures appear in order
  calls <- c()
  is <- c()
  for(i in 1:n){
    tmp <- affycomp.figure.calls(l[[i]]$what)
    calls <- c(calls,tmp)
    is <- c(is,rep(i,length(tmp)))
  }
  is <- is[order(calls)]
  calls <- sort(calls)

  for(i in seq(along=calls)){
    which.plot <- which.plot+1
    if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)
      par(ask=TRUE)
    do.call(calls[i],list(l[[is[i]]]))
    par(ask=FALSE)
  }
}


affycomp.figure.calls <- function(what){  
  args <- c("MA","Dilution","Dilution","Signal",
            "Dilution","FC","FC2","FC","FC","SD",
            "MA2","SpikeInSD","LS",
            "MA2","MA2","MA2")
  fignames <- c("1","2","3","4a","4b","5a","5b","6a","6b","7",
                "1b","2b","4c",
                "5c","5d","5e")
  
  paste("affycomp.figure",fignames[args%in%what],sep="")
}  
  


###MvA plot
affycomp.figure1 <- function(l,main="Figure 1",xlim=NULL,ylim=NULL){
  x <- l$a
  y <- l$m
  if(is.null(ylim)) ylim <- range(y,na.rm=TRUE,finite=TRUE)
  if(is.null(xlim)) xlim <- range(x,na.rm=TRUE,finite=TRUE)
  Index <- match(rownames(l$intended),rownames(x))
  xx <- x[Index,]
  yy <- y[Index,]
  x <- x[-Index,]
  y <- y[-Index,]
  fc <- as.character(l$intended)
  colors <- abs(l$intended)
  colors[colors>11] <- 11
  Colors <- rev((c(topo.colors(8)[1:6],rev(heat.colors(8))[4:8])))

  oo <- sample(1:length(x),12000) #pick a few so the plot isnt too busy
  plot(x[oo],y[oo],pch=".",xlim=xlim,ylim=ylim,main=main,xlab="A",ylab="M",
       las=1)
  o <- abs(y[oo])>1
  points((x[oo])[o],(y[oo])[o],pch=".",col="red")
  o1 <- fc=="Inf"
  o2 <- fc=="-Inf"
  text(xx[!o1 & !o2],yy[!o1 & !o2],fc[!o1 & !o2],col=Colors[colors[!o1 & !o2]])
  # HJ
  if(any(o1))
  	text(xx[o1],yy[o1],expression(infinity),col="black")
  if(any(o2))
  	text(xx[o2],yy[o2],expression(-infinity),col="black")
}


###Variance across replicates
affycomp.figure2 <- function(l,main="Figure 2")
  plot(l$sdplotx,l$sdploty,xlab="log expression",ylab="standard error across replicates",main=main,type="l",lwd=3)

##Sensitivity to amount of RNA
affycomp.figure3 <- function(l,main="Figure 3"){
  x <- l$fc1.25
  y <- l$fc20
  plot(x,y,xlab="Log fold change estimate for 1.25 ug",ylab="Log fold change estimate for 20 ug",main=main,pch=".")
  abline(0,1)
  o <- abs(x-y)>1
  points(x[o],y[o],pch=15,col="orange")
  o <- abs(x-y)>log2(3)
  points(x[o],y[o],pch=16,col="red")
}

##obersved expression v nominal expression
affycomp.figure4a <- function(l,main="Figure 4a",equal.lims=FALSE){
  x <- l$plotx
  y <- l$ploty
  XLIM <- range(x,na.rm=TRUE,finite=TRUE)
  YLIM <- range(y,na.rm=TRUE,finite=TRUE)

  if(equal.lims){
    ranges <- c(diff(XLIM),diff(YLIM))
    slope <- max(ranges)/min(ranges)
    if(diff(XLIM)<diff(YLIM))
      XLIM <- mean(XLIM) + slope*(XLIM-mean(XLIM))
    else
      YLIM <- mean(YLIM) + slope*(YLIM-mean(YLIM))
  }

  cols <- as.numeric(as.factor(names(x)))
  plot(x,y,col=cols,xlab="Nominal concentration (in picoMolar)", ylab="Observed expression",main=main,pch=16,xlim=XLIM,ylim=YLIM)
  lines(l$linex,l$liney,lwd=3)
  if(equal.lims) lines(XLIM,YLIM,lwd=2,lty=2)
}

##slope vs concentration
affycomp.figure4b <- function(l,main="Figure 4b"){
  plot(l$slopeplotx,l$slopeploty,xlab="log expression",ylab="slopes of log expression vs log concentration regresssion",main=main,pch=".")
  abline(h=1)
}

###ROC all FC
affycomp.figure5a <- function(l,main="Figure 5a",maxfp=100){
  x <- l$fp
  y <- l$tp
  Index <- x<=maxfp
  y <- y[Index]
  x <- x[Index]
  plot(x,y,xlab="False Positives",ylab="True Positives",main=main,type="l",lwd=3)
}

##ROC FC=2
affycomp.figure5b <- function(l,main="Figure 5b",maxfp=100)
  affycomp.figure5a(l,main=main,maxfp=maxfp)

##observed FC v intendended
affycomp.figure6a <-function(l,main="Figure 6a",xlim=NULL,ylim=NULL){
  x <- l$intended.log.ratios
  y <- l$observed.log.ratios
  if(is.null(ylim)) ylim <- range(y,na.rm=TRUE,finite=TRUE)
  if(is.null(xlim)) xlim <- range(x,na.rm=TRUE,finite=TRUE)
  matplot(x, y, xlab="Nominal log ratio", ylab="observed log ratio",
          main=main, xlim=xlim, ylim=ylim)
  tmp <- tapply(y,x,mean)
  lines(as.numeric(names(tmp)),tmp,lwd=3)
   hs <- l$quantiles
  N <- length(hs)/2
  abline(h=hs,lty=c(N:1,1:N))
}

##observed FC v intendended low conc
affycomp.figure6b <-function(l,main="Figure 6b",xlim=NULL,ylim=NULL){
  x <- l$intended.log.ratios
  y <- l$observed.log.ratios
  Names <- rep(colnames(x),rep(nrow(x),ncol(x)))[l$index.low.signal]
  x <- as.vector(x)[l$index.low.signal]
  y <- as.vector(y)[l$index.low.signal]
  cols <- as.numeric(as.factor(Names))
  if(is.null(ylim)) ylim <- range(y,na.rm=TRUE,finite=TRUE)
  if(is.null(xlim)) xlim <- range(x,na.rm=TRUE,finite=TRUE)
  plot(x,y,col=cols,xlab="Nominal log ratio",ylab="observed log ratio",main=main,xlim=xlim,ylim=ylim)
  tmp <- tapply(y,x,mean)
  lines(as.numeric(names(tmp)),tmp,lwd=3)
  hs <- l$quantiles
  N <- length(hs)/2
  abline(h=hs,lty=c(N:1,1:N))
}

###SD assessment
affycomp.figure7 <- function(l,main="Figure 7"){
  require(splines,quietly = TRUE)
  x <- l$average.log.expression;y <- l$log.ratio
  o <- sample(1:length(x),5000)
  x <- x[o];y <- y[o]
  oo <- which(!is.na(x) & !abs(x)==Inf)
  x <- x[oo]
  y <- y[oo]
  smooth1 <- lm(y~ns(x,7))
  x1 <- sort(x)[seq(1,length(x),length=100)]
  y1 <- smooth1$fitted[order(x)][seq(1,length(x),length=100)]

  set.seed(1)
  o <- sample(1:length(x),5000)
  plot(x[o],y[o],pch=".",xlab="Log expression",ylab="Log (Nominal SD/Observed SD)",main=main)
  lines(x1,y1,lwd=3,col="red")
  abline(h=0)
}
                               

