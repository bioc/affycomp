affycomp.compfigs <- function(l,method.names=NULL,
                              figure1.xlim=c(-4,15),figure1.ylim=c(-10,12),
                              figure6a.xlim=c(-12,12),figure6a.ylim=c(-12,12),
                              figure6b.xlim=c(-3,3),figure6b.ylim=c(-6,6)){
  ##if spikein or all take out the what. otherwise make list
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  n <- length(l[[1]])
  
  scn <- prod(par("mfrow"))
  ask <- dev.interactive()
  which.plot <- 0

  ##next block is just so that figures appear in order
  fororder <- c()
  calls <- c()
  is <- c()
  for(i in 1:n){
    tmp <- affycomp.compfigs.calls(l[[1]][[i]]$what)
    calls <- c(calls,tmp)
    fororder <- c(fororder,affycomp.figure.calls(l[[1]][[i]]$what))
    is <- c(is,rep(i,length(tmp)))
  }
  is <- is[order(fororder)]
  calls <- calls[order(fororder)]

  
  for(i in seq(along=calls)){
    if(length(grep("compfig",calls[i]))>0){
      thelist <- vector(mode="list",length=N)
      for(j in 1:N) thelist[[j]] <- l[[j]][[ is[i] ]]
      which.plot <- which.plot+1
      if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)
        par(ask=TRUE)
      do.call(calls[i],list(thelist,method.names=method.names))
      par(ask=FALSE)
    }
    else{
      for(j in 1:N){
        which.plot <- which.plot+1;
        if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask) par(ask=TRUE)
        tmp <- strsplit(calls[i],"\\.")[[1]][2]
        do.call(calls[i],
                list(l[[j]][[ is[i] ]],main=method.names[j],
                     xlim=get(paste(tmp,"xlim",sep=".")),
                     ylim=get(paste(tmp,"ylim",sep="."))))
        par(ask=FALSE)
      }
    }
  }
}
  
affycomp.compfigs.calls <- function(what){  
  args <- c("MA","Dilution","Dilution","Signal",
            "Dilution","FC","FC2","FC","FC","SD")
  fignames <- c("figure1","compfig2","compfig3","compfig4a","compfig4b",
                "compfig5a","compfig5b","figure6a","figure6b","compfig7")
  
  paste("affycomp.",fignames[args%in%what],sep="")
}  

   
affycomp.compfig2 <- function(l,method.names=as.character(1:length(l)),add.legend=TRUE,main="Figure 2"){
  N <- length(l)
  XLIM <- NA
  YLIM <- NA
  for(i in 1:N){
    XLIM <- range(c(XLIM,l[[i]]$sdplotx),finite=TRUE,na.rm=TRUE)
    YLIM <- range(c(YLIM,l[[i]]$sdploty),finite=TRUE,na.rm=TRUE)
  }
  plot(1,1,xlim=XLIM,ylim=YLIM,xlab="log expression",ylab="standard error across replicates",type="n",main=main)
  for(i in 1:N)
    lines(l[[i]]$sdplotx,l[[i]]$sdploty,col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.6*(XLIM[2]-XLIM[1]),
           YLIM[1]+.9*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=2)
}

affycomp.compfig3 <- function(l,method.names=as.character(1:length(l)),main="Figure 3"){
  N <- length(l)
  tmp <- vector(mode="list",length=N)
  for(i in 1:N){
    x <- l[[i]]$fc1.25
    y <- l[[i]]$fc20
    tmp[[i]] <- x-y
  }
  boxplot(tmp,names=method.names,range=0,col=2:(N+1),main=main,ylab="Difference between log fold changes")
}

affycomp.compfig4b <- function(l,method.names=as.character(1:length(l)),add.legend=TRUE,main="Figure 4b"){
  N <- length(l)
  XLIM <- c(NA,1)
  YLIM <- c(NA,1)
  for(i in 1:N){
    XLIM <- range(c(XLIM,l[[i]]$slopesmoothx),finite=TRUE,na.rm=TRUE)
    YLIM <- range(c(YLIM,l[[i]]$slopesmoothy),finite=TRUE,na.rm=TRUE)
  }
  plot(1,1,xlim=XLIM,ylim=YLIM,xlab="log expression",ylab="regression slopes",type="n",main=main)
  for(i in 1:N)
    lines(l[[i]]$slopesmoothx,l[[i]]$slopesmoothy,col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.7*(XLIM[2]-XLIM[1]),
           YLIM[1]+.98*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=2)
  abline(h=1)
}

affycomp.compfig4a <- function(l,method.names=as.character(1:length(l)),add.legend=TRUE,main="Figure 4a"){
  N <- length(l)
  XLIM <- NA
  YLIM <- NA
  for(i in 1:N){
    XLIM <- range(c(XLIM,as.numeric(l[[i]]$linex)),finite=TRUE,na.rm=TRUE)
    YLIM <- range(c(YLIM,as.numeric(l[[i]]$liney)),finite=TRUE,na.rm=TRUE)
  }
  plot(1,1,xlim=XLIM,ylim=YLIM,,xlab="Nominal concentration (in picoMolar)", ylab="Observed expression",type="n",main=main)
  for(i in 1:N)
    lines(l[[i]]$linex,l[[i]]$liney,col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.6*(XLIM[2]-XLIM[1]),
           YLIM[1]+.4*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=2)
}

affycomp.compfig5a <- function(l,method.names=as.character(1:length(l)),add.legend=TRUE,main="Figure 5a",maxfp=100){
  N <- length(l)
  FP <- vector(mode="list",length=N)
  TP <- vector(mode="list",length=N)
  for(i in 1:N){
    x <- l[[i]]$fp
    y <- l[[i]]$tp
    Index <- x<=maxfp
    TP[[i]] <- y[Index]
    FP[[i]] <- x[Index]
  }
  XLIM <- range(unlist(FP))
  YLIM <- range(unlist(TP))
  plot(1,1,xlab="False Positives",ylab="True Positives",type="n",xlim=XLIM,ylim=YLIM,main=main)
  for(i in 1:N)
    lines(FP[[i]],TP[[i]],col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.6*(XLIM[2]-XLIM[1]),
           YLIM[1]+.3*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=2)

}

affycomp.compfig5b <- function(l,method.names=as.character(1:length(l)),add.legend=TRUE,main="Figure 5b",maxfp=100){
  N <- length(l)
  FP <- vector(mode="list",length=N)
  TP <- vector(mode="list",length=N)
  for(i in 1:N){
    x <- l[[i]]$fp
    y <- l[[i]]$tp
    Index <- x<=maxfp
    TP[[i]] <- y[Index]
    FP[[i]] <- x[Index]
  }
  XLIM <- range(unlist(FP))
  YLIM <- range(unlist(TP))
  plot(1,1,xlab="False Positives",ylab="True Positives",type="n",xlim=XLIM,ylim=YLIM,main=main)
  for(i in 1:N)
    lines(FP[[i]],TP[[i]],col=i+1,lty=i,lwd=3)
  if(add.legend)
    legend(XLIM[1]+.6*(XLIM[2]-XLIM[1]),
           YLIM[1]+.6*(YLIM[2]-YLIM[1]),
           method.names,col=2:(N+1),lty=1:N,lwd=2)
}

affycomp.compfig7 <- function(l,method.names=as.character(1:length(l)),
                              main="Figure 7"){
  N <- length(l)
  tmp <- vector(mode="list",length=N)
  for(i in 1:N) tmp[[i]] <- l[[i]]$log.ratio
  boxplot(tmp,names=method.names,range=0,col=2:(N+1),main=main,ylab="Log ratio of nominal SD and observed SD")
}

