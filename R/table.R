affycompTable <- function(...,Table=NULL,
                          assessment.list=NULL,method.names=NULL){
  if(is.null(Table)) tmp <- tableAll(...,assessment.list=assessment.list,method.names=method.names) else tmp <- Table
  tmp <- tmp[c(1,2, 3,4,5, 7,8, 6, 12,13,14, 21,22,23, 15,16,17),]
  tmp <- data.frame(tmp)
  tmp$whatsgood <- c(0,1,1,0,0,1,1,1,1,0,16,1,0,16,0,1,1)
  tmp$Figure <- c("2","2","3","3","3","4a","4a","4b","5a","5a","5a",
                  "5b","5b","5b","6","6a","6b")
  tmp
}
  
tableAll <- function(...,assessment.list=NULL,method.names=NULL){
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
        tmp <- list(tmp)
        names(tmp) <- tmp[[1]]$what
      }
    }
    l[[i]] <- tmp
  }

  auxcalls <- names(l[[1]])
  auxcalls <- auxcalls[!auxcalls%in%"MA"]
  calls <- paste("table",auxcalls,sep="")
  
  results <- c()
  for(i in seq(along=calls)){
    thelist <- vector(mode="list",length=N)
    for(j in 1:N) thelist[[j]] <- l[[j]][[auxcalls[i]]]
    tmp <- do.call(calls[i],list(thelist,method.names=method.names))
    results <- rbind(results,tmp)
  }
  results
}

tableDilution <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,6,length(l))
  colnames(results) <- method.names
  rownames(results) <- c("Median SD",
                         "R2",
                         "1.25v20 corr",
                         "2-fold discrepancy",
                         "3-fold discrepancy",
                         "Median slope")

  for(i in 1:length(l)){
    results[1,i] <- l[[i]]$mediansd
    results[2,i] <- mean(l[[i]]$R2)
    results[3,i] <- l[[i]]$consitency
    results[4,i] <- l[[i]]$two.fold.discrepancy
    results[5,i] <- l[[i]]$three.fold.discrepancy
    results[6,i] <- l[[i]]$medianbeta

  }
  return(results)
}

tableFC <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,9,length(l))
  colnames(results) <- method.names
  rownames(results) <- c("AUC (FP<10)",
                         "AUC (FP<15)",
                         "AUC (FP<25)",
                         "AUC (FP<100)",
                         "AFP, call if fc>2",
                         "ATP, call if fc>2",
                         "IQR",
                         "Obs-intended-fc slope",
                         "Obs-(low)int-fc slope")
  
  for(i in 1:length(l)){
    results[1:4,i] <- unlist(l[[i]]$area)
    results[5,i] <- mean(l[[i]]$fc2[,1])
    results[6,i] <- mean(l[[i]]$fc2[,2])
    results[7,i] <- -diff(l[[i]]$quantiles[4:3])
    results[8,i] <- mean(l[[i]]$slope)
    results[9,i] <- mean(l[[i]]$low.signal.slope)
  }
  return(results)
}

tableFC2 <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,6,length(l))
  colnames(results) <- method.names
  rownames(results) <- c("FC=2, AUC (FP<10)",
                         "FC=2, AUC (FP<15)",
                         "FC=2, AUC (FP<25)",
                         "FC=2, AUC (FP<100)",
                         "FC=2, AFP, call if fc>2",
                         "FC=2, ATP, call if fc>2")
  
  for(i in 1:length(l)){
    results[1:4,i] <- unlist(l[[i]]$area)
    results[5,i] <- mean(l[[i]]$fc2[,1])
    results[6,i] <- mean(l[[i]]$fc2[,2])
  }
  return(results)
}


tableSignal <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,2,length(l))
  colnames(results) <- method.names
  rownames(results) <- c("Signal detect slope",
                         "Signal detect R2")
  for(i in 1:length(l)){
    results[1,i] <- mean(l[[i]]$slope)
    results[2,i] <- l[[i]]$R2
  }
  return(results)
}


tableSD <- function(l,method.names=NULL){
  N <- length(l)
  if(is.null(method.names)) method.names <- 1:N
  results <- matrix(0,2,length(l))
  colnames(results) <- method.names
  rownames(results) <- c("IQR of log ratio",
                         "Correlation")
  for(i in 1:length(l)){
    results[1,i] <- IQR(l[[i]]$log.ratio)
    results[2,i] <- l[[i]]$corr
  }
  return(results)
}

