fp1 <- setClass("fp1", representation(dens="array", diff="data.frame", dat="data.frame"))

plot.fp1 <- function(x, ..., ylab=c("Density","Density difference"), xlim=NULL) {
  options(warn=-1)
  if (is.null(xlim)) xlim <- range(x@dens[[1]]$x)
  plot(1,type="n",
       ylim=c(0,max(unlist(lapply(x@dens, function(X) {max(X$y)})))),
       xlim=xlim,
       bty="n", ylab=ylab[1], ...)
  for (i in 1:length(x@dens)) {
    lines(x@dens[[i]], col=i)
  }
  
  matplot(x@dens[[1]]$x,x@diff, type='l', pch=16,
          lwd=1, bty="n", ylab=ylab[2], xlim=xlim, ...); abline(h=0, lty=3)
}
setMethod("plot","fp1", function(x, ...) plot.fp1(x, ...))

fpPlot <- function(...) {
  plot.fp1(...) # just a wrapper
}

fpGet <- function(dat, n=512, bw='nrd0') {
  # dat: nx2 dataframe or matrix with in col 1: RT; col 2: condition
  if (is.matrix(dat)) dat <- as.data.frame(dat)
  
  rng <- range(dat[[1]])
  dens <- tapply(dat[[1]], dat[[2]], density, from=rng[1], to=rng[2], n=n, bw=bw)
  diff <- NULL
  for (i in 2:length(dens)) {
    for (j in 1:(i-1)) {
      
      # difference (density-based method)
      tmp <- data.frame(dens[[i]]$y - dens[[j]]$y)
      names(tmp) <- paste(i,j,sep='-')
      if (is.null(diff)) diff <- tmp else diff <- cbind(diff, tmp)
    }
  }
  fp1(dens=dens, diff=diff, dat=dat)
}

fpDensDiff <- function(object) {
  if (!is.list(object)) {
    object <- list(object)
  }
  
  .get.diffs <- function(X) {
    lwr <- min(unlist(lapply(X@dens, function(Y) {which.max(diff(Y$y))})))#min(unlist(lapply(X@dens, function(Y) {which.max(Y$y)})))
    upr <- max(unlist(lapply(X@dens, function(Y) {which.min(diff(Y$y))})))
    
    sapply(1:ncol(X@diff), function(i) {
      y <- X@diff[,i][lwr:upr]
      x <- X@dens[[1]]$x[lwr:upr]
      y1 <- y[!is.na(y)&y!=Inf&y!=0]
      x1 <- x[!is.na(y)&y!=Inf&y!=0]
      index <- which.min(abs(y1))
      x1[index]
    })
  }
  
  roots <- lapply(object, .get.diffs)
  root <- array(dim=c(dim(object[[1]]@diff)[2],length(object)))
  for (i in 1:length(roots)) {
    root[,i] <- unlist(roots[[i]])
  }
  root
}

fpAnova <- function(object, stat="BF", na.rm=TRUE) {
  require(BayesFactor)
  bf <- p <- NULL
  tmp <- fpDensDiff(object)
  tmp <- data.frame(x=c(tmp), cross=factor(1:nrow(tmp)),pp=factor(rep(1:ncol(tmp),each=nrow(tmp))))
  # because tmp is a pp x cross array, we need to test across pp whether the ratios 
  # differ.
  if (na.rm) tmp <- tmp[!is.na(tmp$cross),]
  if (stat=="BF"|stat=="both") {
    bf <- anovaBF(x~cross+pp, whichRandom="pp", data=tmp, progress=F)
  }
  if (stat=="p"|stat=="both") {
    p <- summary(aov(x~cross+Error(pp/cross), data=tmp))
  }
  list(BF=bf, p=p)
}

dnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: quantiles
  #mean/sd: vector of 2
  #p: mixture prop
  p*dnorm(x, mean[1], sd[1]) + (1-p)*dnorm(x, mean[2], sd[2])
}

pnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: probabilities
  #mean/sd: vector of 2
  #p: mixture prop
  p*pnorm(x, mean[1], sd[1]) + (1-p)*pnorm(x, mean[2], sd[2])
}

qnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: quantiles
  #mean/sd: vector of 2
  #p: mixture prop
  p*qnorm(x, mean[1], sd[1]) + (1-p)*qnorm(x, mean[2], sd[2])
}

rnormMix <- function(n, mean=c(0,1), sd=c(1,1), p=1) {
  #n number of obs
  #mean/sd: vector of 2
  #p: mixture prop
  ifelse(sample(0:1,n, replace=T, prob=c(p,1-p)),rnorm(n,mean[1], sd[1]), rnorm(n,mean[2], sd[2]))
}
