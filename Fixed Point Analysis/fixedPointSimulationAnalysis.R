rm(list=ls())
graphics.off()

#install.packages("BayesFactor")
library(BayesFactor)

scepath = file.path("C:","Users", "littled", "Dropbox","Work", "2015 Mixture Models", "Fixed Point Analysis")
setwd(scepath)
source(file.path(scepath, "fp.R"))

bandwidth = 20

# Mixture of Serial and Parallel processing
dat <- read.csv("mixedserpar.csv", header=FALSE) # Columnes are RT, condition, subject number
names(dat) = c("rt", "cond", "pp")
dat1 <- dat[dat$pp == 1, ]

## compute the list of fp1 objects
res <- tapply(1:nrow(dat), dat$pp, function(X) {fpGet(dat[X,], 1000, bw=bandwidth)})
res1 <- tapply(1:nrow(dat1), dat1$pp, function(X) {fpGet(dat1[X,], 1000, bw=bandwidth)}) # Just plot one subject as an example
# You can access different "fields" using res[[1]]@dens, @dat, @diff

## get the crossing points and plot for comparison
crosses=fpDensDiff(res)
boxplot(t(crosses), frame.plot=FALSE,xlab="Crossing point", ylab="Condition pair",names=c("1-2","2-3","1-3"), horizontal=TRUE)

## call fpAnova, with stat="both" to do both a Bayesian and a frequentist test
bf <- fpAnova(res, stat="both")
bf$BF
bf$p

## compute one fp object for plotting
fpobject <- fpGet(dat1, 1000, bw=bandwidth)

## plot it
par(mfrow=c(1,2))
plot(fpobject)

## alternatively plot as single figures
#plot.fp1(fpobject)
#fpPlot(fpobject)

########################################
# Inhibitory parallel model
# Mixture of Serial and Parallel processing
idat <- read.csv("inhibitory.csv", header=FALSE) # Columnes are RT, condition, subject number
names(dat) = c("rt", "cond", "pp")
idat1 <- dat[dat$pp == 1, ]

## compute the list of fp1 objects
ires <- tapply(1:nrow(idat), idat$pp, function(X) {fpGet(idat[X,], 1000, bw=bandwidth)})
ires1 <- tapply(1:nrow(idat1), idat1$pp, function(X) {fpGet(idat1[X,], 1000, bw=bandwidth)}) # Just plot one subject as an example
# You can access different "fields" using res[[1]]@dens, @dat, @diff

## get the crossing points and plot for comparison
icrosses=fpDensDiff(ires)
boxplot(t(icrosses), frame.plot=FALSE,xlab="Crossing point", ylab="Condition pair",names=c("1-2","2-3","1-3"), horizontal=TRUE)

## call fpAnova, with stat="both" to do both a Bayesian and a frequentist test
ibf <- fpAnova(ires, stat="both")
ibf$BF
ibf$p

## compute one fp object for plotting
ifpobject <- fpGet(idat1, 1000, bw=bandwidth)

## plot it
par(mfrow=c(1,2))
plot(ifpobject)

## alternatively plot as single figures
#plot.fp1(fpobject)
#fpPlot(fpobject)