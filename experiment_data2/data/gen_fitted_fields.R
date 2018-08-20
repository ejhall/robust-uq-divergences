# Load data from SPE10.
# Fit covariance model for 1-D slice from permeability data top layer, edge 
# (i.e. create empiracle variogram, call RFfit).
# Also fit covariance model to 1-D slide to 'lacking' data. 
# Generate gaussian paths (independently) in 100 batches of 1000 sample paths.
# #

library(R.matlab)
library(RandomFields)
library(ggplot2)
library(reshape2)
library(corpcor)
library(parallel)

# options for RandomFields
RFoptions(seed=NA, spConform=TRUE)

# number RE draws
numQs = 1e5
# number of samples per bin
M <- 1000
# number of bins
Mb <- 100


# Fix: top layer of sequence, 100 m in from edge.
data <- readMat("/Users/erichall/Research/hybrid_div_rpde/code_up_rpde/SPEcase2a/spe10data.mat")
dpath = "/Users/erichall/Research/hybrid_div_rpde/backup_data/experiment3/run20170608_c/"
x = 6
z = 1
# 1-D slice of 
y <- c(1:220)
d <- length(y)
pc <- seq(0.1,0.9,0.1)*220
daty <- data$y
Kx <- data$Kx
estmodel.gauss <- RMgauss(var=NA, scale=NA) + RMtrend(mean=NA) + RMnugget(var=0.05) 

# # nominal model P is a random sample of 70% of data
# yP <- sort(sample(y, pc[7], replace=FALSE, prob=NULL))
# # plot sample points
# plot(yP, data$Kx[x,yP,z], 'l', xlab='x', ylab='a(x)')
# points(yP, rep(0, 154), col="red")
# # fit nominal model
# P <-  RFfit(estmodel.gauss, data=data$Kx[x,yP,z], 
#             lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), 
#                        list("$", var=0.0001, list("RMnugget"))),
#             upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), 
#                        list("$", var=0.5, list("RMnugget"))), loggaus=TRUE)
# SP <- RFcovmatrix(P, data$y)
# tmp <-P@table['ml']
# MP <- rep(tmp[22,], length(data$y))
# SPinv <- pseudoinverse(SP)
yAlt <- setdiff(y, yP)

## Sampling RE
cat("Running parallel to obtain (RE, yQ)")
# Calculate the number of cores
no_cores <- detectCores() %/% 2
# Initiate cluster
cl <- makeCluster(no_cores)
# make global variables accessible
clusterExport(cl, c('estmodel.gauss', 'yAlt', 'yP', 'x', 'daty', 'z', 'SP', 'SPinv', 'MP', 'd', 'Kx', 'RFfit', 'RFcovmatrix'))
# parallel apply and prduce matrix formatted: c(RE, yQ) by numQs
runs <- parSapply(cl, as.character(1:numQs),
                  function(i){
                    yQ <- sort(union(sample(yAlt, 22, replace=FALSE, prob=NULL), yP))
                    Q <- RFfit(estmodel.gauss, data=Kx[x,yQ,z],
                               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0),
                                          list("$", var=0.0001, list("RMnugget"))),
                               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20),
                                          list("$", var=0.5, list("RMnugget"))), loggaus=TRUE)
                    SQ <- RFcovmatrix(Q, daty)
                    tmp <-Q@table['ml']
                    MQ <- rep(tmp[22,], length(daty))
                    r <- 0.5 * (log(det(SP)) - log(det(SQ)) + sum(diag(SPinv %*% SQ)) + t(MQ - MP) %*% SPinv %*% (MQ - MP) - d)
                    c(RE=r[1], y=t(yQ))
                  })

stopCluster(cl)

RE <- runs[1,]
hist(RE, breaks = pretty(RE, n=100), freq = FALSE)
hist(RE, breaks = pretty(RE, n=50), freq = FALSE)

plot(y, Kx[x,y,z], 's', xlab='x', ylab='a(x)', ylim=c(-10, 250))
points(yAlt, rep(-10,66), col='green')
# points(yAltb, rep(-20,66), col='blue')
# points(yAltc, rep(-30,66), col='red')

imax <- which.max(RE)
ymax <- runs[2:177, imax]
Qmax <- RFfit(estmodel.gauss, data=data$Kx[x,ymax,z], 
              lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), 
                         list("$", var=0.0001, list("RMnugget"))),
              upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), 
                         list("$", var=0.5, list("RMnugget"))),
              loggaus=TRUE)

imin <- which.min(RE)
ymin <- runs[2:177, imin]
Qmin <- RFfit(estmodel.gauss, data=data$Kx[x,ymin,z], 
              lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), 
                         list("$", var=0.0001, list("RMnugget"))),
              upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), 
                         list("$", var=0.5, list("RMnugget"))),
              loggaus=TRUE)

imed <- which.min(abs(runs[1,]-median(runs[1,])))
ymed <- runs[2:177, imed]
Qmed <- RFfit(estmodel.gauss, data=data$Kx[x,ymed,z], 
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), 
                          list("$", var=0.0001, list("RMnugget"))),
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), 
                          list("$", var=0.5, list("RMnugget"))),
               loggaus=TRUE)

imean <- which.min(abs(runs[1,]-mean(runs[1,])))
ymean <- runs[2:177, imean]
Qmean <- RFfit(estmodel.gauss, data=data$Kx[x,ymean,z], 
              lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), 
                         list("$", var=0.0001, list("RMnugget"))),
              upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), 
                         list("$", var=0.5, list("RMnugget"))),
              loggaus=TRUE)

writeMat(paste(dpath, "RE.mat", sep=""), RE=t(as.matrix(runs[1,])), yP=t(as.matrix(yP)), 
         yQ=t(as.matrix(runs[2:177,])), 
         imax=t(as.matrix(imax)), imin=t(as.matrix(imin)),
         imed=t(as.matrix(imed)), imean=t(as.matrix(imean)), fixNames=TRUE)

#

locpath = paste(c("P", "Qmax", "Qmin", "Qmed", "Qmean"), '/', sep="")
foreach(i=1:Mb, .packages='RandomFields') %do% { # can RandomFields be called in parallel? Better not to be safe...
  # specify file to save to
  atag <- paste("a", i, ".mat", sep="")
  afile <- paste(paste(dpath, "Gauss_", sep=""), locpath, atag, sep="")
  # simulate M samples from geostatistical model
  z.P <- RFsimulate(P, x=data$y,  n=M)
  z.Qmax <- RFsimulate(Qmax, x=data$y, n=M)
  z.Qmin <- RFsimulate(Qmin, x=data$y, n=M)
  z.Qmed <- RFsimulate(Qmed, x=data$y, n=M)
  z.Qmean <- RFsimulate(Qmean, x=data$y, n=M)
  # save samples
  writeMat(afile[grep('P', afile, value=FALSE)], a=t(as.matrix(exp(z.P@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[grep('Qmax', afile, value=FALSE)], a=t(as.matrix(exp(z.Qmax@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[grep('Qmin', afile, value=FALSE)], a=t(as.matrix(exp(z.Qmin@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[grep('Qmed', afile, value=FALSE)], a=t(as.matrix(exp(z.Qmed@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[grep('Qmean', afile, value=FALSE)], a=t(as.matrix(exp(z.Qmean@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
}