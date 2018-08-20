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

dpath = "/Users/erichall/Research/hybrid_div_rpde/backup_data/experiment3/run20170602/"
# number of samples per bin
M <- 1000
# number of bins
Mb <- 200

RFoptions(seed=NA, spConform=TRUE)
data <- readMat("/Users/erichall/Research/hybrid_div_rpde/code_uq_rpde/SPEcase2a/spe10data.mat")

# Fix: top layer of sequence, 100 m in from edge.
idx = 6
idz = 1
# 1-D slice of 
iy <- c(1:220)
pc <- seq(0.1,0.9,0.1)*220
iy9 <- sort(sample(iy,  pc[9], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 90% of full data
iy8 <- sort(sample(iy9, pc[8], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 80% of full data
iy7 <- sort(sample(iy8, pc[7], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 70% of full data
iy6 <- sort(sample(iy7, pc[6], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 60% of full data
iy5 <- sort(sample(iy6, pc[5], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 50% of full data
iy4 <- sort(sample(iy5, pc[4], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 40% of full data
iy3 <- sort(sample(iy4, pc[3], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 30% of full data
iy2 <- sort(sample(iy3, pc[2], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 20% of full data
iy1 <- sort(sample(iy2, pc[1], replace=FALSE, prob=NULL)) # take a (uniform) random sample of only 10% of full data
writeMat(paste(dpath, "yis.mat", sep=""), iy9=t(iy9), iy8=t(iy8), iy7=t(iy7), iy6=t(iy6), iy5=t(iy5), iy4=t(iy4), 
         iy3=t(iy3), iy2=t(iy2), iy1=t(iy1), fixNames=TRUE)

emp.vario <- RFempiricalvariogram(x=data$celly, data=data$Kx[idx,1:220,idz])

estmodel.gauss <- RMgauss(var=NA, scale=NA) + RMtrend(mean=NA) + RMnugget(var=0.05) 

emp.10 <- RFfit(estmodel.gauss, data=data$Kx[idx,1:220,idz],
                lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
                upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
                loggaus=TRUE)
emp.9 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy9,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.8 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy8,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.7 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy7,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.6 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy6,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.5 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy5,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.4 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy4,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.3 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy3,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.2 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy2,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)
emp.1 <- RFfit(estmodel.gauss, data=data$Kx[idx,iy1,idz],
               lower=list("+", list("$", var=1, scale=1, list("RMgauss")), list("RMtrend", mean=0), list("$", var=0.0001, list("RMnugget"))), 
               upper=list("+", list("$", var=6, scale=50, list("RMgauss")), list("RMtrend", mean=20), list("$", var=0.5, list("RMnugget"))), 
               loggaus=TRUE)

# plot(emp.vario)
# plot(emp.10, ylim=c(0,11), xlim=c(0,90), main="Gauss model [100%]")
# plot(emp.9,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [90%]")
# plot(emp.8,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [80%]")
# plot(emp.7,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [70%]")
# plot(emp.6,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [60%]")
# plot(emp.5,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [50%]")
# plot(emp.4,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [40%]")
# plot(emp.3,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [30%]")
# plot(emp.2,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [20%]")
# plot(emp.1,  ylim=c(0,11), xlim=c(0,90), main="Gauss model [10%]")

locpath = paste("EmpGauss_", seq(1,10,1), '/', sep="")
for (i in 201:500){
  # specify file to save to
  atag = paste("a", i, ".mat", sep="")
  afile = paste(dpath, locpath, atag, sep="")
  # simulate M samples from geostatistical model
  z.10 <- RFsimulate(emp.10, x=data$y,  n=M)
  z.9  <- RFsimulate(emp.9,  x=data$y, n=M)
  z.8  <- RFsimulate(emp.8,  x=data$y, n=M)
  z.7  <- RFsimulate(emp.7,  x=data$y, n=M)
  z.6  <- RFsimulate(emp.6,  x=data$y, n=M)
  z.5  <- RFsimulate(emp.5,  x=data$y, n=M)
  z.4  <- RFsimulate(emp.4,  x=data$y, n=M)
  z.3  <- RFsimulate(emp.3,  x=data$y, n=M)
  z.2  <- RFsimulate(emp.2,  x=data$y, n=M)
  z.1  <- RFsimulate(emp.1,  x=data$y, n=M)
  # save samples
  writeMat(afile[10], a=t(as.matrix(exp(z.10@data))), x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[9],  a=t(as.matrix(exp(z.9@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[8],  a=t(as.matrix(exp(z.8@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[7],  a=t(as.matrix(exp(z.7@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[6],  a=t(as.matrix(exp(z.6@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[5],  a=t(as.matrix(exp(z.5@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[4],  a=t(as.matrix(exp(z.4@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[3],  a=t(as.matrix(exp(z.3@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[2],  a=t(as.matrix(exp(z.2@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  writeMat(afile[1],  a=t(as.matrix(exp(z.1@data))),  x=t(as.matrix(data$y)), fixNames=TRUE)
  if (i==101){
    #elem <- seq(10,2190,20)
    SIGMA10 <- RFcovmatrix(emp.10, data$y)
    SIGMA9  <- RFcovmatrix(emp.9, data$y)
    SIGMA8  <- RFcovmatrix(emp.8, data$y)
    SIGMA7  <- RFcovmatrix(emp.7, data$y)
    SIGMA6  <- RFcovmatrix(emp.6, data$y)
    SIGMA5  <- RFcovmatrix(emp.5, data$y)
    SIGMA4  <- RFcovmatrix(emp.4, data$y)
    SIGMA3  <- RFcovmatrix(emp.3, data$y)
    SIGMA2  <- RFcovmatrix(emp.2, data$y)
    SIGMA1  <- RFcovmatrix(emp.1, data$y)
    tmp <-emp.10@table['ml']
    mu10 <- rep(tmp[22,], length(data$y))
    tmp <-emp.9@table['ml']
    mu9 <- rep(tmp[22,], length(data$y))
    tmp <-emp.8@table['ml']
    mu8 <- rep(tmp[22,], length(data$y))
    tmp <-emp.7@table['ml']
    mu7 <- rep(tmp[22,], length(data$y))
    tmp <-emp.6@table['ml']
    mu6 <- rep(tmp[22,], length(data$y))
    tmp <-emp.5@table['ml']
    mu5 <- rep(tmp[22,], length(data$y))
    tmp <-emp.4@table['ml']
    mu4 <- rep(tmp[22,], length(data$y))
    tmp <-emp.3@table['ml']
    mu3 <- rep(tmp[22,], length(data$y))
    tmp <-emp.2@table['ml']
    mu2 <- rep(tmp[22,], length(data$y))
    tmp <-emp.1@table['ml']
    mu1 <- rep(tmp[22,], length(data$y))
    writeMat(paste(dpath, "params2.mat", sep=""), SIGMA10=as.matrix(SIGMA10),  
             SIGMA9=as.matrix(SIGMA9), SIGMA8=as.matrix(SIGMA8), SIGMA7=as.matrix(SIGMA7), 
             SIGMA6=as.matrix(SIGMA6), SIGMA5=as.matrix(SIGMA5), SIGMA4=as.matrix(SIGMA4), 
             SIGMA3=as.matrix(SIGMA3), SIGMA2=as.matrix(SIGMA2), SIGMA1=as.matrix(SIGMA1), 
             mu10=as.matrix(mu10), mu9=as.matrix(mu9), mu8=as.matrix(mu8), mu7=as.matrix(mu7), 
             mu6=as.matrix(mu6), mu5=as.matrix(mu5), mu4=as.matrix(mu4), mu3=as.matrix(mu3), 
             mu2=as.matrix(mu2), mu1=as.matrix(mu1), fixNames=TRUE)
  }
}
 
   
# a.10 <- data.frame(x=data$y, a=exp(z.10@data))
# a.9 <- data.frame(x=data$y, a=exp(z.9@data))
# a.8 <- data.frame(x=data$y, a=exp(z.8@data))
# a.7 <- data.frame(x=data$y, a=exp(z.7@data))
# a.6 <- data.frame(x=data$y, a=exp(z.6@data))
# a.5 <- data.frame(x=data$y, a=exp(z.5@data))
# a.4 <- data.frame(x=data$y, a=exp(z.4@data))
# a.3 <- data.frame(x=data$y, a=exp(z.3@data))
# a.2 <- data.frame(x=data$y, a=exp(z.2@data))
# a.1 <- data.frame(x=data$y, a=exp(z.1@data))
# 
# m.10 <- melt(a.10, id.vars = "x")
# m.9 <- melt(a.9, id.vars = "x")
# m.8 <- melt(a.8, id.vars = "x")
# m.7 <- melt(a.7, id.vars = "x")
# m.6 <- melt(a.6, id.vars = "x")
# m.5 <- melt(a.5, id.vars = "x")
# m.4 <- melt(a.4, id.vars = "x")
# m.3 <- melt(a.3, id.vars = "x")
# m.2 <- melt(a.2, id.vars = "x")
# m.1 <- melt(a.1, id.vars = "x")
# 
# ggplot(m.10, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [100%], M=4 samples")
# ggplot(m.9, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [90%], M=4 samples")
# ggplot(m.8, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [80%], M=4 samples")
# ggplot(m.7, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [70%], M=4 samples")
# ggplot(m.6, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [60%], M=4 samples")
# ggplot(m.5, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [50%], M=4 samples")
# ggplot(m.4, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [40%], M=4 samples")
# ggplot(m.3, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [30%], M=4 samples")
# ggplot(m.2, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [20%], M=4 samples")
# ggplot(m.1, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   ggtitle("Fitted Gauss model [10%], M=4 samples")