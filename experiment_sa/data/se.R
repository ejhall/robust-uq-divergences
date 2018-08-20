# GENERATE SAMPLE PATHS
# Generate M sample paths of a lognormal process z on [0,1] with a step size dx = 1/(N-1), where log z 
# has squared exponential (SE or Gaussian) covariance k(r) = v * exp(-r^2/s^2), where var=v, scale=s, 
# with mean mu and nugget variance nv.
# Outputs z an MxN matrix where the ith row contains ith path.

library(RandomFields)
library(ggplot2)
library(reshape)
library(R.matlab)

dpath = "/Users/erichall/Documents/Backup_data/uqii_data_2/"

RFoptions(seed=NA, spConform=TRUE)
# number of samples per bin
M <- 1000
# number of bins
Mb <- 100

# fine mesh
N <- 2^10+1
dx <- 1/(N-1)
x <- seq(0.0, 1.0, dx)
Nc <- 2^6+1
dxc <- 1/(Nc-1)
xc <- seq(0.0, 1.0, dxc)

plab <- c("mu", "v", "s", "nv")
# mu <- 1.0           # trend
v <-  4.0             # sill
# s <-  0.05          # correlation length
# nv <- 0.05          # nugget effect

for (mu in seq(0.8,1.3,.1)){
  for (nv in seq(0.005,0.1,0.01)){
    for(s in seq(0.005,0.1,0.01)){
      pval <- c(mu, v, s, nv)
      # define the geostatistical model for the squared exponential with trend and nugget
      model <- RMgauss(var=v, scale=s) + RMtrend(mean=mu) + RMnugget(var=nv) 
      for (i in 1:Mb){
        # specify file to save to
        atag = paste("a", i, ".mat", sep="")
        locpath = paste("SEsim_", paste(plab,"=",pval,sep="",collapse = "_"), "/", sep="")
        afile = paste(dpath, locpath, atag, sep="")
        # simulate M samples from geostatistical model
        z <- RFsimulate(model=model, x=x, n=M)
        # data <- data.frame(x=x, a=exp(z@data))
        # different view of data for plotting
        # m <- melt(data, id.vars = "x")
        # save samples 
        writeMat(afile, a=t(as.matrix(exp(z@data))), x=t(as.matrix(x)), fixNames=TRUE)
        # if first time samples are being generated for these parameters, then also save various params:
        if(i==1){
          Nc <- 2^6+1
          dxc <- 1/(Nc-1)
          xc <- seq(0.0, 1.0, dxc)
          SIGMA6 <- RFcovmatrix(model, xc)
          Nc <- 2^7+1
          dxc <- 1/(Nc-1)
          xc <- seq(0.0, 1.0, dxc)
          SIGMA7 <- RFcovmatrix(model, xc)
          Nc <- 2^8+1
          dxc <- 1/(Nc-1)
          xc <- seq(0.0, 1.0, dxc)
          SIGMA8 <- RFcovmatrix(model, xc)
          Nc <- 2^9+1
          dxc <- 1/(Nc-1)
          xc <- seq(0.0, 1.0, dxc)
          SIGMA9 <- RFcovmatrix(model, xc)
          Nc <- 2^10+1
          dxc <- 1/(Nc-1)
          xc <- seq(0.0, 1.0, dxc)
          SIGMA0 <- RFcovmatrix(model, xc)
          #perturbations in sigma [sigma = sqrt(v)]
          modeldsig <- RMgauss(var=2*sqrt(v), scale=s)
          SIGMAdsig <- RFcovmatrix(modeldsig, x)
          modeldv <- RMgauss(var=1, scale=s)
          SIGMAdv <- RFcovmatrix(modeldv, x)
          #perturbations in ell [ell = s]
          f <- function(i, j) abs(i - j)**2
          Z <- outer(x, x, f)
          SIGMAdell <- Z*SIGMA0/s**3
          #perturbations in tau [tau = sqrt(nv)]
          modeldtau <- RMnugget(var=2*sqrt(nv))
          SIGMAdtau <- RFcovmatrix(modeldtau, x)
          modeldnv <- RMnugget(var=1)
          SIGMAdnv <- RFcovmatrix(modeldnv, x)
          # write covariance matrices
          paramfile = paste(dpath, locpath, "params.mat", sep="")
          sigsfile = paste(dpath, locpath, "sigs.mat", sep="")
          dsigsfile = paste(dpath, locpath, "dsigs.mat", sep="")
          writeMat(paramfile, SIGMA0=SIGMA0, sig=sqrt(v), v=v, ell=s, s=s, mu=mu, 
                   nv=nv, tau=sqrt(nv), numPerBin=M, nf=N, dxf=dx, fixNames=TRUE)
          writeMat(dsigsfile, SIGMA0=SIGMA0, SIGMAdsig=SIGMAdsig, SIGMAdell=SIGMAdell, 
                   SIGMAdtau=SIGMAdtau, SIGMAdv=SIGMAdv, SIGMAdnv=SIGMAdnv, fixNames=TRUE)
          writeMat(sigsfile, SIGMA6=SIGMA6, SIGMA7=SIGMA7, SIGMA8=SIGMA8, SIGMA9=SIGMA9, 
                   SIGMA0=SIGMA0, fixNames=TRUE)
        }
      }
    }
  }
}
  
# plot
# ggplot(m, aes(x=x, y=value, colour=variable)) + geom_line() +
#   xlab("x") + ylab("conductivity a(x)") + theme(legend.position="none") +
#   #ggtitle(paste("SE covariance with sigma=", sqrt(v), ", ell=", s, ", tau=", round(sqrt(nv),4), ", mu=", mu, sep=""))
#   ggtitle(bquote(Covariance~with~hyperparameters~mu = .(mu)))
#plot(model, xlim=c(0,1), fct.type="Variogram") #additional option dim=2
#Cc <- chol2inv(chol(cmat))
#dc <- diag(Cc)