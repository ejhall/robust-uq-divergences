library(ggplot2)
library(reshape2)
library(R.matlab)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

dat <- readMat('FIM_screening.mat')
r = c(0.005, 0.045, 0.085)
J = c("$\\theta_1 = \\mu$  ", "$\\theta_2 = \\sigma^2$  ", "$\\theta_3 = \\ell$  ", "$\\theta_4 = \\tau^2$")
J2 = c("$\\mu$  ", "$\\sigma^2$  ", "$\\ell$  ", "$\\tau^2$")
ls = rep(r,3)
ts = sort(rep(r,3))
#df <- data.frame(ell=ls, tau=ts, sensitivity=dat$SW) # normalized
df <- data.frame(ell=ls, tau=ts, sensitivity=dat$RW) # not normalized
names(df)[3:6] <- c("Dmu", "Dsigma", "Dell", "Dtau")

# take logarithmic derivative:
df['Dmu'] <- df['Dmu'] * 0.8
df['Dsigma'] <- df['Dsigma'] * 4.0
df['Dell'] <- df['Dell'] * df['ell']
df['Dtau'] <- df['Dtau'] * df['tau']

m <- melt(df, id.vars=c("ell", "tau"), measure.vars=c("Dmu", "Dsigma", "Dell", "Dtau"), variable.name="index")
m$ell <- factor(m$ell, levels=r, labels=paste(paste("$\\ell ", r, sep="="),"$", sep=""))
m$tau <- factor(m$tau, levels=r, labels=paste(paste("$\\tau^2 ", r, sep="="),"$", sep=""))

myColors <- brewer.pal(4,"OrRd")
names(myColors) <- levels(m$index)
colScale <- scale_colour_manual(name="Screening index $J(i,i)$ for $\\theta_i$:", values=myColors)
filScale <- scale_fill_manual(name="Screening index $J(i,i)$ for $\\theta_i$:", values=myColors, labels=J) 

tikz(file = "FIG3_v3.tex", standAlone=TRUE, width = 6, height = 4)
ggplot(m, aes(index, value)) + geom_col(aes(fill=index), color='black') + 
  facet_grid(ell ~ tau)  +
  scale_x_discrete(labels=J2) +
  theme_bw() + filScale +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio = 0.5,
        legend.key.size = unit(0.5,"cm"), legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8),
        legend.text.align=0, legend.position = "bottom", legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8)) 
endoffile <- dev.off() 

axis.text.x=element_blank(),
axis.ticks.x=element_blank(),

