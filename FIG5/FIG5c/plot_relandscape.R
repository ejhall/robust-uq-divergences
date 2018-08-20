library(R.matlab)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

#setwd("~/Research/hybrid_div_rpde/code_uq_rpde/experiment3_data_idea1/idea1_Rplots/")
dat <- readMat('s_nv_re.mat')
dat <- data.frame(dat)
vars = c("$\\epsilon(\\ell)$", "$\\epsilon(\\tau^2)$", "$\bar{\\mathcal{R}}$")
names(dat)[1] <- paste("s")
names(dat)[2] <- paste("nv")
names(dat)[3] <- paste("re")

pt <- data.frame(s=0.0,nv=0.0,val="Nominal")

tikz(file = "FIG4c.tex", standAlone=TRUE, width=2.5, height=2.5)
ggplot(dat, aes(nv,s)) + geom_raster(aes(fill=re), interpolate = TRUE, hjust=0.5, vjust=0.5) + 
  scale_fill_distiller(palette = "Spectral", name="$\\bar{\\mathcal{R}}$") +  
  coord_cartesian(xlim = c(-0.04,0.05), ylim = c(0.0,0.09)) +
  geom_point(data=pt, aes(shape=val), size=2) + 
  geom_text(data=pt, aes(label=val), hjust=0, vjust=0.5, nudge_x = .004) +
  xlab("$\\epsilon(\\tau^2)$") + ylab("$\\epsilon(\\ell)$") +
  theme_bw() + guides(shape=FALSE) +
  theme(panel.grid = element_blank(), legend.text.align=1,
        legend.key.width = unit(0.3,"in"), legend.key.height = unit(.1, "in"), 
        legend.margin = margin(c(0,0,-0.1,0), unit="in"),
        legend.title=element_text(size=10), aspect.ratio = 1, legend.position = "top", legend.justification = 1)
endoffile <- dev.off() 
