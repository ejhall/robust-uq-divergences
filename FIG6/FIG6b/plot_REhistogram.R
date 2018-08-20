library(ggplot2)
library(RColorBrewer)
library(R.matlab)
library(scales)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[12pt]{article}")

dat <- readMat('REdat_outv6_R.mat')
dat <- data.frame(dat)
names(dat)[1] <- paste("nom")
names(dat)[2] <- paste("RE")

nomlab = c("$P_{1}$", "$P_{2}$", "$P_{3}$")
dat$nom <- factor(dat$nom, levels=c(1,2,3), labels=nomlab)
legtxt <- "Nominal"

myColors <- brewer.pal(6,"Set1")
myColors <- myColors[3:5]
names(myColors) <- levels(dat$nom)
colScale <- scale_colour_manual(name=legtxt, values=myColors, labels=nomlab)
filScale <- scale_fill_manual(name=legtxt, values=myColors, labels=nomlab) 

tikz(file = "FIG5b.tex", standAlone=TRUE, width=4, height=3)
ggplot(dat, aes(RE)) + 
  geom_freqpoly(aes(color=nom, linetype=nom), bins=50, alpha=0.8, size=1.2) + #bins or binwidth
  #can switch to density by adding y=..density.., to aes for reqpoly
  scale_linetype_discrete(name=legtxt, labels=nomlab) + colScale + theme_bw() + 
  #scale_y_continuous(labels=scales::scientific_format(digits=1)) +
  scale_y_continuous(breaks=c(0,20000,40000), labels=c("$0$", "$2 \\times 10^4$", "$4 \\times 10^4$")) +
  xlab("relative entropy") + ylab("frequency") +
  theme(legend.text.align=0, legend.position = c(0.735,0.9), aspect.ratio=.5,
        legend.direction = "horizontal",
        legend.key.size = unit(0.5,"cm"),
        legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8), 
      legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
      axis.text = element_text(size=8))
endoffile <- dev.off() 
