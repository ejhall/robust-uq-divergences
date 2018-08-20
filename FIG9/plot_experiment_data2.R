library(R.matlab)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

dat <- readMat('20170731_idea2_outv6.mat')
dat <- data.frame(dat)
names(dat)[1] <- paste("nom")
names(dat)[2] <- paste("categ")
names(dat)[3] <- paste("obsv")
names(dat)[4] <- paste("alt")
names(dat)[5] <- paste("value")

nomlab = c("$P_{1}$", "$P_{2}$", "$P_{3}$")
obsvlab = c("$g_1(\\bar{u})$", "$g_2(\\bar{u})$", "$g_3(\\bar{u})$")
altlab = c("$Q_{\\mathrm{max}}$","$Q_{\\mathrm{min}}$","$Q_{\\mathrm{mean}}$","$Q_{\\mathrm{med}}$")
categlab = c("$\\Xi_{+}$", "$\\mathcal{E}$", "$\\Xi_{-}$")

dat$nom <- factor(dat$nom, levels=c(1,2,3), labels=nomlab)
dat$categ <- factor(dat$categ, levels=c(2,1,3), labels=categlab)
dat$obsv <- factor(dat$obsv, labels= obsvlab)
dat$alt <- factor(dat$alt, labels=altlab)

myColors <- brewer.pal(3,"Set1")
myColors <- myColors[c(2,1,3)]
names(myColors) <- levels(dat$categ)
colScale <- scale_colour_manual(name="Quantity", values=myColors, labels=categlab)
filScale <- scale_fill_manual(name="Quantity", values=myColors, labels=categlab) 

stat1 <- ddply(dat, .(nom, categ, obsv, alt), summarize, value=mean(value))
tikz(file = "FIG8v2.tex", standAlone=TRUE, width=6, height=6)
ggplot(dat, aes(alt, value)) + 
  geom_line(dat=stat1, aes(y=value, group=categ, linetype=categ, color=categ), alpha=0.6, size=1.2) + 
  geom_boxplot(aes(color=categ), outlier.size = 1, outlier.shape = 1, outlier.alpha = 0.5, position=position_dodge(0), 
               show.legend = FALSE) + 
  facet_grid(obsv ~ nom, scales="free_y") + xlab("alternative models") +
  colScale + theme_bw() + scale_linetype(name="Quantity", labels=categlab) +
  ylab("weak error") + 
  theme(aspect.ratio=0.5, 
        legend.text.align=0, legend.position = c(0.1615,0.268), legend.direction = "horizontal",
        legend.key.size = unit(0.39,"cm"),
        legend.margin = margin(c(0.03,0.03,0.03,0.03), unit="cm"),
        legend.title=element_text(size=8), 
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8), axis.title=element_text(size=8)) 
endoffile <- dev.off() 
