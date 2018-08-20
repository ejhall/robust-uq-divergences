library(R.matlab)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

#setwd("~/Research/hybrid_div_rpde/code_uq_rpde/experiment3_data_idea1/idea1_Rplots")
dat <- readMat('20170731_m5e5_mskp1e3_Pmod5_g124_outv6_R.mat')
dat <- data.frame(dat)
names(dat)[1] <- paste("categ")
names(dat)[2] <- paste("obsv")
names(dat)[3] <- paste("alt")
names(dat)[4] <- paste("value")

obsvlab = c("$g_1(\\bar{u})$", "$g_2(\\bar{u})$", "$g_3(\\bar{u})$")
altlab = c("$Q_{10}$","$Q_{20}$","$Q_{30}$","$Q_{40}$","$Q_{60}$",
           "$Q_{70}$","$Q_{80}$","$Q_{90}$","$Q_{100}$")

categlab = c("$\\Xi_{+}$", "$\\mathcal{E}$", "$\\Xi_{-}$")
dat$categ <- factor(dat$categ, levels=c(2,1,3), labels=c("XIp", "WE", "XIm"))
dat$obsv <- factor(dat$obsv, labels= obsvlab)
dat$alt <- factor(dat$alt, labels=altlab)

myColors <- brewer.pal(3,"Set1")
myColors <- myColors[c(2,1,3)]
names(myColors) <- levels(dat$categ)
colScale <- scale_colour_manual(name="Quantity", values=myColors, labels=categlab)
filScale <- scale_fill_manual(name="Quantity", values=myColors, labels=categlab) 

dat1 <- subset(dat, obsv!=paste(obsvlab[3])) #
stat1 <- ddply(dat1, .(categ, obsv, alt), summarize, value=mean(value))

tikz(file = "FIG6v2.tex", standAlone=TRUE, width = 6, height = 3)
ggplot(dat1, aes(alt, value)) + 
  geom_line(dat=stat1, aes(y=value, group=categ, linetype=categ, color=categ), alpha=0.6, size=1.2) + 
  geom_boxplot(aes(color=categ), outlier.size = 1, outlier.shape = 1, outlier.alpha = 0.5, position=position_dodge(0), 
               show.legend = FALSE) + 
  facet_grid(.~ obsv, scales="free_y") + xlab("alternative models") +
  colScale + theme_bw() + scale_linetype(name="Quantity", labels=categlab) +
  ylab("weak error") + xlab("alternative models") +
  theme(legend.direction = "horizontal",
        legend.key.size = unit(0.5,"cm"), legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8),
        aspect.ratio=0.5, legend.text.align=0, legend.position = c(0.825,0.875),
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        axis.text = element_text(size=8), axis.title=element_text(size=8)) # 
endoffile <- dev.off() 

dat2 <- subset(dat, obsv==paste(obsvlab[3]))
stat2 <- ddply(dat2, .(categ, obsv, alt), summarize, value=mean(value))
tikz(file = "FIG7v2.tex", standAlone=TRUE, width = 3, height = 3)
ggplot(dat2, aes(alt, value)) +  
  geom_line(dat=stat2, aes(y=value, group=categ, linetype=categ, color=categ), alpha=0.6, size=1.2) +
  # geom_line(stat="smooth", method = "loess", se = FALSE, 
  #           linetype='dashed', alpha=0.5, data=stat2, aes(group=categ, color=categ)) + 
  geom_boxplot(aes(color=categ), outlier.size = 1, outlier.shape = 1, outlier.alpha = 0.5, position=position_dodge(0),
               show.legend = FALSE) + 
  facet_grid(.~obsv, scales="free_y") + xlab("alternative models") +
  colScale + theme_bw() + scale_linetype(name="Quantity", labels=categlab) +
  ylab("weak error") + xlab("alternative models") +
  theme(legend.direction = "horizontal",
        legend.key.size = unit(0.5,"cm"), legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8),
        aspect.ratio=0.5, legend.text.align=0, legend.position = c(0.635,0.875),
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        axis.text = element_text(size=8), axis.title=element_text(size=8))
endoffile <- dev.off() 
