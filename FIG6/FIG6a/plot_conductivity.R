library(ggplot2)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

load("conduct_avbvc.R")

conduct <- cbind(categ=rep(0,220), y, val=Kx[x,y,z])
excP1 <- cbind(categ=rep(1,66), y=yAlta, val=rep(-20,66))
excP2 <- cbind(categ=rep(2,66), y=yAltb, val=rep(-35,66))
excP3 <- cbind(categ=rep(3,66), y=yAltc, val=rep(-50,66))

categlab <- c("conductivity", "$P_{1}$", "$P_{2}$", "$P_{3}$")
dat <- data.frame(rbind(conduct, excP1, excP2, excP3))
dat$categ <- factor(dat$categ, levels=c(0,1,2,3), labels=categlab)
legtxt <- "Data excluded"

myColors <- brewer.pal(6,"Set1")
myColors <- myColors[3:5]
names(myColors) <- levels(dat$categ)[2:4]
colScale <- scale_colour_manual(name=legtxt, values=myColors, labels=categlab[2:4])
filScale <- scale_fill_manual(name=legtxt, values=myColors, labels=categlab[2:4]) 

formatMeter<- function(x){ x*3.048 }

dat_cond <- subset(dat, categ==paste(categlab[1]))
dat_exc <- subset(dat, categ!=paste(categlab[1]))
tikz(file = "FIG5av2.tex", standAlone=TRUE, width=4, height=3)
ggplot(dat, aes(y)) +
  geom_hline(yintercept = 0, color="red", size=0.5) +
  geom_step(data=dat_cond,aes(y=val, group=categ), size=1) +
  geom_point(data=dat_exc, aes(y=val, color=categ, shape=categ), size=1) +
  colScale + theme_bw() + 
  scale_x_continuous(labels=formatMeter) +
  scale_y_continuous(breaks = seq(0,200,50)) +
  xlab("spatial $x$ ($\\mathrm{m}$)") + ylab("conductivity $a(x)$") + 
  theme(legend.text.align=0, legend.position = c(0.325,0.9), 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5,"cm"),
        legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8), 
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"), 
        aspect.ratio=0.5, 
        axis.text = element_text(size=8)) +
  scale_shape(name=legtxt, labels=categlab[2:4], solid=FALSE) 
endoffile <- dev.off() 
