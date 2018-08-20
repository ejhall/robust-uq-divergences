library(R.matlab)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{article}")

perturblab <- c('$\\epsilon(\\ell)$, perturbation in $\\ell$', '$\\epsilon(\\tau^2)$, perturbation in $\\tau^2$')
categlab <- c("$\\mathcal{E}$", "$\\hat{\\Delta}$", "$\\Xi_{\\pm}$", "$B_{\\pm}$", "$C_{\\pm}$")
observlab <- c('$g_1(\\bar{u})$', '$g_2(\\bar{u})$', '$g_3(\\bar{u})$')
quantdisplaylab <- c("$\\mathcal{E}$", "$\\hat{\\Delta}$", "$\\Xi_{+}$", "$\\Xi_{-}$", 
                     "$B_{+}$", "$B_{-}$", "$C_{+}$", "$C_{-}$")

dat1 <- readMat('20170728_mod1a.mat')
df1 <- data.frame(perturblab[1], dat1[[5]])
dat2 <- readMat('20170728_mod1b.mat')
df2 <- data.frame(perturblab[2], dat2[[5]])

names(df1) <- c('perturb', 'quant', 'observ', 'epsilon', 're', 'cstar', 'value')
names(df2) <- c('perturb', 'quant', 'observ', 'epsilon', 're', 'cstar', 'value')
df <- rbind(df1, df2)
quantlab <- unlist(dat1[[4]])
df$quant <- factor(df$quant, labels=quantlab)
df$observ <- factor(df$observ, labels=observlab)
df <- ddply(df, .(perturb, quant, observ, epsilon, re, cstar, value), function(x) 
  c(pm={if(grepl("p", x$quant)){"p"} else if(grepl("m", x$quant)){"m"} else{"c"}},
  categ={if(grepl("WE", x$quant)){"WE"} else if(grepl("SI", x$quant)){"SI"} 
    else if(grepl("XI", x$quant)){"XI"} else if(grepl("BI", x$quant)){"BI"} 
    else if(grepl("CI", x$quant)){"CI"} else {"OO"}}))
df$categ <- factor(df$categ, levels=c("WE","SI","XI","BI","CI"))
rm(df1,df2,dat1,dat2)

dfstat <- ddply(df, .(perturb, quant, categ, pm, observ, epsilon), function(x) 
  c(value=mean(x$value), errmin=mean(x$value) - 2*sd(x$value), errmax=mean(x$value) + 2*sd(x$value), 
    re=mean(x$re), cstar=mean(x$cstar)))
dfstat <- subset(dfstat, categ!="SI")
dfstat$value <- as.numeric(dfstat$value)
dfstat$re <- as.numeric(dfstat$re)
dfstat$cstar <- as.numeric(dfstat$cstar)
dfstat$errmax <- as.numeric(dfstat$errmax)
dfstat$errmin <- as.numeric(dfstat$errmin)
dfstat$categ <- factor(dfstat$categ, levels=c("WE", "XI","BI","CI"))
myColors <- brewer.pal(4,"Set1")
names(myColors) <- levels(dfstat$categ)

tikz(file = "FIG9av2.tex", standAlone=TRUE, width = 3, height = 2)
ggplot(subset(dfstat, observ==observlab[3] & perturb==perturblab[1]), aes(re)) +
  geom_line(aes(y=value, group=quant, color=categ, linetype=categ), size=1.2) +
  geom_point(aes(y=value, color=categ, shape=categ), size=2, alpha=0.9) +
  geom_ribbon(aes(ymin=errmin, ymax=errmax, group=quant, fill=categ), alpha=0.2) +
  facet_grid(observ ~ perturb, scales="free") +
  ylab("weak error") + 
  xlab("$\\bar{\\mathcal{R}}$") + #scale_x_log10() +
  scale_colour_manual(name="Quantity", values=myColors, labels=categlab[c(1,3:5)]) +
  scale_fill_manual(name="Quantity", values=myColors, labels=categlab[c(1,3:5)])  +
  scale_shape(name="Quantity", labels=categlab[c(1,3:5)], solid=FALSE) +
  theme_bw() +
  theme(legend.position="none", aspect.ratio = 0.5,
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8), axis.title=element_text(size=8))
endoffile <- dev.off()

tikz(file = "FIG9bv2.tex", standAlone=TRUE, width = 3, height = 2)
ggplot(subset(dfstat, observ==observlab[3] & perturb==perturblab[2]), aes(re)) +
  geom_line(aes(y=value, group=quant, color=categ, linetype=categ), size=1.2) +
  geom_point(aes(y=value, color=categ, shape=categ), size=2, alpha=0.9) +
  geom_ribbon(aes(ymin=errmin, ymax=errmax, group=quant, fill=categ), alpha=0.2) +
  facet_grid(observ ~ perturb, scales="free") +
  ylab("weak error") + 
  xlab("$\\bar{\\mathcal{R}}$") + #scale_x_log10() +
  scale_colour_manual(name="Quantity", values=myColors, labels=categlab[c(1,3:5)]) +
  scale_fill_manual(name="Quantity", values=myColors, labels=categlab[c(1,3:5)])  +
  scale_shape(name="Quantity", labels=categlab[c(1,3:5)], solid=FALSE) +
  scale_linetype(name="Quantity", labels=categlab[c(1,3:5)]) + 
  theme_bw() + ylim(-0.55,0.4) +
  theme(aspect.ratio = 0.5, legend.direction = "horizontal",
        legend.position = c(0.5,0.13), legend.key.size = unit(0.5,"cm"), 
        legend.margin = margin(c(0.01,0.01,0.01,0.01), unit="cm"),
        legend.title=element_text(size=8), legend.text.align=0,
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8), axis.title=element_text(size=8))
endoffile <- dev.off()
