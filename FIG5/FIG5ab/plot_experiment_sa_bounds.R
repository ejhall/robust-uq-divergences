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

# plot of re vs quantity
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

theta <- dfstat$perturb
theta <- ifelse(theta==perturblab[1], 0.005, 0.045)
dfstat$value <- theta * dfstat$value / dfstat$epsilon
dfstat$errmin <- theta * dfstat$errmin / dfstat$epsilon
dfstat$errmax <- theta * dfstat$errmax / dfstat$epsilon
thetalab = c("$\\hat{\\Delta}$", "$\\frac{\\theta}{\\epsilon} \\, \\Xi_{\\pm}$")
legstyle <- theme(legend.direction = "horizontal", legend.position = "bottom", 
                  legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"))

tikz(file = "FIG4a.tex", standAlone=TRUE, width=2.25, height=3.4)
ggplot(subset(dfstat, categ %in% c("WE", "XI") & perturb==perturblab[1]), aes(epsilon)) + 
  geom_line(aes(y=value, group=quant, color=categ, linetype=categ), size=1.2) + 
  geom_point(aes(y=value, color=categ, shape=categ), size=2, alpha=0.9) +
  geom_ribbon(aes(ymin=errmin, ymax=errmax, group=quant, fill=categ), alpha=0.2) +
  facet_grid(observ ~ perturb, scales="free_x") +
  scale_x_log10(breaks=seq(0.01,0.05,0.02)) + xlab("$\\epsilon$") +
  scale_colour_manual(name="Quantity", values=myColors, labels=thetalab) + 
  scale_fill_manual(name="Quantity", values=myColors, labels=thetalab)  + 
  scale_linetype(name="Quantity", labels=thetalab) + 
  scale_shape(name="Quantity", labels=thetalab, solid=FALSE) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), aspect.ratio = 0.5) +
  theme(legend.text.align=0, legend.position = c(0.58,0.6), legend.direction = "horizontal",
        legend.key.size = unit(0.5,"cm"), legend.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
        legend.title=element_text(size=8),
        legend.background = element_rect(color = "black", fill = "white", size=0.2, linetype = "solid"),
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8))
endoffile <- dev.off() 

tikz(file = "FIG4b.tex", standAlone=TRUE, width=2.25, height=3.4)
ggplot(subset(dfstat, categ %in% c("WE", "XI") & perturb==perturblab[2]), aes(epsilon)) + 
  geom_line(aes(y=value, group=quant, color=categ, linetype=categ), size=1.2) + 
  geom_point(aes(y=value, color=categ, shape=categ), size=2, alpha=0.9) +
  geom_ribbon(aes(ymin=errmin, ymax=errmax, group=quant, fill=categ), alpha=0.2) +
  facet_grid(observ ~ perturb, scales="free_x") +
  scale_x_log10(breaks=seq(0.01,0.05,0.02)) + xlab("$\\epsilon$") +
  scale_colour_manual(name="Quantity", values=myColors, labels=thetalab) + 
  scale_fill_manual(name="Quantity", values=myColors, labels=thetalab)  + 
  scale_linetype(name="Quantity", labels=thetalab) + 
  scale_shape(name="Quantity", labels=thetalab, solid=FALSE) +
  theme_bw() +
  theme(axis.title.y=element_blank(), legend.position="none", aspect.ratio = 0.5,
        strip.text.x = element_text(size=8, margin = margin(c(0.1,0,.1,0), unit="cm")),
        strip.text.y = element_text(size=8, margin = margin(c(0,.1,0,.1), unit="cm")),
        axis.text = element_text(size=8))
endoffile <- dev.off() 
