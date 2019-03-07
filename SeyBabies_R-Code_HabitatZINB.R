#################### Prepare dataset

setwd("C:/Working Directory")
babyhabitat <- read.csv("SeyBabies_Data_Habitat1.csv")

#################### Load packages

library(dplyr)
library(lattice)
source("HighstatLibV10.R")
library(Matrix)
library(lme4)
library(nlme)
library(mgcv)
library(pscl)
library(MASS)
library(visreg)
library(ggplot2)
library(gridExtra)
library(coefplot)
library(glmmTMB)
source("predict.zeroinfl.R")

babyhabitat$island <- as.factor(babyhabitat$island)
babyhabitat$site <- as.character(babyhabitat$site)
babyhabitat$habitat <- as.factor(babyhabitat$habitat)
str(babyhabitat)

#Mean-center covariates

babyhabitat$sandrubble.sc <- scale(babyhabitat$SANDRUBBLE)
babyhabitat$cca.sc <- scale(babyhabitat$CCA)
babyhabitat$macroalgae.sc <- scale(babyhabitat$MACROALGAE)
babyhabitat$herbivores.sc <- scale(babyhabitat$herbivores)
babyhabitat$urchins.sc <- scale(babyhabitat$urchins)
babyhabitat$complexity.sc <- scale(babyhabitat$complexity)
babyhabitat$habitatpat.sc <- scale(babyhabitat$habitatpat)
babyhabitat$habitatgra.sc <- scale(babyhabitat$habitatgra)
babyhabitat$habitatcar.sc <- scale(babyhabitat$habitatcar)


head(babyhabitat)

#################### Data exploration


# Outliers

par(mfrow = c(1, 2))
boxplot(babyhabitat$juvtotal, 
        main = "Total Count")
dotchart(babyhabitat$juvtotal, 
         xlab = "Range of data", 
         ylab = "Order of the data")


# Collinearity X
MyVar <- c("sandrubble.sc", "cca.sc", "macroalgae.sc", "herbivores.sc", 
           "urchins.sc", "complexity.sc", "habitatgra.sc")
pairs(babyhabitat[,MyVar], 
      lower.panel = panel.cor)
corvif(babyhabitat[,MyVar])
#with dummy variables for factor 'habitat'
# all VIF < 3


# Relationship X & Y

Myxyplot(babyhabitat,MyVar, "juvtotal")


# Interactions

Myxyplot(babyhabitat,MyVar, "complexity.sc")
Myxyplot(babyhabitat,MyVar, "macroalgae.sc")
Myxyplot(babyhabitat,MyVar, "herbivores.sc")
Myxyplot(babyhabitat,MyVar, "urchins.sc")
Myxyplot(babyhabitat,MyVar, "sandrubble.sc")
Myxyplot(babyhabitat,MyVar, "cca.sc")
Myxyplot(babyhabitat,MyVar, "habitatgra.sc")

# Dummy variables for reef types
babyhabitat$habitatgra <- as.factor(babyhabitat$habitatgra)
babyhabitat$habitatcar <- as.factor(babyhabitat$habitatcar)
babyhabitat$habitatpat <- as.factor(babyhabitat$habitatpat)

# Zero inflation
sum(babyhabitat$juvtotal == 0)
100 * sum(babyhabitat$juvtotal == 0) / nrow(babyhabitat)

#################### Comparison of available zero-indlated model packages

library(sjPlot)
library(here)
library(tidyverse)

df<-babyhabitat

library(glmmTMB)

## Fit TMB with nbinom
m1<-glmmTMB(juvtotal ~ sandrubble.sc + macroalgae.sc + herbivores.sc + cca.sc +
              complexity.sc + habitatgra.sc + complexity.sc:habitatgra.sc +  
              macroalgae.sc:complexity.sc + macroalgae.sc:herbivores.sc + (1 | site),
            ziformula= ~.,
            data=df,
            family=nbinom1)
summary(m1)

## Fit TMB with nbinom2
m2<-glmmTMB(juvtotal ~ sandrubble.sc + macroalgae.sc + herbivores.sc + cca.sc +
              complexity.sc + habitatgra.sc + complexity.sc:habitatgra.sc +  
              macroalgae.sc:complexity.sc + macroalgae.sc:herbivores.sc + (1 | site),
            ziformula= ~.,
            data=df,
            family=nbinom2)
summary(m2)


AIC(m1, m2) ## m2 better

par(mfrow=c(2,2), mar=c(2,2,2,2))
hist(resid(m2))
plot(fitted(m2), resid(m2))
plot(fitted(m2), df$juvtotal)

## Fit original model
library(pscl)
m3<-zeroinfl(juvtotal ~ sandrubble.sc + macroalgae.sc + herbivores.sc + cca.sc +
               complexity.sc + habitatgra.sc + complexity.sc:habitatgra.sc +  
               macroalgae.sc:complexity.sc + macroalgae.sc:herbivores.sc,
             data=df,
             dist = 'negbin')
summary(m3)

AIC(m1, m2, m3) ## m2 better

par(mfrow=c(2,2), mar=c(2,2,2,2))
m.focal<-m3
hist(resid(m.focal))
plot(fitted(m.focal), resid(m.focal))
plot(predict(m.focal), df$juvtotal)
#visreg::visreg(m.focal)

## compare coefficients from m3 and m2

coef2<-data.frame(est=summary(m2)$coef$cond[,'Estimate'], model = 'Count')
coef2<-rbind(coef2, data.frame(est=summary(m2)$coef$zi[,'Estimate'], model = 'Zero-inflated'))
coef2$covariate<-rep(c('Intercept', 'Sand & Rubble', 'Macroalgae', 'Herbivores', 'CCA', 'Complexity', 'Granitic reef type', 'Complexity * reef type', 'Macroalgae * complexity', 'Macroalgae * herbivores'), 
                     times=2)
coef2$se<-c(summary(m2)$coef$cond[,'Std. Error'], summary(m2)$coef$zero[,'Std. Error'])
coef2$method<-'TMB'

coef3<-data.frame(est=summary(m3)$coef$count[,'Estimate'], model = 'Count')
coef3<-coef3[!rownames(coef3) == 'Log(theta)',] ## drop dispersion estimates
coef3<-rbind(coef3, data.frame(est=summary(m3)$coef$zero[,'Estimate'], model = 'Zero-inflated'))
coef3$covariate<-rep(c('Intercept', 'Sand & Rubble', 'Macroalgae', 'Herbivores', 'CCA', 'Complexity', 'Granitic reef type', 'Complexity * reef type', 'Macroalgae * complexity', 'Macroalgae * herbivores'), times=2)
coef3$se<-c(summary(m3)$coef$count[1:10,'Std. Error'], summary(m3)$coef$zero[,'Std. Error'])
coef3$method<-'PSCL'

coef<-rbind(coef2, coef3)

## separate intercept from covariates in plots to see standard errors

g1<-ggplot(coef[!coef$covariate=='Intercept',], aes(covariate, est, col=method)) + 
  geom_pointrange(aes(ymin = est - 2*se, ymax= est + 2*se), position=position_dodge(width=0.4)) +
  coord_flip() + facet_wrap(~model) + geom_hline(yintercept=0, linetype=2) +
  labs( x = '', y = '') +
  theme(legend.position = 'none',
        legend.title=element_blank())

g2<-ggplot(coef[coef$covariate=='Intercept',], aes(covariate, est, col=method)) + 
  geom_pointrange(aes(ymin = est - 2*se, ymax= est + 2*se), position=position_dodge(width=0.4)) +
  coord_flip() + facet_wrap(~model, scales='free') + #geom_hline(yintercept=0, linetype=2) +
  labs( x = '', y = 'Coefficient estimate') +
  theme(legend.position = 'bottom')

cowplot::plot_grid(g1, g2, nrow=2, align = 'v', rel_heights=c(1, 0.5))

# use pscl due to more conservative errors and ability to extract Pearson residuals

#################### Start of analysis

# Create formulae

formula2 <- juvtotal ~ sandrubble.sc + macroalgae.sc + herbivores.sc + cca.sc +
  complexity.sc + habitatgra.sc + complexity.sc:habitatgra.sc +  
  macroalgae.sc:complexity.sc + macroalgae.sc:herbivores.sc

modelZINB2 <- zeroinfl(formula = formula2, dist = "negbin", link = "logit",
                       data = babyhabitat)


step(modelZINB2)

################### Model validation

E1 <- resid(modelZINB2, type = "pearson")
N  <- nrow(babyhabitat)
p  <- length(coef(modelZINB2))
sum(E1^2) / (N - p)

plot(predict(modelZINB2), resid(modelZINB2))
plot(resid(modelZINB2))

summary(modelZINB2)

E2 <- resid(modelZINB2, type = "pearson")
F2 <- fitted(modelZINB2)
op <- par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals"
plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
abline(h = 0, lty = 2)


#################### Model plotting

# ggplot theme

theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

theme_set(theme_sleek())


#################### Effect sizes

# Extract coefficients and standard errors from ZERO model summary
coefs = as.data.frame(summary(modelZINB2)$coefficients$zero[,1:2])
names(coefs)[2] = "se" 
coefs$vars = rownames(coefs)
coefs$Estimate
coefs<-coefs[!coefs$vars %in% c('(Intercept)', 'Log(theta)'),]
coefs

# Coefficient plot
g1 <- ggplot(coefs, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="#ccd5dd") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se), 
                lwd=1, colour="#f44323", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#2ca1db", width=0) +
  geom_point(size=4, pch=21, fill="#112047") +
  theme(axis.text.x = element_text(angle=1)) + 
  coord_flip() + 
  scale_y_reverse() + 
  labs(y = 'Coefficient estimate', x = '') +
  scale_x_discrete(limits= c('complexity.sc:habitatgra.sc', 'macroalgae.sc:complexity.sc', 'macroalgae.sc:herbivores.sc', 'complexity.sc', 'cca.sc', 'macroalgae.sc', 'habitatgra.sc', 'sandrubble.sc', 'herbivores.sc'), 
                   labels = c('Granitic reef * complexity', 'Macroalgae * complexity', 'Macroalgae * herbivores', 'Complexity', 'CCA', 'Macroalgae', 'Granitic reef', 'Sand & rubble', 'Herbivores'))

# Extract coefficients and standard errors from COUNT model summary
coefs = as.data.frame(summary(modelZINB2)$coefficients$count[,1:2])
names(coefs)[2] = "se" 
coefs$vars = rownames(coefs)
coefs$Estimate


coefs<-coefs[!coefs$vars %in% c('(Intercept)', 'Log(theta)'),]

# Coefficient plot
g2 <- ggplot(coefs, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="#C4CFD0") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se), 
                lwd=1, colour="#F24D29", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#1DACE8", width=0) +
  geom_point(size=4, pch=21, fill="#1C366B") +
  theme(axis.text.x = element_text(angle=0)) + coord_flip() +
  labs(y = 'Coefficient estimate', x = '') +
  scale_x_discrete(limits= c('macroalgae.sc:herbivores.sc', 'complexity.sc:habitatgra.sc', 'macroalgae.sc:complexity.sc', 'herbivores.sc', 'habitatgra.sc', 'cca.sc', 'complexity.sc', 'sandrubble.sc', 'macroalgae.sc'), 
                   labels = c('Macroalgae * herbivores', 'Granitic reef * complexity', 'Macroalgae * complexity', 'Herbivores', 'Granitic reef', 'CCA', 'Complexity', 'Sand & rubble', 'Macroalgae'))


cowplot::plot_grid(g1, g2, nrow=1, align = 'v', rel_heights=c(1, 0.5))

#################### Model trends

summary(modelZINB2)

#Partial residual plots

visreg(modelZINB2, 'herbivores.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'cca.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'complexity.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'sandrubble.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'macroalgae.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'habitatgra.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'complexity.sc', by = 'habitatgra.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'macroalgae.sc', by = 'complexity.sc', gg = TRUE, overlay = TRUE)
visreg(modelZINB2, 'macroalgae.sc', by = 'herbivores.sc', gg = TRUE, overlay = TRUE)

# Setup master prediction dataframe. All effects = 0
master<-data.frame(macroalgae.sc = rep(0,30),
                   complexity.sc = rep(0,30),
                   sandrubble.sc = rep(0,30),
                   herbivores.sc = rep(0,30),
                   cca.sc = rep(0,30),
                   macroalgae.sc = rep(0,30),
                   habitatgra.sc = rep(0,30))

# Figure 3
# Macroalgae as focal covariate
pred1<-master
pred1$macroalgae.sc<-seq(min(babyhabitat$macroalgae.sc),max(babyhabitat$macroalgae.sc), length.out=30)
pred1$macroalgae<-seq(min(babyhabitat$MACROALGAE),max(babyhabitat$MACROALGAE), length.out=30)
pred1$estimate<-predict(modelZINB2, newdata=pred1, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred1, MC = 1000, se = TRUE, type='response')[[2]]
pred1$upper<-cis$upper
pred1$lower<-cis$lower

# Sand 
pred2<-master
pred2$sandrubble.sc<-seq(min(babyhabitat$sandrubble.sc),max(babyhabitat$sandrubble.sc), length.out=30)
pred2$sandrubble<-seq(min(babyhabitat$SANDRUBBLE),max(babyhabitat$SANDRUBBLE), length.out=30)
pred2$estimate<-predict(modelZINB2, newdata=pred2, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred2, MC = 1000, se = TRUE, type='response')[[2]]
pred2$upper<-cis$upper
pred2$lower<-cis$lower

# CCA
pred3 <- master
pred3$cca.sc<-seq(min(babyhabitat$cca.sc),max(babyhabitat$cca.sc), length.out=30)
pred3$cca<-seq(min(babyhabitat$CCA),max(babyhabitat$CCA), length.out=30)
pred3$estimate<-predict(modelZINB2, newdata=pred3, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred3, MC = 1000, se = TRUE, type='response')[[2]]
pred3$upper<-cis$upper
pred3$lower<-cis$lower

# Herbivores
pred4 <- master
pred4$herbivores.sc<-seq(min(babyhabitat$herbivores.sc),max(babyhabitat$herbivores.sc), length.out=30)
pred4$herbivores<-seq(min(babyhabitat$herbivores),max(babyhabitat$herbivores), length.out=30)
pred4$estimate<-predict(modelZINB2, newdata=pred4, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred4, MC = 1000, se = TRUE, type='response')[[2]]
pred4$upper<-cis$upper
pred4$lower<-cis$lower

# Complexity
pred5 <- master
pred5$complexity.sc<-seq(min(babyhabitat$complexity.sc),max(babyhabitat$complexity.sc), length.out=30)
pred5$complexity<-seq(min(babyhabitat$complexity),max(babyhabitat$complexity), length.out=30)
pred5$estimate<-predict(modelZINB2, newdata=pred5, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred5, MC = 1000, se = TRUE, type='response')[[2]]
pred5$upper<-cis$upper
pred5$lower<-cis$lower

# Granite
pred6 <- master
pred6$habitatgra.sc<-seq(min(babyhabitat$habitatgra.sc),max(babyhabitat$habitatgra.sc), length.out=30)
pred6$estimate<-predict(modelZINB2, newdata=pred6, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred6, MC = 1000, se = TRUE, type='response')[[2]]
pred6$upper<-cis$upper
pred6$lower<-cis$lower

# Final figures of focal covariates with all other effects held to zero

p1 <- ggplot(pred1, aes(x = macroalgae, y = estimate)) +
  geom_line() + labs(x = 'Macroalgae cover in %', y = expression(paste("Predicted juvenile coral count * ", m^-2)))
p1 <- p1 + geom_ribbon(data = pred1, aes(ymin = lower, ymax = upper),
                 alpha = 0.1)

p2 <- ggplot(pred2, aes(x = sandrubble, y = estimate)) +
  geom_line() + labs(x = 'Sand & rubble cover in %', y = '')
p2 <- p2 + geom_ribbon(data = pred2, aes(ymin = lower, ymax = upper),
                 alpha = 0.1)

p3 <- ggplot(pred3, aes(x = cca, y = estimate)) +
  geom_line() + labs(x = 'CCA cover in %', y = '')
p3 <- p3 + geom_ribbon(data = pred3, aes(ymin = lower, ymax = upper),
                 alpha = 0.1)

p4 <- ggplot(pred4, aes(x = herbivores, y = estimate)) +
  geom_line() + labs(x = expression(paste("Herbivore biomass in kg * ", ha^-1)), y = '')
p4 <- p4 + geom_ribbon(data = pred4, aes(ymin = lower, ymax = upper),
                 alpha = 0.1)

fig3 <- list(p1, p2, p3, p4)
marrangeGrob(fig3A, nrow=1, ncol=4, top='')

# Figure 4
# Macroalgae * complexity interactions
pred7<-master
pred7$macroalgae.sc<-seq(min(babyhabitat$macroalgae.sc),max(babyhabitat$macroalgae.sc), length.out=30)
pred7$macroalgae<-seq(min(babyhabitat$MACROALGAE),max(babyhabitat$MACROALGAE), length.out=30)
pred7<-do.call('rbind', replicate(3, pred7, simplify=FALSE))
pred7$complexity.sc<-rep(c(0, 2, 3.5), each = 30)
pred7$complexity<-rep(c(0, 2, 3.5), each = 30)

pred7$estimate<-predict(modelZINB2, newdata=pred7, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred7, MC = 1000, se = TRUE, type='response')[[2]]
pred7$upper<-cis$upper	
pred7$lower<-cis$lower

p7 <- ggplot(pred7, aes(macroalgae, estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.1) + 
  labs(x = 'Macroalgae cover in %', y = expression(paste("Predicted juvenile coral count * ", m^-2))) + 
  facet_wrap(~complexity)

# Macroalgae * herbivore interaction
pred8<-master
pred8$macroalgae.sc<-seq(min(babyhabitat$macroalgae.sc),max(babyhabitat$macroalgae.sc), length.out=30)
pred8$macroalgae<-seq(min(babyhabitat$MACROALGAE),max(babyhabitat$MACROALGAE), length.out=30)

pred8<-do.call('rbind', replicate(3, pred8, simplify=FALSE))
quants<-quantile(babyhabitat$herbivores.sc)
pred8$herbivores.sc<-rep(c(quants[2], quants[3], quants[4]), each = 30)
pred8$herbivores<-rep(c(100, 300, 500), each = 30)

pred8$estimate<-predict(modelZINB2, newdata=pred8, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred8, MC = 1000, se = TRUE, type='response')[[2]]
pred8$upper<-cis$upper	
pred8$lower<-cis$lower

p8 <- ggplot(pred8, aes(macroalgae, estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.1) + 
  facet_wrap(~herbivores) + 
  labs(x = 'Macroalgae cover in %', y = '')

fig4 <- list(p7, p8)
marrangeGrob(fig4, nrow=1, ncol=2, top='')

# Figure 5
# Complexity * reef type interaction
pred9<-master
pred9$complexity.sc<-seq(min(babyhabitat$complexity.sc),max(babyhabitat$complexity.sc), length.out=30)
pred9$complexity<-seq(min(babyhabitat$complexity),max(babyhabitat$complexity), length.out=30)

pred9<-do.call('rbind', replicate(2, pred9, simplify=FALSE))
pred9$habitatgra.sc<-c(rep(min(babyhabitat$habitatgra.sc),30),
                       rep(max(babyhabitat$habitatgra.sc),30))

pred9$estimate<-predict(modelZINB2, newdata=pred9, type='response')
cis<-predict.zeroinfl(modelZINB2, newdata = pred9, MC = 1000, se = TRUE, type='response')[[2]]
pred9$upper<-cis$upper	
pred9$lower<-cis$lower

p9 <- ggplot(pred9, aes(complexity, estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.1) + 
  facet_wrap(~habitatgra.sc) + 
  labs(x = 'Complexity', y = expression(paste("Predicted juvenile coral count * ", m^-2)))
p9