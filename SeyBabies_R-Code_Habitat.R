#################### Prepare dataset

babyhabitat <- read.csv("SeyBabies_Data_Habitat.csv")

names(babyhabitat)
babyhabitat$island <- as.factor(babyhabitat$island)
babyhabitat$site <- as.factor(babyhabitat$site)
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

MyVar <- c("sandrubble.sc", "cca.sc", "macroalgae.sc", "herbivores.sc", 
           "urchins.sc", "complexity.sc", "site", "habitat")

Mydotplot(babyhabitat[,MyVar])

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

#################### Data exploration


# Outliers

par(mfrow = c(1, 2))
boxplot(babyhabitat$juvtotal, 
        main = "Total Count")
dotchart(babyhabitat$juvtotal, 
         xlab = "Range of data", 
         ylab = "Order of the data")

# Identify outliers

par(mfrow = c(1, 1))
plot(x = babyhabitat$complexity.sc, 
     y = babyhabitat$juvtotal)
identify(x = babyhabitat$complexity.sc, 
         y = babyhabitat$juvtotal)

# Remove outliers

babyhabitat <- babyhabitat[-58, ]
dim(babyhabitat)
babyhabitat <- babyhabitat[-48, ]
dim(babyhabitat)
babyhabitat <- babyhabitat[-45, ]
dim(babyhabitat)

plot(x = babyhabitat$complexity.sc, 
     y = babyhabitat$juvtotal)

# Covariates

MyVar <- c("sandrubble.sc", "cca.sc", "macroalgae.sc", "herbivores.sc", 
           "urchins.sc", "complexity.sc", "site","island", "habitat")
Mydotplot(babyhabitat[,MyVar])
#drop island - captured by site

MyVar <- c("sandrubble.sc", "cca.sc", "macroalgae.sc", "herbivores.sc", 
           "urchins.sc", "complexity.sc", "site", "habitat")
Mydotplot(babyhabitat[,MyVar])


# Collinearity X

pairs(babyhabitat[,MyVar], 
      lower.panel = panel.cor)
corvif(babyhabitat[,MyVar])
# drop site - VIF = 7.56
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


# Zero inflation
sum(babyhabitat$juvtotal == 0)
100 * sum(babyhabitat$juvtotal == 0) / nrow(babyhabitat)

#################### Start of analysis

# Create formulae

formula1 <- juvtotal ~ sandrubble.sc + cca.sc + macroalgae.sc + urchins.sc +
  habitatgra.sc:complexity.sc + macroalgae.sc:complexity.sc + 
  herbivores.sc:macroalgae.sc

# Fit ZIP
modelZIP1 <- zeroinfl(formula = formula1, dist = "poisson", link = "logit",
                     data = babyhabitat)

E1 <- resid(modelZIP1, type = "pearson")
N  <- nrow(babyhabitat)
p  <- length(coef(modelZIP1))
sum(E1^2) / (N - p)

# Fit ZINB
modelZINB1 <- zeroinfl(formula = formula1, dist = "negbin", link = "logit",
                      data = babyhabitat)

E1 <- resid(modelZINB1, type = "pearson")
N  <- nrow(babyhabitat)
p  <- length(coef(modelZINB1))
sum(E1^2) / (N - p)

#################### Model selection

step(modelZINB1)

formula2 <- juvtotal ~ sandrubble.sc + cca.sc + macroalgae.sc + 
  habitatgra.sc:complexity.sc + macroalgae.sc:complexity.sc + 
  herbivores.sc:macroalgae.sc

modelZINB2 <- zeroinfl(formula = formula2, dist = "negbin", link = "logit",
                       data = babyhabitat)

step(modelZINB2)

AIC(modelZINB1, modelZINB2)

E1 <- resid(modelZINB2, type = "pearson")
N  <- nrow(babyhabitat)
p  <- length(coef(modelZINB1))
sum(E1^2) / (N - p)

summary(modelZINB2)
coef(modelZINB2)

################### Model validation

plot(resid(modelZINB2))
plot(babyhabitat$sandrubble.sc ~ resid(modelZINB2))
plot(babyhabitat$cca.sc ~ resid(modelZINB2))
plot(babyhabitat$macroalgae.sc ~ resid(modelZINB2))
plot(babyhabitat$complexity.sc ~ resid(modelZINB2))
plot(babyhabitat$herbivores.sc ~ resid(modelZINB2))


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

ggplot(coefs, aes(vars, Estimate)) + 
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
  scale_x_discrete(limits= c('habitatgra.sc:complexity.sc', 'macroalgae.sc:herbivores.sc', 'macroalgae.sc:complexity.sc', 'cca.sc', 'sandrubble.sc', 'macroalgae.sc'), 
                   labels = c('Granitic reef * complexity', 'Macroalgae * herbivores', 'Macroalgae * complexity', 'CCA', 'Sand & rubble', 'Macroalgae'))

# Extract coefficients and standard errors from COUNT model summary
coefs = as.data.frame(summary(modelZINB2)$coefficients$count[,1:2])
names(coefs)[2] = "se" 
coefs$vars = rownames(coefs)
coefs$Estimate


coefs<-coefs[!coefs$vars %in% c('(Intercept)', 'Log(theta)'),]

# Coefficient plot
ggplot(coefs, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="#C4CFD0") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se), 
                lwd=1, colour="#F24D29", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#1DACE8", width=0) +
  geom_point(size=4, pch=21, fill="#1C366B") +
  theme(axis.text.x = element_text(angle=0)) + coord_flip() +
  labs(y = 'Coefficient estimate', x = '') +
  scale_x_discrete(limits= c('habitatgra.sc:complexity.sc', 'macroalgae.sc:herbivores.sc', 'macroalgae.sc:complexity.sc', 'cca.sc', 'sandrubble.sc', 'macroalgae.sc'), 
                   labels = c('Granitic reef * complexity', 'Macroalgae * herbivores', 'Macroalgae * complexity', 'CCA', 'Sand & rubble', 'Macroalgae'))


#################### Model trends

#sand & rubble

plottrend.sand <- visreg(modelZINB2, 'sandrubble.sc', gg = TRUE, 
                         xlab='Cover in %', 
                         ylab='')
sr.axis<-data.frame(sr = round(seq(min(babyhabitat$SANDRUBBLE), 
                                   max(babyhabitat$SANDRUBBLE), length.out=4),0),
                    sr.sc=   seq(min(babyhabitat$sandrubble.sc), 
                                 max(babyhabitat$sandrubble.sc), length.out=4))
g1 <- plottrend.sand + ggtitle('Sand & rubble') +
  scale_x_continuous(breaks = sr.axis$sr.sc, labels=sr.axis$sr)

#macroalgae

plottrend.ma <- visreg(modelZINB2, 'macroalgae.sc', gg = TRUE, 
                         xlab='', 
                         ylab='Juvenile coral density per quadrat')
ma.axis<-data.frame(ma = round(seq(min(babyhabitat$MACROALGAE), 
                                   max(babyhabitat$MACROALGAE), length.out=4),0),
                    ma.sc=   seq(min(babyhabitat$macroalgae.sc), 
                                 max(babyhabitat$macroalgae.sc), length.out=4))
g2 <- plottrend.ma + ggtitle('Macroalgae') +
  scale_x_continuous(breaks = ma.axis$ma.sc, labels=ma.axis$ma)

#CCA

plottrend.cca <- visreg(modelZINB2, 'cca.sc', gg = TRUE,
                         xlab='', 
                         ylab='')
cca.axis<-data.frame(cca = round(seq(min(babyhabitat$CCA), 
                                   max(babyhabitat$CCA), length.out=4),0),
                    cca.sc=   seq(min(babyhabitat$cca.sc), 
                                 max(babyhabitat$cca.sc), length.out=4))
g3 <- plottrend.cca + ggtitle('CCA') +
  scale_x_continuous(breaks = cca.axis$cca.sc, labels=cca.axis$cca)

g <- list(g2, g1, g3)
marrangeGrob(g, nrow=1, ncol=3, top='')


#macroalgae : complexity

plottrend.macomp <- visreg(modelZINB2, 'macroalgae.sc', by = 'complexity.sc', gg = TRUE,
                          xlab='', 
                          ylab='Juvenile coral density per quadrat')
ma.axis<-data.frame(ma = round(seq(min(babyhabitat$MACROALGAE), 
                                   max(babyhabitat$MACROALGAE), length.out=4),0),
                    ma.sc=   seq(min(babyhabitat$macroalgae.sc), 
                                 max(babyhabitat$macroalgae.sc), length.out=4))
g4 <- plottrend.macomp + scale_x_continuous(breaks = ma.axis$ma.sc, labels=ma.axis$ma)

#macroalgae : herbivores

plottrend.maherb <- visreg(modelZINB2, 'macroalgae.sc', by = 'herbivores.sc', gg = TRUE,
                           xlab='Macroalgae cover in %', 
                           ylab='Juvenile coral density per quadrat')
ma.axis<-data.frame(ma = round(seq(min(babyhabitat$MACROALGAE), 
                                   max(babyhabitat$MACROALGAE), length.out=4),0),
                    ma.sc=   seq(min(babyhabitat$macroalgae.sc), 
                                 max(babyhabitat$macroalgae.sc), length.out=4))
g5 <- plottrend.maherb + scale_x_continuous(breaks = ma.axis$ma.sc, labels=ma.axis$ma)
    
g2 <- list(g4, g5)
marrangeGrob(g2, nrow=2, ncol=1, top='')

#complexity : habitat                           
                           
plottrend.comphabitat <- visreg(modelZINB2, 'complexity.sc', by = 'habitatgra.sc', gg = TRUE,
       xlab='Structural complexity', 
       ylab='Juvenile coral density per quadrat')
comp.axis<-data.frame(comp = round(seq(min(babyhabitat$complexity), 
                                   max(babyhabitat$complexity), length.out=4),0),
                    comp.sc=   seq(min(babyhabitat$complexity.sc), 
                                 max(babyhabitat$complexity.sc), length.out=4))
habitat.labs <- c('Carbonate & patch', 'Granite')
plottrend.comphabitat + scale_x_continuous(breaks = comp.axis$comp.sc, labels=comp.axis$comp)