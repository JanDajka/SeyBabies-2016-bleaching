#################### Prepare dataset

babytemp3 <- read.csv("SeyBabies_Data_Temporal.csv")

babytemp3$year <- as.factor(babytemp3$year)
babytemp3$site <- as.factor(babytemp3$site)
str(babytemp3)

#################### Load packages

library(lattice)
library(nlme)
library(lme4)
library(visreg)
library(MASS)
source("HighstatLibV10.R")
library(ggplot2)
library(tidyverse)
library(gridExtra)

#################### Check response variable

hist(babytemp3$totaldens, breaks = 30)

# All corals
sum(babytemp3$totaldens == 0)
100 * sum(babytemp3$totaldens == 0) / nrow(babytemp3)

# Acropora
sum(babytemp3$acrodens == 0)
100 * sum(babytemp3$acrodens == 0) / nrow(babytemp3)

#Favites
sum(babytemp3$favidens == 0)
100 * sum(babytemp3$favidens == 0) / nrow(babytemp3)

#Porites
sum(babytemp3$poridens == 0)
100 * sum(babytemp3$poridens == 0) / nrow(babytemp3)

par(mfrow = c(2,2))
dotchart(babytemp3$totaldens, main = "All corals")
dotchart(babytemp3$acrodens, main = "Acropora")
dotchart(babytemp3$favidens, main = "Favites")
dotchart(babytemp3$poridens, main = "Porites")

#Outliers

par(mfrow = c(1,1))
plot(x = babytemp3$acrodens, 
     y = babytemp3$totaldens)
identify(x = babytemp3$acrodens, 
         y = babytemp3$totaldens)

#Remove outliers
babytemp3 <- babytemp3[-5, ]
dim(babytemp3)
babytemp3 <- babytemp3[-22, ]
dim(babytemp3)
babytemp3 <- babytemp3[-34, ]
dim(babytemp3)

plot(x = babytemp3$acrodens, 
     y = babytemp3$totaldens)

par(mfrow = c(2, 2))
boxplot(babytemp3$totaldens ~ babytemp3$year, main = "All corals")
boxplot(babytemp3$acrodens ~ babytemp3$year, main = "Acropora")
boxplot(babytemp3$favidens ~ babytemp3$year, main = "Favites")
boxplot(babytemp3$poridens ~ babytemp3$year, main = "Porites")

par(mfrow = c(1, 1))

#################### Check covariates

MyVar <- c("year", "site")
pairs(babytemp3[,MyVar], 
      lower.panel = panel.cor)

#################### Fit model for all corals

LM1_total <- lm(totaldens ~ year, data = babytemp3)

#Test for overdispersion
E1 <- resid(LM1_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(LM1_total))
sum(E1^2) / (N - p)
#strongly Overdispersed

#Account for overdispersion with negative bionomial distribution
GLM3_total <- glm.nb(totaldens ~ year, link = identity, data = babytemp3)
warnings(GLM3_total) 
# due to fitting negative binomial distribution, 
# no other distribution we tested could address observer bias of different years 
# AND overdispersion

#Test for overdispersion
E1 <- resid(GLM3_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLM3_total))
sum(E1^2) / (N - p)
#no dispersion

#Try random effects

GLS1_total <- gls(totaldens ~ year, data = babytemp3)

summary(GLS1_total)
visreg(GLS1_total)
plot(GLS1_total)

LME1_total <- lme(totaldens ~ year, random =  ~1 | site, data = babytemp3)

E1 <- resid(LME1_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(LME1_total))
sum(E1^2) / (N - p)

anova(GLS1_total, LME1_total)

summary(LME1_total)
plot(LME1_total)
visreg(LME1_total)

GLMM2_total <- glmer.nb(totaldens ~ year + (1 | site), 
                        data = babytemp3)
E1 <- resid(GLMM2_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLMM2_total))
sum(E1^2) / (N - p)

AIC(GLMM2_total, LME1_total, GLM3_total)

# Better stick with fixed effects

#Investigate model fit
visreg(GLM3_total)
plot(GLM3_total)

#Model validation all corals
E3 <- resid(GLM3_total, type = "pearson")
F3 <- fitted(GLM3_total)
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

plot(x = babytemp3$totaldens, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

boxplot(E3 ~ year, data = babytemp3)

#################### Model result

summary(GLM3_total)
TukeyHSD(aov(GLM3_total))

#################### Fit model for Acropora

GLM4_acro <- glm.nb(acrodens ~ year, link = log, data = babytemp3)
E1 <- resid(GLM4_acro, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLM4_acro))
sum(E1^2) / (N - p)

plot(GLM4_acro)
summary(GLM4_acro)
TukeyHSD(aov(GLM4_acro))

#Model validation Acropora
E3 <- resid(GLM4_acro, type = "pearson")
F3 <- fitted(GLM4_acro)
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

plot(x = babytemp3$acrodens, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

boxplot(E3 ~ year, data = babytemp3)

#################### Fit model for Favites

GLM4_favi <- glm.nb(favidens ~ year, link = log, data = babytemp3)
E1 <- resid(GLM4_favi, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLM4_favi))
sum(E1^2) / (N - p)

plot(GLM4_favi)
summary(GLM4_favi)
TukeyHSD(aov(GLM4_favi))

#Model validation Favites
E3 <- resid(GLM4_favi, type = "pearson")
F3 <- fitted(GLM4_favi)
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

plot(x = babytemp3$favidens, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

boxplot(E3 ~ year, data = babytemp3)

#################### Fit model for Porites

GLM4_pori <- glm.nb(poridens ~ year, link = log, data = babytemp3)
E1 <- resid(GLM4_pori, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLM4_pori))
sum(E1^2) / (N - p)

plot(GLM3_pori)
summary(GLM4_pori)
TukeyHSD(aov(GLM4_pori))

#Model validation Porites
E3 <- resid(GLM4_pori, type = "pearson")
F3 <- fitted(GLM4_pori)
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

plot(x = babytemp3$poridens, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

boxplot(E3 ~ year, data = babytemp3)

#################### Model interpretation

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

#################### All 4 model fits

visreg(GLM3_total, gg = TRUE) + 
  labs(title = "All corals", x = '', y = 'Predicted model values')
visreg(GLM4_acro, gg = TRUE) + 
  labs(title = "Acropora", x = '', y = 'Predicted model values')
visreg(GLM4_favi, gg = TRUE) + 
  labs(title = "Favites", x = '', y = 'Predicted model values')
visreg(GLM4_pori, gg = TRUE) + 
  labs(title = "Porites", x = '', y = 'Predicted model values')

#################### Model visualisation

# turning into long format, keeping total density only
temp <- gather(babytemp3, taxa, cover, -year, -site, -area) %>%
  filter(taxa %in% c('totaldens'))

# now estimate mean and se across sites
se<-function(x) sd(x)/sqrt(length(x))
temp<-temp %>% group_by(year, taxa) %>%
  summarise(mean = mean(cover), se = se(cover))

temp$year<-as.factor(temp$year)
temp

# Basic line plot with points
p1 <- ggplot(data=temp, aes(x=year, y=mean, group=1)) +
  geom_line(linetype = "dashed", color = "#f44323")+
  geom_point(size = 3, color="#112047")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.02, color ="#2ca1db") +
  annotate(geom="text", x=0.5, y=15, label="A", color = "#636363") +
  labs(x = '', y = expression(paste("Juvenile coral density * ", m^-2)))

# turning into long format, keeping total density only
temp <- gather(babytemp3, taxa, cover, -year, -site, -area) %>%
  filter(taxa %in% c('acrodens', 'favidens', 'poridens'))
temp$year<-as.factor(temp$year)

str(temp)

p2 <- ggplot(data=temp, aes(x=taxa, y=cover, fill=year)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Set2") + 
  annotate(geom="text", x=0.5, y=5.7, label="B", color = "#636363") +
  labs(x = '', y = expression(paste("Juvenile coral density * ", m^-2))) + 
  scale_x_discrete(labels = c(expression(paste(italic("Acropora"))), 
                              expression(paste(italic("Favites"))), 
                               expression(paste(italic("Porites")))))

p <- list(p1, p2)
marrangeGrob(p, nrow=2, ncol=1, top='')
