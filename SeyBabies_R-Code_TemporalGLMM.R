#################### Prepare dataset

setwd("C:/Working directory")

babytemp3 <- read.csv("SeyBabies_Data_Temporal3.csv")

#################### Load packages & data housekeeping

library(lattice)
library(nlme)
library(lme4)
library(MASS)
library(ggplot2)
library(tidyverse)
library(emmeans)
library(gridExtra)

babytemp3<-babytemp3 %>% group_by(year, area, site) %>%
  summarise(total = mean(total), acropora = mean(acropora), 
            favites = mean(favites), porites = mean(porites))
babytemp3$total<-round(babytemp3$total/babytemp3$area)
babytemp3$acropora<-round(babytemp3$acropora/babytemp3$area)
babytemp3$favites<-round(babytemp3$favites/babytemp3$area)
babytemp3$porites<-round(babytemp3$porites/babytemp3$area)

babytemp3$year <- as.factor(babytemp3$year)
babytemp3$site <- as.factor(babytemp3$site)
str(babytemp3)


#################### Check response variable

hist(babytemp3$total, breaks = 30)

# All corals
sum(babytemp3$total == 0)
100 * sum(babytemp3$total == 0) / nrow(babytemp3)

# Acropora
sum(babytemp3$acropora == 0)
100 * sum(babytemp3$acropora == 0) / nrow(babytemp3)

#Favites
sum(babytemp3$favites == 0)
100 * sum(babytemp3$favites == 0) / nrow(babytemp3)

#Porites
sum(babytemp3$porites == 0)
100 * sum(babytemp3$porites == 0) / nrow(babytemp3)

#Outliers

par(mfrow = c(2,2))
dotchart(babytemp3$total, main = "All corals")
dotchart(babytemp3$acropora, main = "Acropora")
dotchart(babytemp3$favites, main = "Favites")
dotchart(babytemp3$porites, main = "Porites")

#site as random effect

par(mfrow = c(1,1))
plot(babytemp3$total ~ babytemp3$site)

#year as fixed effect

par(mfrow = c(2, 2))
boxplot(babytemp3$total ~ babytemp3$year, main = "All corals")
boxplot(babytemp3$acropora ~ babytemp3$year, main = "Acropora")
boxplot(babytemp3$favites ~ babytemp3$year, main = "Favites")
boxplot(babytemp3$porites ~ babytemp3$year, main = "Porites")

par(mfrow = c(1, 1))

#################### Fit model for all corals

#Linear model with offset to account for different quadrat sizes between years

LM1_total <- lm(total ~ year + offset(area), data = babytemp3)

#Test for overdispersion
E1 <- resid(LM1_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(LM1_total))
sum(E1^2) / (N - p)
#Overdispersed

#Random effects

LME1_total <- lme(total ~ year + offset(area), random =  ~1 | site, 
                  data = babytemp3)

E1 <- resid(LME1_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(LME1_total))
sum(E1^2) / (N - p)
#Underdispersed



#Try GLMM
GLMM1_total <- glmer(total ~ year + (1 | site) + offset(area), 
                     family = poisson, data = babytemp3)
E1 <- resid(GLMM1_total, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLMM1_total))
sum(E1^2) / (N - p)
#indicates GLMM

AIC(GLMM1_total, LME1_total, LM1_total)

#Investigate GLMM and LME fit
par(mfrow = c(3, 2))
plot(resid(GLMM1_total), babytemp3$total)
plot(resid(LME1_total), babytemp3$total)
plot(fitted(GLMM1_total), resid(GLMM1_total))
plot(fitted(LME1_total), resid(LME1_total))
hist(resid(LME1_total))
hist(resid(GLMM1_total))

#use GLMM

#################### All corals model result

summary(GLMM1_total)
emmeans(GLMM1_total, list(pairwise ~ year), adjust = "tukey")


#################### Fit model for Acropora

GLMM1_acro <- glmer(acropora ~ year + (1 | site) + offset(area), 
                   family = poisson, data = babytemp3)
E1 <- resid(GLMM1_acro, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLMM1_acro))
sum(E1^2) / (N - p)

hist(resid(GLMM1_acro))
plot(fitted(GLMM1_acro), resid(GLMM1_acro))

#################### Acropora model result

summary(GLMM1_acro)
emmeans(GLMM1_acro, list(pairwise ~ year), adjust = "tukey")

#################### Fit model for Favites
LME1_favi <- lme(favites ~ year + offset(area), random =  ~1 | site, 
                 data = babytemp3)

E1 <- resid(LME1_favi, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(LME1_favi))
sum(E1^2) / (N - p)

hist(resid(LME1_favi))
plot(fitted(LME1_favi), resid(LME1_favi))

#################### Favites model result

summary(LME1_favi)
emmeans(LME1_favi, list(pairwise ~ year), adjust = "tukey")

#################### Fit model for Porites

GLMM1_pori <- glmer(porites ~ year + (1 | site) + offset(area), 
                    family = poisson, data = babytemp3)
E1 <- resid(GLMM1_pori, type = "pearson")
N  <- nrow(babytemp3)
p  <- length(coef(GLMM1_pori))
sum(E1^2) / (N - p)

hist(resid(GLMM1_pori))
plot(fitted(GLMM1_pori), resid(GLMM1_pori))

summary(GLMM1_pori)
emmeans(GLMM1_pori, list(pairwise ~ year), adjust = "tukey")

#################### Observed values visualisation

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

# turning into long format, keeping total density only
temp <- gather(babytemp3, taxa, cover, -year, -site, -area) %>%
  filter(taxa %in% c('total'))

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
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.05, color ="#2ca1db") +
  annotate(geom="text", x=0.5, y=14, label="A", color = "#636363") +
  labs(x = '', y = expression(paste("Mean juvenile coral density * ", m^-2)))

# turning into long format, keeping densities for selected taxa
temp <- gather(babytemp3, taxa, cover, -year, -site, -area) %>%
  filter(taxa %in% c('acropora', 'favites', 'porites'))
temp <- temp %>% group_by(taxa, year, site) %>%
  summarise(cover = mean(cover)) 
temp$year<-as.factor(temp$year)
temp$cover <- as.numeric(temp$cover)
str(temp)
head(temp)

p2 <- ggplot(data=temp, aes(x=taxa, y=cover, fill=year)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Set2") + 
  annotate(geom="text", x=0.5, y=12, label="B", color = "#636363") +
  labs(x = '', y = expression(paste("Juvenile coral density * ", m^-2))) + 
  scale_x_discrete(labels = c(expression(paste(italic("Acropora"))), 
                              expression(paste(italic("Favites"))), 
                               expression(paste(italic("Porites")))))

p <- list(p1, p2)
marrangeGrob(p, nrow=2, ncol=1, top='')
