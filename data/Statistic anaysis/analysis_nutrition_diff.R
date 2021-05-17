#this scrips measures the mean diff on nutrition between condition each paired days
library(knitr)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(readr)
library(knitr)
library(lme4)
library(dplyr)
library(lmerTest)
library(emmeans)


#Does not contain FP 116 because of not delivered data.

#summarise macro
#Added sup_pro as whey protein isolate, 50 grams each group, each timepoint,
#90 grams glucose and 360 kcal.

#This script examine the mean difference on nutrition intake between groups each
#pried days.

fp.weight <-read_excel("data/results/ribose_dxa.xlsx") %>%
        select(weight, subject) %>%
        print()

nutha <- read_excel("data/results/Ribose_nutrition_result.xlsx")%>%
        select(timepoint, subject, group, calories, fat, carbohydrates, protein, sup_pro, sup_gluc, kcal_glu) %>%
        mutate(timepoint = if_else(timepoint %in% c("T1","T2"), 
                                   "1",
                                   if_else(timepoint %in% c("3", "4"),
                                           "2",
                                           if_else(timepoint %in% c("5", "6"),
                                                   "3",
                                                   if_else(timepoint %in% c("7", "8"),
                                                           "4",
                                                           if_else(timepoint %in% c("9", "10"),
                                                                   "5",
                                                                   if_else(timepoint %in% c("T3", "T4"),
                                                                           "6",
                                                                           "na"))))))) %>%
        group_by(timepoint, subject, group) %>%
        inner_join(fp.weight) %>%
        summarise(protein = sum(protein + sup_pro),
                  fat = sum(fat),
                  calories = sum(calories + kcal_glu),
                  carbohydrates = sum(carbohydrates + sup_gluc),
                  weight = sum(weight)) %>%
        mutate(proprkg = (protein/weight)) %>%
  print()
        
        
        
#change data
#Protein

pro_change <- nutha %>%
        select(timepoint, subject, group, protein,) %>%
        print()


prolm <- lmer(protein ~ timepoint + timepoint:group + (1|subject), data = pro_change)

plot(prolm)

summary(prolm)

#No difference on protein between conditions

#Fat
fat_change <- nutha %>%
        select(timepoint, subject, group, fat,) %>%
        print()

fatlm <- lmer(fat  ~ timepoint + timepoint:group +(1|subject), data = fat_change)

plot(fatlm)

summary(fatlm)

#No difference on fat between conditions

#Carbohydrates
carb_change <- nutha %>%
        select(timepoint, subject, group, carbohydrates,) %>%
        print()

carblm <- lmer(carbohydrates  ~ timepoint + timepoint:group + (1|subject), data = carb_change)

plot(carblm)

summary(carblm)
#No difference on fat between conditions


#calories
cal_change <- nutha %>%
        select(timepoint, subject, group, calories,) %>%
        print()



callm <- lmer(calories  ~ 0 + timepoint + timepoint:group + (1|subject), data = cal_change)

plot(callm)

summary(callm)

#There is different intake on total kcal each dat, but no difference between conditions

#Protein pr bodyweigh


prowe <- nutha %>%
        select(timepoint, subject, group, proprkg) %>%
        print()

prokglm <- lmer(proprkg  ~ 0 + timepoint + timepoint:group + (1|subject), data = prowe)

plot(prokglm)

summary(prokglm)

#There is a difference on total protein intake pr/kg each paired day, but
#with no diffence between conditions
