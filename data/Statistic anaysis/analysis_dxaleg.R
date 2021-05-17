#DXA pr. leg.
#subject 104, 110 and 111 is not in table<- did not complete intervention

library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(broom)
library(knitr)
library(emmeans)

#Load dxa data

#Makge grams into kg
dxaleg <-  read_excel("data/results/ribose_dxaprleg.xlsx") %>%
        select(subject, sex, leg, leanmassg, fatmassg, totalmasskg) %>%
        #make grams into kg
        mutate(leanmass = leanmassg/1000,
               fatmass = fatmassg/1000)%>%
        select(subject, sex, leg, leanmass, fatmass, totalmasskg) %>%
        print()
#Make analysis bewteen groups. 
dxa.stat <- dxaleg %>%
        select(subject, leg, leanmass, fatmass, totalmasskg) %>%
        print()

#lean mass between legs, within subject.
#no difference bewtween legs

leanmass.lm <- lmer(leanmass ~ leg + leg:subject + (1|subject), data = dxa.stat)

plot(leanmass.lm)

summary(leanmass.lm)

emmeans(leanmass.lm, specs = ~ "leg") %>%
        confint()

# fat mass between legs, within subject
#no difference bewtween legs

fat.lm <- lmer(fatmass ~ leg + leg:subject +(1|subject), data = dxa.stat)

plot(fat.lm)

summary(fat.lm)

emmeans(fat.lm, specs = ~ "leg") %>%
        confint()

#Totalmass between legs, within subject.
#no difference bewtween legs

total.lm <- lmer(totalmasskg ~leg + leg:subject + (1|subject), data = dxa.stat)

plot(total.lm)

summary(total.lm)

supsup <- emmeans(total.lm, specs = ~ "leg") %>%
        confint() %>%
        print()


#Table showing leanmass, fatmass and total mass in kg pr. leg
dxalegtable <- dxaleg %>%
        select(leg, leanmass, fatmass, totalmasskg) %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= leanmass:totalmasskg) %>%
        group_by(leg, variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(leg, variable, stat) %>%
        pivot_wider(names_from = leg,
                    values_from = stat) %>%
        print()

#table in markdown
dxalegtable %>%
        kable()


saveRDS(dxalegtable, file = "./data/derdata/analysis_dxa/dxalegtable.RDS")
