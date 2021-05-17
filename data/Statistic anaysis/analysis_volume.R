library(knitr)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(readr)
library(tinytex)
library(tidyr)
library(knitr)
library(broom)
library(arsenal) 
library(lme4)
library(dplyr)
library(lmerTest)
library(emmeans)

# to combine rpe values and session score
rpe <- readRDS(file = "./data/derdata/analysis_training/rpe.RDS") %>%
        print()

scejoin <- readRDS("./data/derdata/analysis_training/scjoin.RDS") %>%
        print()


# sort training data
tr.volume <- read_excel("data/results/ribose_volume.xlsx") %>%
        select(subject, timepoint, lp.volume, ke.volume, supplement) %>%
        inner_join(rpe) %>%
        inner_join(scejoin) %>%
        mutate(timepoint = if_else(timepoint %in% c("T1","T2"), 
                                   "Day 1",
                                   if_else(timepoint %in% c("D3", "D4"),
                                           "Day 2",
                                           if_else(timepoint %in% c("D5", "D6"),
                                                   "Day 3",
                                                   if_else(timepoint %in% c("D7", "D8"),
                                                           "Day 4",
                                                           if_else(timepoint %in% c("D9", "D10"),
                                                                   "Day 5",
                                                                   if_else(timepoint %in% c("T3", "T4"),
                                                                           "Day 6",
                                                                           "na"))))))) %>%
        group_by(timepoint, group, subject) %>%
        summarise(lp = sum(lp.volume),
                ke = sum(ke.volume),
                rpe =sum(rpe, na.rm =TRUE),
                sc =sum(`session score`, na.rm = TRUE)) %>%
        pivot_wider(names_from = group,
                    values_from = c(lp, ke, rpe, sc)) %>%
        print()

#training data table with rpe and sc        
trvol.table <- tr.volume %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= lp_glucose:sc_placebo) %>%
        group_by(timepoint, variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        ungroup() %>%
        separate(variable, into = c("variable", "group")) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(timepoint, variable, group, stat,) %>%
        pivot_wider(names_from = variable,
                    values_from = stat) %>%
        print()

saveRDS(trvol.table, file = "./data/derdata/analysis_training/traing_table.RDS")

vol.change



tot.volchange <- lmerTest::lmer(rpe ~ group + (1|subject), 
                          data = trvol.table)

        

Tot.volume <- ribose_volume <- read_excel("data/results/ribose_volume.xlsx") %>%
        #select variables
        select(subject, timepoint, leg, tot.volume, group) %>%
        group_by(group) %>%
        #change in %
        mutate(tot.change = (tot.volume / mean(tot.volume, na.rm = TRUE)) * 100) %>%
        select(subject, group, tot.change, timepoint) %>%
        #how to get change from a baseline to a post (?)
        mutate(timepoint = if_else(timepoint %in% c("T1","T2"), 
                                   "Baseline",
                                   if_else(timepoint %in% c("T3", "T4"),
                                   "post",
                                   ""))) %>%
        
        print()


totvolume.change <- lmer(Tot.volume ~ group + (1|subject), 
                      data = Tot.volume)


plot(volume.change)

summary(volume.change)