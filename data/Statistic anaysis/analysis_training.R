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



training <- read_excel("data/results/ribose_training.xlsx") %>%
        select(subject, timepoint, leg, exercise, set, repetitions, load, group) %>%
        mutate(timepoint = if_else(timepoint %in% c("T1","T2"), 
                                   "Day 1",
                                   if_else(timepoint %in% c("3", "4"),
                                           "Day 2",
                                           if_else(timepoint %in% c("5", "6"),
                                                   "Day 3",
                                                   if_else(timepoint %in% c("7", "8"),
                                                           "Day 4",
                                                           if_else(timepoint %in% c("9", "10"),
                                                                   "Day 5",
                                                                   if_else(timepoint %in% c("T3", "T4"),
                                                                           "Day 6",
                                                                           "na"))))))) %>%
        group_by(exercise) %>%
        pivot_wider(names_from = exercise,
                    values_from = load) %>%
        print()
        group_by(group) %>%
        pivot_wider(names_from = group,
                    values_from = load) %>%
        print()
