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

isok60 <-  read_excel("data/results/Ribose_isokinetic.xlsx") %>%
        select(subject, timepoint, test, peak, leg, testpoint, group) %>%
        mutate(timepoint = if_else(timepoint %in% c("-7", "-5", "-1"),
                                   "baseline",
                                   if_else(timepoint %in% c("4", "5"),
                                           "pre",
                                           if_else(timepoint %in% c("8", "9"),
                                                   "midt",
                                                   if_else(timepoint %in% c("t3", "t4"),
                                                           "post",
                                                           ""))))) %>%
        filter(test %in% c("isok60")) %>%
        group_by(timepoint, group) %>%
        select(timepoint, test, peak, group,) %>%
        print()
        summarise(m = mean(isok60)) %>%
        print()


isok60change <- lmer(change ~ group + supplement + (1|subject), 
                     data = isok60)