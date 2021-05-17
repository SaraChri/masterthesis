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


onerm <- read_excel("data/results/ribose_1rm.xlsx", na="NA") %>%
        #selec needed variables
        select(subject, timepoint, exercise, leg, load) %>%
        #fixed names to get mean 
        mutate(timepoint = gsub("-", "", timepoint),
               timepoint = if_else(timepoint %in% c("7", "5"), 
                                   "1RM",
                                   "post")) %>%
        #mean pr exercise in kg
        group_by(exercise, timepoint, subject) %>%
        summarise(load = mean (load, na.rm=TRUE)) %>%
        pivot_wider(names_from = "exercise",
                    values_from = "load") %>%
        print()

