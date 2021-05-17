
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(broom)
library(knitr)
library(emmeans)
library(knitr)

# descriptive table of subject
#all variables is in kg excepts for age and height.

#make grams into kg
table.dxa <- read_excel("data/results/ribose_dxa.xlsx") %>%
        select(sex, subject, age, height, weight, leanmasskg, fatmasskg, ffmg) %>%
        #make grams into kg
        mutate(ffkg = ffmg/1000) %>%
        select(age, height, weight, leanmasskg, fatmasskg, ffkg) %>%
        print()

# dxa table on total group.
dxa.table2 <- table.dxa %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= age:ffkg) %>%
        select(variable, values) %>%
        group_by(variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(variable, stat) %>%
        mutate(variable = if_else(variable %in% c("age"),
                                  "Age",
                                  if_else(variable %in% c("weight"),
                                          "Weight",
                                          if_else(variable %in% c("ffkg"),
                                                  "Fat free mass (kg)",
                                                  if_else(variable %in% c("height"),
                                                          "Height (cm)",
                                                          if_else(variable %in% c("leanmasskg"),
                                                                  "Lean mass (kg)",
                                                                  if_else(variable %in% c("fatmasskg"),
                                                                          "Fatt mass (kg)",
                                                                          ""))))))) %>%
        print()

dxa.table2 %>%
        kable()

saveRDS(dxa.table2, file = "./data/derdata/analysis_dxa/dxatable.RDS")
        
              
#2. take out kg of subject and merge it with nutrition to protein pr kg pr group

 
saveRDS(fp.weight, file = "./data/derdata/analysis_dxa/fp_weight.RDS")

