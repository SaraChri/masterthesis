#This script contains script for sorting nutrition data and script for 
#table

#Load packages
library(knitr)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(readr)
library(tinytex)
library(tidyr)
library(flextable)


#nutritiontable with pro.pr kg bodyweight pr. group

#laod dxa data and take out weight and subject to join with nutrition data
fp.weight <-read_excel("data/results/ribose_dxa.xlsx") %>%
        select(weight, subject) %>%
        print()

#nutrition data.
#the data is in grams and 
# added 44 grams protein as sup_protein for ingestion of whey protein isolat,
#And calculated in to the total amount of proteins ingested pr day.
#added 90g of glucose in to carbohydrate 
#added 360 kcal from glucose intake and calculated in total kcal pr day
#one participant did not deliver data on nutrition

#Load nutrition data.
nut.weight <- read_excel("data/results/Ribose_nutrition_result.xlsx")%>%
        select(timepoint, meal, subject, protein, fat, calories,
               carbohydrates, group, sup_pro, sup_gluc, kcal_glu) %>%
        #changed names to compare groups each paired days
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
        #take in weight
        inner_join(fp.weight) %>%
        #summarise each variable to get total pr timepoint.
        #Whey protein supplement is added to total protein ingestion
        summarise(protein = sum(protein + sup_pro),
                  fat = sum(fat),
                  calories = sum(calories + kcal_glu),
                  carbohydrates = sum(carbohydrates + sup_gluc),
                  weight = sum(weight)) %>%
        #summarise to get protein pr weight.
        mutate(proprkg = (protein/weight)) %>%
        pivot_wider(names_from = group,
                    values_from = c(fat, protein, calories, carbohydrates, weight, proprkg)) %>%
        print()

#nutrition table with pro/weight. The table shows mean and sd of total intake 
#in grams pr. group, pr. timepoint.
#How much kcal on glucose supplement?

nuttable <- nut.weight %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= fat_glucose:proprkg_placebo) %>%
        group_by(timepoint, variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        ungroup() %>%
        separate(variable, into = c("variable", "group")) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(timepoint, variable, group, stat,) %>%
        pivot_wider(names_from = variable,
                    values_from = stat) %>%
        select(timepoint, group, calories, carbohydrates, fat, protein,
               proprkg) %>%
        #changed names of group for simplicity 
        mutate(group = if_else(group %in% c("glucose"),
                               "G",
                               if_else(group %in% c("placebo"),
                                       "P",
                                       ""))) %>%
        print()

#Check the table
nuttable %>%
        kable() 

#saved shortcut for table
saveRDS(nuttable, file = "./data/derdata/analysis_nutrition/nuttable.RDS")

##---analysis---##






