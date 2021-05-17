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

#RPE

rpe <- read_excel("data/results/ribose_rpe.xlsx") %>%
        select(subject, timepoint, rpe, group) %>%
        print()

rpe.change <- lmer(rpe ~ timepoint + group + (1|subject), 
                   data = rpe)

summary(rpe.change)

plot(rpe.change)


saveRDS(rpe, file = "./data/derdata/analysis_training/rpe.RDS")

#session score

sc <- read_excel("data/results/ribose_rpe.xlsx") %>%
        select(subject, timepoint, group, `session score`) %>%
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
        #group_by(timepoint, subject, group) %>%
        print()

plot.score <- sc %>%
        select(timepoint, subject, group, `session score`) %>%
        ungroup %>%
        print()


sc.change <- lmer(`session score` ~ timepoint + timepoint:group +(1|subject),
                  data = plot.score)

plot(sc.change)

summary(sc.change)

saveRDS(sc, file = "./data/derdata/analysis_training/scanalysis.RDS")

### join sc score to training table

scjoin <- read_excel("data/results/ribose_rpe.xlsx") %>%
        select(subject, timepoint, group, `session score`) %>%
        print()

saveRDS(scjoin, file = "./data/derdata/analysis_training/scjoin.RDS")


##Table seccionscore only

sctable <- sc %>%
        select(subject, timepoint, group, `session score`) %>%
        group_by(timepoint, group) %>%
        summarise(m = mean(`session score`),
                  s = sd(`session score`)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(timepoint, group, stat) %>%
        mutate(group = if_else(group %in% c("glucose"),
                               "G",
                               if_else(group %in% c("placebo"),
                                       "P",
                                       ""))) %>%
        print()

sctable %>%
        kable() %>%
        kable_styling(bootstrap_options = c("Striped", "hover", "condensed"), full_width = FALSE) %>%
        print()

saveRDS(sctable, file = "./data/derdata/analysis_training/sctable.RDS")
