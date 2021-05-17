#### Blood glucose change analysis

## Author: Kristian Lian
## Project: Ribose

## This script analyses and plots the difference in blood glucose between supplement groups, with x axis
# scaling (using sample_time instead of the mutated factor "time" as in the original analysis)

## Data
# Subject
# Timepoints: T1 (baseline, before supplements and RT), T3 (posttest leg #1), and T4 (posttest leg #2)
# Sample_time: Minutes following protein ingestion (0-270)
# Time: time of day
# Glu: blood glucose, mmol/L
# Lak: blood lactat, mmol/L
# Supplement: glucose or placebo

# Packages
library(tidyverse)
library(readxl)
library(nlme)
library(lme4)
library(knitr)
library(broom)
library(emmeans)
library(dplyr)
library(ggplot2)
library(ggpubr)
## Data handling

gluc.dat <- read_excel("data/results/ribose_gluctest.xlsx") %>%
        print()

gluc.dat$sample_time <- as.character(gluc.dat$sample_timepoint)

glu.dat <- gluc.dat %>%
        mutate(time = if_else(sample_time == "0",
                              "baseline",
                              if_else(sample_time == "45",
                                      "min45",
                                      if_else(sample_time == "90",
                                              "min90",
                                              if_else(sample_time == "120",
                                                      "min120",
                                                      if_else(sample_time == "135",
                                                              "min135",
                                                              if_else(sample_time == "150",
                                                                      "min150",
                                                                      if_else(sample_time == "270",
                                                                              "min270", sample_time))))))),
               time = factor(time, levels = c("baseline", "min45", "min90", "min120", "min135", "min150", "min270"))) %>%
        mutate(glu = as.numeric(glu),
               lak = as.numeric(lak)) %>%
        print()

#total change


glumean <- glu.dat %>%
        select(subject, time, glu, supplement) %>%
        group_by(time, supplement) %>%
        summarise(glumean = mean(glu, na.rm = TRUE),
                  glusd = sd(glu, na.rm = TRUE)) %>%
        print()


gluchangeplot <- ggplot(glumean, aes(group = supplement, y = glumean, x = time)) +
        geom_line(position = pos, aes(linetype = supplement)) +
        labs(x = "Minutes", y = "Blood glucose levels \n(mmol/L)\n", fill = "Supplement" ) +
        scale_x_discrete(labels=c("baseline" = "0", "min45" = "45", "min90" = "90",
                                   "min120" = "120", "min135" = "135", "min150" = "150",
                                   "min270" = "270"),
                         expand = c(0,0)) +
        geom_errorbar(aes(x = time, ymin = glumean - glusd, ymax = glumean + glusd),
                          width = 0.3) +
        geom_point(shape = 16, size = 3, aes(color=supplement)) +
        theme_classic()

saveRDS(gluchangeplot, "./data/derdata/analysis_training/glutotalplot.RDS")
## Change data

change_dat.glu <- glu.dat %>%
        dplyr::select(subject, time, glu, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(glu = mean(glu, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = glu) %>%
        ungroup() %>%
        mutate(change.45 = log(min45)-log(baseline),
               change.90 = log(min90)-log(baseline),
               change.120 = log(min120)-log(baseline),
               change.135 = log(min135)-log(baseline),
               change.150 = log(min150)-log(baseline),
               change.270 = log(min270)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.45, change.90, change.120, change.135, 
               change.150, change.270) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.45:change.270)) %>%
        print()

## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

glu.m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.glu)

plot(glu.m1)

summary(glu.m1)

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are fold change values (changeble in the mutate function)


gluconfint.m1 <- confint(emmeans(glu.m1, specs = ~"supplement|time")) %>%
        data.frame()

percentchange <- gluconfint.m1 %>%
        mutate(gluchange = (exp(emmean) -1)*100) %>%
        print()

glucong <- round(gluconfint.m1$upper.CL, 2)

glusd <- round(gluconfint.m1$df, 2)



#----------------------Script under has not been used-------------------------


## Emmeans figure
pos <- position_dodge(width = 0.2)
#fold change figure
#not used
glufoldchangeplot <- gluconfint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before = 2) %>%
        mutate(time.c = gsub("change.", "", time), 
               time.c = as.numeric(time.c)) %>%
        ggplot(aes(time.c, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      size = 0.5,
                      position = pos) +
        geom_line(position = pos,aes(linetype = supplement)) +
        geom_point(shape = 17, position = pos, size = 3, aes(color = supplement)) +
        labs(x = "Minutes", y = "Blood glucose levels \n(fold change mmol/L)\n") +
        scale_x_continuous(limits = c(0, 300), breaks = c(0, 45, 90, 120, 135, 150, 270),
                         expand = expansion(0)) +
        theme_classic() 

saveRDS(glufoldchangeplot, "./data/derdata/analysis_training/glufoldplot.RDS")


## Table
tab1 <- cbind(coef(summary(glu.m1)), data.frame(confint(glu.m1))[1:13, ]) 

# make figures into 1
#not used
glucosefigure <- ggarrange(gluchangeplot, glufoldchangeplot,
                       labels = c("A", "B"),
                       ncol = 1, nrow = 2,
                       align = c("hv"),
                       common.legend = TRUE, legend = "top")

glucosefigure1 <- annotate_figure(glucosefigure,
                                   top = text_grob("Blood glucose", face = "bold"),
                                   fig.lab = "Figur 3", fig.lab.face = "bold")

saveRDS(glucosefigure1, "./data/derdata/analysis_training/glucosefigure_ab.RDS")
