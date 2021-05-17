#### Humac analysis pre vs. post fifth RT session


## Author: Kristian Lian and Sara Moen

# Purpose: This script plots mean torque per supplement (both through intervention and pre vs. post) results from the ribose project, 
# and analyses the data per test (isometric, isokinetic 60, isokinetic 240) in a linear model.

## Time-points
# D-1: Baseline, before any supplementation or training
# D4, D5, D8 and D9: Day 4, 5, 8 and 9 of the intervention, humac testing of the leg that performed
# RT the preceding day
# T3: Post testing leg #1 (leg that started the intervention). Leg #1 is tested four times at T3/T4:
# Test 1 leg 1: 1.5hrs after protein ingestion, 45min before RT (T3)
# Test 2 leg 1: 30min after RT (T3)
# Test 3 leg 1: 2hrs after RT (T3)
# Test 4 leg 1: ~23hrs after RT (T4)
# Test 1 serve as a post test for the 5 RT sessions and pre test before the sixth session, test 2,
# 3, and 4 serve as post test following sixth session
# T4 and 13 follow the same design for leg #2

## Data
# Date of testing
# Subject
# Test type: isok.60 (isokinetic 60), isok.240 (isokinetic 240), isom (isometric)
# Peak.torque: Highest peak torque from each test
# Leg: left or right leg
# Supplement: glucose or placebo

# Packages
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(broom)
library(knitr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(scales)


## Data handling

humac <- read_excel("data/results/ribose.humac.xlsx") %>%
        mutate(time = if_else(timepoint == "D-1", 
                              "baseline", 
                              if_else(timepoint %in% c("D4", "D5"), 
                                      "test1", 
                                      if_else(timepoint %in% c("D8", "D9"), 
                                              "test2", 
                                              if_else(timepoint %in% c("T3", "T4") & acute %in% c("rest", "post30min", "post2h"),
                                                      "test3",
                                                      if_else(acute == "post23h", "test4", timepoint)))))) %>%
        mutate(time = factor(time, levels = c("baseline", "test1", "test2", "test3", "test4")), 
               acute = factor(acute, levels = c("rest", "post30min", "post2h", "post23h")), 
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        print()

rest.dat <- humac %>%
        filter(acute == "rest" ) %>%
        print()

# Figure, exploratory

pos <- position_dodge(width = 0.2)

rest.dat <- humac %>%
        filter(acute == "rest" ) %>%
        print()


## Change data

# Isometric
isom.dat <- rest.dat %>%
        filter(test == "isom") %>%
        print()

change_dat <- isom.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

# Isok.60

isok60.dat <- rest.dat %>%
        filter(test == "isok.60") %>%
        print()

change_dat2 <- isok60.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

## Isok.240

isok240.dat <- rest.dat %>%
        filter(test == "isok.240") %>%
        print()

change_dat3 <- isok240.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

#Average total peak data

#Isom
isom.peak<- humac %>%
        filter(test == "isom") %>%
        print()

isom.peak %>%
        select(subject, time, test, peak.torque, supplement) %>%
        group_by(time, supplement) %>%
        summarise(m =mean(peak.torque),
               s =sd(peak.torque)) %>%
        print()


## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

# Mean of all subjects 

# Isometric
isok.m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat)
plot(isok.m1)

summary(isok.m1)

# Ism.60
ism.m2 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat2)
plot(ism.m2)

summary(ism.m2)


# Isok.240
ism.m3 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat3)
plot(ism.m3)

summary(ism.m3)

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

# Isometric

confint.m1 <- confint(emmeans(isok.m1, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.60

confint.m2 <- confint(emmeans(ism.m2, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.240

confint.m3 <- confint(emmeans(ism.m3, specs = ~"supplement|time")) %>%
        data.frame()

#emmeans % change

#isometric
m1change <- confint.m1 %>%
        mutate(change.m1 = (exp(emmean) -1)*100) %>%
        print()
#isok.60
m2change <- confint.m2 %>%
        mutate(change.m2 = (exp(emmean) -1)*100) %>%
        print()

#isok.240
m3change <- confint.m3 %>%
        mutate(change.m3 = (exp(emmean) -1)*100) %>%
        print()


## Emmeans figures

# Isom
m1.plot <- m1change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 16, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isometric \n(%log change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))

# Isok 60

m2.plot <- m2change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 16, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isokinetic 60 \n(%log change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))


# Isok 240

m3.plot <- m3change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 16, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isokinetic 240 \n(%log change)\n", fill = "Supplement",
             labels=c("-0.15" ="15 %", "-0.10" = "10 %", "-0.05" = "5 %",
                      "0" = "0 %", "0.05" = "5 %", "0.10" = "10 %")) +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))

#take all plots into 1 plot, names alphabetically
humacfold <- ggarrange(m2.plot, m3.plot, m1.plot,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2,
                    align = c("hv"),
                    common.legend = TRUE, legend = "top")

humactimefigure <- annotate_figure(humacfold,
                                 top = text_grob("Isokinetic and isometric time effect", face = "bold"),
                                 fig.lab = "Figure 5", fig.lab.face = "bold")


saveRDS(humactimefigure, "./data/derdata/analysis_training/humactimefig.RDS")

