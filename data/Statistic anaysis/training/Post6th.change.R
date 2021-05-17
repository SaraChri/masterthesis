#### Humac analysis pre vs. post sixth RT session


## Author: Kristian Lian ans Sara Moen

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


## Change data

# Isometric
isom.dat <- humac %>%
        filter(test == "isom",
               time %in% c("test3", "test4")) %>%
        print()

change_dat4 <- isom.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

# Isok.60

isok60.dat <- humac %>%
        filter(test == "isok.60",
               time %in% c("test3", "test4")) %>%
        print()

change_dat5 <- isok60.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

## Isok.240

isok240.dat <- humac %>%
        filter(test == "isok.240",
               time %in% c("test3", "test4")) %>%
        print()

change_dat6 <- isok240.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()


## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

# Mean of all subjects 

# Isometric
m4 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat4)
plot(m4)

summary(m4)

# Isok.60
m5 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat5)
plot(m5)

summary(m5)

# Isok.240
m6 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat6)
plot(m6)

summary(m6)

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

# Isometric

confint.m4 <- confint(emmeans(m4, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Isok.60

confint.m5 <- confint(emmeans(m5, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.240

confint.m6 <- confint(emmeans(m6, specs = ~"supplement|time")) %>%
        data.frame()

## Emmeans figure 

pos <- position_dodge(width = 0.2)

## Emmeans % change

#isometric
m4change <- confint.m4 %>%
        mutate(change.m4 = (exp(emmean) -1)*100) %>%
        print()
#isok.60
m5change <- confint.m5 %>%
        mutate(change.m5 = (exp(emmean) -1)*100) %>%
        print()

#isok.240
m6change <- confint.m6 %>%
        mutate(change.m6 = (exp(emmean) -1)*100) %>%
        print()

# Isom

m4.plot <- m4change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement,)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 18, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "30min", "change.3" = "2hr",
                                  "change.4" = "23hr")) +
        labs(x = "", y = "Isometric \n(%log change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Isok 60

m5.plot <- m5change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 18, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "30min", "change.3" = "2hr",
                                  "change.4" = "23hr")) +
        labs(x = "", y = "Isokinetic 60 \n(%log change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Isok 240

m6.plot <- m6change %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype =  supplement)) +
        geom_point(shape = 18, position = pos, size = 3, aes(color = supplement)) +
        scale_x_discrete(labels=c("change.1" = "B", "change.2" = "30min", "change.3" = "2hr",
                                  "change.4" = "23hr")) +
        labs(x = "", y = "Isokinetic 240 \n(%log change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


#take all plots into 1 plot, names alphabetically
isopost <- ggarrange(m5.plot, m6.plot, m4.plot,
                       labels = c("A", "B", "C"),
                       ncol = 2, nrow = 2,
                     align = c("hv"),
                       common.legend = TRUE, legend = "top")

isopostfigure <- annotate_figure(isopost,
                top = text_grob("Isokinetic and isometric post test", face = "bold"),
                fig.lab = "Figure 6", fig.lab.face = "bold")


saveRDS(isopostfigure, "./data/derdata/analysis_training/isopostchange.RDS")
