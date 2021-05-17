### Total volume change, both absolute values and fold-change

# Packages
library(readxl)
library(tidyverse)
library(knitr)
library(lme4)
library(broom)
library(emmeans)
library(dabestr)
library(magrittr)
library(ggpubr)


# Data
tot.vol <- read_excel("data/results/ribose_volume.xlsx")  %>%
        # take out variables in interest, to find change of  total volume
        select(subject, timepoint, tot.volume, supplement) %>%
        print()

# Data handling
#pair so glucose and place is on the same timepoint
tot.volh <- tot.vol %>%
        mutate(time = if_else(timepoint %in% c("T1", "T2"),
                              "baseline",
                              if_else(timepoint %in% c("D3", "D4"),
                                      "session2",
                                      if_else(timepoint %in% c("D5", "D6"),
                                              "session3",
                                              if_else(timepoint %in% c("D7", "D8"),
                                                      "session4",
                                                      if_else(timepoint %in% c("D9", "D10"),
                                                              "session5",
                                                              if_else(timepoint %in% c("T3", "T4"),
                                                                      "session6", timepoint))))))) %>%
        mutate(time = factor(time, levels = c("baseline", "session1", "session2", "session3", "session4", "session5", "session6")),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
print()


## Absolute data for volume change pr. group pr. time in mean and sd
vol.exp <- tot.volh %>%
        select(subject, supplement, time, tot.volume) %>%
        group_by(time, supplement) %>%
        summarise(mean.vol = mean(tot.volume),
                  sd.vol = sd(tot.volume)) %>%
        arrange(supplement) %>%
        print()

# fix position
pos <- position_dodge(width = 0.2)

#line plot of mean total volume each group + 
#use this one, showing total volume change in each group.
tot.volineplot <- ggplot(vol.exp, aes(group = supplement, y = mean.vol, x = time)) +
        geom_line(aes (color=supplement)) + geom_point(aes(color=supplement)) +
        scale_color_manual(values = c("turquoise3", "coral1")) +
        labs(x = "", y = "Training volume \n(kg)\n ", fill = "Supplement" ) + 
        scale_x_discrete(labels=c("baseline" = "B", "session2" = "2", "session3" = "3",
                                  "session4" = "4", "session5" = "5", "session6" = "6")) +
        scale_y_continuous(limits = c(2500, 10000)) +
        geom_errorbar(aes(x = time, ymin = mean.vol - sd.vol, ymax = mean.vol + sd.vol, color = factor(supplement),
                          width = 0.2)) +
        theme_classic()

saveRDS(tot.volineplot, "./data/derdata/analysis_training/totvol_lineplot.RDS")

# Fold-change
totchange_dat <- tot.volh %>%
        dplyr::select(subject, time, tot.volume, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(tot.volume = mean(tot.volume, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = tot.volume) %>%
        # print()
        
        ungroup() %>%
        mutate(change.2 = log(session2)-log(baseline),
               change.3 = log(session3)-log(baseline),
               change.4 = log(session4)-log(baseline),
               change.5 = log(session5)-log(baseline),
               change.6 = log(session6)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4, change.5, change.6) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.6)) %>%
        print()


# Fold-change model
totvol.m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                            data = totchange_dat)
plot(totvol.m1)

summary(totvol.m1)

totvolconfint.m1 <- confint(emmeans(totvol.m1, specs = ~"supplement|time")) %>%
        data.frame()

#change to %
m1totvolchange <- totvolconfint.m1 %>%
        mutate(totvolconfint.m1 = (exp(emmean) -1)*100) %>%
        print()


