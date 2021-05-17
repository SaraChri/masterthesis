#load packages
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




#protein results
#UBF
#MURF
#only pre-post
proteinanalyse <- read_excel("data/results/Ribose_proteinanalyse.xlsx")
        
        dat <- read_excel("data/results/Ribose_proteinanalyse.xlsx") %>%
                # Group data set per row to allow rowwise average
                rowwise() %>%
                # Calculate average of total protein and background
                mutate(total.protein = mean(c(total.protein1, total.protein2)), 
                       background = mean(c(background1, background2))) %>%
                ungroup() %>%
                
                
                
                #  On UNGROUPED data, remove background from protein 
                select(gel, well, sample,subject, supplement,timepoint, total.protein, background, ubf, murf) %>%
                mutate(total.protein = total.protein - background) %>%
                
                # Group by gel and calculate relative signal for bakground and 
                # each target protein
                group_by(gel) %>%
                mutate(total.protein = total.protein / mean(total.protein, na.rm = TRUE), 
                       ubf = ubf / mean(ubf, na.rm = TRUE), 
                       murf = murf / mean(murf, na.rm = TRUE), 
                       ubf = ubf / total.protein, 
                       murf = murf / total.protein) %>%
                
                # Average over multiple samples (UBF) 
                group_by(sample, supplement, subject, timepoint) %>%
                summarise(ubf = mean(ubf, na.rm = TRUE), 
                          murf = mean(murf, na.rm = TRUE)) %>%
                # Fix timepoint info:
                # remove rna1/2 info,
                # change T1 and T2 to pre, and T3/T4 to post
                mutate(timepoint = gsub("rna1", "", timepoint), 
                       timepoint = gsub("rna2", "", timepoint), 
                       timepoint = if_else(timepoint %in% c("T1","T2"), 
                                           "pre", 
                                           "post")) %>%
                ungroup() %>%
                select(-sample) %>%
                
                pivot_wider(names_from = timepoint, 
                            values_from = c(ubf, murf)) %>%
                
                # Add change scores on the log scale
                # Change order of glucose/placebo
                mutate(ubf_change = log(ubf_post) - log(ubf_pre), 
                       murf_change = log(murf_post) - log(murf_pre), 
                       supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        arrange(subject) %>%
                print()
        
        ### Exploratory data 
        
        dat.ubf <- dat %>%
                filter(ubf_change < -1) %>%
                print()
        
        ggplot(dat.ubf, aes(supplement, ubf_change)) + geom_point()
        
        
        
        
        
        
        
        ubf.m1 <- lmer(ubf_change ~ ubf_pre + supplement + (1|subject), 
                       data = dat)
        
        plot(ubf.m1)
        
        summary(ubf.m1)
        
        
        emmeans(ubf.m1, specs = ~ "supplement") %>%
                confint()
        
        
        
        murf.m1 <- lmer(murf_change ~ murf_pre + supplement + (1|subject), 
                        data = dat)
        
        plot(murf.m1)
        
        summary(murf.m1)
        
        
        emmeans(murf.m1, specs = ~ "supplement") %>%
                confint()

dat.murf <- dat %>%
        filter(murf_change < -1) %>%
        print()

plot.murf <- ggplot(dat.murf, aes(supplement, murf_change)) +
        geom_point()
        