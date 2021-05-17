##### PROTEIN MURF ANALYSIS ###############


library(tidyverse)
library(readxl)
library(lmerTest)
library(emmeans)
library(pbkrtest)
library(ggplot2)
library(ggpubr)
library(dabestr)
library(pbkrtest)
library(scales)

dat2 <- read_excel("./data/results/Ribose_proteinanalyse.xlsx") %>%
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
        print()

#plot pre-post change

murfplot <- dat2 %>%
        select(subject, supplement, murf_pre, murf_post) %>%
        pivot_longer(names_to = "time",
                     values_to = "murf",
                     cols = (murf_pre:murf_post)) %>%
        mutate(time = gsub("murf_", "", time),
               time = factor(time, levels = c("pre", "post"))) %>%
        group_by(time, supplement) %>%
        print()

murfm2 <- lmerTest::lmer(murf ~ time + supplement + supplement:time + (1|subject),
                     data = murfplot)
plot(murfm2)

summary(murfm2)

confint.murf<- confint(emmeans(murfm2, specs = ~"time|supplement")) %>%
        data.frame() %>% 
        print()

#Mean protein conten
murfemmeans <- confint.murf %>%
  mutate(murfchange = (exp(emmean) -1)*100) %>%
  print()

#fixed position
pos <- position_dodge(width = 0.2)

# total MuRF1 protein concentrtion between groups from pre-post
murfm2.plot <- confint.murf %>%
        ggplot(aes(time, emmean, group = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos, aes(linetype = supplement)) +
        geom_point(shape = 16, position = pos, size = 3, aes(color = supplement)) +
  scale_color_manual(values = c("turquoise3", "coral1")) +
        labs(x = "Time", y = "Murf \n(Total)\n", fill = "Supplement") +
        theme_classic()

### logfold change and summary between conditions in post test
murfchange3 <- dat2 %>% 
        print()

murf.m4 <- lmer(murf_change ~ murf_pre + supplement + (1|subject), 
                data = murfchange3)
plot(murf.m4)

summary(murf.m4)

#get put confints
murf_ci <- cbind(coef(summary(murf.m4)),
      confint(murf.m4)[c(3:5),])

#pvaluse
murf_pval <- data.frame(murf_ci)[3,5]
murf_pval <- scales::pvalue(murf_pval)

confint.murfchange <- confint(emmeans(murf.m4, specs = ~"supplement")) %>%
        data.frame() %>%
        print()

percentchange <- confint.murfchange %>%
  mutate(murfchange = (exp(emmean) -1)*100) %>%
  print()



### fold change change between groups
#Glucose has a fold change on 22% and glucose decreases -5%

murfcahngeplotten <- confint.murfchange%>%
        ggplot(aes(supplement, emmean, group = supplement)) + 
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_point(shape = 17, position = pos, size = 3, aes(color = supplement)) + 
  scale_color_manual(values = c("turquoise3", "coral1")) +
        labs(x = "Post", y = "MuRF1 \n(fol change)\n", fill = "Supplement") +
        theme_classic()
        
#plot on ci between groups.
#plot log change, gluce 26% lower than placebo
murf_ciplot <- murf_ci %>%
        data.frame() %>%
        mutate(coef = rownames(.)) %>%
        filter(coef == "supplementglucose") %>%
        ggplot(aes(x= 1, y = Estimate)) +
        geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error,
                           width = 0.2,)) +
                geom_point(shape = 18, size = 4, color = "yellow2") +
        scale_y_continuous(position = "right") + 
        scale_x_continuous(limits = (1), breaks = c(1)) + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        labs(x = "", y = "MuRF1 % diff", position = "right") +
        theme_classic()
#---place figures in one---#

Murffinal <- ggarrange(murfm2.plot, 
                       ggarrange(murfcahngeplotten, murf_ciplot,
                                 ncol = 2, labels = c("", ""),
                                 font.label = list(),
                                 legend = FALSE),
                       nrow = 2,
                       labels = "",
                       font.label = list(),
                       legend = "right")
                       
Murftwofigure <- annotate_figure(Murffinal,
                                   top = text_grob("Murf", face = "bold"),
                                   fig.lab = "Figure 7", fig.lab.face = "bold")

saveRDS(Murftwofigure , "./data/derdata/analysis_murf/murffigures.RDS")



