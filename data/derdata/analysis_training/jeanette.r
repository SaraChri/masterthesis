a<<### Styrke Kneekstensjon

# Packages
library(tidyverse)
library(readxl)



# Read data
meta <- read_excel("./data/strength-data.xlsx", sheet = "metadata") %>%
        mutate(subject = as.character(subject)) %>%
        print()

str_ke <- read_excel("./data/strength-data.xlsx") %>%
        separate(sample, into = c("subject", "leg", "time"), sep = "_") %>%
        separate(time, into = c("time", "rep"), sep = "-") %>%
        select(subject:rep, knee_ex) %>%
        
        group_by(subject, leg, time) %>%
        summarise(knee_ex = mean(knee_ex, na.rm = TRUE)) %>%
        inner_join(meta) %>%
        filter(reps != "CTRL") %>%
        select(subject, sex, bmi, leg, time, reps, knee_ex) %>%
        
        ungroup() %>%
        filter(time %in% c("T1", "T2", "T3", "T4")) %>%
        group_by(reps, time) %>%
        summarise(m= mean(knee_ex, na.rm = TRUE), 
                  s= sd(knee_ex, na.rm = TRUE)) %>%
        print()

knee_change <- read_excel("./data/strength-data.xlsx") %>%
        separate(sample, into = c("subject", "leg", "time"), sep = "_") %>%
        separate(time, into = c("time", "rep"), sep = "-") %>%
        select(subject:rep, knee_ex) %>%
        
        
        group_by(subject, leg, time) %>%
        summarise(knee_ex = mean(knee_ex, na.rm = TRUE)) %>%
        inner_join(meta) %>%
        filter(reps != "CTRL") %>%
        select(subject, sex, bmi, leg, time, reps, knee_ex) %>%
        pivot_wider(names_from = time,
                    values_from = knee_ex) %>%
        rowwise() %>% 
        
        mutate(baseline = mean(c(T1, T2), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(ch.t3 = log(T3) - log(baseline), # økning fra baseline til t3
               ch.t4 = log(T4) - log(baseline), # økning fra baseline til t4
               baseline = baseline - mean(baseline, na.rm = TRUE), 
               bmi = bmi - mean(bmi, na.rm = TRUE)) %>%
        print()


#-1)*100) %>%
        

#summarise(m = mean(endring, na.rm = TRUE),
 #         s = sd(endring, na.rm = TRUE )) %>%


### Create model ##########

# load lmerTest
library(lmerTest)

kn_m1 <- lmer(ch.t4 ~ baseline  + reps + (1|subject), 
           data = knee_change) 

plot(kn_m1)

summary(kn_m1)


kn_m2 <- lmer(ch.t4 ~ baseline  + sex + reps + (1|subject), 
              data = knee_change) 

plot(kn_m2)

summary(kn_m2)

kn_m3 <- lmer(ch.t4 ~ baseline  + sex + bmi* reps + (1|subject), 
              data = knee_change) 

plot(kn_m3)

summary(kn_m3)



library(emmeans)
emmeans(kn_m1, specs = ~ "reps")
emmeans(kn_m2, specs = ~ "reps")
emmeans(kn_m3, specs = ~ "reps") 
      #  data.frame() %>%
       # ggplot(aes(reps, emmean)) +
#        geom_point() +
 #       geom_smooth(method = "lm")


#CI intervall, spesifisering av hvilken rader man vil ha
kn_es <- cbind(coef(summary(kn_m1)),
      data.frame(confint(kn_m1))[c(3,4,5), ]) 
kn_est_plot <- kn_es %>%
        mutate(reps="RM30",
               cf= row.names(.),
              Estimate=(exp(Estimate)-1)*100,
              upper=(exp(X97.5..)-1)*100,
              lower=(exp(X2.5..)-1)*100)%>%
        filter(cf=="repsRM30") %>%
        add_row(reps="RM10") %>%
ggplot(aes(reps, Estimate)) +
        geom_hline(yintercept = 0,
                   lty=2,
                   color="gray80")+ #linje 

        geom_point() +
        geom_errorbar(aes(ymin=lower,ymax=upper),
                      width=0.2) +
        scale_y_continuous(limits = c(-12, 5)) +
        theme_classic() 



library(ggplot2)






#em means, estimering i per gruppe


#knee_change %>%
 #       two.group.unpaired <- knee_change %>%
  #      select(subject, reps, ch.t4) %>%
   #     print()
# dabest(timepoint, reps,
  #     idx = c( "reps", "10RM",
   #             "reps", "30RM"), 
    #   paired = FALSE) %>% 
     #   print()

library(dabestr)

sum_knee_change <- knee_change %>%
        group_by(reps)%>% # gruppering 
        summarise(ch.t4 = mean(ch.t4, na.rm= TRUE)) %>% #gjennomsnittet av t4 per gruppe
        mutate(ch.t4=(exp(ch.t4)-1)*100) %>%
        print()

in_data_plot <- ggplot(knee_change, aes(reps, (exp(ch.t4)-1)*100, group = subject)) + #aes, hva som skal være på x og y, gruppering 
        geom_point() + # bare punkter
        geom_line()+ # linjer som kobler sammen samme person
        geom_point(data = filter(sum_knee_change, #gjennomsnittspunkt for 10 RM
                                 reps== "RM10"), 
                   aes(group=NULL), 
                   position= position_nudge(x =-0.2), 
                  
                   shape= 21, fill= "blue", size= 3) + #gjennomsnitt punkt per gruppe
                  geom_point(data = filter(sum_knee_change, # gjennomsnitt punkt for 30RM
                                 reps== "RM30"), 
                   aes(group=NULL), 
                   position= position_nudge(x = 0.2), 
                   
                   shape= 21, fill= "blue", size= 3) +
                theme_classic()+
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank()) 


        library(cowplot)
plot_grid(in_data_plot, kn_est_plot, nrow = 2, 
          align = "v")
        
?position_nudge
        
        
        
        
   #     geom_point(position = position_jitter(width = 0.2))+
   #     ggtitle("Knee-extension") +
    #    labs(x= "Tid", 
     #        y = " Prostentvis Endring, baseline til T4(mean\U00B1 95% CI")+
      #  theme(legend.position= "none") 

knee_change%>%
        
        select(subject, sex, bmi, reps, ch.t4) %>%
        pivot_wider(names_from = reps, values_from = ch.t4) %>%
        mutate(df = RM30 - RM10) %>%
        ggplot(aes(bmi, df, color = sex)) + geom_point() +
        geom_smooth(method = "lm")

             
        