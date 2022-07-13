# Set working directory
setwd("/home/dima/Desktop/R/FCS")

library(tidyverse)
library(RColorBrewer)
library(flowCore)
library(flowAI)
library(ggcyto)
library(listviewer)
library(flowAI)
library(gridExtra)
library(CytoExploreR)
library(CytoExploreRData)
library(ggpubr)

my_theme <- theme_bw() +
  theme(axis.text.x = element_text(size=14, angle=0), 
        strip.text.x = element_text(size = 14, angle = 90),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_line(colour = "gray")
        )   

# Reading the file
myfile <- list.files(path="/home/dima/Desktop/R/FCS/", pattern=".fcs")
my_fcs <- flowCore::read.flowSet(myfile, path="/home/dima/Desktop/R/FCS")

lgcl <- logicleTransform( w = 0.5, t= 10000, m =5)
trans <- transformList(colnames(my_fcs[,5:12]), lgcl)
after <- transform(my_fcs, trans)
autoplot(after[[1]])
###stop here to look at autoplot


PIvsFSCA <- cbind(after[[1]]@exprs[,"FSC-A"], 
                  after[[1]]@exprs[,"FL2-A"]) %>% as_tibble()

colnames(PIvsFSCA) <- c("FSC_A", "PE_A")

#take Log Scale for FSC_A
PIvsFSCA <- PIvsFSCA %>% mutate(log_fsc = log2(FSC_A))
max_fsc <- max(PIvsFSCA$log_fsc)
min_fsc <- min(PIvsFSCA$log_fsc)

##
threshold <- 6

#Calculate survival in moving window
stat_av <- PIvsFSCA %>% 
  mutate(fsc_window = round(log_fsc, 1)) %>%
  group_by(fsc_window) %>%
  summarise(dead = sum(PE_A > threshold), n = n()) %>%
  mutate(survival = dead/n) %>%
  as_tibble() %>%
  dplyr::filter(n > 100)


###Plotting
fig1 <- PIvsFSCA %>% ggplot(aes(x = log_fsc, y = PE_A)) + 
  geom_point(alpha = 0.05) + 
  geom_hline(yintercept = threshold, lty = 3, color = "red") +
  xlab("Log2(FSC-A)") +
  ylab("PI fluorescence, a.u.") +
  my_theme

fig2 <- stat_av %>% ggplot(aes(x = fsc_window, y = survival)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  ylim(0,1) +
  xlim(min_fsc, max_fsc) +
  xlab("Log2(FSC-A)") +
  ylab("Survival") +
  stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.8, aes(label = ..rr.label..)) +
  my_theme

ggarrange(fig1, fig2, ncol = 1, align = "v")
  
lm(stat_av$survival ~ stat_av$fsc_window)