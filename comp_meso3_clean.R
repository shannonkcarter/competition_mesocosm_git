#############################################################
#############################################################
###---------------------PREP WORKSPACE--------------------###
#############################################################
#############################################################

## Set working directory and clear memory
setwd("C:\\Users\\Shannon\\Google Drive\\CompetitionMesocosm_S18_cleaning")
rm(list = ls(all = T))

## Load required packages
library(lme4)
library(car)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(multcomp)
library(wesanderson)

## Load universal plotting elements
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.minor = element_blank(), 
                 panel.grid.major = element_blank(),
                 axis.text  = element_text(size = rel(1.3),colour = "black"),
                 axis.title = element_text(size = rel(1.5)),
                 axis.line  = element_line(colour = "black"))

#############################################################
###---------------------DATA PROCESSING-------------------###
#############################################################

###---LOAD AND PROCESS RAW DATA------------------------------

## Metamorph timing and size data
metamorphs_raw <- read.csv("raw_metamorphs.csv", header = T)

## Format date variable
metamorphs_raw$date <- as.Date(strptime(metamorphs_raw$date,"%d-%b"))
metamorphs_raw$doy  <- as.numeric(strftime(metamorphs_raw$date, format = "%j"))

## Add time to emergence for each frog
metamorphs_raw$timetoemerge <- metamorphs_raw$doy - 105

write.csv(metamorphs_raw, "individual_results.csv")

###---MAKE METAMORPH ACCUMULATION DATASET--------------------

## Add a column of ones to sum as a counting mechanism
ind_results$ones <- 1

## Function to replace NA masses with global mean mass for biomass counts
biomassxNA <- function(x) {
  x[which(is.na(x))] = mean(ind_results$mass_mg, na.rm = T)
  return(cumsum(x))
}

## For each tank, add a cumulative metamorph number and biomass
meta_accum <- ind_results %>%
  group_by(tank) %>%
  mutate(cumulative_metas = cumsum(ones),
            cumulative_biom  = biomassxNA(mass_mg))
meta_accum <- select(meta_accum, date, tank, block, trt, order, sync, ones,
                     mass_mg, doy, timetoemerge, cumulative_metas, cumulative_biom)

## One line per metamorph collected
write.csv(meta_accum, "metamorph_accumulation.csv")

###---FINAL SURVIVORS & BIOMASS BY TREATMENT---------------------------

## Summarize to show final (i.e., max) metamorphs and biomass in each tank
tank_results <- meta_accum %>%
  group_by(tank) %>%
  summarize(meta = max(cumulative_metas),
            biom = max(cumulative_biom),
            mass = mean(mass_mg, na.rm = T),
            emer = mean(timetoemerge),
            emsd = sd(timetoemerge))

## Add tank treatment info back in
treatments <- read.csv("raw_tanks.csv", header = T)
tank_results <- inner_join(tank_results, treatments, by = 'tank')

## Select relevant columns and make 1 row per tank
tank_results <- unique(select(tank_results, tank, block, trt, order, sync, 
                           meta, biom, mass, emer, emsd))

## Add rana survival
teardown <- read.csv("raw_teardown.csv", header = T)   # numbers of all organisms removed at the end
teardown_rana <- select(teardown, tank, rana_surv)     # select just relevant variable to make the join simpler
tank_results <- inner_join(tank_results, teardown_rana, by = "tank")
tank_results$rana <- tank_results$rana_surv/30

## Add proportion survivors (45 tads introduced per tank)
tank_results$surv <- as.numeric(tank_results$meta/45)

## Reorder the variables in an intuitive way and drop unnecessary
tank_results <- select(tank_results, tank, block, trt, order, sync,
                       surv, biom, mass, emer, emsd, rana)

## One line per tank with each response variable
write.csv(tank_results, "tank_results.csv")

###---MEAN AND SE HELPER FUNCTION----------------------------

## if needing to re-run this, do so with no packages loaded, because major interference w/ tidyverse
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  #detach(package:dplyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###---GET MEANS AND SES FOR ALL VARIABLES---------------------

## SURVIVAL
means_surv <- summarySE(data = tank_results, measurevar = "surv", 
                        groupvars = c("order", "sync"))
write.csv(means_surv, "means_surv.csv")

## BIOMASS
means_biom <- summarySE(data = tank_results, measurevar = "biom", 
                        groupvars = c("order", "sync"))
write.csv(means_biom, "means_biom.csv")

## PER CAPITA MASS
means_mass <- summarySE(data = ind_results, measurevar = "mass_mg", 
                        groupvars = c("order", "sync"), na.rm = T)
write.csv(means_mass, "means_mass.csv")

## EMERGENCE TIME
means_emer <- summarySE(data = ind_results, measurevar = "timetoemerge", 
                        groupvars = c("order", "sync"))
write.csv(means_emer, "means_emer.csv")

## EMERGENCE S.D.
means_emsd <- summarySE(data = tank_results, measurevar = 'emsd', 
                        groupvars = c('order', 'sync'))
write.csv(means_emsd, 'means_emsd.csv')

## RANA SURVIVAL
means_rana <- summarySE(data = subset(tank_results, subset = (order != 'cont')), 
                        measurevar = 'rana', groupvars = c('order', 'sync'))
write.csv(means_rana, 'means_rana.csv')
###---ADD CONTROL ADJ MEANS TO TANK_RESULTS------------------

## First, subset trt_means, which already has control means
trt_means_control <- subset(trt_means, subset = (order == 'cont'))
trt_means_control <- select(trt_means_control, sync, variable, mean)

## Now, make it wide format to match tank_results
trt_means_control <- trt_means_control %>%
  spread(key = variable, value = mean) %>%
  rename(surv_cont = surv, biom_cont = biom, mass_cont = mass, emer_cont = emer, emsd_cont = emsd)

## Now add control means to each tank
tanks_cont <- merge(x = treatments, y = trt_means_control, by = 'sync', all = T)

## And now match tanks to tank_results
tank_results_contadj <- merge(x = tanks_cont, 
                              y = subset(tank_results, subset = (order != 'cont')),
                              by = 'tank', all.y = T)

## Now, subtract controls from means to get relative means
tank_results_contadj <- tank_results_contadj %>%
  mutate(surv_adj = surv - surv_cont,
         biom_adj = biom - biom_cont,
         mass_adj = mass - mass_cont,
         emer_adj = emer - emer_cont,
         emsd_adj = emsd - emsd_cont) %>%
  select(tank, block.x, trt.x, order.x, sync.x, surv_adj, biom_adj, emer_adj, mass_adj, emsd_adj) %>%
  rename(block = block.x, trt = trt.x, order = order.x, sync = sync.x)

## Finally, write a csv or add to existing survival csv
write.csv(tank_results_contadj, 'tank_results_contadj.csv')

#############################################################
###-----------------------LOAD DATA-----------------------###
#############################################################

###---LOAD DATA----------------------------------------------

## Treatment assignments for each tank
treatments <- read.csv("raw_tanks.csv", header = T)

## Tidy-ed raw metamorph data - 1 row per metamorph collected
ind_results <- read.csv("individual_results.csv", header = T)
ind_results$date <- as.Date(ind_results$date)

## Metamorph accumulation data (by tank and day)
meta_accum <- read.csv("metamorph_accumulation.csv", header = T)
meta_accum$date <- as.Date(meta_accum$date)
#meta_accum_alldates <- read.csv("metamorph_accumulation_alldates.csv")

## Tank results (1 line per tank for 5 key response variables)
tank_results <- read.csv("tank_results.csv", header = T)

## Tank results adjusted to appropriate mean (1 line per tank for 5 key response variables)
tank_results_contadj <- read.csv('tank_results_contadj.csv', header = T)

## Means and SEs for key variables
trt_means <- read.csv("trt_means.csv", header = T)  # made by copy-pasting above dfs together manually

###---REDUNDANT VARIABLE FORMATTING--------------------------

## Put factor levels in intuitive order for all dfs
ind_results$sync  <- factor(ind_results$sync, levels = c("low", "med", "high"))
ind_results$order <- factor(ind_results$order, levels = c('early', 'same', 'late', 'cont'))

meta_accum$sync  <- factor(meta_accum$sync, levels = c("low", "med", "high"))
meta_accum$order <- factor(meta_accum$order, levels = c('early', 'same', 'late', 'cont'))

tank_results$sync  <- factor(tank_results$sync, levels = c("low", "med", "high"))
tank_results$order <- factor(tank_results$order, levels = c('early', 'same', 'late', 'cont'))

trt_means$sync  <- factor(trt_means$sync, levels = c("low", "med", "high"))
trt_means$order <- factor(trt_means$order, levels = c('early', 'same', 'late', 'cont'))

#############################################################
#############################################################
###-------------------------PLOTS-------------------------###
#############################################################
#############################################################

###---SURVIVAL (HYLA) ----------------------------------------

### CHOSEN MODEL ###
msurv_rel <- lmer(data = tank_results_contadj, surv_adj ~ sync * order + (1 | block))
Anova(msurv_rel)

### RAW MEANS ###
surv_means <- ggplot(subset(trt_means, subset = (variable == 'surv' & order != 'cont')),
                     aes(x = order, y = mean, group = sync, fill = sync, color = sync, shape = sync)) + mytheme +
  geom_hline(yintercept = 0.58, size = 1, linetype = 'dashed', color = "#fbb973") +
  geom_hline(yintercept = 0.69, size = 1, linetype = 'dashed', color = "#f96968") +
  geom_hline(yintercept = 0.62, size = 1, linetype = 'dashed', color = "#4b1d1d") +
  #geom_point(data = survivors_nocontrol, aes(x = order, y = prop_surv, group = sync, fill = sync, shape = sync),
  #           position = position_dodge(width = 0.1)) + 
  geom_line(size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  ylab('proportion survival') +
  xlab('mean Hyla hatching (relative to Rana)') +
  labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony', color = "hatching\nsynchrony") +
  theme(legend.position = c(0.1, 0.22),
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))  

### MEANS RELATIVE TO CONTROL ###
surv_means_rel <- ggplot(subset(trt_means, subset = (variable == 'surv' & order != 'cont')), 
                         aes(x = order, y = adj, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed', color = "black") +
  geom_line(size = 1, position = position_dodge(width = 0.1)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  ylab(expression(atop("proportion survival", paste(" (", Delta, " control)")))) +
  xlab(NULL) +
  labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = c(0.1, 0.22),
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_shape_manual(values = c(21, 24, 23))
surv_means_rel

surv_means_rel1 <- ggplot(subset(trt_means, subset = (variable == 'surv' & order != 'cont')), 
                         aes(x = order, y = adj, shape = sync, fill = sync, color = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed', color = "black") +
  geom_line(size = 0.5, alpha = 0.8, position = position_dodge(width = 0.1)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_point(size = 5, color = "black", position = position_dodge(width = 0.1)) +
  ylab(expression(atop("proportion survival", paste(" (", Delta, " control)")))) +
  xlab('mean Hyla hatching (relative to Rana)') +
  #xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony', color = 'hatching\nsynchrony') +
  theme(legend.position = c(0.05, 0.22),
    #legend.position = c(0.2, 0.22),
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  #scale_fill_brewer(palette = 'Dark2') +
  #scale_color_brewer(palette = 'Dark2') +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  #scale_color_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23))
surv_means_rel1 

###---BIOMASS EXPORT (HYLA)----------------------------------

### CHOSEN MODEL ###
mbiom_rel <- lmer(data = tank_results_contadj, biom_adj ~ sync * order + (1 | block))
Anova(mbiom_rel)

### RAW MEANS ###
biom_means <- ggplot(subset(trt_means, subset = (variable == 'biom' & order != 'cont')), 
                     aes(x = order, y = mean, color = sync, shape = sync, fill = sync, group = sync)) + mytheme +
  #geom_point(data = subset(survivors, subset = (order != 'cont')), aes(x = order, y = biomass, fill = sync, shape = sync),
  #           position = position_dodge(width = 0.1)) + 
  geom_hline(yintercept = 6127, size = 1, linetype = 'longdash', color = "#fbb973") +
  geom_hline(yintercept = 6103, size = 1, linetype = 'dashed', color = "#f96968") +
  geom_hline(yintercept = 5035, size = 1, linetype = 'dashed', color = "#4b1d1d") +
  geom_line(size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  ylab("total biomass export (mg)") + 
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

### MEANS RELATIVE TO CONTROL ###
biom_means_rel <- ggplot(subset(trt_means, subset = (variable == 'biom' & order != 'cont')), 
                         aes(x = order, y = adj/1000, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 0, size = 0.7, linetype = 'dashed', color = "black") +
  geom_line(size = 0.5, alpha = 0.8, position = position_dodge(width = 0.1), aes(color = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(color = sync, ymin = (adj - se)/1000, ymax = (adj + se)/1000, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  ylab(expression(paste("total biomass export (g, ", Delta, " control)"))) +
  ylab(expression(atop("total biomass export", paste(" (g, ", Delta, " control)")))) +
  #xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

###---PER CAPITA MASS (HYLA)---------------------------------

### CHOSEN MODEL ###
mmass_rel <- lmer(data = tank_results_contadj_27, mass_adj ~ sync * order + (1 | block))
Anova(mmass_rel)

### RAW MEANS ###
mass_means <- ggplot(subset(trt_means, subset = (variable == 'mass' & order != 'cont')), 
                     aes(x = order, y = log(mean), color = sync, group = sync, fill = sync, shape = sync)) + mytheme +
  geom_hline(yintercept = log(236.5), size = 1, linetype = 'dashed', color = "#fbb973") +
  geom_hline(yintercept = log(197.4), size = 1, linetype = 'dashed', color = "#f96968") +
  geom_hline(yintercept = log(177.2), size = 1, linetype = 'dashed', color = "#4b1d1d") +
  geom_line(size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = log(mean - se), ymax = log(mean + se), width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  ylab("log per capita mass (mg)") + 
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

### MEANS RELATIVE TO CONTROL ###
mass_means_rel <- ggplot(subset(trt_means, subset = (variable == 'mass' & order != 'cont')), 
                         aes(x = order, y = adj, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed', color = "black") +
  geom_line(size = 0.5, alpha = 0.8, position = position_dodge(width = 0.1), aes(color = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(color = sync, ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  #ylab(expression(paste("per capita mass (mg, ", Delta, " control)"))) +
  ylab(expression(atop("per capita mass", paste(" (mg, ", Delta, " control)")))) +
  #xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

###---EMERGENCE MEAN (HYLA)----------------------------------

### CHOSEN MODEL ###
memer_rel <- lmer(data = tank_results_contadj, emer_adj ~ sync * order + (1 | block))
Anova(memer_rel)

### RAW MEANS ###
emer_means <- ggplot(subset(trt_means, subset = (variable == 'emer' & order != 'cont')), 
                     aes(x = order, y = mean, color = sync, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 33.4, size = 1, linetype = 'dashed', color = "#fbb973") +
  geom_hline(yintercept = 43.8, size = 1, linetype = 'dashed', color = "#f96968") +
  geom_hline(yintercept = 55.8, size = 1, linetype = 'dashed', color = "#4b1d1d") +
  geom_line(size = 0.7, position = position_dodge(width = 0.1)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  ylab("days to emergence") + 
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

### MEANS RELATIVE TO CONTROL ###
emer_means_rel <- ggplot(subset(trt_means, subset = (variable == 'emer' & order != 'cont')), 
                         aes(x = order, y = adj, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed', color = "black") +
  geom_line(size = 0.5, alpha = 0.8, position = position_dodge(width = 0.1), aes(color = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(color = sync, ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  #ylab(expression(paste("mean days to emergence (", Delta, " control)"))) +
  ylab(expression(atop("mean days to emergence", paste(" (", Delta, " control)")))) +
  #xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_shape_manual(values = c(21, 24, 23))
emer_means_rel

###---EMERGENCE S.D. (HYLA)----------------------------------

## CHOSEN MODEL
library(nlme)
memersd_rel4 <- lme(data = surv_contadj,
                    emersd_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | sync)))

### RAW MEANS ###
emsd_means <- ggplot(subset(trt_means, subset = (variable == 'emsd' & order != 'cont')), 
                        aes(x = order, y = mean, color = sync, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_hline(yintercept = 10.3, size = 1, linetype = 'dashed', color = "#fbb973") +
  geom_hline(yintercept = 16.88, size = 1, linetype = 'dashed', color = "#f96968") +
  geom_hline(yintercept = 15.11, size = 1, linetype = 'dashed', color = "#4b1d1d") +
  geom_line(size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  ylab("s.d., days to emergence") + 
  xlab('mean Hyla hatching (relative to Rana)') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23))

### MEANS RELATIVE TO CONTROLS ###
emsd_means_rel <- ggplot(subset(trt_means, subset = (variable == 'emsd' & order != 'cont')), 
                            aes(x = order, y = adj, shape = sync, fill = sync, group = sync)) + #mytheme +
  geom_hline(yintercept = 0, size = 1, linetype = 'dashed', color = "black") +
  geom_line(size = 0.5, alpha = 0.8, position = position_dodge(width = 0.1), 
            aes(color = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(color = sync, ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  #ylab(expression(paste("s.d., days to emergence (", Delta, " control)"))) +
  ylab(expression(atop("s.d., days to emergence", paste(" (", Delta, " control)")))) +
  xlab('mean Hyla hatching (relative to Rana)') +
  #xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13),
        axis.text  = element_text(size = 12)) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_shape_manual(values = c(21, 24, 23))
emsd_means_rel

###---RANA SURVIVAL------------------------------------------

### CHOSEN MODEL ###
mrana1 <- glmer(data = subset(tank_results, subset = (order != 'cont')), 
                cbind(rana * 30, 30 - rana * 30) ~ order * sync + (1 | block), family = 'binomial')
Anova(mrana1)

### RAW MEANS ###
rana_means <- ggplot(subset(trt_means, subset = (variable == 'rana')),
                     aes(x = order, y = mean, group = sync, fill = sync, shape = sync)) + mytheme +
  geom_line(size = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = (mean - se), ymax = (mean + se), width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  xlab(expression(paste('mean ',  italic(' H. versicolor'),' hatching (relative to',  italic(' R. sphenocephala'),")"))) +
  labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  ylab(expression(paste('proportion survival ',  italic(' R. sphenocephala'),))) +
  theme(legend.position = c(0.1, 0.82),
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23))  

### CORRELATION BETWEEN RANA AND HYLA SURVIVAL ###
# no correlation = grounds to distrust rana survival data
hyla_rana <- ggplot(tank_results, aes(x = surv, y = rana)) + theme_bw() +
  geom_point() + 
  geom_smooth(method = 'lm') +
  ylab('proportion R. sphenocephala survival') + xlab('proprotion H. versicolor survival')
hyla_rana
m_hr <- lm(data = subset(tank_results, subset = (order != 'cont')),
               surv ~ rana)

###---CONTROLS FOR REL PLOT----------------------------------

## SEPARATING RESPONSE VARIABLES
trt_means_control <- subset(trt_means, subset = (order == 'cont'))

surv_cont <- ggplot(subset(trt_means_control, subset = (variable == 'surv')),
                    aes(x = order, y = mean, shape = sync, fill = sync, color = sync, group = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.8),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0.4, 0.9) +
  labs(color = 'hatching\nsynchrony', fill = 'hatching\nsynchrony', shape = 'hatching\nsynchrony') +
  xlab(NULL) + ylab('proportion survival')

biom_cont <- ggplot(subset(trt_means_control, subset = (variable == 'biom')),
                    aes(x = order, y = mean, shape = sync, fill = sync, color = sync, group = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.8),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(3000, 7000) + 
  labs(color = 'hatching\nsynchrony', fill = 'hatching\nsynchrony', shape = 'hatching\nsynchrony') +
  xlab(NULL) + ylab('total biomass export (mg)')

mass_cont <- ggplot(subset(trt_means_control, subset = (variable == 'mass')),
                    aes(x = order, y = mean, shape = sync, fill = sync, color = sync, group = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.8),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(color = 'hatching\nsynchrony', fill = 'hatching\nsynchrony', shape = 'hatching\nsynchrony') +
  xlab(NULL) + ylab('per capita mass (mg)')

emer_cont <- ggplot(subset(trt_means_control, subset = (variable == 'emer')),
                    aes(x = order, y = mean, shape = sync, fill = sync, color = sync, group = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.8),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(color = 'hatching\nsynchrony', fill = 'hatching\nsynchrony', shape = 'hatching\nsynchrony') +
  xlab(NULL) + ylab('mean days to emergence')

emsd_cont <- ggplot(subset(trt_means_control, subset = (variable == 'emsd')),
                    aes(x = order, y = mean, shape = sync, fill = sync, color = sync, group = sync)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.8),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  theme(axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(color = 'hatching\nsynchrony', fill = 'hatching\nsynchrony', shape = 'hatching\nsynchrony') +
  xlab(NULL) + ylab('s.d., days to emergence')

## ALL TOGETHER

trt_means1 <- trt_means
trt_means1$variable <- factor(trt_means$variable, levels = c('surv', 'biom', 'mass', 'emer', 'emsd'))
levels(trt_means1$variable) <- c('survival', 'biomass', 'per capita mass', 'mean emergence', 's.d. emergence')

control_means <- ggplot(subset(trt_means1, subset = (order == 'cont')),
                        aes(x = sync, y = mean, shape = sync, fill = sync, group = sync)) + mytheme +
  geom_line(size = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.2),
                aes(ymin = mean - se, ymax = mean + se, width = 0, color = sync)) +
  geom_point(size = 5, position = position_dodge(width = 0.2)) +
  xlab(NULL) + ylab('control means') +
  theme(legend.position = c(0.7, 0.25),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.background = element_rect(fill = 'white', color = 'white'),
        strip.text = element_text(size = 12, face = 'bold')) +  
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  #scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) + 
  scale_shape_manual(values = c(21, 24, 23)) +
  labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  facet_wrap(~ variable, scales = 'free', labeller = label_wrap_gen(width = 12, multi_line = T))
control_means

###---EMERGENCE DISTRIBUTIONS---------------------------------

## LAVA PLOT SPLIT BY TREATMENTS, COLOR = SYNC, WRAP = ORDER
high_sync <- subset(ind_results, subset = (sync == 'high'))
low_sync  <- subset(ind_results, subset = (sync == 'low'))

emer_lava_sync <- ggplot(high_sync, aes(date, group = as.factor(block))) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               fill = "#4b1d1d", colour = "#4b1d1d",
               geom = "ribbon", position = "identity") +
  facet_grid(order ~. , scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text  = element_text(size = 12),
        axis.line.x  = element_line(size = 1),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15),
        axis.title   = element_text(size = 18)) + 
  ylab("density of emerging individuals") + xlab("date") +
  stat_density(data = low_sync, size = 1.25, alpha = 0.75, adjust = 1,
               aes(date, color = as.factor(block), fill = as.factor(block), ymax = 0, ymin = -..density..),
               fill = "#fbb973", colour = "#fbb973",
               geom = "ribbon", position = "identity") +
  geom_hline(yintercept = 0, size = 1)
emer_lava_sync

## POOLING REPLICATES, WES ANDERSON COLOR SCHEME
emer_lava_sync1 <- ggplot(high_sync, aes(date)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               fill = "#4b1d1d", colour = "#4b1d1d",
               geom = "ribbon", position = "identity") +
  facet_grid(order ~. , scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text  = element_text(size = 12),
        axis.line.x  = element_line(size = 1),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15),
        axis.title   = element_text(size = 18)) + 
  ylab("density of emerging individuals") + xlab("date") +
  stat_density(data = low_sync, size = 1.25, alpha = 0.75, adjust = 1,
               aes(date, color = as.factor(block), fill = as.factor(block), ymax = 0, ymin = -..density..),
               fill = "#fbb973", colour = "#fbb973",
               geom = "ribbon", position = "identity") +
  geom_hline(yintercept = 0, size = 1)
emer_lava_sync1

## BOXPLOT
fig4 <- ggplot(ind_results, aes(x = order, y = date, color = sync, fill = sync)) +
  geom_boxplot(size = 1, alpha = 0.65, position = position_dodge(0.85)) +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  scale_fill_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  labs(fill = 'hatching\nsynchrony', color = 'hatching\nsynchrony') +
  ylab('date of emergence') + xlab('mean arrival (relative to Rana)') + 
  coord_flip()

###---TREATMENT GRID-----------------------------------------
###---EMPTY GRID FOR TRT PRESENTATION------------------------

## Empty plot for showing factorial treatment design
## Could consider putting in toy data, or just add distributions in powerpoint/inkscape

#means_compiled$order <- factor(means_compiled$order, levels = c('cont', 'early', 'same', 'late'))
View(trt_means)
trt_grid <- ggplot(trt_means, x = order, y = se) + theme_bw() +
  facet_grid(sync ~ order) +
  xlab("date") + ylab('number hatching') + 
  theme(strip.text = element_text(size = 16),
        axis.title = element_text(size = 18))
trt_grid

####---MS FIGURES---####

## SUPPLEMENT - TREATMENT MEANS, ABSOLUTE 
figS1 <- plot_grid(surv_means, biom_means, mass_means, emer_means, emsd_means, #rana_means,
          nrow = 2, labels = c("A", "B", "C", "D", "E"))#, "F"))
#tiff("figS1.tiff", height = 22, width = 28, units = 'cm', res = 1200)
#plot(figS1)
#dev.off()

## FIGURE 2 - CONTROL MEANS
fig2 <- plot_grid(surv_cont, biom_cont, mass_cont, emer_cont, emsd_cont, 
          nrow = 1, labels = c('A', 'B', 'C', 'D', 'E'),
          rel_widths = c(1, 1, 1, 1, 1.5))
fig2
#tiff("fig2.tiff", height = 7, width = 28, units = 'cm', res = 1200)
#plot(fig2)
#dev.off()

## FIGURE 3 - TREATMENT MEANS, RELATIVE TO CONTROLS
fig3 <- plot_grid(surv_means_rel1, biom_means_rel, mass_means_rel, emer_means_rel, emsd_means_rel,
          nrow = 2, labels = c('A', 'B', 'C', 'D', 'E'))
fig3
#tiff("fig3.tiff", height = 22, width = 28, units = 'cm', res = 1200)
#plot(fig3)
#dev.off()

## FIGURE 4 - 
fig4
#tiff("fig4.tiff", height = 13, width = 15, units = 'cm', res = 1200)
#plot(fig4)
#dev.off()

#############################################################
#############################################################
###-----------------------ANALYSIS------------------------###
#############################################################
#############################################################

###---HELPFUL RESOURCES--------------------------------------

## Modelling 101
# https://www.ssc.wisc.edu/sscc/pubs/MM/MM_Models.html
# https://arxiv.org/ftp/arxiv/papers/1308/1308.5499.pdf

## Mixed models in lme4
# http://lme4.r-forge.r-project.org/book/Ch1.pdf

## Basic guide to model formulation and selection
# https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html

## Plot diagnostics for glms
# https://stats.stackexchange.com/questions/121490/interpretation-of-plot-glm-model
# https://www.ssc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html

## Continuous dependent variable with ordinal independent variable
# https://stats.stackexchange.com/questions/33413/continuous-dependent-variable-with-ordinal-independent-variable

## Comparison of different glm/lm packages and functions
# http://glmm.wikidot.com/pkg-comparison

## Comparing model using numeric vs. unordered factor vs. ordered factor independent variable
# https://stackoverflow.com/questions/25735636/interpretation-of-ordered-and-non-ordered-factors-vs-numerical-predictors-in-m


#############################################################
###--------------------ABSOLUTE MODELS--------------------###
#############################################################

###---SURVIVAL (HYLA)---------------------------------------- 

### CHECKING ASSUMPTIONS ###

# is data normally distributed?
qqp(tank_results$surv, 'norm') # looks pretty normal, though it's customary to use binomial error structure

# are variances equal across treatments? yes
ggplot(tank_results, aes(x = order, y = surv)) +
  geom_boxplot() +
  facet_grid(~ sync)

### COARSE MODELS ###

# one explanatory factor at a time. bad fit
msurv0 <- lmer(data = tank_results, 
               surv ~ order + (1|block), REML = F)

# dropping controls, interaction, block as RE
msurv1 <- lmer(data = subset(tank_results, subset = (order != 'cont')), 
               surv ~ sync * order + (1|block), REML = F)

# dropping controls, no interaction, block as RE
msurv2 <- lmer(data = subset(tank_results, subset = (order != 'cont')),
               surv ~ sync + order + (1|block), REML = F)

# dropping controls, no interaction, block as FE
msurv3 <- lm(data = subset(tank_results, subset = (order != 'cont')), 
             surv ~ sync + order + block)

# with controls, interaction, block as RE
msurv4 <- lmer(data = tank_results, 
               surv ~ sync * order + (1|block), REML = F)

### MODEL SELECTION ###

# comparing m1, m2, m3. Do m4 separately because it uses a different dataset
anova(msurv1, msurv2, msurv3)   # msurv1 is preferred, lowest AIC
Anova(msurv1)

# comparing m1 to m4- almost identical.
summary(msurv4)  # adv of msurv4 is inclusion of controls
Anova(msurv4)    # adv of msurv1 is option for ordinal order and sync

### MODEL DIAGNOSTICS ###

# do this for all candidate models to compare fit
plot(msurv4)   # looks like maybe less variance on the right, but is this just because there are fewer observations there?
plot(ranef(msurv4))
coefficients(msurv4)
res_msurv4 <- residuals(msurv4)
plot(res_msurv4)
qqnorm(res_msurv4)
qqline(res_msurv4)
plot(predict(msurv4), tank_results$surv, xlim = c(0, 1), ylim = c(0, 1)) # predicted vs. actual
abline(a = 0, b = 1)  # actual vs. predicted looks good, plot by trt to check for systematic misses

### FINETUNING SELECTED MODEL ###
### LEAVING THIS IS TO SHOW FULL MODEL SELECTION PROCESS ###
### BUT CHOSE NOT TO USE ANY ORDINAL/NUMERIC FACTORS ###

# using numeric codes for synchrony and it totally changes the results. why so different?
#msurv4a <- lmer(data = survivors_nocontrol, prop_surv ~ sync_numeric * order_numeric + (1|block), na.action = na.omit,  REML = F)

# using factored numeric codes or sync & order -- same as msurv4, as expected
#msurv4b <- lmer(data = survivors_nocontrol, prop_surv ~ as.factor(sync_numeric) * as.factor(order_numeric) + (1|block), na.action = na.omit,  REML = F)

# sync as an ordered factor, keeping controls (so order stays a factor)
#msurv4c <- lmer(data = survivors, prop_surv ~ ordered(sync) * order + (1|block), na.action = na.omit, REML = F)

# both response variables as ordered factors - same Anova, but coefs are complicated. 
#msurv4d <- lmer(data = survivors_nocontrol, prop_surv ~ as.ordered(sync) * as.ordered(order) + (1|block), na.action = na.omit, REML = F)

### LAST MODEL TO COMPARE- BINOMIAL ERROR STRUCTURE ###

# basically the same fit to data, going for msurv5 because it's customary for survival data to be binomial
msurv5 <- glmer(data = subset(tank_results, subset = (order != 'cont')), 
                cbind(surv*45, 45 - surv*45) ~ (order) * (sync) + (1 | block), family = 'binomial')

###---BIOMASS EXPORT (HYLA)----------------------------------

### CHECKING ASSUMPTIONS ###

# is data normally distributed?
qqp(tank_results$biom, 'norm')
qqp(log(tank_results$biom), 'norm') # both look ok, one outlier on log. check fit of model both ways

# are variances equal across treatments? appears not across order
ggplot(tank_results, aes(x = order, y = biom)) +
  geom_boxplot() +
  facet_grid(~ sync)

### COARSE MODELS ###

# with controls, interaction, block as RE 
mbiom1 <- lmer(data = tank_results, 
               biom ~ sync * order + (1 | block), REML = F)

# without controls, interaction, block as RE - same formulation as survival
mbiom4 <- lmer(data = subset(tank_results, subset = (order != 'cont')), 
               biom ~ (sync) * (order)  + (1 | block), REML = F)
mbiom4_log <- lmer(data = subset(tank_results, subset = (order != 'cont')),
                   log(biom) ~ (sync) * (order)  + (1 | block), na.action = na.omit, REML = F)
anova(mbiom4, mbiom4_log)

### MODEL DIAGNOSTICS ###
## Compare mbiom4_log with mbiom1
plot(mbiom4_log)
res_mbiom4_log <- residuals(mbiom4_log)
plot(res_mbiom4_log)
qqnorm(res_mbiom4_log)
qqline(res_mbiom4_log)

###---PER CAPITA MASS (HYLA)---------------------------------

### CHECKING ASSUMPTIONS ###

# is data normally distributed? no, lognormal
ind_resultsxNA <- na.omit(ind_results)
qqp(ind_resultsxNA$mass_mg, 'norm')       # bad
qqp(log(ind_resultsxNA$mass_mg), 'norm')  # much better

# are variances equal across treatments? yes
ggplot(ind_resultsxNA, aes(x = order, y = mass_mg)) +
  geom_boxplot() +
  facet_grid(~ sync)

### COARSE MODELS ###

# interaction, tank as RE
# just taking out controls (m1 vs m5) has a huge impact on all 3 p-values-- shouldn't it mostly just affect synchrony?
mmass1 <- lmer(data = ind_resultsxNA, 
               log(mass_mg) ~ sync * order + (1 | tank), REML = F)

# interaction, block as FE, tank as RE
mmass2 <- lmer(data = ind_resultsxNA, 
               log(mass_mg) ~ sync * order + block + (1 | tank), REML = F)

# interaction, block and tank as independent REs
mmass3 <- lmer(data = ind_resultsxNA, 
               log(mass_mg) ~ sync * order + (1 | tank) + (1 | block), REML = F)

# interaction, no controls
mmass5 <- lmer(data = subset(ind_resultsxNA, subset = (order != 'cont')), 
               log(mass_mg) ~ sync * order + (1 | tank), REML = F)

### MODEL SELECTION ###

# comparing m1-3 
anova(mmass1, mmass2, mmass3)   # mmass3 is preferred, but not convinced these can be compared...

### MODEL DIAGNOSTICS ###

# do this for all candidate models to compare fit
plot(mmass5) 
res_mmass5 <- residuals(mmass5)
plot(res_mmass5)   # residuals have some pattern
qqnorm(res_mmass5)
qqline(res_mmass5) # amazing
# too many individuals to have a good fit-- model tank means instead

### MODELLING TANK MEANS INSTEAD OF INDIVIDUALS ###
qqp(tank_results$mass, 'norm')    # ok, but should remove outlier
tank_results27 <- subset(tank_results, subset = (tank != 27))
qqp(tank_results27$mass, 'norm') # much better

# does fewer survivors mean larger individuals?
plot(tank_results$surv, tank_results$mass) # no

# modelling tank mean masses
mmass6 <- lmer(data = subset(tank_results27, subset = (order != 'cont')), 
               mass ~ sync * order + (1 | block), REML = F)

# do this for all candidate models to compare fit
plot(mmass6)
res_mmass6 <- residuals(mmass6)
plot(res_mmass6)
qqnorm(res_mmass6)
qqline(res_mmass6)

###---EMERGENCE MEAN (HYLA)--------------------------------

### CHECKING ASSUMPTIONS ###

# is data normally distributed?
qqp(ind_results$timetoemerge, 'norm')      # definitely not normal 
qqp(log(ind_results$timetoemerge), 'norm') # or lognormal

# are variances equal across treatments? pretty close
ggplot(ind_results, aes(x = order, y = timetoemerge)) +
  geom_boxplot() +
  facet_grid(~ sync)

### COARSE MODELS ###

### Using individual time to emerge data as response ###
memer1 <- lmer(data = ind_results, 
               timetoemerge ~ sync * order + (1 | tank), REML = F)

### MODEL DIAGNOSTICS ###
# not a good fit, look to model tank means instead
plot(memer1)
res_memer1 <- residuals(memer1) # residuals have a definite pattern
plot(res_memer1)
qqnorm(res_memer1)
qqline(res_memer1)

### MODELLING TANK MEANS INSTEAD OF INDIVIDUALS ###

# is data normally distributed?
qqp(tank_results$emer, 'norm')

# are variances equal across treatments? not really
ggplot(tank_results, aes(x = order, y = emer)) +
  geom_boxplot() +
  facet_grid(~ sync)

# modelling tank emergence means
memer2 <- lmer(data = subset(tank_results, subset = (order != 'cont')), 
               emer ~ sync * order + (1 | block), REML = F)

# assessing fit
plot(memer2)
res_memer2 <- residuals(memer2)
plot(res_memer2)
qqnorm(res_memer2)
qqline(res_memer2)

###---EMERGENCE S.D. (HYLA)----------------------------------

### CHECKING ASSUMPTIONS ###

# is data normally distributed? meh
qqp(tank_results$emsd, 'norm')

# are variances equal across treatments? no
ggplot(tank_results, aes(x = order, y = emsd)) +
  geom_boxplot() +
  facet_wrap(~ sync)

### COARSE MODELS ###

memsd1 <- lmer(data = subset(tank_results, subset = (order != 'cont')),
               emsd ~ order * sync + (1 | block), REML = F)

### MODEL DIAGNOSTICS ###
plot(memer_sd1) 
plot(ranef(memer_sd1))
res_memer_sd1 <- residuals(memer_sd1)
plot(res_memer_sd1)  
qqnorm(res_memer_sd1)
qqline(res_memer_sd1) # doesn't look amazing
plot(predict(memer_sd1), survivors_nocontrol$emer_sd, xlim = c(10, 60), ylim = c(10, 60))
abline(a = 0, b = 1)  

## Accounting for unequal variances in both factors, then 1 at a time
library(nlme)
memsd2 <- lme(data = subset(tank_results, subset = (order != 'cont')),
                 emsd ~ order * sync, 
                 random = ~1 | block, 
                 weights = varComb(varIdent(form = ~1 | order), 
                                   varIdent(form = ~1 | sync)))

memsd3 <- lme(data = subset(tank_results, subset = (order != 'cont')),
                 emsd ~ order * sync, 
                 random = ~1 | block, 
                 weights = varComb(varIdent(form = ~1 | order)))

memsd4 <- lme(data = subset(tank_results, subset = (order != 'cont')),
                 emsd ~ order * sync, 
                 random = ~1 | block, 
                 weights = varComb(varIdent(form = ~1 | sync)))

# AICs very close but memer_sd4 preferred, only unequal var in sync
anova(memsd2, memsd3, memsd4) 
summary(memer_sd4)

## Don't know if this is correct way to pull out test statistics...
anova(memsd4)

## Checking fit
plot(memsd4)
res_memsd4 <- residuals(memsd4)
plot(res_memsd4)
qqnorm(res_memsd4)
qqline(res_memsd4)

###---RANA SURVIVAL------------------------------------------

### CHECKING ASSUMPTIONS ###

# is data normally distributed? no
qqp(tank_results$rana, 'norm') 

# are variances equal across treatments? no
ggplot(tank_results, aes(x = order, y = rana)) +
  geom_boxplot() +
  facet_grid(~ sync)

### COARSE MODELS ###

## Matching hyla survival formulation
mrana1 <- glmer(data = subset(tank_results, subset = (order != 'cont')), 
                cbind(rana * 30, 30 - rana * 30) ~ order * sync + (1 | block), family = 'binomial')
Anova(mrana1)

## Pretty good fit for Rana
plot(mrana1)
res_mrana1 <- residuals(mrana1)
plot(res_mrana1)
qqnorm(res_mrana1)
qqline(res_mrana1)

## Accounting for unequal variances in both factors, then 1 at a time-- cannot specify binomial family
## Not significantly changing things, so sticking with mrana1
library(nlme)
mrana2 <- lme(data = subset(tank_results, subset = (order != 'cont')),
              cbind(rana * 30, 30 - rana * 30) ~ order * sync, 
              random = ~1 | block, 
              weights = varComb(varIdent(form = ~1 | order), 
                                varIdent(form = ~1 | sync)))

mrana3 <- lme(data = subset(tank_results, subset = (order != 'cont')),
              cbind(rana * 30, 30 - rana * 30) ~ order * sync, 
              random = ~1 | block, 
              weights = varComb(varIdent(form = ~1 | order)))

mrana4 <- lme(data = subset(tank_results, subset = (order != 'cont')),
              cbind(rana * 30, 30 - rana * 30) ~ order * sync, 
              random = ~1 | block, 
              weights = varComb(varIdent(form = ~1 | sync)))

# AICs very close but mrana3 preferred, only unequal var in order
anova(mrana2, mrana3, mrana4) 
summary(mrana3)
Anova(mrana3)

## Checking fit
plot(mrana3)
res_mrana3 <- residuals(mrana3)
plot(res_mrana3)
qqnorm(res_mrana3)
qqline(res_mrana3)


#############################################################
###--------------------RELATIVE MODELS--------------------###
#############################################################

###---SURVIVAL (HYLA)----------------------------------------

# is data normally distributed? yes
qqp(tank_results_contadj$surv_adj, 'norm')

# model
msurv_rel <- lmer(data = tank_results_contadj, surv_adj ~ sync * order + (1 | block))
Anova(msurv_rel)

# diagnostics
summary(msurv_rel)
plot(msurv_rel)   
res_msurv_rel <- residuals(msurv_rel)
plot(res_msurv_rel)
qqnorm(res_msurv_rel)
qqline(res_msurv_rel)
plot(predict(msurv_rel), tank_results_contadj$surv_adj, xlim = c(-0.5, 0.1), ylim = c(-0.5, 0.1)) # predicted vs. actual
abline(a = 0, b = 1)

###---BIOMASS EXPORT (HYLA)----------------------------------

# is data normally distributed? yes
qqp(tank_results_contadj$biom_adj, 'norm')

# model
mbiom_rel <- lmer(data = tank_results_contadj, biom_adj ~ sync * order + (1 | block))
Anova(mbiom_rel)

# diagnostics
summary(mbiom_rel)
plot(mbiom_rel)   
res_mbiom_rel <- residuals(mbiom_rel)
plot(res_mbiom_rel)
qqnorm(res_mbiom_rel)
qqline(res_mbiom_rel)
plot(predict(mbiom_rel), tank_results_contadj$biom_adj, xlim = c(-4000, 2000), ylim = c(-4000, 2000)) # predicted vs. actual
abline(a = 0, b = 1)

###---PER CAPITA MASS (HYLA)---------------------------------

# is data normally distributed? w/o t27, yes
tank_results_contadj_27 <- subset(tank_results_contadj, subset = (tank !=27))
qqp(tank_results_contadj_27$mass_adj, 'norm') # looks pretty normal, but could remove T27

# model
mmass_rel <- lmer(data = tank_results_contadj_27, mass_adj ~ sync * order + (1 | block))
Anova(mmass_rel)

# diagnostics
summary(mmass_rel)
plot(mmass_rel)   
res_mmass_rel <- residuals(mmass_rel)
plot(res_mmass_rel)
qqnorm(res_mmass_rel)
qqline(res_mmass_rel)
plot(predict(mmass_rel), surv_contadj_27$mass_adj, xlim = c(-60, 60), ylim = c(-60, 60)) # predicted vs. actual
abline(a = 0, b = 1)

###---EMERGENCE MEAN (HYLA)----------------------------------

# is data normally distributed? yes
qqp(tank_results_contadj$emer_adj, 'norm')

# model
memer_rel <- lmer(data = tank_results_contadj, emer_adj ~ sync * order + (1 | block))
Anova(memer_rel)

# diagnostics
summary(memer_rel)
plot(memer_rel)   
res_memer_rel <- residuals(memer_rel)
plot(res_memer_rel)
qqnorm(res_memer_rel)
qqline(res_memer_rel)
plot(predict(memer_rel), tank_results_contadj$emer_adj, xlim = c(-20, 80), ylim = c(-20, 80)) # predicted vs. actual
abline(a = 0, b = 1)

###---EMERGENCE S.D. (HLYA)----------------------------------

# is data normally distributed? close, some deviation at tail
qqp(tank_results_contadj$emsd_adj, 'norm')

# variances not equal across treatments? no
ggplot(tank_results_contadj, aes(x = order, y = emsd_adj)) +
  geom_boxplot() +
  facet_wrap(~ sync)

# model - first, see model fit without accounting for unequal variances
memsd_rel <- lmer(data = tank_results_contadj, emsd_adj ~ sync * order + (1 | block))
Anova(memsd_rel)

# diagnostics - not great, so let's account for unequal variances
summary(memsd_rel)
plot(memsd_rel)   
res_memsd_rel <- residuals(memsd_rel)
plot(res_memsd_rel)
qqnorm(res_memsd_rel)
qqline(res_memsd_rel)
plot(predict(memsd_rel), tank_results_contadj$emsd_adj, xlim = c(0, 40), ylim = c(0, 40)) # predicted vs. actual
abline(a = 0, b = 1)

### Now, accounting for unequal variances ###
library(nlme) # can't be done in lme4 package

# first in both factors
memsd_rel2 <- lme(data = tank_results_contadj,
                    emsd_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | order), 
                                      varIdent(form = ~ 1 | sync)))

# just order
memsd_rel3 <- lme(data = tank_results_contadj,
                    emsd_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | order)))

# just synchrony
memsd_rel4 <- lme(data = tank_results_contadj,
                    emsd_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | sync)))

# AICs very close but memercv_rel4 preferred, only unequal var in sync
anova(memsd_rel2, memsd_rel3, memsd_rel4) 
Anova(memsd_rel4) # need to confirm this is the best way to pull out test statistics

# diagnostics
plot(memsd_rel4)
res_memsd_rel4 <- residuals(memsd_rel4)
plot(res_memsd_rel4)
qqnorm(res_memsd_rel4)
qqline(res_memsd_rel4)
plot(predict(memsd_rel4), tank_results_contadj$emsd_adj, xlim = c(0, 40), ylim = c(0, 40))
abline(a = 0, b = 1) # not horrible


#############################################################
###------------------PREDICTxACTUAL PLOTS-----------------###
#############################################################

## These need to be reformulated with new data structures, but useful for assessing fit of models

###---ALL MODELS---------------------------------------------
tank_results_nocont <- subset(tank_results, subset = (order != 'cont'))
msurv_rel <- lmer(data = tank_results_contadj, surv_adj ~ sync * order + (1 | block))
mbiom_rel <- lmer(data = tank_results_contadj, biom_adj ~ sync * order + (1 | block))
tank_results_contadj_27 <- subset(tank_results_contadj, subset = (tank !=27))
mmass_rel <- lmer(data = tank_results_contadj_27, mass_adj ~ sync * order + (1 | block))
memer_rel <- lmer(data = tank_results_contadj, emer_adj ~ sync * order + (1 | block))
memsd_rel4 <- lme(data = tank_results_contadj, emsd_adj ~ order * sync, random = ~ 1 | block, weights = varComb(varIdent(form = ~ 1 | sync)))
mrana1 <- glmer(data = subset(tank_results, subset = (order != 'cont')), cbind(rana * 30, 30 - rana * 30) ~ order * sync + (1 | block), family = 'binomial')

###---SURVIVAL (HYLA)----------------------------------------

# chose to plot with raw data, and also don't think this is the right 
# way to plot model outputs, but keeping in case there's something useful
### ADD PREDICTED VALUE TO DF FOR EACH TANK ###
tank_results_nocont$predict_msurv <- predict(msurv_rel, type = 'response')

ggplot(tank_results_nocont, aes(x = surv, y = predict_msurv)) +
  geom_point()

### TRYING TO EXTRACT COEFFICIENTS, BUT NO LUCK ###
#surv_coef <- as.data.frame(summary(msurv)$coef)
#colnames(surv_coef) <- c('mean', 'se', 'z', 't')
#surv_coef$mean_corr <- c(0.559, 0.752, 0.593, 0.493, 0.411, 0.578, 0.333, 0.366, 0.426)
#surv_coef$sync  <- factor(surv_coef$sync, levels = c("low", "med", "high"))
#surv_coef$order <- factor(surv_coef$order, levels = c('early', 'same', 'late', 'cont'))

### NOW SUMMARIZE BY TREATMENT ###
predict_msurv_summary <- summarySE(data = survivors_nocontrol, measurevar = "predict_msurv", groupvars = c("order", "sync"))
predict_msurv_summary$sync  <- factor(predict_msurv_summary$sync, levels = c("low", "med", "high"))
predict_msurv_summary$order <- factor(predict_msurv_summary$order, levels = c('early', 'same', 'late', 'cont'))

### CHECK PREDICTED VALUES AGAINST DATA ###
surv_predictxactual <- ggplot(survivors_nocontrol, 
                              aes(x = order, y = prop_surv)) + #mytheme +
  geom_point(size = 5, color = 'black', alpha = 0.5) +
  geom_point(aes(x = order, y = predict_msurv),
             size = 5, color = 'springgreen4', alpha = 0.5) +
  ylab("number of survivors") + xlab('mean arrival (relative to Rana)') +
  facet_grid(. ~ sync)
surv_predictxactual

### PLOT PREDICTED MEANS (PAPER FIGURE) ###
# with means from raw data and se from model 
surv_predict_means <- ggplot(surv_coef,
                             aes(x = order, y = mean_corr, group = sync, fill = sync, shape = sync)) + mytheme +
  geom_hline(yintercept = 0.58, size = 1, linetype = 'dashed', color = "#bdc9e1") +
  geom_hline(yintercept = 0.69, size = 1, linetype = 'dashed', color = "#67a9cf") +
  geom_hline(yintercept = 0.62, size = 1, linetype = 'dashed', color = "#016c59") +
  geom_point(data = subset(survivors, subset = (order != 'cont')), aes(x = order, y = prop_surv, group = sync, fill = sync, shape = sync),
             position = position_dodge(width = 0.1)) + 
  geom_line(size = 1, position = position_dodge(width = 0.1)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(ymin = mean_corr - se, ymax = mean_corr + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  ylab("proportion survival") + xlab(NULL) + labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = c(0.1, 0.22),
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23))  
surv_predict_means 

# with means and se calculated from model predicted values
# error bars unreasonably small
surv_predict_means2 <- ggplot(predict_msurv_summary,
                              aes(x = order, y = predict_msurv, group = sync, fill = sync, shape = sync)) + mytheme +
  geom_hline(yintercept = 0.58, size = 1, linetype = 'dashed', color = "#bdc9e1") +
  geom_hline(yintercept = 0.69, size = 1, linetype = 'dashed', color = "#67a9cf") +
  geom_hline(yintercept = 0.62, size = 1, linetype = 'dashed', color = "#016c59") +
  geom_point(data = survivors_nocontrol, aes(x = order, y = prop_surv, group = sync, fill = sync, shape = sync),
             position = position_dodge(width = 0.1)) + 
  geom_line(size = 1, position = position_dodge(width = 0.1)) +
  geom_errorbar(size = 1, position = position_dodge(width = 0.1),
                aes(ymin = predict_msurv - se, ymax = predict_msurv + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.1)) +
  ylab("proportion survival") + xlab(NULL) + labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = c(0.1, 0.22),
        axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10)) +  
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23))  

###---BIOMASS (HYLA)----------------------------------------

### MAKE A DATAFRAME OF PREDICTION FOR EACH TANK ###
survivors_nocontrol$predict_mbiom <- predict(mbiom_rel)

### CHECK PREDICTED VALUES AGAINST DATA ###
biom_predictxactual <- ggplot(survivors_nocontrol, 
                              aes(x = order, y = biomass)) + 
  geom_point(size = 5, color = 'black', alpha = 0.5) +
  geom_point(aes(x = order, y = predict_mbiom),
             size = 5, color = 'springgreen4', alpha = 0.5) +
  ylab("total biomass export (mg)") + xlab('mean arrival (relative to Rana)') +
  facet_grid(. ~ sync)
biom_predictxactual
###---PER CAPITA MASS (HYLA)---------------------------------

### ADD PREDICTED MASS FOR EACH INDIVIDUAL TO DF ###
metamorphsxNA_nocont$log_predict_mass <- predict(mmass)
### NOW SUMMARIZE BY TREATMENT ###
predict_mmass_summary <- summarySE(data = metamorphsxNA_nocont, measurevar = "log_predict_mass", groupvars = c("order", "sync"))
predict_mmass_summary$sync  <- factor(predict_mmass_summary$sync, levels = c("low", "med", "high"))
predict_mmass_summary$order <- factor(predict_mmass_summary$order, levels = c('early', 'same', 'late', 'cont'))

### MODEL WITH TANK MEANS INSTEAD OF INDIVIDUALS
survivors_nocontrol_27 <- subset (survivors_nocontrol, subset = (tank != 27))
mmass6 <- lmer(data = survivors_nocontrol_27, mass_mean ~ sync * order + (1 | block), na.action = na.omit, REML = F)
survivors_nocontrol_27$predict_mmass <- predict(mmass6)
### NOW SUMMARIZE BY TREATMENT ###
predict_mmass6_summary <- summarySE(data = metamorphsxNA_nocont, measurevar = "log_predict_mass", groupvars = c("order", "sync"))
predict_mmass6_summary$sync  <- factor(predict_mmass_summary$sync, levels = c("low", "med", "high"))
predict_mmass6_summary$order <- factor(predict_mmass_summary$order, levels = c('early', 'same', 'late', 'cont'))

### CHECK PREDICTED VALUES AGAINST DATA ###
mass_predictxactual <- ggplot(survivors_nocontrol_27, 
                              aes(x = order, y = mass_mean)) +
  geom_point(size = 2, color = 'black') +
  geom_boxplot(aes(x = order, y = predict_mmass),
               size = 2, color = 'springgreen4') +
  ylab("per capita mass (mg)") + xlab('mean arrival (relative to Rana)') +
  facet_grid(. ~ sync)
mass_predictxactual

###---EMERGENCE MEAN (HYLA)----------------------------------

survivors_nocontrol$predict_memer <- predict(memer)

### CHECK PREDICTED VALUES AGAINST DATA ###
emer_predictxactual <- ggplot(survivors_nocontrol, 
                              aes(x = order, y = emer_mean)) + #mytheme +
  #geom_point(size = 5, color = 'black', alpha = 0.5) +
  geom_boxplot(size = 2, color = 'black') +
  #geom_point(data = predict_memer, aes(x = order, y = predicted_mean),
  #           size = 5, color = 'springgreen4', alpha = 0.5) +
  geom_boxplot(aes(x = order, y = predict_memer),
               size = 2, color = 'springgreen4', alpha = 0.7) +
  ylab("cumulative metamorphs") + xlab('mean arrival (relative to Rana)') +
  facet_grid(. ~ sync)
emer_predictxactual # sexy af




#############################################################
#############################################################
###----------------------EXTRA BITS-----------------------###
#############################################################
#############################################################

###---DIGGING INTO RANA SURVIVAL-----------------------------
## Looking at Rana sizes
rana_size  <- read.csv("raw_ranasize.csv", header = T)
treatments <- read.csv("raw_tanks.csv", header = T)
rana_size <- left_join(rana_size, treatments, "tank")
rana_svl <- subset(rana_size, subset = (measurement == 'svl'))
rana_svl$sync  <- factor(rana_svl$sync, levels = c("low", "med", "high"))
rana_svl$order <- factor(rana_svl$order, levels = c('early', 'same', 'late', 'cont'))

# if rana data is reliable, we expect 'late' trts to be bigger bc introduced earlier
# otherwise, might be indication that some rana died at metamorphosis and our survival data is bad
# late are a little larger, but not totally convincing and not ordered
ggplot(rana_svl, aes(x = order, y = value)) + theme_bw() +
  geom_boxplot(size = 2) 

## Was there a bias in how many HV remained across treatments? no
teardown <- read.csv('raw_teardown.csv', header = T)
teardown_hyla <- ggplot(teardown, aes(x = order, y = hyla_tads)) + theme_bw() +
  geom_boxplot(size = 2) +
  facet_wrap(~ sync) +
  ylab("number of hyla at end of exp") + xlab("mean") 

m_hylatads <- lmer(data = teardown,
                    hyla_tads ~ order * sync + (1 | block))
Anova(m_hylatads)

# no correlation between the number of rana and hyla remaining
ggplot(teardown, aes(x = hyla_tads, y = rana_tads)) + 
  geom_point()
###---PLOTTING TADPOLE SIZES ACROSS INTRODUCTIONS------------

## Measurements of a subset of organisms added to mesocosms at each introduction
introsizes <- read.csv("raw_introsizes.csv", header = T)

## Were sizes standardized across introductions?
introsize_plot <- ggplot(introsizes, aes(x = as.factor(intro), y = hw, color = sp)) + mytheme +
  geom_boxplot(size = 2) +
  #geom_point(alpha = 0.5) +
  ylim(0, 4)+
  ylab("headwidth (mm)") + xlab("introduction number")
introsize_plot

## Verifying measurements-- are length and hw fairly correlated?
l_hw <- ggplot(introsizes, aes(x = hw, y = l)) + mytheme +
  geom_point() +
  facet_wrap(~sp, scale = 'free') +
  geom_smooth(method = 'lm') +
  ylab("length (mm)") + xlab("headwidth (mm)")
l_hw

###---TEST SIZE DIFFERENCES ACROSS INTROS WITH LM-------------

## HYLA
hyla_m <- lm(hw ~ intro, data = subset(introsizes, subset = (sp == 'hyla')))
summary(hyla_m)

hyla_m2 <- lm(l ~ intro, data = subset(introsizes, subset = (sp == 'hyla')))
summary(hyla_m2)

# RANA
rana_m <- lm(hw ~ intro, data = subset(introsizes, subset = (sp == 'rana')))
summary(rana_m)

rana_m2 <- lm(l ~ intro, data = subset(introsizes, subset = (sp == 'rana')))
summary(rana_m2)
###---VOLKER'S METHOD FOR MEANS/SES--------------------------

### VOLKER'S METHOD FOR GETTING MEANS/SES ###
# make a treatment combination matrix
trt <- expand.grid(sync = levels(survivors_nocontrol$sync), order = levels(survivors_nocontrol$order))

#then use your model (GLMMODEL) to get new output
# can't get this line to work. how to group predictions by trt?
newdatatrt <- cbind(trt, predict(msurv)) #, level = 0, se.fit = T)) # get predicted values

newdatatrt$prop_surv <- (newdatatrt$fit)
newdataS1$LCI <- (newdataS1$fit - (1.95 * newdataS1$se.fit)) # lower confidence interval
newdataS1$HCI <- (newdataS1$fit + (1.95 * newdataS1$se.fit)) # upper confidence interval
ggplot(survivors_nocontrol, aes(x = order, y = prop_surv, color = sync)) + 
  geom_point() + 
  stat_summary(fun.data = mean_se) + 
  geom_ribbon(data = newdatatrt, aes(ymin = LCI, ymax = HCI, fill = sync), alpha = 0.2) +
  geom_line(data = newdatatrt, aes(x = order, y = prop_surv, color = sync), size = 2) + 
  theme_bw()

###---MODELLING EMER C.V. (ABANDONED FOR S.D.)---------------

qqp(surv_contadj$emercv_adj, 'norm') # looks pretty normal

# variances not equal across treatments
ggplot(surv_contadj, aes(x = order, y = emercv_adj)) +
  geom_boxplot() +
  facet_wrap(~ sync)

# first, see model fit without accounting for unequal variances
memercv_rel <- lmer(data = surv_contadj, emercv_adj ~ sync * order + (1 | block))
Anova(memercv_rel)

summary(memercv_rel)
plot(memercv_rel)   
res_memercv_rel <- residuals(memercv_rel)
plot(res_memercv_rel)
qqnorm(res_memercv_rel)
qqline(res_memercv_rel)
plot(predict(memercv_rel), surv_contadj$emercv_adj, xlim = c(0, 0.4), ylim = c(0, 0.4)) # predicted vs. actual
abline(a = 0, b = 1)

## Now, accounting for unequal variances
library(nlme) # can't be done in lme4 package

# First in both factors
memercv_rel2 <- lme(data = surv_contadj,
                    emercv_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | order), 
                                      varIdent(form = ~ 1 | sync)))

# Just order
memercv_rel3 <- lme(data = surv_contadj,
                    emercv_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | order)))

# Just synchrony
memercv_rel4 <- lme(data = surv_contadj,
                    emercv_adj ~ order * sync, 
                    random = ~ 1 | block, 
                    weights = varComb(varIdent(form = ~ 1 | sync)))

# AICs very close but memercv_rel4 preferred, only unequal var in sync
anova(memercv_rel2, memercv_rel3, memercv_rel4) 
summary(memercv_rel4)

## Pull out test statistics
anova(memercv_rel4)

## Checking fit
plot(memercv_rel4)
plot(ranef(memercv_rel4))
res_memercv_rel4 <- residuals(memercv_rel4)
plot(res_memercv_rel4)
qqnorm(res_memercv_rel4)
qqline(res_memercv_rel4)
plot(predict(memercv_rel4), surv_contadj$emercv_adj, xlim = c(-0.2, 0.4), ylim = c(-0.2, 0.4))
abline(a = 0, b = 1) # not horrible



#############################################################
###--------------------PROGRESS PLOTS---------------------###
#############################################################

## Some of these need to be adjusted to new df names
###---PRELIMINARY METAMORPH ANALYSIS: BLOCKS & TOTALS--------

## Cumulative metamorphs collected, across treatments
metamorph_accumulation <- ggplot(ind_results, aes(x = date)) + mytheme +
  stat_bin(aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  #geom_smooth(aes(y = ..y..)) +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20)) +
  ylab("cumulative metamorphs collected") 

## Cumulative metamorphs collected by block-- big block effect
## Most of block effect caused in first couple weeks-- first low synchrony tads didn't take off?
metamorph_accumulation_block <- ggplot(ind_results, aes(x = date, color = as.factor(block))) + mytheme +
  stat_bin(data = subset(ind_results, block == 1), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  stat_bin(data = subset(ind_results, block == 2), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  stat_bin(data = subset(ind_results, block == 3), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  stat_bin(data = subset(ind_results, block == 4), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  stat_bin(data = subset(ind_results, block == 5), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  stat_bin(data = subset(ind_results, block == 6), aes(y = cumsum(..count..)), geom = "step", size = 2, binwidth = 1) +
  ylab("cumulative count")

## Accumulation by tank -- need to make new df for this to work. not necessary, except for appendix maybe
metamorph_accumulation_tank <- ggplot(meta_accum_alldates, aes(x = as.Date(date), y = cumulative_metas, color = as.factor(tank))) + mytheme +
  geom_point(alpha = 0.25) + 
  geom_smooth(size = 2, se = F) +
  theme(legend.position = 'none') +
  ylab("cumulative metamorphs collected") + xlab("date")

## Total number of metamorphs so far coming from each treatment
e <- ggplot(ind_results, aes(x = sync)) + 
  geom_histogram(stat = 'count') + theme_bw() +
  facet_grid(. ~ order) +
  xlab("synchrony") + ylab("number of metamorphs")

## Showing variation in metamorph collection within treatment across replicated
f <- ggplot(ind_results, aes(x = as.factor(block))) + 
  geom_histogram(stat = 'count') + #mytheme +
  facet_grid(sync ~ order) +
  xlab("block") + ylab("number of metamorphs")

###---MAKE METAMORPH ACCUMULATION DATASET, ALL DATES---------

### NOT WORKING NOW, BUT NECESSARY FOR FOLLOWING PLOTS ###
## Drop some columns so it's easier to work with
meta_accum_cut <- select(meta_accum, date, tank, ones, mass_mg)

## Make a dataframe with all date and tank combos
date_tank <- expand.grid(date = seq(as.Date("2018-05-07"), as.Date("2018-09-14"), by = "days"), tank = c(1:60))

### THIS IS NOT RIGHT ###
## Merge metamorph accumulation with expanded date/tank df
meta_accum_alldates <- merge(meta_accum_cut, date_tank, by.x = c('date', 'tank'), all.x = T, all.y = T)

## Missing values mean 0 metamorphs collected
meta_accum_alldates$ones[is.na(meta_accum_alldates$ones)] <- 0
meta_accum_alldates$mass_mg[is.na(meta_accum_alldates$mass_mg)] <- 0

## Sum multiple entries in a given date and tank for total found that day
meta_accum_alldates <- meta_accum_alldates %>%
  group_by(date, tank) %>%
  mutate(daily_metas = sum(ones),
         daily_biom  = sum(mass_mg))

## Add a cumulative column
meta_accum_alldates <- meta_accum_alldates %>%
  group_by(tank) %>%
  mutate(cumulative_metas = cumsum(daily_metas),
         cumulative_biom  = cumsum(daily_biom))

## Add back in treatment info for each tank
treatments <- read.csv("tanks.csv", header = T)
meta_accum_alldates <- meta_accum_alldates %>%
  left_join(treatments, by = 'tank')

## One line per tank per day
write.csv(meta_accum_alldates, "metamorph_accumulation_alldates.csv")

###---METAMORPH ACCUMULATION CURVES--------------------------

## Accumulation of metamorphs over time, tanks separated
## Can also make this with replicates pooled, using accum_trt in "data processing"
metamorph_accumulation_trt <- ggplot(meta_accum_alldates, aes(x = as.Date(date), y = cumulative_metas/45, color = as.factor(block))) + mytheme +
  geom_point(alpha = 0.15) + 
  geom_smooth(size = 2, se = F) +
  #geom_smooth(size = 2, se = F,
  #            method = 'glm', method.args = list(family = 'binomial')) + 
  facet_grid(sync ~ order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        legend.position = 'none') +
  ylab("proportion emergence") + xlab("date")

## Highlighting tanks with Pantala invasion
meta_accum_pantala <- subset(meta_accum_alldates, subset = (tank == 9:10 | tank == 30 | tank ==58))
metamorph_accumulation_trt_pantala <- ggplot(meta_accum_alldates, aes(x = as.Date(date), y = cumulative_metas/45, color = as.factor(block))) + mytheme +
  geom_point(alpha = 0.15) + 
  geom_smooth(size = 2, se = F, alpha = 0.5) +
  #geom_smooth(size = 2, se = F,
  #            method = 'glm', method.args = list(family = 'binomial')) + 
  geom_point(data = meta_accum_pantala, aes(x = as.Date(date), y = cumulative_metas/45),
             alpha = 0.5, color = 'red') + 
  geom_smooth(data = meta_accum_pantala, aes(x = as.Date(date), y = cumulative_metas/45),
              color = 'red', size = 2, se = F) +
  facet_grid(sync ~ order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        legend.position = 'none') +
  ylab("proportion emergence") + xlab("date")

## Biomass export over time, tanks separated
biomass_accumulation_trt <- ggplot(meta_accum_alldates, aes(x = as.Date(date), y = cumulative_biom, color = as.factor(block))) + mytheme +
  geom_point(alpha = 0.5) +
  geom_smooth(size = 2, se = F) +
  #method = 'lm', formula = y ~ poly(x, 2)) + 
  facet_grid(sync ~ order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        legend.position = 'none') +
  ylab("cumulative biomass export") + xlab("date")


###---TREATMENT EFFECTS: SURVIVAL--------------------

## BOXPLOT
surv_box <- ggplot(tank_results, aes(x = order, y = surv, color = sync)) + mytheme +
  geom_boxplot(size = 2) +
  ylab("proportion survivors") + xlab('mean arrival (relative to Rana)')

## MEAN +/- 95% CI
surv_mean <- ggplot(subset(trt_means, subset = (variable == 'surv')), 
                    aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  ylab("number of survivors") + xlab('mean arrival (relative to Rana)')

## MEAN WITH CONFIDENCE BANDS, CONTINUOUS MEAN (NO CONTROLS)
surv_abs_nocont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'surv')), 
                          aes(x = order_numeric, y = mean, color = sync)) + mytheme +
  geom_point(size = 4) + 
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = mean - se, ymax = mean + se, fill = sync, color = NULL)) +
  ylab(expression(paste("number of survivors"))) +
  xlab('mean arrival (relative to Rana)') 

## MEAN +/- 95% CI, RELATIVE TO CONTROLS
surv_cont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'surv')), 
                    aes(x = order_numeric, y = adj, color = sync)) + mytheme +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  geom_point(size = 4) + #, position = position_dodge(width = 0.5)) +
  #geom_errorbar(size = 3, position = position_dodge(width = 0.5),
  #              aes(ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_line(size = 2) + #, position = position_dodge(width = 0.5)) +
  geom_ribbon(aes(ymin=adj - se, ymax=adj + se, fill = sync, color = NULL), alpha=0.4) +
  ylab(expression(paste("number of survivors (", Delta, " control)"))) +
  xlab('mean arrival (relative to Rana)') 

## VIOLIN PLOT, ABSOLUTE
surv_violin <- ggplot(survivors, aes(x = order, y = no_metas, fill = sync, color = sync)) + mytheme +
  geom_violin(size = 1.5, draw_quantiles = 0.5, position = position_dodge(width = 0.75), alpha = 0.5) + 
  geom_point(position = position_dodge(width = 0.75)) +
  ylab('number of survivors') + xlab('mean arrival (relative to Rana)')

###---TREATMENT EFFECTS: BIOMASS-----------------------------

## BOXPLOT
biom_box <- ggplot(survivors, aes(x = order, color = sync, y = biomass)) + mytheme + 
  geom_boxplot(size = 2) +
  #geom_point(position = position_dodge(width = 0.5)) +
  ylab('biomass export (mg)') + xlab('mean arrival (relative to Rana)')

## MEAN +/- 95% CI
biom_mean <- ggplot(subset(means_compiled, subset = (variable == 'biom')), 
                    aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  ylab("total biomass export (mg)") + xlab('mean arrival (relative to Rana)')

## MEAN WITH CONFIDENCE BANDS, CONTINUOUS MEAN (NO CONTROLS)
biom_abs_nocont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'biom')), 
                          aes(x = order_numeric, y = mean, color = sync)) + mytheme +
  geom_point(size = 4) + 
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = mean - se, ymax = mean + se, fill = sync, color = NULL)) +
  ylab("total biomass export (mg)") +
  xlab('mean arrival (relative to Rana)') 

## MEAN +/- 95% CI, RELATIVE TO CONTROLS
biom_cont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'biom')), 
                    aes(x = order_numeric, y = adj, color = sync)) + mytheme +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  geom_point(size = 4) + #, position = position_dodge(width = 0.5)) +
  #geom_errorbar(size = 3, position = position_dodge(width = 0.5),
  #              aes(ymin = adj - ci, ymax = adj + ci, width = 0)) +
  geom_line(size = 2) + #, position = position_dodge(width = 0.5)) +
  geom_ribbon(aes(ymin = adj - se, ymax = adj + se, fill = sync, color = NULL), alpha=0.4) +
  ylab(expression(paste("total biomass export (", Delta, " control)"))) +
  xlab('mean arrival (relative to Rana)') 

## VIOLIN PLOT, ABSOLUTE
biom_violin <- ggplot(survivors, aes(x = order, y = biomass, fill = sync, color = sync)) + mytheme +
  geom_violin(size = 1.5, draw_quantiles = 0.5, position = position_dodge(width = 0.75), alpha = 0.5) + 
  geom_point(position = position_dodge(width = 0.75)) +
  ylab('total biomass export (mg)') + xlab('mean arrival (relative to Rana)')

###---TREATMENT EFFECTS: PER CAPITA MASS---------------------------

## BOXPLOT
mass_box <- ggplot(metamorphs, aes(x = order, y = mass_mg, color = sync)) + mytheme +
  geom_boxplot(size = 2) +
  ylab('per capita mass (mg)') + xlab('mean arrival (relative to Rana)')

## MEAN +/- 95% CI
mass_mean <- ggplot(subset(means_compiled, subset = (variable == 'mass')), 
                    aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_line(size = 2) +
  #geom_point(data = metamorphs, aes(x = order, y = mass_mg, color = sync), alpha = 0.25, position = 'jitter') +
  ylab("per capita mass") + xlab('mean arrival (relative to Rana)')

## MEAN WITH CONFIDENCE BANDS, CONTINUOUS MEAN (NO CONTROLS)
mass_abs_nocont <- ggplot(subset(means_compiled, subset = (variable == 'mass')), 
                          aes(x = order_numeric, y = mean, color = sync)) + mytheme +
  geom_point(size = 4) +
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = mean - se, ymax = mean + se, fill = sync, color = NULL)) +
  ylab("per capita mass (mg)") +
  xlab('mean arrival (relative to Rana)') 

## MEAN +/- 95% CI, RELATIVE TO CONTROLS
mass_cont <- ggplot(subset(means_compiled, subset = (variable == 'mass')), 
                    aes(x = order_numeric, y = adj, color = sync)) + mytheme +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  geom_point(size = 4) +
  #geom_errorbar(size = 2, position = position_dodge(width = 0.5),
  #              aes(ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = adj - se, ymax = adj + se, fill = sync, color = NULL), alpha=0.4) +
  ylab(expression(paste("per capita mass (", Delta, " control)"))) +
  xlab('mean arrival (relative to Rana)') 

## VIOLIN PLOT, ABSOLUTE
## is this right? treatment differences are really obscured here
mass_violin <- ggplot(metamorphs, aes(x = order, y = mass_mg, fill = sync, color = sync)) + mytheme +
  geom_violin(size = 1.5, draw_quantiles = 0.5, position = position_dodge(width = 0.75), alpha = 0.5) + 
  geom_point(position = position_dodge(width = 0.75)) +
  ylab('per capita mass (mg)') + xlab('mean arrival (relative to Rana)')

## LAVA PLOT
mass_lava <- ggplot(subset(metamorphs, subset = (order != 'cont')),
                    aes(mass_mg)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1,
               aes(ymax = ..density..,  ymin = -..density..),
               fill = "grey37", colour = "grey27",
               geom = "ribbon", position = "identity") +
  facet_grid(sync ~ order)

###---TREATMENT EFFECTS: EMERGENCE TIME--------------------------

## BOXPLOT
emer_box <- ggplot(metamorphs, aes(x = order, y = timetoemerge, color = sync)) + mytheme +
  geom_boxplot(size = 2) +
  ylab('days to emergence') + xlab('mean arrival (relative to Rana)')

## MEAN +/- 95% CI
emer_mean <- ggplot(subset(means_compiled, subset = (variable == 'emer')), 
                    aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  ylab("days to emergence") + xlab('mean arrival (relative to Rana)')

## MEAN WITH CONFIDENCE BANDS, CONTINUOUS MEAN (NO CONTROLS)
emer_abs_nocont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'emer')), 
                          aes(x = order_numeric, y = mean, color = sync)) + mytheme +
  geom_point(size = 4) + 
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = mean - se, ymax = mean + se, fill = sync, color = NULL)) +
  ylab("time to emergence (days)") +
  xlab('mean arrival (relative to Rana)') 

## MEAN +/- 95% CI, RELATIVE TO CONTROLS
emer_cont <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == 'emer')), 
                    aes(x = order_numeric, y = adj, color = sync)) + mytheme +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  geom_point(size = 4) + #, position = position_dodge(width = 0.5)) +
  #geom_errorbar(size = 3, position = position_dodge(width = 0.5),
  #              aes(ymin = adj - se, ymax = adj + se, width = 0)) +
  geom_line(size = 2) + #, position = position_dodge(width = 0.5)) +
  geom_ribbon(aes(ymin = adj - se, ymax = adj + se, fill = sync, color = NULL), alpha=0.4) +
  ylab(expression(paste("time to emergence (", Delta, " control)"))) +
  xlab('mean arrival (relative to Rana)') 

## VIOLIN PLOT, ABSOLUTE
emer_violin <- ggplot(metamorphs, aes(x = order, y = timetoemerge, fill = sync, color = sync)) + mytheme +
  geom_violin(size = 1.5, draw_quantiles = 0.5, position = position_dodge(width = 0.75), alpha = 0.5) +
  ylab("time to emergence (days)") + xlab('mean arrival (relative to Rana)')

## LAVA PLOT
emer_lava <- ggplot(metamorphs, aes(date, color = sync, fill = sync)) + mytheme +
  stat_density(size = 1.25, alpha = 0.25, adjust = 1,
               aes(ymax = ..density..,  ymin = -..density..),
               geom = "ribbon", position = 'identity') +
  facet_grid(. ~ order)

## LAVA PLOT SPLIT BY TREATMENTS, COLOR = ORDER, WRAP = SYNC
early <- subset(metamorphs, subset = (order == 'early'))
late  <- subset(metamorphs, subset = (order == 'late'))

emer_lava_order <- ggplot(early, aes(date)) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               fill = "tomato1", colour = "tomato2",
               geom = "ribbon", position = "identity") +
  facet_grid(sync ~. , switch = "y", scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_blank(),
        axis.line.x  = element_line(size = 1),
        strip.text.y = element_text(size = 15, angle = 180)) +
  ylab(NULL) + xlab("Date") +
  stat_density(data = late, size = 1.25, alpha = 0.75, adjust = 1,
               aes(date, ymax = 0, ymin = -..density..),
               fill = "navyblue", colour = "navyblue",
               geom = "ribbon", position = "identity")

## LAVA PLOT SPLIT BY TREATMENTS, COLOR = SYNC, WRAP = ORDER
high_sync <- subset(ind_results, subset = (sync == 'high'))
low_sync  <- subset(ind_results, subset = (sync == 'low'))

emer_lava_sync <- ggplot(high_sync, aes(date, group = as.factor(block))) + mytheme +
  stat_density(size = 1.25, alpha = 0.75, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               fill = "#016c59", colour = "#016c59",
               geom = "ribbon", position = "identity") +
  facet_grid(order ~. , scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text  = element_text(size = 12),
        axis.line.x  = element_line(size = 1),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15),
        axis.title   = element_text(size = 18)) + 
  ylab("density of emerging individuals") + xlab("date") +
  stat_density(data = low_sync, size = 1.25, alpha = 0.75, adjust = 1,
               aes(date, color = as.factor(block), fill = as.factor(block), ymax = 0, ymin = -..density..),
               fill = "#99d8c9", colour = "#99d8c9",
               geom = "ribbon", position = "identity") +
  geom_hline(yintercept = 0, size = 1)
emer_lava_sync

emer_lava_sync2 <- ggplot(high_sync, aes(date)) + mytheme +
  stat_density(size = 1.25, alpha = 0.85, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               fill = "#bdc9e1", colour = "#bdc9e1",
               geom = "ribbon", position = "identity") +
  facet_grid(order ~. , scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text.y  = element_text(size = 12),
        axis.line.x  = element_line(size = 1),
        strip.text.y = element_text(size = 15)) + #, angle = 180)) +
  ylab("density") + xlab("date") +
  stat_density(data = low_sync, size = 1.25, alpha = 0.5, adjust = 1,
               aes(date, ymax = ..density.., ymin = 0),
               fill = "#016c59", colour = "#016c59",
               geom = "ribbon", position = "identity")
emer_lava_sync2

# can use this data frame to put the lighter colored distribution on top, but actually looks weird
metamorphs_nomed <- subset(metamorphs, subset = (sync != 'med'))
metamorphs_nomed$sync <- factor(metamorphs_nomed$sync, levels = c("high", "low"))
metamorphs$order <- factor(metamorphs$order, levels = c('cont', 'early', 'same', 'late'))

emer_lava_sync3 <- ggplot(subset(metamorphs, subset = (sync != "med")), aes(x = date, fill = sync, color = sync, alpha = sync)) + mytheme +
  stat_density(size = 1.25, adjust = 1,
               aes(ymax = ..density..,  ymin = 0),
               geom = "ribbon", position = "identity") +
  facet_grid(order ~. , scales = 'free') +
  theme(axis.ticks   = element_blank(), 
        axis.text  = element_text(size = 12),
        axis.line.x  = element_line(size = 1),
        axis.title = element_text(size = 15),
        strip.text.y = element_text(size = 13),
        legend.position = c(0.7, 0.92)) + 
  ylab("density of emerging individuals") + xlab("date") +
  labs(color = "hatching\nsynchrony", fill = 'hatching\nsynchrony', alpha = 'hatching\nsynchrony') +
  scale_fill_manual(values = c("#99d8c9", "#016c59")) +
  scale_color_manual(values = c("#99d8c9", "#016c59"))+
  scale_alpha_manual(values = c(0.85, 0.6)) 
emer_lava_sync3 

###---DISTRIBUTIONS: SURVIVAL--------------------------------

## SURVIVAL, GLOBALLY
surv_dist <- ggplot(survivors, aes(x = no_metas)) + mytheme +
  #geom_histogram(binwidth = 3) + 
  geom_density(size = 1.5, alpha = 0.75, color = 'darkslategray', fill = 'darkslategray') +
  xlab("number of survivors")

## SURVIVAL, BY TREATMENT
surv_dist_trt <- ggplot(subset(survivors, subset = (order != 'cont')), aes(x = no_metas)) + mytheme +
  #geom_histogram(binwidth = 3) +
  geom_density(size = 1.5, alpha = 0.75, color = 'darkslategray', fill = 'darkslategray') +
  facet_grid(sync ~ order) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  xlab("number of survivors")

###---DISTRIBUTIONS: BIOMASS---------------------------------

## BIOMASS, GLOBALLY
biom_dist <- ggplot(survivors, aes(x = biomass)) + mytheme +
  #geom_histogram(binwidth = 3) + 
  geom_density(size = 1.5, alpha = 0.75, color = 'darkslategray', fill = 'darkslategray') +
  xlab("biomass export (mg)")

## BIOMASS, BY TREATMENT
biom_dist_trt <- ggplot(subset(survivors, subset = (order != 'cont')), aes(x = biomass)) + mytheme +
  #geom_histogram(binwidth = 3) +
  geom_density(size = 1.5, alpha = 0.75, color = 'darkslategray', fill = 'darkslategray') +
  facet_grid(sync ~ order) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  xlab("biomass export (mg)")

###---DISTRIBUTIONS: PER CAPITA MASS-------------------------

## MASS, GLOBALLY
mass_dist <- ggplot(metamorphs, aes(x = mass_mg)) + mytheme +
  geom_density(size = 1.5, color = "darkslategray", fill = "darkslategray", alpha = 0.75) +
  xlab('per captia mass (mg)')

## MASS, BY TREATMENT
mass_dist_trt <- ggplot(metamorphs, aes(x = mass_mg)) + mytheme +
  #geom_histogram() + 
  geom_density(size = 1.25, color = "darkslategray", fill = "darkslategray", alpha = 0.75) +
  facet_grid(sync ~ order) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  xlab('per capita mass (mg)')

###---DISTRIBUTIONS: EMERGENCE TIME--------------------------

## EMERGENCE TIMING, GLOBALLY
emer_dist <- ggplot(metamorphs, aes(x = date)) + mytheme +
  geom_density(size = 1.5, color = "darkslategray", fill = "darkslategray", alpha = 0.75) +
  xlab('emergence date')

## EMERGENCE TIMING, BY TREATMENT
emer_dist_trt <- ggplot(ind_results, aes(x = date, color = as.factor(block), fill = as.factor(block))) + mytheme +
  #geom_histogram() +
  geom_density(size = 1.5, alpha = 0.5) + #, color = "darkslategray", fill = "darkslategray") + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = 'none') +
  facet_grid(sync ~ order) +
  xlab("date") + ylab("density of emerging individuals")
emer_dist_trt

## A version with number of metamorphs per day
metamorphs$ones <- 1
accum_tank_day <- ddply(metamorphs, .(tank, doy), transform,
                        metas_today = sum(ones))
accum_tank_day <- select(accum_tank_day, tank, block, trt, order, sync, doy, timetoemerge, metas_today)
accum_tank_day <- unique(accum_tank_day)

## Get mean and se for each day
accum_tank_day_means <- summarySE(data = accum_tank_day, measurevar = "metas_today", groupvars = c("order", "sync", "doy"))
write.csv(accum_tank_day_means, "means_accumulation.csv")
means_accumulation <- read.csv("means_accumulation.csv")

emergence_sum <- ggplot(means_accumulation, aes(x = doy, y = metas_today)) + theme_bw() +
  geom_point() +
  #geom_smooth(span = 0.5, se = F) +
  #geom_ribbon(aes(ymin = 0, ymax = metas_today)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  facet_grid(sync ~ order) +
  xlab("emergence date") +
  ylim(0,4)
emergence_sum

emergence_plot <- ggplot(accum_tank_day, aes(x = doy, y = metas_today, color = as.factor(block))) + theme_bw() +
  geom_point(alpha = 0.5) +
  #geom_line() +
  #geom_ribbon(aes(ymin = 0, ymax = metas_today)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  facet_grid(sync ~ order) +
  xlab("emergence date")
emergence_plot

emergence_sum <- ggplot(means_accumulation, aes(x = doy, y = metas_today)) + theme_bw() +
  geom_point(stat = "summary", fun.y = mean, alpha = 0.5) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  facet_grid(sync ~ order) +
  xlab("emergence date")
emergence_sum

emer_predict <- ggplot(predict_memer, aes(x = date)) + mytheme +
  geom_density(size = 1.5, alpha = 0.75, color = "darkslategray", fill = "darkslategray") + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  facet_grid(sync ~ order) +
  xlab("emergence date")
emer_predict

emer_dist_tank <- ggplot(metamorphs, aes(x = date, group = block, color = as.factor(block), fill = as.factor(block))) + mytheme +
  geom_density(size = 1, alpha = 0.5) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = 'none') +
  facet_grid(sync ~ order) +
  xlab("emergence date")

###---MASS OVER TIME-----------------------------------------

## MASS OVER TIME, GLOBALLY
mass_time <- ggplot(metamorphs, aes(x = date, y = mass_mg)) + mytheme +
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm', size = 2, se = F) +
  ylab("per capita mass (mg)") 

## MASS OVER TIME, BY TREATMENT
mass_time_trt <- ggplot(metamorphs, aes(x = date, y = mass_mg)) + mytheme +
  geom_point(alpha = 0.5, aes(color = as.factor(block))) + 
  geom_smooth(method = 'lm') +
  facet_grid(sync ~ order) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  ylab('per capita mass (mg)')

## MASS OVER TIME, BY TREATMENT, SUMMARIZED
mass_time_trt2 <- ggplot(metamorphs, aes(x = date, y = mass_mg, color = order)) + mytheme +
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm', size = 2, se = F) +
  ylab("per capita mass (mg)") +
  facet_grid(. ~ sync)

###---PROGRESS PLOTS-----------------------------------------

## TREATMENT EFFECTS, BOX PLOTS
surv_box
biom_box
mass_box
emer_box

## TREATMENT EFFECTS, MEANS +/- 95% CI
surv_mean
biom_mean
mass_mean
emer_mean

## TREATMENT EFFECTS, CONTINUOUS, ABSOLUTE, NO CONTROLS
# add control values for reference
surv_abs_nocont
biom_abs_nocont
mass_abs_nocont
emer_abs_nocont

## TREATMENT EFFECTS, RELATIVE TO CONTROLS
surv_cont 
biom_cont
mass_cont
emer_cont

## TREATMENT EFFECTS, VIOLIN PLOTS
surv_violin
biom_violin
mass_violin
emer_violin

## DISTRIBUTIONS, GLOBAL
surv_dist
biom_dist
mass_dist
emer_dist

## MASS OVER TIME
mass_time
mass_time_trt
###---LIST OF PLOTS------------------------------------------

## METAMORPH ACCUMULATION
metamorph_accumulation
metamorph_accumulation_block
metamorph_accumulation_trt
metamorph_accumulation_trt_pantala
biomass_accumulation_trt
metamorph_accumulation_tank

## DISTRIBUTIONS, TREATMENTS
surv_dist_trt
biom_dist_trt
mass_dist_trt
emer_dist_trt 
emer_dist_tank

## MISCELLANY
mass_lava        # same as mass_dist_trt, but maybe looks cooler. way to make sync a color, not a facet?
emer_lava        # cool, but not too readable
emer_lava_order  # orange is early, blue is late
emer_lava_sync   # orange is high synchrony, blue is low synchrony

## PREDICTED VS ACTUAL RESULTS (CODE BELOW)
surv_predictxactual  # see fit of model with predicted vs actual values by trt
biom_predictxactual
mass_predictxactual  # seems to over-predict often plus squashes all variation

#############################################################
###------------PLOT VERSIONS FOR PRESENTATIONS------------###
#############################################################

###---TREATMENT EFFECTS: SURVIVAL----------------------------

## CONTROLS ONLY 
surv_mean_cont <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == "surv")), 
                         aes(x = sync, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("number of survivors") + xlab('hatching synchrony') +
  theme(legend.position = 'none') +
  ylim(10, 40)
surv_mean_cont 

## TREATMENTS ONLY
surv_mean_trt <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == "surv")), 
                        aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("number of survivors") + xlab('mean arrival (relative to Rana)') +
  theme(legend.position = 'none') 
surv_mean_trt 

###---TREATMENT EFFECTS: PER CAPITA MASS---------------------

## CONTROLS ONLY
mass_mean_cont <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == "mass")), 
                         aes(x = sync, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("per capita mass at emergence (mg)") + xlab('hatching synchrony') +
  theme(legend.position = 'none')
mass_mean_cont

## TREATMENTS ONLY
mass_mean_trt <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == "mass")), 
                        aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("per capita mass at emergence (mg)") + xlab('mean arrival (relative to Rana)') +
  theme(legend.position = 'none')
mass_mean_trt

###---TREATMENT EFFECTS: EMERGENCE TIME----------------------

## CONTROLS ONLY
emer_mean_cont <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == "emer")), 
                         aes(x = sync, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("days to emergence") + xlab("hatching synchrony") +
  theme(legend.position = 'none') +
  ylim(30, 60)
emer_mean_cont

## TREATMENTS ONLY
emer_mean_trt <- ggplot(subset(means_compiled, subset = (order != 'cont' & variable == "emer")), 
                        aes(x = order, y = mean, color = sync)) + mytheme +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 3, position = position_dodge(width = 0.5),
                aes(ymin = mean - ci, ymax = mean + ci, width = 0)) +
  ylab("days to emergence") + xlab("mean arrival (relative to Rana)") +
  theme(legend.position = 'none')
emer_mean_trt

###---COMMITTEE MEETING FIGURES------------------------------

## Controls only
means_compiled$variable  <- factor(means_compiled$variable, levels = c("surv", "biom", "mass", "emer"))

controls <- ggplot(subset(means_compiled, subset = (order == 'cont')),
                   aes(x = sync, y = mean, fill = sync, group = sync, shape = sync)) + mytheme +
  geom_errorbar(size = 1, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  facet_wrap(.~ variable, nrow = 2, scales = 'free') +
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  xlab('hatching synchrony') + ylab("mean +/- se")
controls

cont_surv <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == 'surv')),
                    aes(x = sync, y = mean, fill = sync, group = sync, shape = sync)) + mytheme +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  ylim(0.2, 1) +
  #labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 18)) +
  xlab(NULL) + ylab("proportional survival")
cont_surv

cont_biom <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == 'biom')),
                    aes(x = sync, y = mean, fill = sync, group = sync, shape = sync)) + mytheme +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  ylim(3000, 8000) +
  #labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 18)) +
  xlab(NULL) + ylab("total biomass export (mg)")
cont_biom

cont_mass <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == 'mass')),
                    aes(x = sync, y = mean, fill = sync, group = sync, shape = sync)) + mytheme +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  #ylim(3000, 8000) +
  #labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 18)) +
  xlab('hatching synchrony') + ylab("per capita mass (mg)")
cont_mass

cont_emer <- ggplot(subset(means_compiled, subset = (order == 'cont' & variable == 'emer')),
                    aes(x = sync, y = mean, fill = sync, group = sync, shape = sync)) + mytheme +
  geom_errorbar(size = 2, position = position_dodge(width = 0.5),
                aes(ymin = mean - se, ymax = mean + se, width = 0)) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#bdc9e1", "#67a9cf", "#016c59")) +
  scale_shape_manual(values = c(21, 24, 23)) +
  #ylim(3000, 8000) +
  #labs(shape = "hatching\nsynchrony", fill = 'hatching\nsynchrony') +
  theme(legend.position = 'none',
        axis.title = element_text(size = 18)) +
  xlab('hatching synchrony') + ylab("days to emergence")
cont_emer
grid.arrange(cont_surv, cont_biom, cont_mass, cont_emer)
