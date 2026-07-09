# Alice Palmer
#11/28/22

#read in data
farmData <- read.csv("2012_fig5_farm_phenotype.csv")

#need to calculate stem density
#h2o_g is stem volume
#stem_g is dry weight

library(dplyr)
library(emmeans)
library(glmmTMB)
library(car)
library(ggplot2)
library(stringr)
library(forcats)
library(patchwork)

farmData <- farmData %>% 
  mutate(stem_density = stem_g/h2o_g)

#need total branches
#opp_br are opposite branches <-two branches for each opp_br
#alt_br are alternating branches

farmData <- farmData %>% 
  mutate(total_br = alt_br+(2*opp_br))

#Is there a difference between winteri and annuus in these traits? ####

#Anova of stem density ####
#fit
density_all <- glm(stem_density ~ sp + sp:pop, data=farmData)

#assumptions
plot(density_all, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(density_all) #species and sp:pop are both significant
#sp: p = 0.0001034
#sp:pop: p = 0.0007738
#sp:pop is the effect of population once we've taken species into account

#plot
stemDensityPlot1 <- ggplot(data=farmData, mapping=aes(x=pop, y=stem_density, color=sp))+
  geom_boxplot(outliers=F)+
  geom_jitter()
stemDensityPlot1


#recode populations to differentiate allopatric and sympatric annuus in metapop column
HW <- c("Boyd", "HWb", "HWE")
HA_s <- c("Acad", "Canal", "DRey")
HA_a <- c("LA", "Man")

farmData <- farmData %>% 
  mutate(metapop = ifelse(pop %in% HW, "HW",
                          ifelse(pop %in% HA_s, "HA_s",
                                 ifelse(pop %in% HA_a, "HA_a",
                                        NA))) )

#reorder to group by metapop
farmData <- farmData %>% 
  mutate(pop = fct_relevel(pop, c("Acad", "Canal", "DRey", "Boyd", "HWE", "HWb", "Man", "LA")))

stemDensityPlot <- ggplot(data=farmData,
                           aes(x=pop, y=stem_density, color=metapop))+
  geom_boxplot(outliers=F)+
  geom_jitter(width=0.3, size=0.3, alpha=0.8)+
  scale_color_manual(values = c("HW" = "#225ea8", 
                                "HA_s" ="#a1d99b",
                                "HA_a"="#31a354")) + #change colors
  theme_light(base_size = 6)+ #change theme
  xlab("")+ #need to add axis/facet/key labels
  ylab("Stem Density")+
  scale_x_discrete(labels = c('Academy','Canal','Del Rey', 'Boyd', 'Valley', 'Yokuts', 'Manteca', 'Mountain'))+ #change population labels to full names
  theme(legend.position = "none")+ #will add for right-most figure
  theme(axis.text.x=element_text(angle = 25, hjust = 0.82))+
  ggtitle("D")#add letter title for multipanel figure
stemDensityPlot

farmData %>% 
  filter(!is.na(stem_density)) %>% 
  group_by(sp) %>% 
  summarise(meanStemDensity = mean(stem_density))

#Anova of days to flowering ####
#fit
fdays_all <- glm(fdays ~ sp + sp:pop, data=farmData)

#assumptions
plot(fdays_all, which=c(1,2,4,5))
#looks good

#ANOVA
Anova(fdays_all) #species and sp:pop are both significant
#sp: p = 1.482e-10
#sp:pop: p = 0.03629

farmData %>% 
  filter(!is.na(fdays)) %>% 
  group_by(sp) %>% 
  summarise(meanDays = mean(fdays))

fDaysPlot <- ggplot(data=farmData,
                          aes(x=pop, y=fdays, color=metapop))+
  geom_boxplot(outliers=F)+
  geom_jitter(width=0.3, size=0.3, alpha=0.8)+
  scale_color_manual(values = c("HW" = "#225ea8", #change colors
                                "HA_s" ="#a1d99b",
                                "HA_a"="#31a354"),
                     name = "", #remove legend title
                     labels = c("Allopatric H. annuus", "Sympatric H. annuus", "H. winteri")) + 
  theme_light(base_size = 6)+ #change theme
  xlab("")+ #need to add axis/facet/key labels
  ylab("Days to Flowering")+
  scale_x_discrete(labels = c('Academy','Canal','Del Rey', 'Boyd', 'Valley', 'Yokuts', 'Manteca', 'Mountain'))+ #change population labels to full names
  theme(legend.position = "right")+ #will add for right-most figure
  theme(axis.text.x=element_text(angle = 25, hjust = 0.82))+
  ggtitle("E")#add letter title for multipanel figure
fDaysPlot


#Anova of height ####
#fit
height_all <- glm(height ~ sp + sp:pop, data=farmData)

#assumptions
plot(height_all, which=c(1,2,4,5))
#looks good

#ANOVA
Anova(height_all) #species is significant. sp:pop is not
#sp: p = 3.183e-06
#sp:pop: p = 0.1196

farmData %>% 
  filter(!is.na(height)) %>% 
  group_by(sp) %>% 
  summarise(meanDays = mean(height))

heightPlot <- ggplot(data=farmData,
                          aes(x=pop, y=height, color=metapop))+
  geom_boxplot(outliers=F)+
  geom_jitter(width=0.3, size=0.3, alpha=0.8)+
  scale_color_manual(values = c("HW" = "#225ea8", 
                                "HA_s" ="#a1d99b",
                                "HA_a"="#31a354")) + #change colors
  theme_light(base_size = 6)+ #change theme
  xlab("")+ #need to add axis/facet/key labels
  ylab("Height")+
  scale_x_discrete(labels = c('Academy','Canal','Del Rey', 'Boyd', 'Valley', 'Yokuts', 'Manteca', 'Mountain'))+ #change population labels to full names
  theme(legend.position = "none")+ #will add for right-most figure
  theme(axis.text.x=element_text(angle = 25, hjust = 0.82))+
  ggtitle("A")#add letter title for multipanel figure
heightPlot

#Anova of circumference ####
#fit
circum_all <- glm(stem_circum ~ sp + sp:pop, data=farmData)

#assumptions
plot(circum_all, which=c(1,2,4,5))
#looks good

#ANOVA
Anova(circum_all) #neither is significant
#sp: p = 0.06249
#sp:pop: p = 0.13459

circumPlot <- ggplot(data=farmData,
                          aes(x=pop, y=stem_circum, color=metapop))+
  geom_boxplot(outliers=F)+
  geom_jitter(width=0.3, size=0.3, alpha=0.8)+
  scale_color_manual(values = c("HW" = "#225ea8", 
                                "HA_s" ="#a1d99b",
                                "HA_a"="#31a354")) + #change colors
  theme_light(base_size = 6)+ #change theme
  xlab("")+ #need to add axis/facet/key labels
  ylab("Stem Circumference")+
  scale_x_discrete(labels = c('Academy','Canal','Del Rey', 'Boyd', 'Valley', 'Yokuts', 'Manteca', 'Mountain'))+ #change population labels to full names
  theme(legend.position = "none")+ #will add for right-most figure
  theme(axis.text.x=element_text(angle = 25, hjust = 0.82))+
  ggtitle("B")#add letter title for multipanel figure
circumPlot

#Anova of branches ####
#fit
br_all <- glm(total_br ~ sp + sp:pop, data=farmData)

#assumptions
plot(br_all, which=c(1,2,4,5))
#looks good

#ANOVA
Anova(br_all) #species is significant. sp:pop is not
#sp: p = 9.677e-11
#sp:pop: p = 0.6954

farmData %>% 
  filter(!is.na(total_br)) %>% 
  group_by(sp) %>% 
  summarise(meanDays = mean(total_br))

branchesPlot <- ggplot(data=farmData,
                          aes(x=pop, y=total_br, color=metapop))+
  geom_boxplot(outliers=F)+
  geom_jitter(width=0.3, size=0.3, alpha=0.8)+
  scale_color_manual(values = c("HW" = "#225ea8", 
                                "HA_s" ="#a1d99b",
                                "HA_a"="#31a354")) + #change colors
  theme_light(base_size = 6)+ #change theme
  xlab("")+ #need to add axis/facet/key labels
  ylab("Total Branches")+
  scale_x_discrete(labels = c('Academy','Canal','Del Rey', 'Boyd', 'Valley', 'Yokuts', 'Manteca', 'Mountain'))+ #change population labels to full names
  theme(legend.position = "none")+ #will add for right-most figure
  theme(axis.text.x=element_text(angle = 25, hjust = 0.82))+
  ggtitle("C")#add letter title for multipanel figure

branchesPlot

#there is a difference between winteri and annuus in stem density, days until flowering, total number of branches, and height. there is no difference in stem circumference. 

#add plots together

#Is there a difference between Boyd and the other winteri populations in these traits? ####

#remove annuus populations from dataset
winteriFarmData <- farmData %>% 
  filter(sp=="HW")

#Anova of stem density winteri ####
#fit
density_winteri <- glm(stem_density ~ pop, data=winteriFarmData)

#assumptions
plot(density_winteri, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(density_winteri) #pop is significant

#pop: p = 0.006614

#Tukey's HSD
contrast(emmeans(density_winteri, ~pop), method="tukey")
#HWb-Boyd p = 0.0120 *
#HWE-Boyd p = 0.9901
#HWE-HWb p = 0.0234 *
#HWb is the odd one out, not Boyd

#Anova of days-to-flowering winteri ####
#fit
fdays_winteri <- glm(fdays ~ pop, data=winteriFarmData)

#assumptions
plot(fdays_winteri, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(fdays_winteri) #pop is significant
#pop: p = 0.009932

#Tukey's HSD
contrast(emmeans(fdays_winteri, ~pop), method="tukey")
#HWb-Boyd p = 0.1638 
#HWE-Boyd p = 0.0075 *
#HWE-HWb p = 0.4771 
#Boyd is different from HWE but not HWb

winteriFarmData %>% 
  filter(!is.na(fdays)) %>% 
  group_by(pop) %>% 
  summarise(meanDays = mean(fdays))

#Anova of height winteri ####
#fit
height_winteri <- glm(height ~ pop, data=winteriFarmData)

#assumptions
plot(height_winteri, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(height_winteri) #pop is significant
#pop: p = 0.03793

#Tukey's HSD
contrast(emmeans(height_winteri, ~pop), method="tukey")

#HWb-Boyd p = 0.9993
#HWE-Boyd p = 0.0648
#HWE-HWb p = 0.0668 

#no post hoc comparisions are significant

#Anova of stem circumference winteri ####
#fit
circum_winteri <- glm(stem_circum ~ pop, data=winteriFarmData)

#assumptions
plot(circum_winteri, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(circum_winteri) #pop is not significant
#pop: p = 0.0787

#Anova of branches winteri ####
#fit
br_winteri <- glm(total_br ~ pop, data=winteriFarmData)

#assumptions
plot(br_winteri, which=c(1,2,4,5))
#looks okay

#ANOVA
Anova(br_winteri) #pop is not significant
#pop: p = 0.543

#Overall results: ####

#There is a difference between winteri and annuus in stem density, days until flowering, total number of branches, and height. There is no difference in stem circumference. 

#HWb, but not Boyd, differs from the other winteri populations in stem density. 
#Boyd is different from HWE, but not HWb, in days-til-flowering
#No post-hoc comparisons were significant for height, and stem circumference and branches were not significantly different when only winteri populations were compared.

#Therefore, stem density, days until flowering, height, and total number of branches are traits which differ between winteri and annuus
#Boyd was not significantly different from the rest of the winteri populations in any of these traits

  