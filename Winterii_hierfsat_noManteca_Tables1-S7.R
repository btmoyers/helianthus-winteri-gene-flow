#Author: Alice & Brook
#Date: 1/13/23 updated 3/2/26

#THIS IS WITH ONLY MAN1, MAN2, AND MAN3 REMOVED

library(tidyverse)
library(hierfstat)
library(adegenet)

#import VCF with heirfstat's read.VCF function
sunflowers_vcf <- read.VCF("~/Documents/Alice_winterii/filtered_hw_no_manteca.recode.vcf", BiAllelic = T, convert.chr = FALSE)

#import sample metadata
popfile <- read.table("~/Documents/Alice_winterii/HW_populations.txt", header = T) 

#define dataset with only loci
loci <- as.matrix(sunflowers_vcf)
hw_loci <- as.matrix(sunflowers_vcf[popfile$species=="HW",])
ha_loci <- as.matrix(sunflowers_vcf[popfile$species=="HA",])

# F stats for dosage (population-level)
pop_dos <- fs.dosage(dos = loci, pop=popfile$population)
pop_dos$Fi
pop_dos$FsM
pop_dos$Fst2x2

#Academy       Boyd      Canal     DelRey   Manteca  Mountain    Valley    Yokuts
#Academy          NA 0.08819353 0.06215245 0.06787791 0.1033789 0.1428047 0.1249145 0.1317715
#Boyd     0.08819353         NA 0.08064630 0.11818692 0.1382013 0.1636330 0.1283641 0.1315419
#Canal    0.06215245 0.08064630         NA 0.08417891 0.1123918 0.1505258 0.1187420 0.1227066
#DelRey   0.06787791 0.11818692 0.08417891         NA 0.1200466 0.1763354 0.1508581 0.1514555
#Manteca  0.10337886 0.13820129 0.11239175 0.12004656        NA 0.1246247 0.1788024 0.1782634
#Mountain 0.14280472 0.16363296 0.15052581 0.17633539 0.1246247        NA 0.2186619 0.2125167
#Valley   0.12491452 0.12836405 0.11874202 0.15085806 0.1788024 0.2186619        NA 0.0559838
#Yokuts   0.13177151 0.13154194 0.12270658 0.15145554 0.1782634 0.2125167 0.0559838        NA

pop_dos$Fs

#Academy      Boyd     Canal    DelRey    Manteca   Mountain    Valley    Yokuts      All
#Fis 0.1495887 0.2112510 0.1895570 0.1472200  0.2716172 0.31360789 0.1606337 0.1688843 0.201545
#Fst 0.1052245 0.1637493 0.1450014 0.1782208 -0.1178780 0.03979645 0.2662897 0.2770037 0.132176

# F stats for dosage (species-level)
spp_dos <- fs.dosage(dos = loci, pop=popfile$species)
spp_dos$Fi
spp_dos$FsM
spp_dos$Fst2x2
spp_dos$Fs

#HA        HW       All
#Fis  0.32111793 0.2462226 0.2836702
#Fst -0.02287187 0.1788921 0.0780101

pi.dosage(dos = hw_loci) #837.5608
theta.Watt.dosage(dos = hw_loci) #679.2634
TajimaD.dosage(dos = hw_loci) #0.8336905

pi.dosage(dos = ha_loci) #1052.552
theta.Watt.dosage(dos = ha_loci) #1018.537
TajimaD.dosage(dos = ha_loci) #0.1147344

# You can see the pop structure here, along with a few weirdos
beta <- beta.dosage(dos = as.matrix(sunflowers_vcf))

image(z = beta, main="Kinship and Inbreeding (diagonal) \n in California Helianthus",xlab="",ylab="")
# Aca4 looks weird, has the most missing data. Also Man 4, which could be due to introgression??

## below not correct for dosage data, need to convert to genotypes if we want to use
## could just do str replacement (0 --> 00, 1 --> 01, 2 --> 11) 
# this annoyingly is not working
loci_dip <- as.data.frame(loci) %>% 
  replace_values(., "0" ~ "00", "1" ~ "01", "2" ~ "11")

#define levels
pop <- as.data.frame(popfile$population)
spp <- as.data.frame(popfile$species)

#estimate hierarchical f stats
spp_wc <- varcomp.glob(levels=cbind(spp,pop), loci)
spp_wc$overall
spp_wc$F

# more pop level variance than species level, but negative for individual?
#hmmm, is the data being interpreted correctly by these functions?

pg <- basic.stats(cbind(spp,loci)) #heterozygousities and Fsts and etc
pg$overall
#     Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.2540  0.1847  0.1870  0.0023  0.1893  0.0046  0.0124  0.0245 -0.3757  0.0057 

pg$pop.freq

w <- wc(cbind(spp,loci))

w$FST #weirdly quite different from the one from basic.states
#0.144706

#test winteri and annuus

#add sp column
#library(dplyr)
#sunflowers 
#test.between(loci, test=sunflowers$species, nperm=1000, rand.unit = sample)

#convert to structure file (for use later)
sunflowers$indv <- row.names(sunflowers) #add column with individual names
library(dplyr)
sunflowers <- sunflowers %>% #move indv column so it's first
  select(indv, everything())
write.struct(sunflowers,ilab=sunflowers$indv,pop=sunflowers$pop,fname="filtered-winteri-no-Manteca.str")
