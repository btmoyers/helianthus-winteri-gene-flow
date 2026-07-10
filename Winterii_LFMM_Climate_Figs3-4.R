setwd("~/Documents/Alice_winterii/")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("LEA")
#BiocManager::install("qvalue")

library(tidyverse)
library(patchwork)

#map stuff
library(terra) 
library(tidyterra)
library(geodata)
library(sf) #yes

#genotype stuff
library(pegas)
library(adegenet)
library(LEA)

library(lfmm)     # Used to run LFMM
library(vegan)    # Used to run PCA & RDA
library(qvalue)   # Used to post-process LFMM output

library(pcadapt)

#First climate analysis
# download climate variables as a SpatRaster (only actually downloads the first time)
r <- worldclim_global(var="bio", res=5) #fun, discovered that the worldclim_country function for the USA does not download data appropriately
r
crs(r)

# set a SpatExtent
s <- ext(x = c(-124,-115,33,39))

# crop the SpatRaster by the SpatExtent
r <- crop(x = r, y = s)

# also download a vector of CA, as an sf object
ca <- USAboundaries::us_states(resolution = "low", states = "California")

# make it an sp object
ca.sp <- sf::as_Spatial(ca)

# load pop metadata
pop_file <- read.table("HW_populations.txt", header = TRUE)

# creating coordinates
coords <- data.frame(x=pop_file$long,y=pop_file$lat) 
coords_uniq <- pop_file %>% 
  select(species, population, speciescolor, popcolor, lat, long) %>% 
  distinct

# Adding coordinates(points) to United States Bioclim data 
values <- terra::extract(x = r, coords)

# create a data frame for spatial coordinates and bioclimm variables
env <- cbind(pop_file[,1:4], coords, values)
head(env)

# Plot United States Bioclim variable 12 (annual precip) with spatial coordinates
# working on this
autoplot(r[[12]]) +
  geom_spatvector(data=ca.sp) +
  theme_bw()

# works fine
plot(r[[18]]) 
plot(ca.sp, add = T)
terra::points(cbind(coords_uniq$long,coords_uniq$lat), col = coords_uniq$speciescolor, pch = 20, cex = .75)

#------------------------------------------------------------------------
# Import/export genomic data in VCF format using LEA package
LEA::vcf2lfmm("filtered_hw_no_manteca.recode.vcf")

# rejects the sites (115) with more than 2 alleles
# creates a .geno, .lfmm, .removed (filtered sites for extra alleles), and .vcfsnp ("SNP info")

#------------------------------------------------------------------------
# The function snmf() can be run on the data with missing genotypes as follows. The completion of the genotypic matrix is based on estimated ancestry coefficients and ancestral genotype frequencies

#project.missing = snmf("~/Documents/Alice_winterii/filtered_hw_no_manteca.recode.geno", K = 1:10,
#                       entropy = TRUE,
#                       repetitions = 10,
#                       iterations = 1000,
#                       project = "new")

## To load the project, use:
project.missing = load.snmfProject("filtered_hw_no_manteca.recode.snmfProject")

##To remove the project, use:
#  remove.snmfProject("filtered_hw_no_manteca.recode.snmfProject")

# plot cross-entropy criterion for all runs in the snmf project (4), Figure S1
plot(project.missing, col = "blue", pch = 19, cex = 1.2)

# select the run with the lowest cross-entropy value for the best K (4)
best = which.min(cross.entropy(project.missing, K = 4))

my.colors <- c("tomato", "lightblue",
              "olivedrab", "gold")

barchart(project.missing, K = 4, run = best,
        border = NA, space = 0,
        col = my.colors,
        xlab = "Individuals",
        ylab = "Ancestry proportions",
        sort.by.Q = F,
        lab = pop_file$sample) -> bp

axis(1, at = 1:length(bp$order), 
     labels = pop_file$population[bp$order], las = 3, 
     cex.axis = .4)

#okay let's save the best run

#K4Q <- Q(project.missing, K = 4, run = best)
#K4Q <- cbind(popfile$population,as.data.frame(K4Q))
#write.csv(K4Q, file = "artesinal-admixture.csv")


# The snmf project data can be used to impute the missing data as follows
#LEA::impute(project.missing, "filtered_hw_no_manteca.recode.geno", method = 'mode', K = 4, run = best)
#----------------------------------------------------------------------------------

# Need to reduce the sampling to one per population?
pred <- env %>% 
  select(-sample, -species, -metapopulation, -x, -y, -ID) %>% 
  distinct()

pred_mat <- as.matrix(pred[,-1],dimnames = list(pred[,1],colnames(pred[,-1])))

pred.pca <- vegan::pca(pred_mat, scale=T)
plot(pred.pca) 
#PC1 differentiates Mountain (7) from 1:4 (Acad, Boyd, Canal, DelRey)
#PC2 differentiates Manteca & DelRey from Yokuts & Valley (which have the same values)
plot(pred.pca, choices =c(1,3))
#PC3 differentiates Boyd and Mountain (with Acad & Canal close) from Manteca (and to a lesser degree Yokuts & Valley)

summary(pred.pca)$cont #first three PCs are relevant with eigenvalue > 1
x <- data.frame(summary(pred.pca)$cont)
#write.table(x, file = "BIOCLIM_PCs.txt", sep = "\t")

screeplot(pred.pca, main = "Screeplot of HW-HA Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
# Or just the first two PCs by the broken stick method

## correlations between the PC axis and predictors:
scores(pred.pca, choices=1:3, display=c("species", scaling=0), digits=3)

predload <- round(scores(pred.pca, choices=1:3, display=c("species", scaling=0), digits=3))
predload
#write.table(predload, file = "BIOCLIM_loading.txt", sep = "\t")

pred.PC123 <- scores(pred.pca, choices=1:3, display=c("sites", scaling=0))

pred <- cbind(pred,pred.PC123)

pop <- c("Academy","Boyd","Canal","DelRey","Yokuts & Valley","Valley","Mountain","Manteca")
Taxon <- c("H. annuus","H. winteri","H. annuus","H. annuus","H. winteri","H. winteri","H. annuus","H. annuus")

envpca <- ggplot(pred, aes(x=PC1, y=PC2, label = pop)) +
  geom_text(check_overlap = T,vjust = 0, nudge_y = -.18, color = c("#31a354","#225ea8","#31a354","#31a354","#225ea8","#225ea8","#31a354","#31a354"), size = 2.3) +
  geom_point(color = c("#31a354","#225ea8","#31a354","#31a354","#225ea8","#225ea8","#31a354","#31a354"),cex = 1.2) +
  labs(x = "PC1 (48.3 PVE)", y = "PC2 (34.0 PVE)") +
  xlim(-1.3,3.25) +
  theme_bw()

envpca

bc12 <- ggplot(pred, aes(x=Taxon, y = wc2.1_5m_bio_12, color = Taxon)) +
  #geom_jitter(show.legend = F) +
  geom_boxplot(show.legend = F) +#, fill = "#00000000") +
  scale_color_manual(values = c("H. winteri" = "#225ea8",
                                "H. annuus" ="#31a354")) +
  labs(y = "Annual precipiation (mm)") +
  theme_bw()

bc13 <- ggplot(pred, aes(x=Taxon, y = wc2.1_5m_bio_13, color = Taxon)) +
  #geom_jitter(show.legend = F) +
  geom_boxplot(show.legend = F) +
  scale_color_manual(values = c("H. winteri" = "#225ea8",
                                "H. annuus" ="#31a354")) +
  labs(y = "Precipitation of Wettest Month") +
  theme_bw()

bc16 <- ggplot(pred, aes(x=Taxon, y = wc2.1_5m_bio_16, color = Taxon)) +
  #geom_jitter(show.legend = F) +
  geom_boxplot(show.legend = F) +
  scale_color_manual(values = c("H. winteri" = "#225ea8",
                                "H. annuus" ="#31a354")) +
  labs(y = "Precipitation of Wettest Quarter") +
  theme_bw()

bc19 <- ggplot(pred, aes(x=Taxon, y = wc2.1_5m_bio_19, color = Taxon)) +
  #geom_jitter(show.legend = F) +
  geom_boxplot(show.legend = F) +
  scale_color_manual(values = c("H. winteri" = "#225ea8",
                                "H. annuus" ="#31a354")) +
  labs(y = "Precipitation of Coldest Quarter") +
  theme_bw()

#all four look so similar, let's just do the annual precip

envpca + bc12 + plot_layout(widths = c(3, 1)) + plot_annotation(tag_levels = 'A')

ggsave("bioclim_pca_Fig4.pdf", device = "pdf", units = "in", width = 6.5, height = 3)
ggsave("bioclim_pca_Fig4.png", device = "png", units = "in", width = 6.5, height = 3)

#Okay, need to expand the population-level PCs back out to the individual dataset
pop.n <- pop_file %>% 
  summarise(n = n(), .by = population) %>% 
  select(n)

pred.PC123 <- data.frame(pred.PC123)

# check that the function is accurate
cbind(pop_file$population,rep(row.names(pred.PC123), times = pop.n$n))

rep_pred.PC123 <- pred.PC123[rep(row.names(pred.PC123), times = pop.n$n), ]

rownames(rep_pred.PC123) <- 1:nrow(rep_pred.PC123)

#------------------------------------------------------------------------
# Running LFMM (Univariate GEA)
# Run lfmm at K=4

# Import imputation results
dat.imp <- read.lfmm("filtered_hw_no_manteca.recode.lfmm_imputed.lfmm")

# Genetic scree plot
gen.pca <- vegan::pca(dat.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")
# K = 3 is best by the broken stick method, but I think 4 is okay

K <- 4

hw.lfmm <- lfmm_ridge(Y=dat.imp, X=rep_pred.PC123$PC2, K=K) ## change K as you see fit

hw.pv <- lfmm_test(Y=dat.imp, X=rep_pred.PC123$PC2, lfmm=hw.lfmm, calibrate="gif")
names(hw.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
hw.pv$gif     #3.86124   Maybe too high?

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(hw.pv$pvalue, main="Unadjusted p-values")        
hist(hw.pv$calibrated.pvalue, main="GIF-adjusted p-values")  #actually it is pretty flat    
#------------------------------------------------------------------------
# How to manually adjust p-values:
# Let's change the GIF and readjust the p-values:
zscore <- hw.pv$score[,1]   # zscores for first predictor, we only have one in our case...
gif <- hw.pv$gif       ## default GIF for this predictor

# Adjusted to 0.5
new.gif <- 3              ## choose your new GIF

adj.pv1 <- pchisq(zscore^2/new.gif, df=1, lower = FALSE)

# hist(argo.pv$pvalue[,1], main="Unadjusted p-values")        
# hist(argo.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=3.0)")
hist(adj.pv1, main="REadjusted p-values (GIF=2.0)")

# the standard adjustment seems better

#------------------------------------------------------------------------
# convert adjusted p values to q values
hw.qv <- qvalue(hw.pv$calibrated.pvalue)$qvalues

length(which(hw.qv < 0.05)) ## how many SNPs have an FDR < 5%? #25

# Using K=4 and default GIF calculated from lfmmm, and an FDR threshold of 0.05, we only detected 25 candidate SNPs under selection in response to our PC2 environmental predictor.

## identify which SNPs these are
snps <- read.table("filtered_hw_no_manteca.recode.vcfsnp")
colnames(snps) <- c("chr","pos","locus","a1","a2","V6","V7","V8","V9")
snps <- cbind(snps, hw.qv)

hw.FDR.PC2 <- snps[which(hw.qv < 0.05),]
write_tsv(hw.FDR.PC2, file = "lfmm.txt")

## add in pcadapt
hw.geno <- read.pcadapt("./filtered_hw_no_manteca.recode.lfmm", type = "lfmm")

hw.pcadapt <- pcadapt(input = hw.geno, K = 10)

plot(hw.pcadapt, option = "screeplot") # a K of 3 is reasonable, although the screeplot continues to decay gradually to 8, could also argue for K = 4

pve <- hw.pcadapt$singular.values^2
pve #[1] 0.10676865 0.05816172 0.04675001 0.03370765 

plot(hw.pcadapt, option = "scores", pop = pop_file$population)
plot(hw.pcadapt, option = "scores", i = 3, j = 4, pop = pop_file$population) # K = 4 does clump populations a little
plot(hw.pcadapt, option = "scores", i = 5, j = 6, pop = pop_file$population) # starts to get less clumpy here

# okay, keep K = 4 for this analysis, maybe use K = 4 above too??
hw.pcadapt <- pcadapt(input = hw.geno, K = 4)

summary(hw.pcadapt)
plot(hw.pcadapt, option="manhattan") #maybe a peak or two, unclear how these are ordered
plot(hw.pcadapt, option="qqplot") #comes way off but looks better than their tutorial
plot(hw.pcadapt, option="stat.distribution") #good fit
hist(hw.pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange") # looks good, so I'm okay with the qqplot

qval <- qvalue(hw.pcadapt$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers) #190 loci
pca.out <- snps[outliers,]

snp_pc <- get.pc(hw.pcadapt, outliers) #some SNPs are associated with PC4, so it is worth keeping

pca.out <- cbind(pca.out,qval[outliers],snp_pc)

write_tsv(pca.out, file = "pcadapt.txt")

# homebrew manhattans
snps <- cbind(snps, hw.pcadapt$pvalues, qval)
chroms <- c("HA412HOCHR01", "HA412HOCHR02", "HA412HOCHR03", "HA412HOCHR04", "HA412HOCHR05", "HA412HOCHR06", "HA412HOCHR07", "HA412HOCHR08", "HA412HOCHR09", "HA412HOCHR10", "HA412HOCHR11", "HA412HOCHR12", "HA412HOCHR13", "HA412HOCHR14", "HA412HOCHR15", "HA412HOCHR16", "HA412HOCHR17")

snps17 <- snps %>% 
  filter(chr %in% chroms) %>% 
  mutate(`-log10p` = -log10(`hw.pcadapt$pvalues`))

ggplot(snps17, aes(x = pos, y = -log10(hw.qv))) +
  facet_wrap(~chr) +
  geom_point() +
  geom_abline(slope = 0, intercept = -log10(0.05), color = "red") +
  theme_bw()

ggplot(snps17, aes(x = pos, y = -log10(qval))) +
  facet_wrap(~chr) +
  geom_point() +
  geom_abline(slope = 0, intercept = -log10(0.05), color = "red") +
  theme_bw()

# look at overlap
overlap <- pca.out %>% 
  dplyr::filter(locus %in% hw.FDR.PC2$locus)
dim(overlap) #11 loci, all associated with PC1 & 3

#write.table(overlap, file = "pcadapt-lfmm-overlap.txt")

pca.out %>% 
  group_by(PC) %>% 
  tally

#1     1    47
#2     2   125
#3     3    11
#4     4     7

# make set of just the PC1,3,4 and LFMM markers
pca.out %>% 
  filter(PC != 2)


# make figure 3 for the paper

pc.scores <- as.data.frame(hw.pcadapt$scores)
colnames(pc.scores) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
pc.scores <- cbind(pop_file,pc.scores)

pop.scores <- pc.scores %>% 
  group_by(population) %>% 
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2),
            PC3 = mean(PC3),
            PC4 = mean(PC4))

library(ggrepel)

pca12 <- ggplot(data = pc.scores, aes(x = PC1, y = PC2, color = species)) +
  geom_point(show.legend=F, cex = 1) +
  geom_label_repel(data = pop.scores, aes(label = population), alpha = 0.75, size = 2,
            color = c("#31a354","#225ea8","#31a354","#31a354","#31a354","#31a354","#225ea8","#225ea8")) +
  scale_color_manual(values = c("HW" = "#225ea8",
                                "HA" ="#31a354")) +
  labs(x = "PC1 (10.7 PVE)", y = "PC2 (5.8 PVE)") +
  theme_bw()

pca34 <- ggplot(data = pc.scores, aes(x = PC3, y = PC4, color = species)) +
  geom_point(show.legend=F, cex = 1) +
  geom_label_repel(data = pop.scores, aes(label = population), alpha = 0.75, size = 2,
                   color = c("#31a354","#225ea8","#31a354","#31a354","#31a354","#31a354","#225ea8","#225ea8")) +
  scale_color_manual(values = c("HW" = "#225ea8",
                                "HA" ="#31a354")) +
  labs(x = "PC3 (4.7 PVE)", y = "PC4 (3.4 PVE)") +
  theme_bw()

pcas <- pca12 + pca34

# now the ancestry
K4Q.art <- read.csv(file = "artesinal-admixture-rearrange.csv")
ind <- 1:78
K4Q.art <- cbind(ind,K4Q.art)
anc <- K4Q.art %>% 
  pivot_longer(cols=c(-ind,-pop), names_to = "K", values_to = "Prop")

adplot <- ggplot(anc, aes(x=ind,y=Prop, fill = K)) +
  geom_col(show.legend = F) +
  scale_fill_manual(values=c("#4daf4a","#984ea3","#377eb8","#e41a1c")) +
  scale_x_continuous(n.breaks = 78) +
  coord_cartesian(ylim = c(0,1),expand = FALSE, clip = "off") +
  labs(x = "Population", y = "Ancestry Proportion") +
  annotate(geom = "text", x = seq_len(length(ind)), y = -0.13, label = K4Q.art$pop, size = 1.75, angle = 90, color = c(rep("#31a354",47),rep("#225ea8", 31))) +
  #annotate(geom = "text", x = 78/2, y = -0.3, label = "Population", size = 4) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.margin = margin(1, 1, 2, 1, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

adplot / pcas + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1, 1.3)) 

ggsave("paper_Fig3.pdf", device = "pdf", units = "in", width = 6.5, height = 5.5)
ggsave("paper_Fig3.png", device = "png", units = "in", width = 6.5, height = 5.5)

#------------------------------------------------------------------------------
