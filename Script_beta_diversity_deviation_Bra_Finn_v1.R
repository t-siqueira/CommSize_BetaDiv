### Community size affects the signals of niche selection and ecological drift 
# on beta diversity in stream organisms: data and script 

## This code calculates 4 types of beta-diversity metric for each of 100 
# watersheds (5 streams) in Brazil and Finland.
# Beta diversity: Sorensen and Bray-Curtis dissimilarity between all 
# pairs 
# Beta deviation from null models: Raup-Crick (vegan version) and 
# Bray-Curtis Beta-deviation (based on the scripts by Chris Catano and 
# Jonathan Myers)

## These beta diversity metrics are modelled against community size, 
# environmental heterogeneity and spatial extent.

# Prepared by Tadeu Siqueira (tsiqueira.bio@gmail.com) in December/2018
# Revised by T. Siqueira in March/2019
#===============================================================================

## Required packages
library(vegan)
library(reshape2)
library(adespatial)

library(ggplot2)
library(grid)
library(ggpubr)
library(gridExtra)


## Helper functions
# This contains the following functions:
# sp.by.site.function; null.model; effect.sizes
source("helper_functions.R") 

### Data preparation ===========================================================
## Load your data
# Data must be a site (rows) by species (columns) matrix
# Some columns may be grouping ID variables (e.g., watersheds, treatments)

bioraw <- read.csv("bio_20_bra.csv", T)

# For the case of our original study on community size-beta diversity,
# we ran analyses separated for tropical and boreal data
# Load boreal data
bioraw_f <- read.csv("bio_20_fin.csv", T) 

# Check local community abudances
sort(rowSums(bioraw[,-c(1,2)]))

# I'm removing sites with less than 10 individuals
good_sites <- rowSums(bioraw[,-c(1,2)])>10
table(good_sites)
bioraw <- bioraw[good_sites,]
dim(bioraw)

# Finland
good_sites_f <- rowSums(bioraw_f[,-c(1,2)])>10
table(good_sites_f)
bioraw_f <- bioraw_f[good_sites_f,]
dim(bioraw_f)

## Create a list in which each element is a compositional matrix (sitexspecies)
# Each element represents a metacommunity replicate
# Here, each of the 20 watersheds is a replicate 
biolist <- split(bioraw[,-c(1,2)], bioraw$catchment)
str(biolist)

# Check if you need to remove any replicate due to insufficient samples
# Here the minimim number of localities is 3
lapply(biolist, FUN=nrow)

# Finland
biolist_f <- split(bioraw_f[,-c(1,2)], bioraw_f$catchment)
str(biolist_f)
lapply(biolist_f, FUN=nrow)

## Keep only the taxa that occur within your metacommunity (if you have)
biolist2 <- lapply(biolist, FUN= function(x) x[,colSums(x)>0])
lapply(biolist2, FUN=nrow)

biolist2_f <- lapply(biolist_f, FUN= function(x) x[,colSums(x)>0])
lapply(biolist2_f, FUN=nrow)


## Estimate community size =====================================================
# Here, because sampling design/effort was standardized
# I am using the mean number of individuals per watershed
# Each watershed has from 3 to 5 streams

local_abund <- rowSums(bioraw[,-c(1,2)])
mean_abund_catch <- aggregate(local_abund, by=list(bioraw$catchment), FUN=mean) 

local_abund_f <- rowSums(bioraw_f[,-c(1,2)])
mean_abund.catch_f <- aggregate(local_abund_f, by=list(bioraw_f$catchment), 
                                FUN=mean)

### Observed beta diversity ====================================================
## Calculate obs. b-diversity for each of the metacommunities
# You can decide which dissimilarity metric you want

# Bray-Curtis
# Estimate the dissimilarity between all pairs of streams in each metacom
beta_obs_bray <- lapply(lapply(biolist2, vegdist, method = "bray"), c)

# Average these pairwise values -- one value per metacomm (watershed)
# Here we have 20 watersheds -- so, 20 b-diversity values
mean_obs_bray <- unlist(lapply(beta_obs_bray, mean))

# Sorensen
beta_obs_sor <- lapply(lapply(biolist2, vegdist, method = "bray", binary=T), c)
mean-obs_sor <- unlist(lapply(beta_obs_sor, mean))

# Finn
beta_obs_bray_f <- lapply(lapply(biolist2_f, vegdist, method = "bray"), c)
mean.obs.bray_f <- unlist(lapply(beta_obs_bray_f, mean))

beta_obs_sor_f <- lapply(lapply(biolist2_f, vegdist, method = "bray", 
                                binary=T), c)
mean_obs_sor_f <- unlist(lapply(beta_obs_sor_f, mean))


### Beta diversity based on null models=========================================
## 1) Raup-Crick
# First I will the Raup–Crick null model to determine how much 
# variation in composition was due only to sampling effects on 
# a-diversity (richness) that are neutral with respect to species identity 
# (Chase et al. 2011)

# Because I am using all occuring species to define the species pool,
# I will calculte the dissimilarity between all pairs using the bioraw data
# I decide to use a fixed gamma diversity (100 streams) because my metacomms
# were small (3 to 5 streams)

beta_raup <- raupcrick (decostand(bioraw[,-c(1,2)], "pa"), 
                         nsimul = 9999, chase = F)

# Scale to vary between -1 and 1
beta_raup <- (beta_raup-.5)*2

# I needed to reorganize the output in beta_raup as it has the dissimilarities
# among all pairs of sites - and I only want those among sites within the same
# metacommunity. So I will create a list, where each element has the 
# dissimilarity values among pairs of streams within a watershed

beta_raup_list <- list()

for (j in 1:nlevels(as.factor(bioraw$catchment))){
  nstream <- table(bioraw$catchment)[j]
  soma <- sum(table(bioraw$catchment)[1:j])-nstream
  beta_raup_list[[j]] <- as.matrix(beta_raup)[(soma+1):(soma+nstream), 
                                                (soma+1):(soma+nstream)]
}

# Now estimate the mean dissimilarity within each watershed
mean_null_raup <- unlist(lapply(beta_raup_list, 
                                 FUN= function(x) mean(as.dist(x))))


# Finland
beta_raup_f <- raupcrick (decostand(bioraw_f[,-c(1,2)], "pa"), 
                        nsimul = 9999, chase = F)

beta_raup_f <- (beta_raup_f-.5)*2

beta_raup_list_f <- list()

for (j in 1:nlevels(as.factor(bioraw_f$catchment))){
  nstream <- table(bioraw.f$catchment)[j]
  soma <- sum(table(bioraw.f$catchment)[1:j])-nstream
  beta_raup_list_f[[j]] <- as.matrix(beta_raup_f)[(soma+1):(soma+nstream), 
                                              (soma+1):(soma+nstream)]
}

mean_null_raup_f <- unlist(lapply(beta_raup_list_f, 
                                FUN= function(x) mean(as.dist(x))))


## 2) Beta deviation from Bray-Curtis ==========================================
# This one needs a specific data format and a number of helper functions

# This is from the helper functions provided by Chris Catano: 
# "MAKE SURE DATA IS IN THE SAME FORMAT AS "indsbyplot" (an array where each 
# replicate - whether it be a different site, treatment, etc.- is an element &
# is in the format shown for "inds.expanded")

# > str(indsbyplot)
# List of 3
# $ 1:'data.frame':	12323 obs. of  4 variables:
# ..$ sp       : int [1:12323] 1 3 6 7 9 10 10 11 13 15 ...
# ..$ quadrat  : int [1:12323] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ plot     : int [1:12323] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ abundance: int [1:12323] 1 1 1 1 1 2 2 1 1 1 ..."

## My data
bio_longformat <- read.csv("bio_20_bra_mod.csv", header=T)
bio_longformat_f <- read.csv("bio_20_fin_mod.csv", header=T)

str(bio_longformat)
bio_longformat <- bio_longformat[good_sites,]
bio_longformat <- bio_longformat[,colSums(bio_longformat[,-c(1,2)])>0]

table(bio_longformat$catchment)

# Finn
bio_longformat_f <- bio_longformat_f[good_sites_f,]
bio_longformat_f <- bio_longformat_f[,colSums(bio_longformat_f[,-c(1,2)])>0]

## Because I am using all sites as the regional pool
# I will treat watersheds as a grouping variable with one level
bio_longformat$stream <- bioraw$stream

bio2 <- melt(bio_longformat, id.vars=c("stream", "catchment"), 
             variable.name = "sp", value.name = "abundance")

bio3 <- data.frame("sp"=as.integer(bio2$sp), "quadrat"=bio2$stream, 
                   "plot"=bio2$catchment, "abundance"=bio2$abundance)
bio3[,"plot"] <- 1 # because my whole data set is the species pool

# Make each individual a row, and then create an array where each element is
# a data.frame corresponding to the whole region

inds_expanded <- bio3[rep(seq.int(1, nrow(bio3)), bio3$abundance), ]
indsbyplot <- split(inds_expanded, inds_expanded$plot)

# Finn
bio_longformat_f$stream <- bioraw_f$stream

bio2_f <- melt(bio_longformat_f, id.vars=c("stream", "catchment"), 
             variable.name = "sp", value.name = "abundance")

bio3_f <- data.frame("sp"=as.integer(bio2_f$sp), "quadrat"=bio2_f$stream, 
                   "plot"=bio2_f$catchment, "abundance"=bio2_f$abundance)
bio3_f[,"plot"] <- 1 # because my whole data set is the species pool

inds_expanded_f <- bio3_f[rep(seq.int(1, nrow(bio3_f)), bio3_f$abundance), ]
indsbyplot_f <- split(inds_expanded_f, inds_expanded_f$plot)

## Create a species x site matrix for each plot (watersheds), where rows are 
# quadrats (streams) and columns are species.
# Here we need to use function "sp.by.site.function" (helper function)
# to a species x sample matrix 

spbysite <- lapply(indsbyplot, sp.by.site.function, species.col = "sp", 
                   site.col = "quadrat")

spbysite_f <- lapply(indsbyplot_f, sp.by.site.function, species.col = "sp", 
                   site.col = "quadrat")

## Calculate mean pair-wise dissimilarity among samples within each replicate. 
# In this example we use Bray-Curtis, but you can use any measure  
# that is most appropriate given the purpose of your particular analysis.
beta_bray_all <- lapply(lapply(spbysite, vegdist, method = "bray"), c)

beta_bray_all_f <- lapply(lapply(spbysite_f, vegdist, method = "bray"), c)

## Run the function "null.model" and calculate the mean and standard deviation
# of the expected (simulated) dissimilarities. 
# "k" is the number of iterations
# Please check the helper function for details about "null.model"

null_models_bray <- lapply(indsbyplot, null.model, k = 9999, species.col = "sp", 
                           site.col = "quadrat")

# Organize observed and simulated dissimilarity values in a list
bray_all_results <- mapply(list, beta_bray_all, null_models_bray, SIMPLIFY = F)

# Finn
null_models_bray_f <- lapply(indsbyplot_f, null.model, k = 9999, 
                             species.col = "sp", site.col = "quadrat")
bray_all_results_f <- mapply(list, beta_bray_all_f, null_models_bray_f, 
                             SIMPLIFY = F)

## Calculate effect sizes and standardized effect sizes 
# with helper functions
# This will result in a list with the following elements:
# ..$ obs.bray
# ..$ sim.bray.mean
# ..$ effect.size
# ..$ std.effect.size

bray_effect_size <- lapply(bray_all_results, FUN=effect.sizes)
str(bray_effect_size)

bray_effect_size_f <- lapply(bray_all_results_f, FUN=effect.sizes)

## Now reorganize these values per watershed as described above for Raup-Crick
# I only load mefa here because its function melt is also part of reshape2
library(mefa)

bray_dist_eff <- vec2dist (bray_effect_size$`1`$std.effect.size, 
                           nrow(bio_longformat))

beta_bray_list <- list()

for (j in 1:nlevels(as.factor(bioraw$catchment))){
  nstream <- table(bioraw$catchment)[j]
  soma <- sum(table(bioraw$catchment)[1:j])-nstream
  beta_bray_list[[j]] <- as.matrix(bray_dist_eff)[(soma+1):(soma+nstream), 
                                                  (soma+1):(soma+nstream)]
}

mean_null_bray <- unlist(lapply(beta_bray_list, 
                                FUN= function(x) mean(as.dist(x))))

# Finn
bray_dist_eff_f <- vec2dist (bray_effect_size_f$`1`$std.effect.size, 
                           nrow(bio_longformat_f))

beta_bray_list_f <- list()

for (j in 1:nlevels(as.factor(bioraw_f$catchment))){
  nstream <- table(bioraw.f$catchment)[j]
  soma <- sum(table(bioraw.f$catchment)[1:j])-nstream
  beta_bray_list_f[[j]] <- as.matrix(bray_dist_eff_f)[(soma+1):(soma+nstream), 
                                                  (soma+1):(soma+nstream)]
}

mean_null_bray_f <- unlist(lapply(beta_bray_list_f, 
                                FUN= function(x) mean(as.dist(x))))


## Env and spatial predictors ==================================================
# Load your env data
# In-stream variables
env <- read.csv("env_20_bra.csv", T)
env_f <- read.csv("env_20_fin.csv", T)

env <- env[good_sites,]
env_f <- env_f[good_sites_f,]

# Catchment scale variables
land400 <- read.csv("land_20_bra.csv", T)
land400_f <- read.csv("land_20_fin.csv", T)

land400 <- land400[good_sites,]
land400_f <- land400_f[good_sites_f,]

names(env)
names(land400)

names(env_f)
names(land400_f)

# Exclude correlated variables
all_env <- cbind(env[,-c(18,19)], land400[,-c(1:5)]) 
all_env_f <- cbind(env_f[,-c(5,21)], land400_f[,-c(1:6)]) 

# Standardize your data if needed
all_env_stnd <- decostand(all_env[,-c(1,2)], "stand")
all_env_stnd_f <- decostand(all_env_f[,-c(1,2)], "stand")

# Organize your env data as a list, just like you did with compositional data  
envlist <- split(all_env_stnd, bioraw$catchment)
env_dists_list <- lapply(envlist, vegdist, method = "eucl")

envlist_f <- split(all_env_stnd_f, bioraw_f$catchment)
env_dists_list_f <- lapply(envlist_f, vegdist, method = "eucl")

# Prepare the spatial data
spa <- read.csv("spa_20.csv", T)
spa_f <- read.csv("spa_20_fin.csv", T)

spa <- spa[good_sites,]
spa_f <- spa_f[good_sites_f,]

# Organize your spa data as a list, just like you did with compositional data
spalist <- split(spa[,-c(1,2)], bioraw$catchment)
spa_obs_eucl <- lapply(spalist, vegdist, method = "eucl")

spalist_f <- split(spa_f[,-c(1,2)], bioraw_f$catchment)
spa_obs_eucl_f <- lapply(spalist_f, vegdist, method = "eucl")

# Use Marti Anderson's betadisper ==============================================
# to estimate mean env dispersion per watershed

env_het <- unlist(
  lapply(
    lapply(
      env_dists_list, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# Post-hoc analysis; among watersheds
permutest(
  betadisper(vegdist(all_env_stnd, "euc"), group= bioraw$catchment, type="med"),
  pairwise = F, permutations = 999)

# Do the same with spatial coordinates to estimate spatial extent per watershed
spa_ext <- unlist(
  lapply(
    lapply(
      spa_obs_eucl, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# And with boreal data
env_het_f <- unlist(
  lapply(
    lapply(
      env_dists_list_f, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# Post-hoc analysis; among watersheds
permutest(
  betadisper(vegdist(all_env_stnd_f, "euc"), group= bioraw_f$catchment, 
             type="med"), pairwise = F, permutations = 999)

spa_ext_f <- unlist(
  lapply(
    lapply(
      spa_obs_eucl_f, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))


## OLS models ==================================================================

mod2_obs_sor <- lm(scale(mean_obs_sor) ~ scale(mean_abund_catch[,2]))
summary(mod2_obs_sor)

# Finn
mod2_obs_sor_f <- lm(scale(mean_obs_sor_f) ~ scale(mean_abund_catch_f[,2]))
summary(mod2_obs_sor_f)


# Now with Bray-Curtis dissimilarities
mod2_obs_bray <- lm(scale(mean_obs_bray) ~ scale(mean_abund_catch[,2]))
summary(mod2_obs_bray)

mod2_obs_bray_f <- lm(scale(mean_obs_bray_f) ~ scale(mean_abund_catch_f[,2]))
summary(mod2_obs_bray_f)

### Now use beta deviations 
# Raup-Crick
# Let's include env heterogeneity and spatial extent

pairs(cbind(scale(mean_abund_catch[,2]), 
            scale(env_het), scale(spa_ext)))

cor(cbind(scale(mean_abund_catch[,2]), 
          scale(env_het), scale(spa_ext)))

mod1_null_raup <- lm(scale(mean_null_raup) ~ scale(mean_abund_catch[,2]) + 
                       scale(env_het) + scale(spa_ext))
summary(mod1_null_raup)

mod1_null_raup_f <- lm(scale(mean_null_raup_f) ~ scale(mean_abund_catch_f[,2]) + 
                       scale(env_het_f) + scale(spa_ext_f))
summary(mod1_null_raup_f)


# Bray-Curtis beta deviations 
mod1_null_bray <- lm(scale(mean_null_bray) ~ scale(mean_abund_catch[,2]) + 
                       scale(env_het) + scale(spa_ext))
summary(mod1_null_bray) 

mod1_null_bray_f <- lm(scale(mean_null_bray_f) ~ scale(mean_abund_catch_f[,2]) + 
                       scale(env_het_f) + scale(spa_ext_f))
summary(mod1_null_bray_f)


### Plots ======================================================================

df_beta_plots <- data.frame(
  bray= c(mean_obs_bray, mean_obs_bray_f),
  sor= c(mean_obs_sor, mean_obs_sor_f),
  raup= c(mean_null_raup, mean_null_raup_f),
  bray.null= c(mean_null_bray, mean_null_bray_f),
  size= c(mean_abund_catch[,2], mean_abund_catch_f[,2]),
  region= rep(c("Tropical", "Boreal"), each=20)
)


## Plots for linear models

plot5 <- ggplot(data = subset(df_beta_plots, region=="Tropical"), 
                aes(x = size, y = raup))  +
  ggtitle("a") +
  geom_point(size=3, color="#D55E00") +
  stat_smooth(method = "lm", formula = y ~ x, se=F, size = 1, color="#D55E00") +
  scale_y_continuous(limits = c(-1, 0.30), 
                     breaks = c(-0.90,  -0.60, -0.30, 0.00, 0.30),
                     labels=c("-0.90",  "-0.60", "-0.30", "0.00", "0.30")) +
  ylab("Beta deviation (Raup-Crick)") +
  scale_x_continuous(name = "") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot6 <- ggplot(data = subset(df_beta_plots, region=="Boreal"), 
                aes(x = size, y = raup)) +
  ggtitle("b") +
  geom_point(size=3, color="#0072B2") +
  scale_y_continuous(name = "", limits = c(-1, 0.30), 
                     breaks = c(-0.90,  -0.60, -0.30, 0.00, 0.30),
                     labels=c("-0.90",  "-0.60", "-0.30", "0.00", "0.30")) +
  scale_x_continuous(name = "") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot7 <- ggplot(data = subset(df_beta_plots, region=="Tropical"), 
                aes(x = size, y = bray.null)) +
  ggtitle("c") +
  scale_y_continuous(limits = c(0,105), 
                     breaks = c(0, 25, 50, 75, 100),
                     labels=c("0.00",  "25.0", "50.0", "75.0", "100.0")) +
  geom_point(size=3, color="#D55E00") +
  stat_smooth(method = "lm", formula = y ~ x, se=F, size = 1, color="#D55E00") +
  ylab("Beta deviation (Bray-Curtis)") +
  xlab("Community size (ind./watershed)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot8 <- ggplot(data = subset(df_beta_plots, region=="Boreal"), 
                aes(x = size, y = bray.null))  +
  ggtitle("d") +
  scale_y_continuous(name = "", limits = c(0,105), 
                     breaks = c(0, 25, 50, 75, 100), 
                     labels=c("0.00",  "25.0", "50.0", "75.0", "100.0")) +
  geom_point(size=3, color="#0072B2") +
  stat_smooth(method = "lm", formula = y ~ x, se=F, size = 1, color="#0072B2") +
  xlab("Community size (ind./watershed)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

grid.newpage()

ggsave(filename = "scats_all.eps",
       plot = grid.arrange(
         arrangeGrob(plot5, plot6, plot7, plot8)
       ), width = 9, height = 9, dpi = 600)

dev.off()


### Plot beta diversity results

plot9 <- ggplot(data = subset(df_beta_plots, region=="Tropical"), 
                aes(x = size, y = sor))  +
  ggtitle("a") +
  geom_point(size=3, color="#D55E00") +
  stat_smooth(method = "lm", formula = y ~ x, se=F, size = 1, color="#D55E00") +
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(name = "") +
  ylab("Beta diversity (Sorensen)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot10 <- ggplot(data = subset(df_beta_plots, region=="Boreal"), 
                 aes(x = size, y = sor)) +
  ggtitle("b") +
  geom_point(size=3, color="#0072B2") +
  scale_y_continuous(name = "", limits = c(0,1)) +
  scale_x_continuous(name = "") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot11 <- ggplot(data = subset(df_beta_plots, region=="Tropical"), 
                 aes(x = size, y = bray)) +
  ggtitle("c") +
  geom_point(size=3, color="#D55E00") +
  scale_y_continuous(limits = c(0,1))+
  stat_smooth(method = "lm", formula = y ~ x, se=F, size = 1, color="#D55E00") +
  ylab("Beta diversity (Bray-Curtis)") +
  xlab("Community size (ind./watershed)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

plot12 <- ggplot(data = subset(df_beta_plots, region=="Boreal"), 
                 aes(x = size, y = bray)) +
  ggtitle("d") + geom_point() + 
  geom_point(size=3, color="#0072B2") +
  scale_y_continuous(name = "", limits = c(0,1))+
  xlab("Community size (ind./watershed)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

grid.newpage()

ggsave(filename = "scats_betadiv_all.eps",
       plot = grid.arrange(
         arrangeGrob(plot9, plot10, plot11, plot12)
       ), width = 9, height = 9, dpi = 600)

dev.off()