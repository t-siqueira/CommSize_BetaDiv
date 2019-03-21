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
#===============================================================================

## Required packages
library(vegan)
library(reshape2)
library(adespatial)

library(ggplot2)
library(grid)
library(ggpubr)
library(gridExtra)


## Help functions
# This contains the following functions:
# sp.by.site.function; null.model; effect.sizes
source("helper_functions.R") 

### Data preparation ===========================================================
## Load your data
# Data must be a site (rows) by species (columns) matrix
# One of the columns may be a grouping variable (e.g., watersheds, threatments)

bioraw <- read.csv("bio_20_bra.csv", T)
#repeat the same analyses with finn data
bioraw.f <- read.csv("bio_20_fin.csv", T) 

str(bioraw)

# Check local community abudances
sort(rowSums(bioraw[,-c(1,2)]))

# I'm removing sites with less than 10 individuals
good.sites <- rowSums(bioraw[,-c(1,2)])>10
table(good.sites)
bioraw <- bioraw[good.sites,]
dim(bioraw)

# Finn
good.sites.f <- rowSums(bioraw.f[,-c(1,2)])>10
table(good.sites.f)
bioraw.f <- bioraw.f[good.sites.f,]
dim(bioraw.f)

## Create a list in which each element is a compositional matrix (sitexspecies)
# representing your metacommunity replicate
# Here, each of the 20 watersheds is a replicate 
biolist <- split(bioraw[,-c(1,2)], bioraw$catchment)
str(biolist)

# Check if you need to remove any replicate due to insufficient samples
lapply(biolist, FUN=nrow)

biolist.f <- split(bioraw.f[,-c(1,2)], bioraw.f$catchment)
str(biolist.f)
lapply(biolist.f, FUN=nrow)

## Keep only the taxa that occur within your grouping variable (if you have)
biolist2 <- lapply(biolist, FUN= function(x) x[,colSums(x)>0])
lapply(biolist2, FUN=nrow)

biolist2.f <- lapply(biolist.f, FUN= function(x) x[,colSums(x)>0])
lapply(biolist2.f, FUN=nrow)


## Estimate community size =====================================================
# Here, because sampling design/effort was standardized
# I am using the mean number of individuals per watershed
local.abund <- rowSums(bioraw[,-c(1,2)])
mean.abund.catch <- aggregate(local.abund, by=list(bioraw$catchment), FUN=mean) 

local.abund.f <- rowSums(bioraw.f[,-c(1,2)])
mean.abund.catch.f <- aggregate(local.abund.f, by=list(bioraw.f$catchment), 
                                FUN=mean)

### Observed beta diversity ====================================================
## Calculate obs. b-diversity for each of the levels of your grouping variable
beta.obs.bray <- lapply(lapply(biolist2, vegdist, method = "bray"), c)
mean.obs.bray <- unlist(lapply(beta.obs.bray, mean))

# Sorensen
beta.obs.sor <- lapply(lapply(biolist2, vegdist, method = "bray", binary=T), c)
mean.obs.sor <- unlist(lapply(beta.obs.sor, mean))

# Finn
beta.obs.bray.f <- lapply(lapply(biolist2.f, vegdist, method = "bray"), c)
mean.obs.bray.f <- unlist(lapply(beta.obs.bray.f, mean))

beta.obs.sor.f <- lapply(lapply(biolist2.f, vegdist, method = "bray", 
                                binary=T), c)
mean.obs.sor.f <- unlist(lapply(beta.obs.sor.f, mean))


### Beta diversity based on null models=========================================
## 1) Raup-Crick
# First I will the Raup–Crick null model to determine how much 
# variation in composition was due only to sampling effects on 
# a-diversity (richness) that are neutral with respect to species identity 
# (Chase et al. 2011)

# Because I am using all occuring species to define the species pool,
# I will calculte the dissimilarity between all pairs using the bioraw data
beta.raup <- raupcrick (decostand(bioraw[,-c(1,2)], "pa"), 
                        nsimul = 9999, chase = F)

# Scale to vary between -1 and 1
beta.raup <- (beta.raup-.5)*2

beta.raup.list <- list()

for (j in 1:nlevels(as.factor(bioraw$catchment))){
  nstream <- table(bioraw$catchment)[j]
  soma <- sum(table(bioraw$catchment)[1:j])-nstream
  beta.raup.list[[j]] <- as.matrix(beta.raup)[(soma+1):(soma+nstream), 
                                              (soma+1):(soma+nstream)]
}

mean.null.raup <- unlist(lapply(beta.raup.list, 
                                FUN= function(x) mean(as.dist(x))))


# Finn
beta.raup.f <- raupcrick (decostand(bioraw.f[,-c(1,2)], "pa"), 
                          nsimul = 9999, chase = F)

beta.raup.f <- (beta.raup.f-.5)*2

beta.raup.list.f <- list()

for (j in 1:nlevels(as.factor(bioraw.f$catchment))){
  nstream <- table(bioraw.f$catchment)[j]
  soma <- sum(table(bioraw.f$catchment)[1:j])-nstream
  beta.raup.list.f[[j]] <- as.matrix(beta.raup.f)[(soma+1):(soma+nstream), 
                                                  (soma+1):(soma+nstream)]
}

mean.null.raup.f <- unlist(lapply(beta.raup.list.f, 
                                  FUN= function(x) mean(as.dist(x))))


## 2) Beta deviation from Bray-Curtis ==========================================
# This one needs a specific data format and a number of functions
# MAKE SURE DATA IS IN THE SAME FORMAT AS "indsbyplot" (an array where each 
# replicate - whether it be a different site, treatment, etc.- is an element &
# is in the format shown for "inds.expanded")

# > str(indsbyplot)
# List of 3
# $ 1:'data.frame':	12323 obs. of  4 variables:
# ..$ sp       : int [1:12323] 1 3 6 7 9 10 10 11 13 15 ...
# ..$ quadrat  : int [1:12323] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ plot     : int [1:12323] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ abundance: int [1:12323] 1 1 1 1 1 2 2 1 1 1 ...

## My data
bio.longformat <- read.csv("bio_20_bra_mod.csv", header=T)
bio.longformat.f <- read.csv("bio_20_fin_mod.csv", header=T)

str(bio.longformat)
bio.longformat <- bio.longformat[good.sites,]
bio.longformat <- bio.longformat[,colSums(bio.longformat[,-c(1,2)])>0]

table(bio.longformat$stream)
table(bio.longformat$catchment)

# Finn
bio.longformat.f <- bio.longformat.f[good.sites.f,]
bio.longformat.f <- bio.longformat.f[,colSums(bio.longformat.f[,-c(1,2)])>0]

## Because I am using all sites as the regional pool
# I will treat watersheds as a grouping variable with one level
bio.longformat$stream <- bioraw$stream

bio2 <- melt(bio.longformat, id.vars=c("stream", "catchment"), 
             variable.name = "sp", value.name = "abundance")

bio3 <- data.frame("sp"=as.integer(bio2$sp), "quadrat"=bio2$stream, 
                   "plot"=bio2$catchment, "abundance"=bio2$abundance)
bio3[,"plot"] <- 1 # because my whole data set is the species pool

inds.expanded <- bio3[rep(seq.int(1, nrow(bio3)), bio3$abundance), ]
indsbyplot <- split(inds.expanded, inds.expanded$plot)

# Finn
bio.longformat.f$stream <- bioraw.f$stream

bio2.f <- melt(bio.longformat.f, id.vars=c("stream", "catchment"), 
               variable.name = "sp", value.name = "abundance")

bio3.f <- data.frame("sp"=as.integer(bio2.f$sp), "quadrat"=bio2.f$stream, 
                     "plot"=bio2.f$catchment, "abundance"=bio2.f$abundance)
bio3.f[,"plot"] <- 1 # because my whole data set is the species pool

inds.expanded.f <- bio3.f[rep(seq.int(1, nrow(bio3.f)), bio3.f$abundance), ]
indsbyplot.f <- split(inds.expanded.f, inds.expanded.f$plot)

## Create a species x site matrix for each plot (watersheds), where rows are 
# quadrats (streams) and columns are species.
spbysite <- lapply(indsbyplot, sp.by.site.function, species.col = "sp", 
                   site.col = "quadrat")

spbysite.f <- lapply(indsbyplot.f, sp.by.site.function, species.col = "sp", 
                     site.col = "quadrat")

## Calculate mean pair-wise dissimilarity among samples within each replicate. 
# In this example we use Bray-Curtis, but you can use any measure  
# that is most appropriate given the purpose of your particular analysis.
beta.bray.all <- lapply(lapply(spbysite, vegdist, method = "bray"), c)

beta.bray.all.f <- lapply(lapply(spbysite.f, vegdist, method = "bray"), c)

## Run the null model and calculate the mean and standard deviation
# of the expected (simulated) dissimilarities. 
# "k" is the number of iterations 
null.models.bray <- lapply(indsbyplot, null.model, k = 9999, species.col = "sp", 
                           site.col = "quadrat")
bray.all.results <- mapply(list, beta.bray.all, null.models.bray, SIMPLIFY = F)

# Finn
null.models.bray.f <- lapply(indsbyplot.f, null.model, k = 9999, 
                             species.col = "sp", site.col = "quadrat")
bray.all.results.f <- mapply(list, beta.bray.all.f, null.models.bray.f, 
                             SIMPLIFY = F)

## Calculate effect sizes
bray.effect.size <- lapply(bray.all.results, FUN=effect.sizes)
str(bray.effect.size)

bray.effect.size.f <- lapply(bray.all.results.f, FUN=effect.sizes)

## Now reorganize these values per watershed ===================================
# I only load mefa here because its function melt is also part of reshape2
library(mefa)

bray.dist.eff <- vec2dist (bray.effect.size$`1`$std.effect.size, 
                           nrow(bio.longformat))

beta.bray.list <- list()

for (j in 1:nlevels(as.factor(bioraw$catchment))){
  nstream <- table(bioraw$catchment)[j]
  soma <- sum(table(bioraw$catchment)[1:j])-nstream
  beta.bray.list[[j]] <- as.matrix(bray.dist.eff)[(soma+1):(soma+nstream), 
                                                  (soma+1):(soma+nstream)]
}

mean.null.bray <- unlist(lapply(beta.bray.list, 
                                FUN= function(x) mean(as.dist(x))))

# Finn
bray.dist.eff.f <- vec2dist (bray.effect.size.f$`1`$std.effect.size, 
                             nrow(bio.longformat.f))

beta.bray.list.f <- list()

for (j in 1:nlevels(as.factor(bioraw.f$catchment))){
  nstream <- table(bioraw.f$catchment)[j]
  soma <- sum(table(bioraw.f$catchment)[1:j])-nstream
  beta.bray.list.f[[j]] <- as.matrix(bray.dist.eff.f)[(soma+1):(soma+nstream), 
                                                      (soma+1):(soma+nstream)]
}

mean.null.bray.f <- unlist(lapply(beta.bray.list.f, 
                                  FUN= function(x) mean(as.dist(x))))


## Env and spatial predictors ==================================================
# Load your env data
# In-stream variables
env <- read.csv("env_20_bra.csv", T)
env.f <- read.csv("env_20_fin.csv", T)

env <- env[good.sites,]
env.f <- env.f[good.sites.f,]

# Catchment scale variables
land400 <- read.csv("land_20_bra.csv", T)
land400.f <- read.csv("land_20_fin.csv", T)

land400 <- land400[good.sites,]
land400.f <- land400.f[good.sites.f,]

names(env)
names(land400)

names(env.f)
names(land400.f)

# Exclude correlated variables
all.env <- cbind(env[,-c(18,19)], land400[,-c(1:5)]) 
all.env.f <- cbind(env.f[,-c(5,21)], land400.f[,-c(1:6)]) 

# Standardize your data if needed
all.env.stnd <- decostand(all.env[,-c(1,2)], "stand")
all.env.stnd.f <- decostand(all.env.f[,-c(1,2)], "stand")

# Organize your env data as a list, just like you did with compositional data  
envlist <- split(all.env.stnd, bioraw$catchment)
env.dists.list <- lapply(envlist, vegdist, method = "eucl")

envlist.f <- split(all.env.stnd.f, bioraw.f$catchment)
env.dists.list.f <- lapply(envlist.f, vegdist, method = "eucl")

# Prepare the spatial data
spa <- read.csv("spa_20.csv", T)
spa.f <- read.csv("spa_20_fin.csv", T)

spa <- spa[good.sites,]
spa.f <- spa.f[good.sites.f,]

spalist <- split(spa[,-c(1,2)], bioraw$catchment)
spa.obs.eucl <- lapply(spalist, vegdist, method = "eucl")

spalist.f <- split(spa.f[,-c(1,2)], bioraw.f$catchment)
spa.obs.eucl.f <- lapply(spalist.f, vegdist, method = "eucl")

# Use Marti Anderson's betadisper ==============================================
# to estimate mean dispersion

env.het <- unlist(
  lapply(
    lapply(
      env.dists.list, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# Post-hoc analysis; among watersheds
permutest(
  betadisper(vegdist(all.env.stnd, "euc"), group= bioraw$catchment, type="med"),
  pairwise = F, permutations = 999)

spa.ext <- unlist(
  lapply(
    lapply(
      spa.obs.eucl, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# Finn
env.het.f <- unlist(
  lapply(
    lapply(
      env.dists.list.f, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))

# Post-hoc analysis; among watersheds
permutest(
  betadisper(vegdist(all.env.stnd.f, "euc"), group= bioraw.f$catchment, 
             type="med"), pairwise = F, permutations = 999)

spa.ext.f <- unlist(
  lapply(
    lapply(
      spa.obs.eucl.f, FUN = function(x) 
        betadisper(x, group=factor(rep("a", nrow(as.matrix(x)))))), 
    FUN = function(x) mean(x$distances)))


## OLS models ==================================================================

mod2.obs.sor <- lm(scale(mean.obs.sor) ~ scale(mean.abund.catch[,2]))
summary(mod2.obs.sor)

# Finn
mod2.obs.sor.f <- lm(scale(mean.obs.sor.f) ~ scale(mean.abund.catch.f[,2]))
summary(mod2.obs.sor.f)


# Now with Bray-Curtis dissimilarities
mod2.obs.bray <- lm(scale(mean.obs.bray) ~ scale(mean.abund.catch[,2]))
summary(mod2.obs.bray)

mod2.obs.bray.f <- lm(scale(mean.obs.bray.f) ~ scale(mean.abund.catch.f[,2]))
summary(mod2.obs.bray.f)

### Now use beta deviations 
# Raup-Crick
# Let's include env heterogeneity and spatial extent

pairs(cbind(scale(mean.abund.catch[,2]), 
            scale(env.het), scale(spa.ext)))

cor(cbind(scale(mean.abund.catch[,2]), 
          scale(env.het), scale(spa.ext)))

mod1.null.raup <- lm(scale(mean.null.raup) ~ scale(mean.abund.catch[,2]) + 
                       scale(env.het) + scale(spa.ext))
summary(mod1.null.raup)

mod1.null.raup.f <- lm(scale(mean.null.raup.f) ~ scale(mean.abund.catch.f[,2]) + 
                         scale(env.het.f) + scale(spa.ext.f))
summary(mod1.null.raup.f)

plot(y=mean.null.raup.f, x=env.het.f)

mod2.null.raup.f <- lm(scale(mean.null.raup.f) ~ scale(env.het.f))
summary(mod2.null.raup.f)


# Bray-Curtis beta deviations 
mod1.null.bray <- lm(scale(mean.null.bray) ~ scale(mean.abund.catch[,2]) + 
                       scale(env.het) + scale(spa.ext))
summary(mod1.null.bray) 

mod1.null.bray.f <- lm(scale(mean.null.bray.f) ~ scale(mean.abund.catch.f[,2]) + 
                         scale(env.het.f) + scale(spa.ext.f))
summary(mod1.null.bray.f)


### Plots ======================================================================

df.beta.plots <- data.frame(
  bray= c(mean.obs.bray, mean.obs.bray.f),
  sor= c(mean.obs.sor, mean.obs.sor.f),
  raup= c(mean.null.raup, mean.null.raup.f),
  bray.null= c(mean.null.bray, mean.null.bray.f),
  size= c(mean.abund.catch[,2], mean.abund.catch.f[,2]),
  region= rep(c("Tropical", "Boreal"), each=20)
)


## Plots for linear models

plot5 <- ggplot(data = subset(df.beta.plots, region=="Tropical"), 
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

plot6 <- ggplot(data = subset(df.beta.plots, region=="Boreal"), 
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

plot7 <- ggplot(data = subset(df.beta.plots, region=="Tropical"), 
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

plot8 <- ggplot(data = subset(df.beta.plots, region=="Boreal"), 
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

plot9 <- ggplot(data = subset(df.beta.plots, region=="Tropical"), 
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

plot10 <- ggplot(data = subset(df.beta.plots, region=="Boreal"), 
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

plot11 <- ggplot(data = subset(df.beta.plots, region=="Tropical"), 
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

plot12 <- ggplot(data = subset(df.beta.plots, region=="Boreal"), 
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



### Alternative metrics of community size ======================================

## Stream width as a proxy for communuty size
stream.width <- aggregate(env$width_m, by=list(bioraw$catchment), FUN=mean) 

plot(x=stream.width[,2], y=mean.obs.bray,pch=20, cex=1.5)
plot(x=stream.width[,2], y=mean.obs.sor, pch=20, cex=1.5)

plot(x=stream.width[,2], y=mean.null.bray, pch=20, cex=1.5)
plot(x=stream.width[,2], y=mean.null.raup, pch=20, cex=1.5)

# Estimate mean (and median) population size

median.pop.size.loc <- lapply(biolist2, FUN= function(x) 
  apply(x, 1, function(y) median(y[y>0])))

mean.pop.size.catch <- unlist(lapply(median.pop.size.loc, FUN=mean))

# Plot beta div against mean(median) pop size
plot(x=mean.pop.size.catch, y=mean.obs.bray, pch=20, cex=1.5)
plot(x=mean.pop.size.catch, y=mean.obs.sor, pch=20, cex=1.5)

plot(x=mean.pop.size.catch, y=mean.null.bray, pch=20, cex=1.5)
plot(x=mean.pop.size.catch, y=mean.null.raup, pch=20, cex=1.5)


