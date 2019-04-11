################################################################################
## Code to run Kraft et al. 2011 null model for beta-diversity
## Script originally created by: J. A. LaManna
## Revised and updated by: Christopher Catano (chcatano@gmail.com)
##
## Reference:
## Kraft et al. (2011) Science, 333(6050) pp. 1755-1758
## DOI: 10.1126/science.1208584
################################################################################

# DESCRIPTION:
# ------------
# This code simulates a dataset with each row representing an individual of a 
# particular species in a particular site/quadrat.
# The code can run these models for multiple datasets at once.  
# In this simulated example, there are three 25 ha forests called "plots".  
# Each forest "plot" contains 625 20x20 m quadrats.  
# The null model then shuffles the locations of these trees (preserving the 
# overall species-abundance distribution in each forest plot and the total 
# number of individuals in each site/quadrat).  

# NOTE: In this example, the site/quadrat is denoted by a number from 1 - 625.  

# Load 'vegan' package to calculate beta-diversity
library(vegan)


#-------------------------------------------------------------------------------
# 1. Simulate data for null model
#-------------------------------------------------------------------------------

# n.plots <- 3				# number of forest plots
# n.sites <- 625			# number of quadrats per forest plot
# n.spp <- 40			  	# number of species in each forest plot

# data <- data.frame("sp" = rep(1:n.spp, times = (n.sites * n.plots)), 
#                   "quadrat" = rep(1:n.sites, each = (n.plots * n.spp)), 
#                   "plot" = rep(rep(1:n.plots, each = n.spp), times = n.sites))

# Create a column of abundances for each species in each plot/quadrat, randomly 
# assigned from a poisson distribution.
# data$abundance <- c(0)
# data$abundance <- rpois(length(data$abundance), 0.5)

# Make each individual a row, and then create an array where each element is
# a data.frame corresponding to each forest plot (3 elements for 3 plots)
# inds.expanded <- data[rep(seq.int(1, nrow(data)), data$abundance), ]
# inds.expanded$value <- c(1)
# indsbyplot <- split(inds.expanded, inds.expanded$plot)


## START HERE WITH YOUR EXISTING DATA. 
## MAKE SURE DATA IS IN THE SAME FORMAT AS "indsbyplot" (an array where each 
## replicate - whether it be a different site, treatment, etc.- is an element &
## is in the format shown for "inds.expanded")
#-------------------------------------------------------------------------------
# 2. Calculate observed beta-diversity 
#-------------------------------------------------------------------------------

# Function to convert 'indsbyplot' to a species x sample matrix 
sp.by.site.function = function(test, species.col = "species", site.col = "site") 
{
  sim_table <- tapply(rep(1, length(test[ , which(colnames(test) == site.col)])), 
                      list(as.vector(test[ , which(colnames(test) == site.col)]), 
                           as.vector(test[ , which(colnames(test) == species.col)])), 
                      sum)
  
  sim_table[is.na(sim_table)] <- 0
  return(sim_table)
}

# Create a species x site matrix for each forest plot, where rows are quadrats
# and columns of species.
# spbysite <- lapply(indsbyplot, sp.by.site.function, species.col = "sp", site.col = "quadrat")

# Calculate median pair-wise dissimilarity among quadrats within each forest plot. 
# In this example we use Bray-Curtis, but you can use any measure of beta-diversity 
# that is most appropriate given the purpose of your particular analysis.
# beta <- lapply(lapply(spbysite, vegdist, method = "bray"), c)
# medianbeta <- lapply(beta, median)


#-------------------------------------------------------------------------------
# 3. Run null model iterations
#-------------------------------------------------------------------------------

# Function to run the null model randomization
null.assembly <- function(test, species.col = "species", site.col = "site") 
{
  a <- sample(test[ , which(colnames(test) == site.col)], 
              length(test[ , which(colnames(test) == site.col)]))
  
  sim_table <- tapply(rep(1, length(a)), 
                      list(a, as.vector(test[ , which(colnames(test) == species.col)])), 
                      sum)
  
  sim_table[is.na(sim_table)] <- 0
  return(sim_table)
}

# Function to run the null model and calculate the mean and standard deviation
# of the expected (simulated) dissimilarities. 
# "k" is the number of iterations (set to 2000 in this example).
null.model <- function(test, k, species.col = "species", site.col = "site") 
{
  null.assemblies.list <- replicate(k, null.assembly(test, species.col = species.col, site.col = site.col), simplify = F)
  
  sim.bray <- lapply(null.assemblies.list, vegdist, method = "bray")
  sim.bray.mean <- apply(simplify2array(sim.bray), 1, mean)
  sim.bray.sd <- apply(simplify2array(sim.bray), 1, sd)
  
  result <- list(ncol(simplify2array(sim.bray)), sim.bray.mean, sim.bray.sd)
  names(result) = c("k", "sim.bray.mean", "sim.bray.sd")
  return(result)
}

# null.models <- lapply(indsbyplot, null.model, k = 1000, species.col = "sp", site.col = "quadrat")
# bray.indsbyburn <- mapply(list, beta, null.models, SIMPLIFY = F)


#-------------------------------------------------------------------------------
# 4. Calculate effect sizes and standardized effect sizes
#-------------------------------------------------------------------------------

effect.sizes <- function(test) 
{
  effect.size <- test[[1]] - test[[2]]$sim.bray.mean
  std.effect.size <- effect.size / test[[2]]$sim.bray.sd
  
  result <- list(test[[1]], test[[2]]$sim.bray.mean, effect.size, std.effect.size)
  names(result) <- c("obs.bray", "sim.bray.mean", "effect.size", "std.effect.size")
  return(result)
}

# bray.effect.size <- lapply(bray.indsbyburn, effect.sizes)

