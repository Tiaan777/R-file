
1.
# 0 Script Info -----------------------------------------------------------

# * Author: BA Grobler & T Strydom
# * Creation Date: 24 Sep 2021
# * Project: Canopy plant composition and structure of Cape subtropical 
             dune thicket are predicted by levels of fire-exposure
# *          gradient
# * Collaborators: T Strydom, BA Grobler, T Kraaij, RM Cowling
# * Aim: To explore the compositional patterns of dune thicket exposed to
# *      different fire frequencies
# * Data file location: Supplemetary files (Raw data S1) (Strydom et al., 2022; PeerJ)

# 1 Initialize: load libraries and data -----------------------------------

library(tidyverse)
library(here)
library(janitor)
library(vegan)
library(pairwiseAdonis)

raw_data <- read.csv(here("1_Data", "DuneThicket_FireExposure_SppCover_Data.csv"),
                      stringsAsFactors = T, header = T, row.names = 1)

spp_data <- raw_data %>% 
  clean_names(case = "sentence") %>%    # Creates 'clean' variable names
  .[ ,3:37]   # Create dataframe for NMDS analysis

#spp_data <- spp_data %>% 
  #mutate(across(everything(), round))

guilds_data <- read.csv(here("1_Data", "DuneThicket_FireExposure_GuildsCover_Data.csv"),
                     stringsAsFactors = T, header = T, row.names = 1) %>% 
  clean_names() %>% 
  dplyr::select(-fire_exposure)


# 2 NMDS analaysis of species abundance data ------------------------------

# 2.1 Create NMDS scree plot function and check no. of 
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "No. of dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

dist_test <- vegdist(spp_data,  method = "bray")

NMDS.scree(dist_test)

# 2.2 NMDS using 2 axes
nmds_spp_k2 <- metaMDS(spp_data, distance = "bray", k = 2, try = 999)
stressplot(nmds_spp_k2)   # Non-metric R-square fit = 0.974; Linear R-sq fit = 0.874
nmds_spp_k2$stress   # Stress = 0.161: Okay, but not great; might beed more axes

ordiplot(nmds_spp_k2, type = "n")   # Set up blank ordination plot
ordihull(nmds_spp_k2, groups = raw_data$fire_exposure, draw = "polygon", col = "grey90", label = F)  # Polygons to indiacte fire exposure groups
orditorp(nmds_spp_k2, display = "species", col = "red", air = 0.01)   # Plot species in 2D space
orditorp(nmds_spp_k2, display = "sites", cex = 1.25, air = 0.01)   # Plot transects in 2D space

guild_vectors <- envfit(nmds_spp_k2, guilds_data, permu = 999)

plot(guild_vectors, p.max = 0.05)

# 2.3 NMDS using 3 axes
nmds_spp_k3 <- metaMDS(spp_data, distance = "bray", k = 3, try = 999)
stressplot(nmds_spp_k3)   # Non-metric R-square fit = 0.986; Linear R-sq fit = 0.909
nmds_spp_k3$stress   # Stress = 0.119: Decent, better than 2D solution

# Plotting axes 1 and 2
ordiplot(nmds_spp_k3, choices = c(1, 2), type = "n")   # Set up blank ordination plot
ordihull(nmds_spp_k3, choices = c(1, 2), groups = raw_data$fire_exposure, draw = "polygon", col = "grey90", label = F)  # Polygons to indiacte fire exposure groups
orditorp(nmds_spp_k3, choices = c(1, 2), display = "species", col = "red", air = 0.01)   # Plot species
orditorp(nmds_spp_k3, choices = c(1, 2), display = "sites", cex = 1.25, air = 0.01)   # Plot transects

# Plotting axes 1 and 3
ordiplot(nmds_spp_k3, choices = c(1, 3), type = "n")   # Set up blank ordination plot
ordihull(nmds_spp_k3, choices = c(1, 3), groups = raw_data$fire_exposure, draw = "polygon", col = "grey90", label = F)  # Polygons to indiacte fire exposure groups
orditorp(nmds_spp_k3, choices = c(1, 3), display = "species", col = "red", air = 0.01)   # Plot species
orditorp(nmds_spp_k3, choices = c(1, 3), display = "sites", cex = 1.25, air = 0.01)   # Plot transects

# Plotting axes 2 and 3
ordiplot(nmds_spp_k3, choices = c(2, 3), type = "n")   # Set up blank ordination plot
ordihull(nmds_spp_k3, choices = c(2, 3), groups = raw_data$fire_exposure, draw = "polygon", col = "grey90", label = F)  # Polygons to indiacte fire exposure groups
orditorp(nmds_spp_k3, choices = c(2, 3), display = "species", col = "red", air = 0.01)   # Plot species
orditorp(nmds_spp_k3, choices = c(2, 3), display = "sites", cex = 1.25, air = 0.01)   # Plot transects


# 3 Cluster analysis of species presence-absence data ---------------------

# calculate Bray-Curtis distance among samples
spp_dist <- vegdist(spp_data, method = "bray")

# cluster communities using Ward's method
spp_clust_ward <- hclust(spp_dist, method = "ward.D")

# plot Ward's cluster diagram
plot(spp_clust_ward, ylab = "Bray-Curtis dissimilarity")

# cluster communities using UPGMA method
spp_clust_avg <- hclust(spp_dist, method = "average")

# plot UPGMA cluster diagram
plot(spp_clust_avg, ylab = "Bray-Curtis dissimilarity")


# 4 Combo: cluster and NMDS analyses ------------------------------------

#par(mfrow = c(1, 2), pty = "s")

# don't plot anything yet
nmds_fig <- ordiplot(nmds_spp_k2, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(nmds_fig, raw_data$Fire.exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(nmds_fig, raw_data$Fire.exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(nmds_fig, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(nmds_fig, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(nmds_fig, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")


# don't plot anything yet
nmds_fig <- ordiplot(nmds_spp_k2, type = "none")

# add confidence ellipses and bars around habitat types and plot species
ordiellipse(nmds_fig, raw_data$Fire.exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(nmds_fig, raw_data$Fire.exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
orditorp(nmds_fig, "species", cex = 0.8, col = "black", air = 0.01)

# overlay the cluster results we calculated earlier
#ordicluster(nmds_spp_k2, spp_clust, col = "gray")


# 5 Statistical tests of similarity ------------------------------

# ANOSIM
summary(anosim(spp_data, grouping = raw_data$fire_exposure, distance = "bray"))
plot(anosim(spp_data, grouping = raw_data$fire_exposure, distance = "bray"))

# ADONIS
adonis(spp_data ~ raw_data$fire_exposure, method = "bray", permutations = 999)
pairwise.adonis(spp_data, raw_data$fire_exposure)
-------------------------------------------------------------------------------------------------------------------

2. 
# 0 Script Info -----------------------------------------------------------

# * Author: BA Grobler & T Strydom
# * Creation Date: 24 Sep 2021
# * Project: Canopy plant composition and strcuture of Cape subtropical 
             dune thicket are predicted by levels of fire-exposure gradient
# * Collaborators: T Strydom, BA Grobler, T Kraaij, RM Cowling
# * Aim: To explore the compositional patterns of dune thicket exposed to
# *      different fire frequencies
# * Data file location: Supplemetary files (Raw data S1) (Strydom et al., 2022; PeerJ)

# 1 Initialize: load libraries and data -----------------------------------

library(tidyverse)
library(here)
library(janitor)
library(vegan)
library(pairwiseAdonis)

raw_data <- read.csv(here("1_Data", "DuneThicket_FireExposure_GuildsCover_Data.csv"),
                     stringsAsFactors = T, header = T, row.names = 1) %>% 
  clean_names()   # Creates 'clean' variable names

spp_data <- raw_data %>% 
  select(-fire_exposure)   # Create dataframe for NMDS analysis

spp_data <- raw_data %>% 
  select(-fire_exposure) %>% 
  mutate(total_cover = rowSums(.)) %>% 
  mutate(hedge_former_prop = hedge_former / total_cover * 100,
            lateral_spreader_prop = lateral_spreader / total_cover * 100,
            vertical_grower_prop = vertical_grower / total_cover * 100) %>% 
  select(contains("_prop"))

row.names(spp_data) <- row.names(raw_data)

# 2 NMDS analaysis of guild cover-abundance data -----------------

# 2.1 NMDS using 2 axes
nmds_spp_k2 <- metaMDS(spp_data, distance = "bray", k = 2, try = 999)
stressplot(nmds_spp_k2)   # Non-metric R-square fit = 1; Linear R-sq fit = 0.999
nmds_spp_k2$stress   # Stress = 0.014: Pretty good, two dimensions should be fine

ordiplot(nmds_spp_k2, type = "n")   # Set up blank ordination plot
ordihull(nmds_spp_k2, groups = raw_data$fire_exposure, draw = "polygon", col = "grey90", label = F)  # Polygons to indiacte fire exposure groups
orditorp(nmds_spp_k2, display = "species", col = "red", air = 0.01)   # Plot species in 2D space
orditorp(nmds_spp_k2, display = "sites", cex = 1.25, air = 0.01)   # Plot transects in 2D space


# 3 Cluster analysis of species presence-absence data ---------------------

# calculate Bray-Curtis distance among samples
spp_dist <- vegdist(spp_data, method = "bray")

# cluster communities using Ward's method
spp_clust_ward <- hclust(spp_dist, method = "ward.D")

# plot Ward's cluster diagram
plot(spp_clust_ward, hang = -1, ylab = "Bray-Curtis dissimilarity")

# cluster communities using UPGMA method
spp_clust_avg <- hclust(spp_dist, method = "average")

# plot UPGMA cluster diagram
plot(spp_clust_avg, ylab = "Bray-Curtis dissimilarity")


# 4 Combo: cluster and NMDS analyses ------------------------------------

# don't plot anything yet
nmds_fig <- ordiplot(nmds_spp_k2, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(nmds_fig, raw_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(nmds_fig, raw_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat
points(nmds_fig, "sites", pch = 15, col = "#66c2a5", select = raw_data$fire_exposure == "Low")
points(nmds_fig, "sites", pch = 16, col = "#8da0cb", select = raw_data$fire_exposure == "Moderate")
points(nmds_fig, "sites", pch = 17, col = "#fc8d62", select = raw_data$fire_exposure == "High")

# Guild names
orditorp(nmds_fig, "species", cex = 0.8, col = "black", air = 0.1)

#points(nmds_fig, "sites", pch = 19, col = cutree(spp_clust_ward, k = 3))
#ordiellipse(nmds_fig, cutree(spp_clust_ward, k = 3), conf = 0.95, label = TRUE)

# overlay the cluster results we calculated earlier
ordicluster(nmds_spp_k2, spp_clust_ward, col = "gray")


# 5 Statistical tests of similarity ------------------------------

# ANOSIM
summary(anosim(spp_data, grouping = raw_data$fire_exposure, distance = "bray"))
plot(anosim(spp_data, grouping = raw_data$fire_exposure, distance = "bray"))

# ADONIS
adonis(spp_data ~ raw_data$fire_exposure, method = "bray", permutations = 999)
pairwise.adonis(spp_data, raw_data$fire_exposure)


3.
-------------------------------------------------------------------------------------------------------------------
# * Author: BA Grobler & T Strydom
# * Creation Date: 21 Aug 2022
# * Project: Composition of subtropical dune thicket along a fire-exposure
# *          gradient
# * Collaborators: T Strydom, BA Grobler, T Kraaij, RM Cowling
# * Aim: To explore the compositional patterns of dune thicket exposed to
# *      different fire frequencies
# * Data file location: Supplemetary files (Raw data S1) (Strydom et al., 2022; PeerJ)

# 1 Initialize: load libraries and set up data ---------------------------------

library(tidyverse)
library(here)
library(janitor)
library(vegan)

raw_data <- read.csv(here("1_Data", "DuneThicket_FireExposure_SppCover_Data.csv"),   # Read in the raw data
                     stringsAsFactors = T, header = T, row.names = 1)

spp_data <- raw_data %>%   # Create species-cover matrix
  clean_names(case = "sentence") %>%
  .[ ,3:37] %>% 
  sqrt(.) %>%   # square-root transformation
  wisconsin(.)  # Wisconson double-standardisation (same as our NMDS)

guilds_data <- read.csv(here("1_Data", "DuneThicket_FireExposure_GuildsCover_Data.csv"),   # Create guilds-cover matrix
                        stringsAsFactors = T, header = T, row.names = 1) %>% 
  clean_names() %>% 
  dplyr::select(-fire_exposure) %>% 
  sqrt(.) %>%   # square-root transformation
  wisconsin(.)  # Wisconson double-standardisation (same as our NMDS)

env_data <- read.csv(here("1_Data", "DuneThicket_Environment.csv"),   # Create environmental matrix
                     stringsAsFactors = T, header = T, row.names = 1)


# 2 Analysis based on species --------------------------------------------------

## db-RDA with all env. factors
dt_rda <- dbrda(spp_data ~ fire_exposure + slope + aspect, data = env_data, dist = "bray")
summary(dt_rda)
plot(dt_rda)

# don't plot anything yet
fig1_rda <- ordiplot(dt_rda, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig1_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig1_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig1_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig1_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig1_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

# fit env. factors and plot centroids (of factor levels)
dt_env_fit <- envfit(dt_rda ~ fire_exposure + slope + aspect, data = env_data, perm = 999, display = "lc")
plot(dt_env_fit, cex = 0.6, col = "black")

# Testing significance of environmental variables
m0 <- dbrda(spp_data ~ 1, data = env_data, dist = "bray")   # Model with no constraints
m1 <- dbrda(spp_data ~ ., data = env_data, dist = "bray")   # Model with all terms

# Set seed to allow replicable permutational results
set.seed(20220824)

# Significance test and marginal effects
anova(m1)
anova(m1, by = "margin", perm = 999)

## Automatic selection of variables by permutation P-values
dt_rda_2 <- ordistep(m0, scope = formula(m1), permutations = 999)
dt_rda_2$anova
summary(dt_rda_2)

## Permutation test for all variables
anova(dt_rda_2)

# Permutation test of "type III" effects - significance when a term is added to
# the model after all other terms
anova(dt_rda_2, by = "margin", perm = 999)
# Terms added sequentially (first to last)
anova(dt_rda_2, by = "term", perm = 999)

# Check variance inflation factor
vif.cca(m1)
vif.cca(dt_rda_2)

## Unconstrained RDA
# don't plot anything yet
figx_rda <- ordiplot(m0, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(figx_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(figx_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(figx_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(figx_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(figx_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")


## RDA using only significant terms
# don't plot anything yet
fig2_rda <- ordiplot(dt_rda_2, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig2_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig2_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig2_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig2_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig2_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

# fit env. factors and plot centroids (of factor levels)
dt_env_fit2 <- envfit(dt_rda_2 ~ fire_exposure + aspect, data = env_data, perm = 999, display = "lc")
plot(dt_env_fit2, cex = 0.6, col = "black")


## Partial db-RDA (aspect partialled out)
dt_rda3 <- dbrda(spp_data ~ fire_exposure + Condition(aspect), data = env_data)
summary(dt_rda3)
anova(dt_rda3, by = "term", perm = 999)
anova(dt_rda3, by = "margin")

# don't plot anything yet
fig3_rda <- ordiplot(dt_rda3, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig3_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig3_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig3_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig3_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig3_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

# fit env. factors and plot centroids (of factor levels)
dt_env_fit3 <- envfit(dt_rda3 ~ fire_exposure, data = env_data, perm = 999, display = "lc")
plot(dt_env_fit3, cex = 0.6, col = "black")


## Statistical tests of community similarity

dist_test <- vegdist(spp_data, method = "bray")   # Distance matrix for tests

# Multivariate homogeneity of groups dispersions (variances)
dt_disp <- betadisper(dist_test, group = raw_data$Fire.exposure)
anova(dt_disp)
plot(dt_disp)

# ANOSIM
summary(anosim(spp_data, grouping = raw_data$Fire.exposure, distance = "bray"))
plot(anosim(spp_data, grouping = raw_data$Fire.exposure, distance = "bray"))

# ADONIS
adonis2(spp_data ~ raw_data$Fire.exposure, method = "bray", permutations = 999)


# 3 Analysis based on guilds ---------------------------------------------------

# Testing significance of environmental variables
m0 <- dbrda(guilds_data ~ 1, data = env_data, dist = "bray")   # Model with no constraints
m1 <- dbrda(guilds_data ~ ., data = env_data, dist = "bray")   # Model with all terms

set.seed(20220824)

anova(m1)   # Check significance of full model
anova(m1, by = "margin", perm = 999)   # Check significance of factor marginal effects

## Automatic selection of variables by permutation P-values
m2 <- ordistep(m0, scope = formula(m1), permutations = 999)
m2$anova
summary(m2)

## Permutation test for all variables
anova(m2)

# Permutation test of "type III" effects - significance when a term is added to
# the model after all other terms
anova(m2, by = "margin", perm = 999)
# Terms added sequentially (first to last)
anova(m2, by = "term", perm = 999)

vif.cca(m1)
vif.cca(m2)

## Unconstrained RDA
# don't plot anything yet
fig0_rda <- ordiplot(m0, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig0_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig0_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig0_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig0_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig0_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

## Constrained RDA - all factors
# don't plot anything yet
fig1_rda <- ordiplot(m1, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig1_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig1_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig1_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig1_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig1_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

# fit env. factors and plot centroids (of factor levels)
dt_env_fit_m1 <- envfit(m1 ~ ., data = env_data, perm = 999, display = "lc")
plot(dt_env_fit_m1, cex = 0.6, col = "black")

## Constrained RDA - significant factors
# don't plot anything yet
fig2_rda <- ordiplot(m2, type = "none")

# add confidence ellipses and bars around habitat types
ordiellipse(fig2_rda, env_data$fire_exposure, conf = 0.95, label = F, lty = "dashed", lwd = 3, col = c("#fc8d62", "#66c2a5", "#8da0cb"))
ordibar(fig2_rda, env_data$fire_exposure, conf = 0.95, label = T, cex = 0.8, col = c("#fc8d62", "#66c2a5", "#8da0cb"))

# plot just the samples, colour by habitat, pch=19 means plot a circle
points(fig2_rda, "sites", pch = 15, col = "#66c2a5", select = raw_data$Fire.exposure == "Low")
points(fig2_rda, "sites", pch = 16, col = "#8da0cb", select = raw_data$Fire.exposure == "Moderate")
points(fig2_rda, "sites", pch = 17, col = "#fc8d62", select = raw_data$Fire.exposure == "High")

# fit env. factors and plot centroids (of factor levels)
dt_env_fit_m2 <- envfit(m2 ~ ., data = env_data, perm = 999, display = "lc")
plot(dt_env_fit_m2, cex = 0.6, col = "black")
