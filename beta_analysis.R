# Packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all
library(microbiome)   # Contains CLR data transformation
library(ggthemes)

# Applying centered log ratio transformation to the ASV data set
S12_CLR <- microbiome::transform(S12_physeq_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S12_NMDS_Eucl <- phyloseq::ordinate(S12_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
theme_set(theme_pander())
plot_ordination(S12_CLR, S12_NMDS_Eucl, color="Season", shape="Station") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))


# Applying centered log ratio transformation to the ASV data set
S16_CLR <- microbiome::transform(S16_physeq_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S16_NMDS_Eucl <- phyloseq::ordinate(S16_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
plot_ordination(S16_CLR, S16_NMDS_Eucl, color="Season", shape="Station") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))



# Applying centered log ratio transformation to the ASV data set
S18_CLR <- microbiome::transform(S18_physeq_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S18_NMDS_Eucl <- phyloseq::ordinate(S18_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
plot_ordination(S18_CLR, S18_NMDS_Eucl, color="Season", shape="Station") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))



# Applying centered log ratio transformation to the ASV data set
COI_CLR <- microbiome::transform(COI_physeq_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
COI_NMDS_Eucl <- phyloseq::ordinate(COI_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
plot_ordination(COI_CLR, COI_NMDS_Eucl, color="Season", shape="Station") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))
