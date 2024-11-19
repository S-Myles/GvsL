# Packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all
library(microbiome)   # Contains CLR data transformation
library(ggthemes)
# Global plot theme setting
theme_set(theme_minimal())

# Applying centered log ratio transformation to the ASV data set
S12_CLR <- microbiome::transform(S12_P_A_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S12_NMDS_Eucl <- phyloseq::ordinate(S12_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
(plot <- plot_ordination(S12_CLR, S12_NMDS_Eucl, color="Station", shape="Season") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(16, 2)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  theme(
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"),
    legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S12_NMDS.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


# Applying centered log ratio transformation to the ASV data set
S16_CLR <- microbiome::transform(S16_P_A_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S16_NMDS_Eucl <- phyloseq::ordinate(S16_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
(plot <- plot_ordination(S16_CLR, S16_NMDS_Eucl, color="Station", shape="Season") +
  geom_point(size=5) +
  scale_shape_manual(values = c(16, 2)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S16_NMDS.png",
  plot = plot, width = 6, height = 6, dpi = 300
)



# Applying centered log ratio transformation to the ASV data set
S18_CLR <- microbiome::transform(S18_P_A_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
S18_NMDS_Eucl <- phyloseq::ordinate(S18_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
(plot <- plot_ordination(S18_CLR, S18_NMDS_Eucl, color="Station", shape="Season") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(16, 2)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S18_NMDS.png",
  plot = plot, width = 6, height = 6, dpi = 300
)



# Applying centered log ratio transformation to the ASV data set
COI_CLR <- microbiome::transform(COI_P_A_filt, "clr")
# Calculating Euclidean distances with a NMDS method for plotting
COI_NMDS_Eucl <- phyloseq::ordinate(COI_CLR, method = "NMDS", distance = "euclidean")
# Plotting 
(plot <- plot_ordination(COI_CLR, COI_NMDS_Eucl, color="Station", shape="Season") + 
  geom_point(size=5) +
  scale_shape_manual(values = c(16, 2)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/COI_NMDS.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
