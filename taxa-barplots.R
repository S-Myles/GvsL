########
# Global Settings
########
# Libraries
library(phyloseq)
library(tidyverse)
# Global plot theme setting
theme_set(theme_bw())
# Define custom labels for facets
custom_labeller <- labeller(
  Depth = c("1" = "1 m", "20" = "20 m", "250" = "250 m"),  # Replace keys with your Depth values
  Station = c("GULD_04" = "The Gully MPA", "LL_07" = "Louisburg Line S7")  # Replace keys with your Station values
)


########
# 12S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S12_P_A_filt_agg <- S12_P_A_filt %>%
  tax_glom(taxrank = "Species") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S12_P_A_filt_agg) %>% 
  mutate(Date = factor(Date, levels = c("2014_S", "2014_F",
                                        "2016_F", "2017_S",
                                        "2017_F", "2018_S")))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Depth ~ Station, labeller = custom_labeller) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      #legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/Fig-barplot/S12-sp-ra_faceted_barplot.png",
  plot = plot,
  width = 8, height = 6,
  dpi = 300
)


########
# 16S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S16_P_A_filt_agg <- S16_P_A_filt %>%
  tax_glom(taxrank = "Phylum") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S16_P_A_filt_agg) %>% 
  mutate(Date = factor(Date, levels = c("2014_S", "2014_F",
                                        "2016_F", "2017_S",
                                        "2017_F", "2018_S")))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Depth ~ Station, labeller = custom_labeller) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      #legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/Fig-barplot/S16-phylum-ra_faceted_barplot.png",
  plot = plot,
  width = 8, height = 6,
  dpi = 300
)



########
# 18S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S18_P_A_filt_agg <- S18_P_A_filt %>%
  tax_glom(taxrank = "Phylum") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S18_P_A_filt_agg) %>% 
  mutate(Date = factor(Date, levels = c("2014_S", "2014_F",
                                        "2016_F", "2017_S",
                                        "2017_F", "2018_S")))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Depth ~ Station, labeller = custom_labeller) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      #legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/Fig-barplot/S18-Phylum-ra_faceted_barplot.png",
  plot = plot,
  width = 8, height = 6,
  dpi = 300
)





########
# COI
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
COI_P_A_filt_agg <- COI_P_A_filt %>%
  tax_glom(taxrank = "Phylum") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(COI_P_A_filt_agg) %>% 
  mutate(Date = factor(Date, levels = c("2014_S", "2014_F",
                                        "2016_F", "2017_S",
                                        "2017_F", "2018_S")))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Depth ~ Station, labeller = custom_labeller) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      #legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/Fig-barplot/COI-Phylum-ra_faceted_barplot.png",
  plot = plot,
  width = 8, height = 6,
  dpi = 300
)
