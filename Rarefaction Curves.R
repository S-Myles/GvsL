library(vegan)
library(ggplot2)
library(viridis)

#Objects created for ease of access
S12_physeq_data
S12_P_A_filt
S12_filt_glom_g
S16_physeq_data
S16_P_A_filt
S16_filt_glom_g
S18_physeq_data
S18_P_A_filt
S18_filt_glom_g
COI_physeq_data
COI_P_A_filt
COI_filt_glom_g


#(merged_data <- merge_phyloseq(S12_physeq_data, S18_physeq_data, COI_physeq_data))

# Random sampling sp accumulation analysis from Vegan
S12_accum_curve <- otu_table(S12_physeq_data) %>% 
  t() %>% 
  specaccum(method = "random")
S16_accum_curve <- otu_table(S16_physeq_data) %>% 
  t() %>% 
  specaccum(method = "random")
S18_accum_curve <- otu_table(S18_physeq_data) %>% 
  t() %>% 
  specaccum(method = "random")
COI_accum_curve <- otu_table(COI_physeq_data) %>% 
  t() %>% 
  specaccum(method = "random")
merged_accum_curve <- otu_table(merged_data) %>% 
  t() %>% 
  specaccum(method = "random")

# Combine results into 1 dataframe
accum_data <- rbind(
  data.frame(SampleNumber = S12_accum_curve$sites, SpeciesDetected = S12_accum_curve$richness, SD = S12_accum_curve$sd, Group = "12S"),
  data.frame(SampleNumber = S16_accum_curve$sites, SpeciesDetected = S16_accum_curve$richness, SD = S16_accum_curve$sd, Group = "16S"),
  data.frame(SampleNumber = S18_accum_curve$sites, SpeciesDetected = S18_accum_curve$richness, SD = S18_accum_curve$sd, Group = "18S"),
  data.frame(SampleNumber = COI_accum_curve$sites, SpeciesDetected = COI_accum_curve$richness, SD = COI_accum_curve$sd, Group = "COI")
  )

# Plot with colors for each group
(plot <- ggplot(accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, fill = Group)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD), alpha = 0.2, color = NA) + 
  scale_color_viridis_d(option = "plasma", end = 0.9) +  # Colorblind-friendly palette
  scale_fill_viridis_d(option = "plasma", end = 0.9) +   # Same palette for fills
  labs(
    title = "ASV Accumulation Curves (Unfiltered)",
    x = "Number of Samples",
    y = "Number of ASVs Detected",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # Bold and center
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "plain"),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ))

# Save
ggsave(
  filename = "outputs/Fig1/sp-accum_unfiltered_data.png",
  plot = plot, width = 6, height = 6, dpi = 300
)







# Bellow shows that there is basically no difference between sites


S12_physeq_data
S12_P_A_filt
S12_filt_glom_g
S16_physeq_data
S16_P_A_filt
S16_filt_glom_g
S18_physeq_data
S18_P_A_filt
S18_filt_glom_g
COI_physeq_data
COI_P_A_filt
COI_filt_glom_g

# Spliting data according to site
gully <- subset_samples(S12_filt_glom_g, Station=="GULD_04")
LL7 <- subset_samples(S12_filt_glom_g, Station=="LL_07")

# 1 curve per site
gully_accum_curve <- otu_table(gully) %>% 
  t() %>% 
  specaccum(method = "random")
ll7_accum_curve <- otu_table(LL7) %>% 
  t() %>% 
  specaccum(method = "random")


site_accum_data <- rbind(
  data.frame(SampleNumber = gully_accum_curve$sites, SpeciesDetected = gully_accum_curve$richness, SD = gully_accum_curve$sd, Group = "Gully"),
  data.frame(SampleNumber = ll7_accum_curve$sites, SpeciesDetected = ll7_accum_curve$richness, SD = ll7_accum_curve$sd, Group = "Louisburg Line")
)


# Plot with different line types
ggplot(site_accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, linetype = Group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD, fill = Group), alpha = 0.2, color = NA) +
  scale_linetype_manual(values = c("solid", "dotted")) +  # Line types for stations
  scale_color_manual(values = c("#440154FF", "#21908CFF")) +  # Colorblind-friendly colors
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +   # Matching fill colors
  labs(
    title = "Species Accumulation Curves by Station",
    x = "Number of Samples",
    y = "Number of Species Detected",
    color = "Group",
    linetype = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

