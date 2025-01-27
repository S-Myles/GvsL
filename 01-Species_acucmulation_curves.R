library(vegan)
library(ggplot2)
library(viridis)


#(merged_data <- merge_phyloseq(S12_physeq_data, S18_physeq_data, COI_physeq_data))

# Random sampling sp accumulation analysis from Vegan
S12_accum_curve <- otu_table(S12_taxfilt_data) %>% 
  t() %>% 
  specaccum(method = "random")
S16_accum_curve <- otu_table(S16_taxfilt_data) %>% 
  t() %>% 
  specaccum(method = "random")
S18_accum_curve <- otu_table(S18_taxfilt_data) %>% 
  t() %>% 
  specaccum(method = "random")
COI_accum_curve <- otu_table(COI_taxfilt_data) %>% 
  t() %>% 
  specaccum(method = "random")
#merged_accum_curve <- otu_table(merged_data) %>% 
#  t() %>% 
#  specaccum(method = "random")

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
    title = "ASV Accumulation Curves",
    x = "Number of Samples",
    y = "Number of ASVs Detected",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Bold and center
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
  filename = "outputs/asv_accum/asv-accum_tax-filt_ASVs.png",
  plot = plot, width = 5, height = 4, dpi = 300
)







# Bellow shows differences between sample attributes



# Spliting data according to site
gully <- subset_samples(COI_taxfilt_data, Station=="GULD_04")
LL7 <- subset_samples(COI_taxfilt_data, Station=="LL_07")

# 1 curve per site
gully_accum_curve <- otu_table(gully) %>% 
  t() %>% 
  specaccum(method = "random")
ll7_accum_curve <- otu_table(LL7) %>% 
  t() %>% 
  specaccum(method = "random")


site_accum_data <- rbind(
  data.frame(SampleNumber = gully_accum_curve$sites, SpeciesDetected = gully_accum_curve$richness, SD = gully_accum_curve$sd, Group = "Gully"),
  data.frame(SampleNumber = ll7_accum_curve$sites, SpeciesDetected = ll7_accum_curve$richness, SD = ll7_accum_curve$sd, Group = "LL07")
)


# Plot with different line types
(plot <- ggplot(site_accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, linetype = Group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD, fill = Group), alpha = 0.2, color = NA) +
  scale_linetype_manual(values = c("solid", "dotted")) +  # Line types for stations
  scale_color_manual(values = c("#007A4D", "#E65100")) +  # Colorblind-friendly colors
  scale_fill_manual(values = c("#007A4D", "#E65100")) +   # Matching fill colors
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, face = 'bold'),
    axis.text.y = element_text(size = 12, face = 'bold'),
    legend.position = 'none'
  ))

# Save
ggsave(
  filename = "outputs/asv_accum/COI-SITES-accum_tax-filt_ASVs.png",
  plot = plot, width = 5, height = 4, dpi = 300
)




# Spliting data according to seasons
spring <- subset_samples(COI_taxfilt_data, Season=="S")
fall <- subset_samples(COI_taxfilt_data, Season=="F")

# 1 curve per site
spring_accum_curve <- otu_table(spring) %>% 
  t() %>% 
  specaccum(method = "random")
fall_accum_curve <- otu_table(fall) %>% 
  t() %>% 
  specaccum(method = "random")


season_accum_data <- rbind(
  data.frame(SampleNumber = spring_accum_curve$sites, SpeciesDetected = spring_accum_curve$richness, SD = spring_accum_curve$sd, Group = "Spring"),
  data.frame(SampleNumber = fall_accum_curve$sites, SpeciesDetected = fall_accum_curve$richness, SD = fall_accum_curve$sd, Group = "Fall")
)


# Plot with different line types
(plot <- ggplot(season_accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, linetype = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD, fill = Group), alpha = 0.2, color = NA) +
    scale_linetype_manual(values = c("solid", "dotted")) +  # Line types for stations
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +  # Colorblind-friendly colors
    scale_fill_manual(values = c("#440154FF", "#21908CFF")) +   # Matching fill colors
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, face = 'bold'),
      axis.text.y = element_text(size = 12, face = 'bold'),
      legend.position = 'none'
    ))

# Save
ggsave(
  filename = "outputs/asv_accum/COI-SEASOIN_accum_tax-filt_ASVs.png",
  plot = plot, width = 5, height = 4, dpi = 300
)




# Spliting data according to site
small <- subset_samples(COI_taxfilt_data, Size.Fraction=="S")
large <- subset_samples(COI_taxfilt_data, Size.Fraction=="L")

# 1 curve per site
small_accum_curve <- otu_table(small) %>% 
  t() %>% 
  specaccum(method = "random")
large_accum_curve <- otu_table(large) %>% 
  t() %>% 
  specaccum(method = "random")


site_accum_data <- rbind(
  data.frame(SampleNumber = small_accum_curve$sites, SpeciesDetected = small_accum_curve$richness, SD = small_accum_curve$sd, Group = "Small"),
  data.frame(SampleNumber = large_accum_curve$sites, SpeciesDetected = large_accum_curve$richness, SD = large_accum_curve$sd, Group = "Large")
)


# Plot with different line types
(plot <- ggplot(site_accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, linetype = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD, fill = Group), alpha = 0.2, color = NA) +
    scale_linetype_manual(values = c("solid", "dotted")) +  # Line types for stations
    scale_color_manual(values = c("#FDD835", "#C2185B")) +  # Colorblind-friendly colors
    scale_fill_manual(values = c("#FDD835", "#C2185B")) +   # Matching fill colors
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, face = 'bold'),
      axis.text.y = element_text(size = 12, face = 'bold'),
      legend.position = 'none'
    ))

# Save
ggsave(
  filename = "outputs/asv_accum/COI-SF-accum_tax-filt_ASVs.png",
  plot = plot, width = 5, height = 4, dpi = 300
)



# Spliting data according to depth
m1 <- subset_samples(COI_taxfilt_data, Depth==1)
m20 <- subset_samples(COI_taxfilt_data, Depth==20)
m250 <- subset_samples(COI_taxfilt_data, Depth==250)

# 1 curve per site
m1_accum_curve <- otu_table(m1) %>% 
  t() %>% 
  specaccum(method = "random")
m20_accum_curve <- otu_table(m20) %>% 
  t() %>% 
  specaccum(method = "random")
m250_accum_curve <- otu_table(m250) %>% 
  t() %>% 
  specaccum(method = "random")


depth_accum_data <- rbind(
  data.frame(SampleNumber = m1_accum_curve$sites, SpeciesDetected = m1_accum_curve$richness, SD = m1_accum_curve$sd, Group = "1 m"),
  data.frame(SampleNumber = m20_accum_curve$sites, SpeciesDetected = m20_accum_curve$richness, SD = m20_accum_curve$sd, Group = "20 m"),
  data.frame(SampleNumber = m250_accum_curve$sites, SpeciesDetected = m250_accum_curve$richness, SD = m250_accum_curve$sd, Group = "250 m")
  )


# Plot with different line types
(plot <- ggplot(depth_accum_data, aes(x = SampleNumber, y = SpeciesDetected, color = Group, linetype = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = SpeciesDetected - SD, ymax = SpeciesDetected + SD, fill = Group), alpha = 0.2, color = NA) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +  # Line types for stations
    scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#08306B")) +  # Colorblind-friendly colors
    scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#08306B")) +   # Matching fill colors
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, face = 'bold'),
      axis.text.y = element_text(size = 12, face = 'bold'),
      legend.position = 'none'
    ))

# Save
ggsave(
  filename = "outputs/asv_accum/COI-depth-accum_tax-filt_ASVs.png",
  plot = plot, width = 5, height = 4, dpi = 300
)
