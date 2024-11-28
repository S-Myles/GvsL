
#########
# 12S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S12_physeq_data)) %>% 
   arrange(desc(sum)))

(S12_rarefied_data <- rarefy_even_depth(S12_physeq_data, sample.size = min(sample_sums(S12_physeq_data)),
                  rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))


# Estimate Richness and find mean
S12_richness <- estimate_richness(S12_rarefied_data, measures = "Chao1")
S12_richness <- as_tibble(sample_data(S12_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S12_richness)
mean(S12_richness$Chao1)

(plot <- ggplot(S12_richness, aes(x = Season, y = Chao1, color = Station)) +
    geom_violin(aes(x = Season, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S12_richness$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S12_richness$Chao1) - 12, mean(S12_richness$Chao1) + 12) + # Set y-axis limits
    scale_color_manual(values = c("#0072B2", "#D55E00")) + # Customize colors
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S12_seasonal_alpha.png",
  plot = plot, width = 6, height = 6, dpi = 300
)



#########
# 16S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S16_physeq_data)) %>% 
   arrange(desc(sum)))

(S16_rarefied_data <- rarefy_even_depth(S16_physeq_data, sample.size = 5000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))



# Estimate Richness and find mean
S16_richness <- estimate_richness(S16_rarefied_data, measures = "Chao1")
S16_richness <- as_tibble(sample_data(S16_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S16_richness)
mean(S16_richness$Chao1)

(plot <- ggplot(S16_richness, aes(x = Season, y = Chao1, color = Station)) +
    geom_violin(aes(x = Season, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S16_richness$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S16_richness$Chao1) - 700, mean(S16_richness$Chao1) + 700) + # Set y-axis limits
    scale_color_manual(values = c("#0072B2", "#D55E00")) + # Customize colors
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S16_seasonal_alpha.png",
  plot = plot, width = 6, height = 6, dpi = 300
)




#########
# 18S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S18_physeq_data)) %>% 
   arrange(desc(sum)))

(S18_rarefied_data <- rarefy_even_depth(S18_physeq_data, sample.size = 4000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))

# Estimate Richness and find mean
S18_richness <- estimate_richness(S18_rarefied_data, measures = "Chao1")
S18_richness <- as_tibble(sample_data(S18_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S18_richness)
mean(S18_richness$Chao1)

(plot <- ggplot(S18_richness, aes(x = Season, y = Chao1, color = Station)) +
    geom_violin(aes(x = Season, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S18_richness$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S18_richness$Chao1) - 275, mean(S18_richness$Chao1) + 275) + # Set y-axis limits
    scale_color_manual(values = c("#0072B2", "#D55E00")) + # Customize colors
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" # Removes the legend completely
    ))

# Save
ggsave(
  filename = "outputs/Fig2/S18_seasonal_alpha.png",
  plot = plot, width = 6, height = 6, dpi = 300
)




#########
# COI
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(COI_physeq_data)) %>% 
   arrange(desc(sum)))

(COI_rarefied_data <- rarefy_even_depth(COI_physeq_data, sample.size = 4000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))

# Estimate Richness and find mean
COI_richness <- estimate_richness(COI_rarefied_data, measures = "Chao1")
COI_richness <- as_tibble(sample_data(COI_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(COI_richness)
mean(COI_richness$Chao1)

(plot <- ggplot(COI_richness, aes(x = Season, y = Chao1, color = Station)) +
  geom_violin(aes(x = Season, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
  geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
  geom_hline(yintercept = mean(COI_richness$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
  ylim(mean(COI_richness$Chao1) - 300, mean(COI_richness$Chao1) + 300) + # Set y-axis limits
  scale_color_manual(values = c("#0072B2", "#D55E00")) + # Customize colors
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18, face = "bold"),
    legend.position = "none" # Removes the legend completely
  ))


# Save
ggsave(
  filename = "outputs/Fig2/COI_seasonal_alpha.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
