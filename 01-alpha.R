library(microbiome)   # Contains alpha diversity metrics calculation function

#########
# 12S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S12_physeq_data)) %>% 
   arrange(desc(sum)))

(S12_rarefied_data <- rarefy_even_depth(S12_physeq_data, sample.size = 5000,
                  rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))


# Estimate alpha metrics and find mean
S12_alpha <- alpha(S12_rarefied_data, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)
S12_alpha <- as_tibble(sample_data(S12_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S12_alpha)
mean(S12_alpha$Chao1)
mean(S12_alpha$Shannon)
mean(S12_alpha$Pielou)

(plot <- ggplot(S12_alpha, aes(x = Size.Fraction, y = Chao1, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S12_alpha$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S12_alpha$Chao1) - 22, mean(S12_alpha$Chao1) + 22) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + 
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S12_Chao1.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S12_alpha, aes(x = Size.Fraction, y = Shannon, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Shannon), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S12_alpha$Shannon), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S12_alpha$Shannon) - 1, mean(S12_alpha$Shannon) + 1) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S12_Shannon.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S12_alpha, aes(x = Size.Fraction, y = Pielou, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Pielou), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S12_alpha$Pielou), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S12_alpha$Pielou) - 0.3, mean(S12_alpha$Pielou) + 0.3) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S12_Pielou.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
###############################################################################




#########
# 16S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S16_physeq_data)) %>% 
   arrange(desc(sum)))

(S16_rarefied_data <- rarefy_even_depth(S16_physeq_data, sample.size = 5000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))

# Estimate alpha metrics and find mean
S16_alpha <- alpha(S16_rarefied_data, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)
S16_alpha <- as_tibble(sample_data(S16_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S16_alpha)
mean(S16_alpha$Chao1)
mean(S16_alpha$Shannon)
mean(S16_alpha$Pielou)


(plot <- ggplot(S16_alpha, aes(x = Size.Fraction, y = Chao1, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S16_alpha$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S16_alpha$Chao1) - 590, mean(S16_alpha$Chao1) + 590) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S16_Chao1.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S16_alpha, aes(x = Size.Fraction, y = Shannon, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Shannon), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S16_alpha$Shannon), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S16_alpha$Shannon) - 1.1, mean(S16_alpha$Shannon) + 1.1) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S16_Shannon.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S16_alpha, aes(x = Size.Fraction, y = Pielou, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Pielou), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S16_alpha$Pielou), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S16_alpha$Pielou) - 0.18, mean(S16_alpha$Pielou) + 0.18) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S16_Pielou.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
###############################################################################





#########
# 18S
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(S18_physeq_data)) %>% 
   arrange(desc(sum)))

(S18_rarefied_data <- rarefy_even_depth(S18_physeq_data, sample.size = 5000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))

# Estimate alpha metrics and find mean
S18_alpha <- alpha(S18_rarefied_data, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)
S18_alpha <- as_tibble(sample_data(S18_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(S18_alpha)
mean(S18_alpha$Chao1)
mean(S18_alpha$Shannon)
mean(S18_alpha$Pielou)


(plot <- ggplot(S18_alpha, aes(x = Size.Fraction, y = Chao1, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S18_alpha$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S18_alpha$Chao1) - 550, mean(S18_alpha$Chao1) + 550) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S18_Chao1.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S18_alpha, aes(x = Size.Fraction, y = Shannon, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Shannon), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S18_alpha$Shannon), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S18_alpha$Shannon) -2.3, mean(S18_alpha$Shannon) + 2.3) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S18_Shannon.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(S18_alpha, aes(x = Size.Fraction, y = Pielou, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Pielou), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(S18_alpha$Pielou), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(S18_alpha$Pielou) - 0.35, mean(S18_alpha$Pielou) + 0.35) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/S18_Pielou.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
###############################################################################



#########
# COI
##########
# Rarefy
(sample_sum_df <- data.frame(sum = sample_sums(COI_physeq_data)) %>% 
   arrange(desc(sum)))

(COI_rarefied_data <- rarefy_even_depth(COI_physeq_data, sample.size = 5000,
                                        rngseed = 199, replace = FALSE, trimOTUs = TRUE, verbose = TRUE))


# Estimate alpha metrics and find mean
COI_alpha <- alpha(COI_rarefied_data, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)
COI_alpha <- as_tibble(sample_data(COI_rarefied_data), rownames = "SampleID") %>% 
  bind_cols(COI_alpha)
mean(COI_alpha$Chao1)
mean(COI_alpha$Shannon)
mean(COI_alpha$Pielou)


(plot <- ggplot(COI_alpha, aes(x = Size.Fraction, y = Chao1, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Chao1), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(COI_alpha$Chao1), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(COI_alpha$Chao1) - 260, mean(COI_alpha$Chao1) + 260) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/COI_Chao1.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(COI_alpha, aes(x = Size.Fraction, y = Shannon, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Shannon), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(COI_alpha$Shannon), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(COI_alpha$Shannon) - 2.2, mean(COI_alpha$Shannon) + 2.2) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/COI_Shannon.png",
  plot = plot, width = 6, height = 6, dpi = 300
)

(plot <- ggplot(COI_alpha, aes(x = Size.Fraction, y = Pielou, color = Season)) +
    geom_violin(aes(x = Size.Fraction, y = Pielou), fill = "gray", alpha = 0.5, color = NA) + # Single violin for all samples per Season
    geom_jitter(size = 5, alpha = 0.7, width = 0.1) +
    geom_hline(yintercept = mean(COI_alpha$Pielou), linetype = "dashed", color = "red", size = 1) + # Add mean line
    ylim(mean(COI_alpha$Pielou) - 0.35, mean(COI_alpha$Pielou) + 0.35) + # Set y-axis limits
    scale_color_manual(values = c("#440154FF", "#21908CFF")) + # Customize colors
    scale_x_discrete(labels = c("Large", "Small")) +
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
  filename = "outputs/alpha/COI_Pielou.png",
  plot = plot, width = 6, height = 6, dpi = 300
)
###############################################################################