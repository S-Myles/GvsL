# Packages
library(microbiome)   # Contains CLR data transformation
library(ggthemes)
library(vegan)
library(compositions) # For CLR transformation

# Global plot theme setting
theme_set(theme_bw())


########
# 12S
########
# -------------- PCA for visualization
########
# CLR Transform
S12_CLR <- microbiome::transform(S12_filt_data, "clr")
otu_clr_t <- t(as.data.frame(otu_table(S12_CLR))) # Transposed for vegan::rda

# Perform PCA using vegan::rda (PCA is applied as the last step in RDA, 
# when there's no metadata, it's only PCA)
S12_PCA <- vegan::rda(otu_clr_t)
# Extract explained variance
explained_variance <- summary(S12_PCA)$cont$importance["Proportion Explained", ] * 100
# Data frame for plotting variance
variance_df <- data.frame(
  PC = paste0("PC", seq_along(explained_variance)),
  Variance = explained_variance) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_along(explained_variance)))) %>% 
  slice(1:20)


# Plot variance explained histogram
(variance_plot <- ggplot(variance_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(
      x = "Principal Components",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12)
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/pca/ONESF-S12_PCA_variance_explained.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

# Plot PCA
# Extract PCA scores for samples
sample_scores <- as.data.frame(scores(S12_PCA, display = "sites")) %>%
  rownames_to_column(var = "SampleID")  # Add sample IDs for merging with metadata
# Extract metadata
metadata <- as.data.frame(sample_data(S12_CLR))
metadata$SampleID <- rownames(metadata)
# Merge PCA scores with metadata
pca_plot_data <- left_join(sample_scores, metadata, by = "SampleID")

# Plotting PCA
(plot <- ggplot(data = pca_plot_data, aes(x = PC1, y = PC2, color = Season, shape = Size.Fraction)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    labs(
      x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
    ) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ))

# Save
ggsave(
  filename = "outputs/beta/pca/ONESF-S12_PCA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


#######
# 12S
# -------------- db-RDA with Aitchison Distance
#######

# Extract CLR-transformed OTU table
otu_clr <- as.data.frame(otu_table(S12_CLR))

# Extract metadata
metadata <- as.data.frame(sample_data(S12_CLR))
# Ensure proper formatting of metadata for RDA
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
metadata <- metadata %>%
  as_tibble() %>%
  mutate(across(c(5, 6, 7, 10, 11), as.factor)) %>%         # Convert selected columns to factors
  column_to_rownames(var = "SampleID") %>%
  select(c(4, 5, 6, 9, 10))                                 # only relevant columns

# Perform db-RDA
rda_result <- rda(t(otu_clr) ~ ., data = metadata)

# Extract RDA scores for samples and biplot arrows
sites_df <- as.data.frame(scores(rda_result, display = "sites")) %>%
  rownames_to_column(var = "SampleID")

arrows_df <- as.data.frame(scores(rda_result, display = "bp")) %>%
  rownames_to_column(var = "labels") %>% 
  mutate(RDA1 = RDA1 * 4, RDA2 = RDA2 * 4)

# Merge metadata with RDA scores
sites_df <- left_join(sites_df, rownames_to_column(metadata, var = "SampleID"), by = "SampleID")

# RDA plot
(plot <- ggplot(data = sites_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes(color = Season, shape = Size.Fraction), size = 4, alpha = 0.8) +  
    # Add segments for metadata arrows
    geom_segment(
      data = arrows_df, 
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text(
      data = arrows_df, 
      aes(x = RDA1, y = RDA2, label = labels),
      size = 5, hjust = 1.2, vjust = 1.2, color = "black"
    ) +                                                  # Vector labels
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" 
    ))

# Save the RDA plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S12_RDA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


# Extract RDA eigenvalues
rda_eigenvalues <- rda_result$CCA$eig  # Eigenvalues for RDA axes
pca_eigenvalues <- rda_result$CA$eig  # Eigenvalues for PCA residual axes
# Combine
eigenvalues <- c(rda_eigenvalues, pca_eigenvalues)
# Calculate variance proportions
variance_proportion <- eigenvalues / sum(eigenvalues) * 100

# Create axis labels
axis_labels <- c(
  paste0("RDA", seq_along(rda_eigenvalues)), 
  paste0("PCA", seq_along(pca_eigenvalues))
)

# Create a data frame for plotting
variance_df <- data.frame(
  Axis = factor(axis_labels, levels = axis_labels),  # Ensure correct order
  Variance = variance_proportion
) %>% 
  slice(1:20)


(variance_plot <- ggplot(variance_df, aes(x = Axis, y = Variance)) +
  geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
  labs(x = "Axes", y = "Variance Explained (%)") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold")
  ))

# Save the plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S12_RDA_PCA_variance_explained_ordered.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

###############################################################################












########
# 16S
########
# -------------- PCA for visualization
########
# CLR Transform
S16_CLR <- microbiome::transform(S16_filt_data, "clr")
otu_clr_t <- t(as.data.frame(otu_table(S16_CLR))) # Transposed for vegan::rda

# Perform PCA using vegan::rda
S16_PCA <- vegan::rda(otu_clr_t)
# Extract explained variance
explained_variance <- summary(S16_PCA)$cont$importance["Proportion Explained", ] * 100
# Data frame for plotting variance
variance_df <- data.frame(
  PC = paste0("PC", seq_along(explained_variance)),
  Variance = explained_variance) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_along(explained_variance)))) %>% 
  slice(1:20)


# Plot variance explained histogram
(variance_plot <- ggplot(variance_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(
      x = "Principal Components",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12)
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/pca/ONESF-S16_PCA_variance_explained.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

# Plot PCA
# Extract PCA scores for samples
sample_scores <- as.data.frame(scores(S16_PCA, display = "sites")) %>%
  rownames_to_column(var = "SampleID")  # Add sample IDs for merging with metadata
# Extract metadata
metadata <- as.data.frame(sample_data(S16_CLR))
metadata$SampleID <- rownames(metadata)
# Merge PCA scores with metadata
pca_plot_data <- left_join(sample_scores, metadata, by = "SampleID")

# Plotting PCA
(plot <- ggplot(data = pca_plot_data, aes(x = PC1, y = PC2, color = Season, shape = Size.Fraction)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    labs(
      x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
    ) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ))

# Save
ggsave(
  filename = "outputs/beta/pca/ONESF-S16_PCA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


#######
# 16S
# -------------- db-RDA with Aitchison Distance
#######

# Extract CLR-transformed OTU table
otu_clr <- as.data.frame(otu_table(S16_CLR))

# Extract metadata
metadata <- as.data.frame(sample_data(S16_CLR))
# Ensure proper formatting of metadata for RDA
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
metadata <- metadata %>%
  as_tibble() %>%
  mutate(across(c(5, 6, 7, 10, 11), as.factor)) %>%         # Convert selected columns to factors
  column_to_rownames(var = "SampleID") %>%
  select(c(4, 5, 6, 9, 10))                     

# Perform db-RDA
rda_result <- rda(t(otu_clr) ~ ., data = metadata)

# Extract RDA scores for samples and biplot arrows
sites_df <- as.data.frame(scores(rda_result, display = "sites")) %>%
  rownames_to_column(var = "SampleID")

arrows_df <- as.data.frame(scores(rda_result, display = "bp")) %>%
  rownames_to_column(var = "labels") %>% 
  mutate(RDA1 = RDA1 * 4, RDA2 = RDA2 * 4)

# Merge metadata with RDA scores
sites_df <- left_join(sites_df, rownames_to_column(metadata, var = "SampleID"), by = "SampleID")

# RDA plot
(plot <- ggplot(data = sites_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes(color = Season, shape = Size.Fraction), size = 4, alpha = 0.8) +  
    # Add segments for metadata arrows
    geom_segment(
      data = arrows_df, 
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text(
      data = arrows_df, 
      aes(x = RDA1, y = RDA2, label = labels),
      size = 5, hjust = 1.2, vjust = 1.2, color = "black"
    ) +                                                  # Vector labels
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" 
    ))

# Save the RDA plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S16_RDA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


# Extract RDA eigenvalues
rda_eigenvalues <- rda_result$CCA$eig  # Eigenvalues for RDA axes
pca_eigenvalues <- rda_result$CA$eig  # Eigenvalues for PCA residual axes
# Combine
eigenvalues <- c(rda_eigenvalues, pca_eigenvalues)
# Calculate variance proportions
variance_proportion <- eigenvalues / sum(eigenvalues) * 100

# Create axis labels
axis_labels <- c(
  paste0("RDA", seq_along(rda_eigenvalues)), 
  paste0("PCA", seq_along(pca_eigenvalues))
)

# Create a data frame for plotting
variance_df <- data.frame(
  Axis = factor(axis_labels, levels = axis_labels),  # Ensure correct order
  Variance = variance_proportion) %>% 
  slice(1:20)


(variance_plot <- ggplot(variance_df, aes(x = Axis, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(x = "Axes", y = "Variance Explained (%)") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16, face = "bold")
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S16_RDA_PCA_variance_explained_ordered.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

###############################################################################








########
# 18S
########
# -------------- PCA for visualization
########
# CLR Transform
S18_CLR <- microbiome::transform(S18_filt_data, "clr")
otu_clr_t <- t(as.data.frame(otu_table(S18_CLR))) # Transposed for vegan::rda

# Perform PCA using vegan::rda
S18_PCA <- vegan::rda(otu_clr_t)
# Extract explained variance
explained_variance <- summary(S18_PCA)$cont$importance["Proportion Explained", ] * 100
# Data frame for plotting variance
variance_df <- data.frame(
  PC = paste0("PC", seq_along(explained_variance)),
  Variance = explained_variance) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_along(explained_variance)))) %>% 
  slice(1:20)


# Plot variance explained histogram
(variance_plot <- ggplot(variance_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(
      x = "Principal Components",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12)
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/pca/ONESF-S18_PCA_variance_explained.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

# Plot PCA
# Extract PCA scores for samples
sample_scores <- as.data.frame(scores(S18_PCA, display = "sites")) %>%
  rownames_to_column(var = "SampleID")  # Add sample IDs for merging with metadata
# Extract metadata
metadata <- as.data.frame(sample_data(S18_CLR))
metadata$SampleID <- rownames(metadata)
# Merge PCA scores with metadata
pca_plot_data <- left_join(sample_scores, metadata, by = "SampleID")

# Plotting PCA
(plot <- ggplot(data = pca_plot_data, aes(x = PC1, y = PC2, color = Season, shape = Size.Fraction)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    labs(
      x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
    ) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ))

# Save
ggsave(
  filename = "outputs/beta/pca/ONESF-S18_PCA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


#######
# 18S
# -------------- db-RDA with Aitchison Distance
#######

# Extract CLR-transformed OTU table
otu_clr <- as.data.frame(otu_table(S18_CLR))

# Extract metadata
metadata <- as.data.frame(sample_data(S18_CLR))
# Ensure proper formatting of metadata for RDA
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
metadata <- metadata %>%
  as_tibble() %>%
  mutate(across(c(5, 6, 7, 10, 11), as.factor)) %>%         # Convert selected columns to factors
  column_to_rownames(var = "SampleID") %>%
  select(4, 5, 6, 9, 10)                                          # Adjust to include only relevant columns

# Perform db-RDA
rda_result <- rda(t(otu_clr) ~ ., data = metadata)

# Extract RDA scores for samples and biplot arrows
sites_df <- as.data.frame(scores(rda_result, display = "sites")) %>%
  rownames_to_column(var = "SampleID")

arrows_df <- as.data.frame(scores(rda_result, display = "bp")) %>%
  rownames_to_column(var = "labels") %>% 
  mutate(RDA1 = RDA1 * 4, RDA2 = RDA2 * 4)

# Merge metadata with RDA scores
sites_df <- left_join(sites_df, rownames_to_column(metadata, var = "SampleID"), by = "SampleID")

# RDA plot
(plot <- ggplot(data = sites_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes(color = Season, shape = Size.Fraction), size = 4, alpha = 0.8) +  
    # Add segments for metadata arrows
    geom_segment(
      data = arrows_df, 
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text(
      data = arrows_df, 
      aes(x = RDA1, y = RDA2, label = labels),
      size = 5, hjust = 1.2, vjust = 1.2, color = "black"
    ) +                                                  # Vector labels
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" 
    ))

# Save the RDA plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S18_RDA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


# Extract RDA eigenvalues
rda_eigenvalues <- rda_result$CCA$eig  # Eigenvalues for RDA axes
pca_eigenvalues <- rda_result$CA$eig  # Eigenvalues for PCA residual axes
# Combine
eigenvalues <- c(rda_eigenvalues, pca_eigenvalues)
# Calculate variance proportions
variance_proportion <- eigenvalues / sum(eigenvalues) * 100

# Create axis labels
axis_labels <- c(
  paste0("RDA", seq_along(rda_eigenvalues)), 
  paste0("PCA", seq_along(pca_eigenvalues))
)

# Create a data frame for plotting
variance_df <- data.frame(
  Axis = factor(axis_labels, levels = axis_labels),  # Ensure correct order
  Variance = variance_proportion) %>% 
  slice(1:20)

(variance_plot <- ggplot(variance_df, aes(x = Axis, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(x = "Axes", y = "Variance Explained (%)") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16, face = "bold")
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/rda/ONESF-S18_RDA_PCA_variance_explained_ordered.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)
###############################################################################






########
# COI
########
# -------------- PCA for visualization
########
# CLR Transform
COI_CLR <- microbiome::transform(COI_filt_data, "clr")
otu_clr_t <- t(as.data.frame(otu_table(COI_CLR))) # Transposed for vegan::rda

# Perform PCA using vegan::rda
COI_PCA <- vegan::rda(otu_clr_t)
# Extract explained variance
explained_variance <- summary(COI_PCA)$cont$importance["Proportion Explained", ] * 100
# Data frame for plotting variance
variance_df <- data.frame(
  PC = paste0("PC", seq_along(explained_variance)),
  Variance = explained_variance) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_along(explained_variance)))) %>% 
  slice(1:20)


# Plot variance explained histogram
(variance_plot <- ggplot(variance_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(
      x = "Principal Components",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12)
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/pca/ONESF-COI_PCA_variance_explained.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)

# Plot PCA
# Extract PCA scores for samples
sample_scores <- as.data.frame(scores(COI_PCA, display = "sites")) %>%
  rownames_to_column(var = "SampleID")  # Add sample IDs for merging with metadata
# Extract metadata
metadata <- as.data.frame(sample_data(COI_CLR))
metadata$SampleID <- rownames(metadata)
# Merge PCA scores with metadata
pca_plot_data <- left_join(sample_scores, metadata, by = "SampleID")

# Plotting PCA
(plot <- ggplot(data = pca_plot_data, aes(x = PC1, y = PC2, color = Season, shape = Size.Fraction)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    labs(
      x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
      y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
    ) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ))

# Save
ggsave(
  filename = "outputs/beta/pca/ONESF-COI_PCA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)



#######
# COI
# -------------- db-RDA with Aitchison Distance
#######

# Extract CLR-transformed OTU table
otu_clr <- as.data.frame(otu_table(COI_CLR))

# Extract metadata
metadata <- as.data.frame(sample_data(COI_CLR))
# Ensure proper formatting of metadata for RDA
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
metadata <- metadata %>%
  as_tibble() %>%
  mutate(across(c(5, 6, 7, 10, 11), as.factor)) %>%         # Convert selected columns to factors
  column_to_rownames(var = "SampleID") %>%
  select(4, 5, 6, 9, 10)                                          # Adjust to include only relevant columns

# Perform db-RDA
rda_result <- rda(t(otu_clr) ~ ., data = metadata)

# Extract RDA scores for samples and biplot arrows
sites_df <- as.data.frame(scores(rda_result, display = "sites")) %>%
  rownames_to_column(var = "SampleID")

arrows_df <- as.data.frame(scores(rda_result, display = "bp")) %>%
  rownames_to_column(var = "labels") %>% 
  mutate(RDA1 = RDA1 * 4, RDA2 = RDA2 * 4)

# Merge metadata with RDA scores
sites_df <- left_join(sites_df, rownames_to_column(metadata, var = "SampleID"), by = "SampleID")

# RDA plot
(plot <- ggplot(data = sites_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes(color = Season, shape = Size.Fraction), size = 4, alpha = 0.8) +  
    # Add segments for metadata arrows
    geom_segment(
      data = arrows_df, 
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    geom_text(
      data = arrows_df, 
      aes(x = RDA1, y = RDA2, label = labels),
      size = 5, hjust = 1.2, vjust = 1.2, color = "black"
    ) +                                                  # Vector labels
    scale_shape_manual(values = c(16, 9)) +
    scale_color_manual(values = c("#440154FF", "#21908CFF")) +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none" 
    ))

# Save the RDA plot
ggsave(
  filename = "outputs/beta/rda/ONESF-COI_RDA.png",
  plot = plot, width = 6, height = 6, dpi = 300
)


# Extract RDA eigenvalues
rda_eigenvalues <- rda_result$CCA$eig  # Eigenvalues for RDA axes
pca_eigenvalues <- rda_result$CA$eig  # Eigenvalues for PCA residual axes
# Combine
eigenvalues <- c(rda_eigenvalues, pca_eigenvalues)
# Calculate variance proportions
variance_proportion <- eigenvalues / sum(eigenvalues) * 100

# Create axis labels
axis_labels <- c(
  paste0("RDA", seq_along(rda_eigenvalues)), 
  paste0("PCA", seq_along(pca_eigenvalues))
)

# Create a data frame for plotting
variance_df <- data.frame(
  Axis = factor(axis_labels, levels = axis_labels),  # Ensure correct order
  Variance = variance_proportion) %>% 
  slice(1:20)


(variance_plot <- ggplot(variance_df, aes(x = Axis, y = Variance)) +
    geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.7) +
    labs(x = "Axes", y = "Variance Explained (%)") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16, face = "bold")
    ))

# Save the plot
ggsave(
  filename = "outputs/beta/rda/ONESF-COI_RDA_PCA_variance_explained_ordered.png",
  plot = variance_plot, width = 6, height = 3.5, dpi = 300
)
