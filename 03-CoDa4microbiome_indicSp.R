# Packages 
library(phyloseq)
library(coda4microbiome)
library(tidyverse)

load("input_objects.RData")
load("results.RData")

# These indicator species analyses are to be done only on 1 size-fraction datasets because
# RDAs showed fundamental impact on the data. 

# Only investigate sample attributes showing a significant association to community organization 
# from rd-RDA anova analysis

# Year should not be a factor to investigate. So only Season, Depth, or Station if significant.

# Depth should be binned to a binary as photic (1 and 20m) or aphotic (250 m)


# So for 
# 12S : none (only year significant with 1 size fraction, station when all samples present)
# 16S : Depth & Season
# 18S : Depth & Season & Station
# COI : Depth & Season 


####
# Data Prep ------------  16S
####
S16_df <- otu_table(S16_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(S16_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
S16_y_season <- metadata[,"Season"] %>% 
  as.factor()
S16_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

S16_df <- column_to_rownames(S16_df, var = "sampleid")



#####
# Cross-sectional (binary) glm penalized net
#####
# Running the model
# Function recognizes the binary outcome so, implements a penalized logistic regression
# First plot generated monitors AUC drop off as feature space is reduced iteratively
# Dashed verticals show highest AUC score and cutoff selection for dimensionality reduction
# This cutoff value is decided based on the lambda argument (default is 1 standard dev)

# 16S Seasonal Model
S16_coda_glmnet_seasonal<-coda_glmnet(x=S16_df,y=S16_y_season, nfolds = 3)

####
# Model evaluation
###
# Extract data from the S16_coda_glmnet_seasonal object
asv_num <- S16_coda_glmnet_seasonal$taxa.num
asv_name <- S16_coda_glmnet_seasonal$taxa.name
coefficients <- S16_coda_glmnet_seasonal$`log-contrast coefficients`

# Combine to a table
S16_seasonal_results <- data.frame(
  ASV_num = asv_num,
  ASV = asv_name,
  Coefficient = coefficients
) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(Coefficient, decreasing = TRUE)]))


(signature_plot <- ggplot(S16_seasonal_results, aes(x = ASV, y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +  # Add reference line at y=0
    labs(x = "Indicator ASV", 
         y = "Log-Contrast Coefficient") +
    scale_y_continuous(breaks = seq(-0.2, 0.6, by = 0.1)) +  # Define the x-axis breaks
    scale_fill_manual(
      values = c("FALSE" = "#440154FF", "TRUE" = "#21908CFF"),  # Assign colors
      labels = c("Fall", "Spring"),  # Legend labels
      name = "Season") +
    theme_minimal() +  # Clean theme
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)) +
    coord_flip()
)

# Save
ggsave(
  filename = "outputs/indicators/S16-coda-seasonal-sp.png",
  plot = signature_plot, width = 8, height = 8, dpi = 300
)

(plot <- S16_coda_glmnet_seasonal$`predictions plot`)



# 16S Depth Model
S16_coda_glmnet_depth<-coda_glmnet(x=S16_df,y=S16_y_depth, nfolds = 3)

####
# Model evaluation
###
# Extract data from the S16_coda_glmnet_seasonal object
asv_num <- S16_coda_glmnet_depth$taxa.num
asv_name <- S16_coda_glmnet_depth$taxa.name
coefficients <- S16_coda_glmnet_depth$`log-contrast coefficients`

# Combine to a table
S16_depth_results <- data.frame(
  ASV_num = asv_num,
  ASV = asv_name,
  Coefficient = coefficients
) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(Coefficient, decreasing = TRUE)]))


(signature_plot <- ggplot(S16_depth_results, aes(x = ASV, y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +  # Add reference line at y=0
    labs(x = "Indicator ASV", 
         y = "Log-Contrast Coefficient") +
    scale_y_continuous(limits = c(-1, 0.5),  # Set the range
                       breaks = seq(-1, 0.5, by = 0.50)) +  # Define the x-axis breaks
    scale_fill_manual(
      values = c("FALSE" = "#A6CEE3", "TRUE" = "#08306B"),  # Assign colors
      labels = c("Photic", "Aphotic"),  # Legend labels
      name = "Sampling Depth") +
    theme_minimal() +  # Clean theme
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)) +
    coord_flip()
)

# Save
ggsave(
  filename = "outputs/indicators/S16-coda-depth-sp.png",
  plot = signature_plot, width = 8, height = 4, dpi = 300
)

(plot <- S16_coda_glmnet_seasonal$`predictions plot`)










####
# Data Prep ------------   18S
####

S18_df <- otu_table(S18_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(S18_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
S18_y_season <- metadata[,"Season"] %>% 
  as.factor()
S18_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

S18_df <- column_to_rownames(S18_df, var = "sampleid")



#####
# Cross-sectional (binary) glm penalized net
#####

# 18S Seasonal Model
S18_coda_glmnet_seasonal<-coda_glmnet(x=S18_df,y=S18_y_season, nfolds = 3)

####
# Model evaluation
###
# Extract data from the S18_coda_glmnet_seasonal object
asv_num <- S18_coda_glmnet_seasonal$taxa.num
asv_name <- S18_coda_glmnet_seasonal$taxa.name
coefficients <- S18_coda_glmnet_seasonal$`log-contrast coefficients`

# Combine to a table
S18_seasonal_results <- data.frame(
  ASV_num = asv_num,
  ASV = asv_name,
  Coefficient = coefficients
) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(Coefficient, decreasing = TRUE)]))


(signature_plot <- ggplot(S18_seasonal_results, aes(x = ASV, y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +  # Add reference line at y=0
    labs(x = "Indicator ASV", 
         y = "Log-Contrast Coefficient") +
    scale_y_continuous(breaks = seq(-0.2, 0.6, by = 0.1)) +  # Define the x-axis breaks
    scale_fill_manual(
      values = c("FALSE" = "#440154FF", "TRUE" = "#21908CFF"),  # Assign colors
      labels = c("Fall", "Spring"),  # Legend labels
      name = "Season") +
    theme_minimal() +  # Clean theme
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)) +
    coord_flip()
)

# Save
ggsave(
  filename = "outputs/indicators/S18-coda-seasonal-sp.png",
  plot = signature_plot, width = 8, height = 6, dpi = 300
)


# 18S Depth Model
S18_coda_glmnet_depth<-coda_glmnet(x=S18_df,y=S18_y_depth, nfolds = 3)

# NO RESULTS HERE EITHER






####
# Data Prep ------------  COI
####
COI_df <- otu_table(COI_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(COI_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
COI_y_season <- metadata[,"Season"] %>% 
  as.factor()
COI_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

COI_df <- column_to_rownames(COI_df, var = "sampleid")


# COI Seasonal Model
COI_coda_glmnet_seasonal<-coda_glmnet(x=COI_df,y=COI_y_season, nfolds = 3)

####
# Model evaluation
###
# Extract data from the COI_coda_glmnet_seasonal object
asv_num <- COI_coda_glmnet_seasonal$taxa.num
asv_name <- COI_coda_glmnet_seasonal$taxa.name
coefficients <- COI_coda_glmnet_seasonal$`log-contrast coefficients`

# Combine to a table
COI_seasonal_results <- data.frame(
  ASV_num = asv_num,
  ASV = asv_name,
  Coefficient = coefficients
) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(Coefficient, decreasing = TRUE)]))


(signature_plot <- ggplot(COI_seasonal_results, aes(x = ASV, y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +  # Add reference line at y=0
    labs(x = "Indicator ASV", 
         y = "Log-Contrast Coefficient") +
    scale_y_continuous(breaks = seq(-0.4, 0.6, by = 0.1)) +  # Define the x-axis breaks
    scale_fill_manual(
      values = c("FALSE" = "#440154FF", "TRUE" = "#21908CFF"),  # Assign colors
      labels = c("Fall", "Spring"),  # Legend labels
      name = "Season") +
    theme_minimal() +  # Clean theme
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)) +
    coord_flip()
)

# Save
ggsave(
  filename = "outputs/indicators/COI-coda-seasonal-sp.png",
  plot = signature_plot, width = 8, height = 6, dpi = 300
)


# COI Depth model
COI_coda_glmnet_depth<-coda_glmnet(x=COI_df,y=COI_y_depth, nfolds = 3)

####
# Model evaluation
###
# Extract data from the S16_coda_glmnet_seasonal object
asv_num <- COI_coda_glmnet_depth$taxa.num
asv_name <- COI_coda_glmnet_depth$taxa.name
coefficients <- COI_coda_glmnet_depth$`log-contrast coefficients`

# Combine to a table
COI_depth_results <- data.frame(
  ASV_num = asv_num,
  ASV = asv_name,
  Coefficient = coefficients
) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(Coefficient, decreasing = TRUE)]))


(signature_plot <- ggplot(COI_depth_results, aes(x = ASV, y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +  # Add reference line at y=0
    labs(x = "Indicator ASV", 
         y = "Log-Contrast Coefficient") +
    scale_y_continuous(limits = c(-0.5, 0.25),  # Set the range
                       breaks = seq(-0.5, 0.25, by = 0.25)) +  # Define the x-axis breaks
    scale_fill_manual(
      values = c("FALSE" = "#A6CEE3", "TRUE" = "#08306B"),  # Assign colors
      labels = c("Photic", "Aphotic"),  # Legend labels
      name = "Sampling Depth") +
    theme_minimal() +  # Clean theme
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)) +
    coord_flip()
)

# Save
ggsave(
  filename = "outputs/indicators/COI-coda-depth-sp.png",
  plot = signature_plot, width = 8, height = 8, dpi = 300
)
