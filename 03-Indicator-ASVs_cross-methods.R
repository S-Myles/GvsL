library(tidyverse)

rda_12S_ASVs <- rda_12S_ASVs %>% rownames_to_column(var= "ASV")
rda_16S_ASVs <- rda_16S_ASVs %>% rownames_to_column(var= "ASV")
rda_18S_ASVs <- rda_18S_ASVs %>% rownames_to_column(var= "ASV")
rda_COI_ASVs <- rda_COI_ASVs %>% rownames_to_column(var= "ASV")


####
# ----------16S
####
####
# Get Taxonomies
####
# Extract taxonomy table from phyloseq object
taxonomy_table <- as.data.frame(tax_table(S16_filt_data))

# Create taxonomy column based on the best available rank
taxonomy_table$Taxonomy <- if_else(
  !is.na(taxonomy_table$Species) & taxonomy_table$Species != "",
  taxonomy_table$Species,  # Use Species if available
  if_else(
    !is.na(taxonomy_table$Genus) & taxonomy_table$Genus != "",
    paste(taxonomy_table$Genus, "sp.", sep = " "),  # Use Genus with " sp."
    if_else(
      !is.na(taxonomy_table$Family) & taxonomy_table$Family != "",
      paste("Family_", taxonomy_table$Family, sep = ""),  # Use Family with "Family_" prefix
      if_else(
        !is.na(taxonomy_table$Order) & taxonomy_table$Order != "",
        paste("Order_", taxonomy_table$Order, sep = ""),  # Use Order with "Order_" prefix
        if_else(
          !is.na(taxonomy_table$Class) & taxonomy_table$Class != "",
          paste("Class_", taxonomy_table$Class, sep = ""),  # Use Class with "Class_" prefix
          NA_character_  # Default to NA if no rank is available
        )
      )
    )
  )
)

# Make ASV IDs the rownames for joining
taxonomy_table <- taxonomy_table %>% rownames_to_column(var = "ASV")

#####
#  Indicator ASVs merging
#####
# Merging 16S Depth indicator ASVs from multiple methods
merged_16S_depth_ASVs <- rda_16S_ASVs %>%
  select(ASV, Depth20, Depth250) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # Select 0.7 permanova scores only
  full_join(indicsp_16S_depth, by = "ASV")

# Add taxonomy
merged_16S_depth_ASVs <- left_join(merged_16S_depth_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, Depth20, Depth250, s.aphotic, s.photic, index, stat, p.value, Dataset, Grouping)

write_csv(merged_16S_depth_ASVs, "outputs/indicators/Merged_indicator_ASVS-16S-Depth.csv")

  


# Merging 16S Seasonal indicator ASVs from multiple methods
merged_16S_season_ASVs <- rda_16S_ASVs %>% 
  select(ASV, SeasonS) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  full_join(indicsp_16S_season, by = "ASV")

# Add taxonomy
merged_16S_season_ASVs <- left_join(merged_16S_season_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, SeasonS, s.F, s.S, index, stat, p.value, Dataset, Grouping)

write_csv(merged_16S_season_ASVs, "outputs/indicators/Merged_indicator_ASVS-16S-Season.csv")




####
# ----------18S
####
####
# Get Taxonomies
####
# Extract taxonomy table from phyloseq object
taxonomy_table <- as.data.frame(tax_table(S18_filt_data))

# Create taxonomy column based on the best available rank
taxonomy_table$Taxonomy <- if_else(
  !is.na(taxonomy_table$Species) & taxonomy_table$Species != "",
  taxonomy_table$Species,  # Use Species if available
  if_else(
    !is.na(taxonomy_table$Genus) & taxonomy_table$Genus != "",
    paste(taxonomy_table$Genus, "sp.", sep = " "),  # Use Genus with " sp."
    if_else(
      !is.na(taxonomy_table$Family) & taxonomy_table$Family != "",
      paste("Family_", taxonomy_table$Family, sep = ""),  # Use Family with "Family_" prefix
      if_else(
        !is.na(taxonomy_table$Order) & taxonomy_table$Order != "",
        paste("Order_", taxonomy_table$Order, sep = ""),  # Use Order with "Order_" prefix
        if_else(
          !is.na(taxonomy_table$Class) & taxonomy_table$Class != "",
          paste("Class_", taxonomy_table$Class, sep = ""),  # Use Class with "Class_" prefix
          NA_character_  # Default to NA if no rank is available
        )
      )
    )
  )
)

# Make ASV IDs the rownames for joining
taxonomy_table <- taxonomy_table %>% rownames_to_column(var = "ASV")


#####
#  Indicator ASVs merging
#####
# Merging 18S Depth indicator ASVs from multiple methods
merged_18S_depth_ASVs <- rda_18S_ASVs %>% 
  select(ASV, Depth20, Depth250) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  full_join(indicsp_18S_depth, by = "ASV")
# Add taxonomy
merged_18S_depth_ASVs <- left_join(merged_18S_depth_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, Depth20, Depth250, s.aphotic, s.photic, index, stat, p.value, Dataset, Grouping)

write_csv(merged_18S_depth_ASVs, "outputs/indicators/Merged_indicator_ASVS-18S-Depth.csv")


# Merging 16S Seasonal indicator ASVs from multiple methods
merged_18S_season_ASVs <- rda_18S_ASVs %>% 
  select(ASV, SeasonS) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  full_join(indicsp_18S_season, by = "ASV")

# Add taxonomy
merged_18S_season_ASVs <- left_join(merged_18S_season_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, SeasonS, s.F, s.S, index, stat, p.value, Dataset, Grouping)

write_csv(merged_18S_season_ASVs, "outputs/indicators/Merged_indicator_ASVS-18S-Season.csv")





####
# ----------COI
####
####
# Get Taxonomies
####
# Extract taxonomy table from phyloseq object
taxonomy_table <- as.data.frame(tax_table(COI_filt_data))

# Create taxonomy column based on the best available rank
taxonomy_table$Taxonomy <- if_else(
  !is.na(taxonomy_table$Species) & taxonomy_table$Species != "",
  taxonomy_table$Species,  # Use Species if available
  if_else(
    !is.na(taxonomy_table$Genus) & taxonomy_table$Genus != "",
    paste(taxonomy_table$Genus, "sp.", sep = " "),  # Use Genus with " sp."
    if_else(
      !is.na(taxonomy_table$Family) & taxonomy_table$Family != "",
      paste("Family_", taxonomy_table$Family, sep = ""),  # Use Family with "Family_" prefix
      if_else(
        !is.na(taxonomy_table$Order) & taxonomy_table$Order != "",
        paste("Order_", taxonomy_table$Order, sep = ""),  # Use Order with "Order_" prefix
        if_else(
          !is.na(taxonomy_table$Class) & taxonomy_table$Class != "",
          paste("Class_", taxonomy_table$Class, sep = ""),  # Use Class with "Class_" prefix
          NA_character_  # Default to NA if no rank is available
        )
      )
    )
  )
)

# Make ASV IDs the rownames for joining
taxonomy_table <- taxonomy_table %>% rownames_to_column(var = "ASV")


#####
#  Indicator ASVs merging
#####
# Merging 18S Depth indicator ASVs from multiple methods
merged_COI_depth_ASVs <- rda_COI_ASVs %>% 
  select(ASV, Depth20, Depth250) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  full_join(indicsp_COI_depth, by = "ASV")
# Add taxonomy
merged_COI_depth_ASVs <- left_join(merged_COI_depth_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, Depth20, Depth250, s.aphotic, s.photic, index, stat, p.value, Dataset, Grouping)

write_csv(merged_COI_depth_ASVs, "outputs/indicators/Merged_indicator_ASVS-COI-Depth.csv")


# Merging 16S Seasonal indicator ASVs from multiple methods
merged_COI_season_ASVs <- rda_COI_ASVs %>% 
  select(ASV, SeasonS) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  full_join(indicsp_COI_season, by = "ASV")
# Add taxonomy
merged_COI_season_ASVs <- left_join(merged_COI_season_ASVs, taxonomy_table, by = "ASV") %>% 
  select(ASV, Taxonomy, SeasonS, s.F, s.S, index, stat, p.value, Dataset, Grouping)

write_csv(merged_COI_season_ASVs, "outputs/indicators/Merged_indicator_ASVS-COI-Season.csv")