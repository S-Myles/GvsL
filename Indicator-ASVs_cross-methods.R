library(tidyverse)

rda_12S_ASVs <- rda_12S_ASVs %>% rownames_to_column(var= "ASV")
rda_16S_ASVs <- rda_16S_ASVs %>% rownames_to_column(var= "ASV")
rda_18S_ASVs <- rda_18S_ASVs %>% rownames_to_column(var= "ASV")
rda_COI_ASVs <- rda_COI_ASVs %>% rownames_to_column(var= "ASV")

####
# ----------16S
####
# Merging 16S Depth indicator ASVs from multiple methods
merged_16S_depth_ASVs <- rda_16S_ASVs %>% 
  select(ASV, Depth20, Depth250) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  merge(indicsp_16S_depth, by = "ASV")

# Merging 16S Seasonal indicator ASVs from multiple methods
merged_16S_season_ASVs <- rda_16S_ASVs %>% 
  select(ASV, SeasonS) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  merge(indicsp_16S_season, by = "ASV")



####
# ----------18S
####
# Merging 18S Depth indicator ASVs from multiple methods
merged_18S_depth_ASVs <- rda_18S_ASVs %>% 
  select(ASV, Depth20, Depth250) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  merge(indicsp_18S_depth, by = "ASV")

# Merging 16S Seasonal indicator ASVs from multiple methods
merged_18S_season_ASVs <- rda_18S_ASVs %>% 
  select(ASV, SeasonS) %>%
  filter(if_any(where(is.numeric), ~ . > 0.7 | . < -0.7)) %>% # select 0.7 permanova scores only
  merge(indicsp_18S_season, by = "ASV")
