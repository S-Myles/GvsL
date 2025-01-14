library(indicspecies)
library(tidyverse)

# These indicator psecies analyses are to be done only on 1 size-fraction datasets because
# RDAs showed fundametal impact on the data. 

# Only investigate sample attributes showing a significant association to community organization 
# from rd-RDA anova analysis

# Year should not be a factor to investigate. So only Season, Depth, or Station if significant.

# Depth should be binned to a binary as photic (1 and 20m) or aphotic (250 m)


# So for 
# 12S : none (only year significant with 1 size fraction, station when all samples present)
# 16S : Depth & Season
# 18S : Depth & Season 
# COI : Depth & Season 


# Here a function to extract only significant results from analyses bellow
extract_multipatt_results <- function(multipatt_result, dataset, grouping_factor) {
  # Extract species names, indicator values, and p-values
  results <- as.data.frame(multipatt_result$sign)
  # Filter to include only ASVs with p-value < 0.05
  results <- results %>%
    rownames_to_column(var = "ASV") %>%
    filter(p.value < 0.05) %>%  # Filter significant ASVs
    mutate(
      Dataset = dataset,
      Grouping = grouping_factor
    )
  
  return(results)
}



####################################################
####    16S Indic Species analysis
######################################################

# Loading the dataset
S16_df <- otu_table(S16_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
# Aligning metadata to ASV table order
metadata <- merge(S16_df, metadata, "sampleid") # Merge by sampleid

# Create sample attribute binary vectors from metadata
S16_y_season <- metadata[,"Season"] %>% 
  as.factor()
S16_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

# Return rownames
S16_df <- column_to_rownames(S16_df, var = "sampleid")



# Indicator species analysis for Season
ind_season <- multipatt(S16_df, S16_y_season, control = how(nperm = 999))
summary(ind_season)

# Indicator species analysis for Depth 
ind_depth <- multipatt(S16_df, S16_y_depth, control = how(nperm = 999))
summary(ind_depth)


# Extracting 16S Results to dataframes
indicsp_16S_season <- extract_multipatt_results(ind_season, "16S", "Season")
indicsp_16S_depth <- extract_multipatt_results(ind_depth, "16S", "Depth")




####################################################
####    18S Indic Species analysis
######################################################

# Loading the dataset
S18_df <- otu_table(S18_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
# Aligning metadata to ASV table order
metadata <- merge(S18_df, metadata, "sampleid") # Merge by sampleid

# Create sample attribute binary vectors from metadata
S18_y_season <- metadata[,"Season"] %>% 
  as.factor()
S18_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

# Return rownames
S18_df <- column_to_rownames(S18_df, var = "sampleid")



# Indicator species analysis for Season
ind_season <- multipatt(S18_df, S18_y_season, control = how(nperm = 999))
summary(ind_season)


# Indicator species analysis for Depth
ind_depth <- multipatt(S18_df, S18_y_depth, control = how(nperm = 999))
summary(ind_depth)


# Extract 18S Results to dataframes
indicsp_18S_season <- extract_multipatt_results(ind_season, "18S", "Season")
indicsp_18S_depth <- extract_multipatt_results(ind_depth, "18S", "Depth")








####################################################
####    COI Indic Species analysis
######################################################

# Loading the dataset
COI_df <- otu_table(COI_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
# Aligning metadata to ASV table order
metadata <- merge(COI_df, metadata, "sampleid") # Merge by sampleid

# Create sample attribute binary vectors from metadata
COI_y_season <- metadata[,"Season"] %>% 
  as.factor()
COI_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

# Return rownames
COI_df <- column_to_rownames(COI_df, var = "sampleid")



# Indicator species analysis for Season
ind_season <- multipatt(COI_df, COI_y_season, control = how(nperm = 999))
summary(ind_season)


# Indicator species analysis for Depth
ind_depth <- multipatt(COI_df, COI_y_depth, control = how(nperm = 999))
summary(ind_depth)



# Extract COI Results to dataframes
indicsp_COI_season <- extract_multipatt_results(ind_season, "COI", "Season")
indicsp_COI_depth <- extract_multipatt_results(ind_depth, "COI", "Depth")


# Optional Cleanup 
rm(ind_depth, ind_season, ind_station, extract_multipatt_results)
rm(list = ls(pattern = "^(COI_|S12_|S16_|S18_)(df|y_depth|y_season|y_station)"))