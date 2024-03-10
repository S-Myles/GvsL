# Packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all

###########
# Load 4 markers metabarcoding data and create phyloseq objects
###########

# Project metadata
metadata <-  read_tsv("./data/metadata.txt") %>%
  column_to_rownames(var = "sampleid")
metadata <- sample_data(metadata) 

####
# 12S MiFish-U2
####
# ASV table (Sample X Taxon, read counts)
S12_ASV_table <- read_tsv("./data/12S_feature-table.tsv") %>%
  column_to_rownames('ASV_ID')     
# Taxonomy reference table
S12_taxonomy <- read_tsv("./data/12S_taxonomy.tsv") %>%
  column_to_rownames('Feature ID') %>%
  as.matrix()

# Respective phyloseq objects
S12_ASV_table <- otu_table(S12_ASV_table, taxa_are_rows = TRUE)         
S12_taxonomy <- tax_table(S12_taxonomy)
# Merging into 1 global object
(S12_physeq_data <-  merge_phyloseq(S12_ASV_table, metadata, S12_taxonomy)) 

# Cleanup if you want
rm(S12_ASV_table, S12_taxonomy)
##################################


####
# 16S V4V5
####
# ASV table (Sample X Taxon, read counts)
S16_ASV_table <- read_tsv("./data/16S_feature-table.tsv") %>%
  column_to_rownames('ASV_ID')      
# Taxonomy reference table
S16_taxonomy <- read_tsv("./data/16S_taxonomy.tsv") %>%
  column_to_rownames('Feature ID') %>%
  as.matrix()

# Respective phyloseq objects
S16_ASV_table <- otu_table(S16_ASV_table, taxa_are_rows = TRUE)         
S16_taxonomy <- tax_table(S16_taxonomy)
# Merging into 1 global object
(S16_physeq_data <-  merge_phyloseq(S16_ASV_table, metadata, S16_taxonomy))

# Cleanup if you want
rm(S16_ASV_table, S16_taxonomy)
##################################


####
# 18S V4
####
# ASV table (Sample X Taxon, read counts)
S18_ASV_table <- read_tsv("./data/18S_feature-table.tsv") %>%
  column_to_rownames('ASV_ID')      
# Taxonomy reference table
S18_taxonomy <- read_tsv("./data/18S_taxonomy.tsv") %>%
  column_to_rownames('Feature ID') %>%
  as.matrix()

# Respective phyloseq objects
S18_ASV_table <- otu_table(S18_ASV_table, taxa_are_rows = TRUE)         
S18_taxonomy <- tax_table(S18_taxonomy)
# Merging into 1 global object
(S18_physeq_data <-  merge_phyloseq(S18_ASV_table, metadata, S18_taxonomy))

# Cleanup if you want
rm(S18_ASV_table, S18_taxonomy)
##################################


####
# COI Leray-XT
####
# ASV table (Sample X Taxon, read counts)
COI_ASV_table <- read_tsv("./data/COI_feature-table.tsv") %>%
  column_to_rownames('ASV_ID')     
# Taxonomy reference table
COI_taxonomy <- read_tsv("./data/COI_taxonomy.tsv") %>%
  column_to_rownames('Feature ID') %>%
  as.matrix()

# Respective phyloseq objects
COI_ASV_table <- otu_table(COI_ASV_table, taxa_are_rows = TRUE)         
COI_taxonomy <- tax_table(COI_taxonomy)
# Merging into 1 global object
(COI_physeq_data <-  merge_phyloseq(COI_ASV_table, metadata, COI_taxonomy)) 

# Cleanup if you want
rm(COI_ASV_table, COI_taxonomy)
###############################################################################







###############################################################################
# Filtering and pre-processing datasets
#########################################

####
# Remove small size fractions
####
(S12_physeq_data <- subset_samples(S12_physeq_data, Size.Fraction=="L"))
(S16_physeq_data <- subset_samples(S16_physeq_data, Size.Fraction=="S"))
(S18_physeq_data <- subset_samples(S18_physeq_data, Size.Fraction=="L"))
(COI_physeq_data <- subset_samples(COI_physeq_data, Size.Fraction=="L"))


####
# Explore and filter by sample sequencing depth
####
# Find outliers
(sample_sum_df <- data.frame(sum = sample_sums(S12_physeq_data)) %>% 
  arrange(desc(sum)))

(sample_sum_df <- data.frame(sum = sample_sums(S16_physeq_data)) %>% 
    arrange(desc(sum)))

(sample_sum_df <- data.frame(sum = sample_sums(S18_physeq_data)) %>% 
    arrange(desc(sum)))

(sample_sum_df <- data.frame(sum = sample_sums(COI_physeq_data)) %>% 
    arrange(desc(sum)))

# Find out more
#out <- boxplot.stats(sample_sum_df$sum)$out
#out_ind <- which(sample_sum_df$sum %in% c(out))
#outliers <- row.names(metadata[out_ind])

# remove samples with arbitrary read depths threshold
# 2000 read depth requirement would limit loss but clean up important bad quality
(S12_physeq_data <- prune_samples(sample_sums(S12_physeq_data) >= 2000 , S12_physeq_data))
(S16_physeq_data <- prune_samples(sample_sums(S16_physeq_data) >= 2000 , S16_physeq_data))
(S18_physeq_data <- prune_samples(sample_sums(S18_physeq_data) >= 2000 , S18_physeq_data))
(COI_physeq_data <- prune_samples(sample_sums(COI_physeq_data) >= 2000 , COI_physeq_data))


####
# ASV abundance filtering
####
(S12_physeq_filt <- filter_taxa(S12_physeq_data, function(x) sum(x) > 100, TRUE))
(S16_physeq_filt <- filter_taxa(S16_physeq_data, function(x) sum(x) > 100, TRUE))
(S18_physeq_filt <- filter_taxa(S18_physeq_data, function(x) sum(x) > 100, TRUE))
(COI_physeq_filt <- filter_taxa(COI_physeq_data, function(x) sum(x) > 100, TRUE))


####
# ASV prevalence filtering
####
(S12_physeq_filt <- filter_taxa(S12_physeq_filt, function(x) sum(x > 1) > (0.1*length(x)), TRUE))
(S16_physeq_filt <- filter_taxa(S16_physeq_filt, function(x) sum(x > 1) > (0.1*length(x)), TRUE))
(S18_physeq_filt <- filter_taxa(S18_physeq_filt, function(x) sum(x > 1) > (0.1*length(x)), TRUE))
(COI_physeq_filt <- filter_taxa(COI_physeq_filt, function(x) sum(x > 1) > (0.1*length(x)), TRUE))
