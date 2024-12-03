# Packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all

###########
# Load 4 markers metabarcoding data and create phyloseq objects
###########

# Project metadata
metadata <-  read_tsv("./data/metadata.txt") %>%
  column_to_rownames(var = "sampleid") %>% 
  mutate(Date = paste(Year, Season, sep = "_"))
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
(S12_physeq_data <- subset_samples(S12_physeq_data, Size.Fraction=="L") %>% 
   filter_taxa(function(x) sum(x > 1) > 0, TRUE))
(S16_physeq_data <- subset_samples(S16_physeq_data, Size.Fraction=="S") %>% 
    filter_taxa(function(x) sum(x > 1) > 0, TRUE))
(S18_physeq_data <- subset_samples(S18_physeq_data, Size.Fraction=="L") %>% 
    filter_taxa(function(x) sum(x > 1) > 0, TRUE))
(COI_physeq_data <- subset_samples(COI_physeq_data, Size.Fraction=="L") %>% 
    filter_taxa(function(x) sum(x > 1) > 0, TRUE))




####
# Explore and filter by sample sequencing depth
####
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
(S12_deep_data <- prune_samples(sample_sums(S12_physeq_data) >= 2000 , S12_physeq_data))
(S16_deep_data <- prune_samples(sample_sums(S16_physeq_data) >= 2000 , S16_physeq_data))
(S18_deep_data <- prune_samples(sample_sums(S18_physeq_data) >= 2000 , S18_physeq_data))
(COI_deep_data <- prune_samples(sample_sums(COI_physeq_data) >= 2000 , COI_physeq_data))




####
# ASV prevalence filtering
####
# Define prevalence of each taxa  (in how many samples did each taxa appear at least once)
prev12S = apply(X = otu_table(S12_deep_data),
                MARGIN = ifelse(taxa_are_rows(S12_deep_data), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
(S12_prev_filt <- filter_taxa(S12_deep_data, function(x) sum(x > 1) > (0.1*length(x)), TRUE))

# Define prevalence of each taxa  (in how many samples did each taxa appear at least once)
prev16S = apply(X = otu_table(S16_deep_data),
                MARGIN = ifelse(taxa_are_rows(S16_deep_data), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
(S16_prev_filt <- filter_taxa(S16_deep_data, function(x) sum(x > 1) > (0.1*length(x)), TRUE))

# Define prevalence of each taxa  (in how many samples did each taxa appear at least once)
prev18S = apply(X = otu_table(S18_deep_data),
                MARGIN = ifelse(taxa_are_rows(S18_deep_data), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
(S18_prev_filt <- filter_taxa(S18_deep_data, function(x) sum(x > 1) > (0.1*length(x)), TRUE))

# Define prevalence of each taxa  (in how many samples did each taxa appear at least once)
prevCOI = apply(X = otu_table(COI_deep_data),
                MARGIN = ifelse(taxa_are_rows(COI_deep_data), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
(COI_prev_filt <- filter_taxa(COI_deep_data, function(x) sum(x > 1) > 1, TRUE))




####
# ASV total abundance filtering
####
(S12_P_A_filt <- filter_taxa(S12_prev_filt, function(x) sum(x) > 20, TRUE))
(S16_P_A_filt <- filter_taxa(S16_prev_filt, function(x) sum(x) > 20, TRUE))
(S18_P_A_filt <- filter_taxa(S18_prev_filt, function(x) sum(x) > 20, TRUE))
(COI_P_A_filt <- filter_taxa(COI_prev_filt, function(x) sum(x) > 20, TRUE))




#######
# Taxonomic agglomeration
#######
# How many genera are present after filtering?
(S12_filt_glom_c = tax_glom(S12_P_A_filt, taxrank = "Class"))
(S16_filt_glom_c = tax_glom(S16_P_A_filt, taxrank = "Class"))
(S18_filt_glom_c = tax_glom(S18_P_A_filt, taxrank = "Class"))
(COI_filt_glom_c = tax_glom(COI_P_A_filt, taxrank = "Class"))

# How many species are present after filtering?
(S12_filt_glom_f = tax_glom(S12_P_A_filt, taxrank = "Family"))
(S16_filt_glom_f = tax_glom(S16_P_A_filt, taxrank = "Family"))
(S18_filt_glom_f = tax_glom(S18_P_A_filt, taxrank = "Family"))
(COI_filt_glom_f = tax_glom(COI_P_A_filt, taxrank = "Family"))

# How many species are present after filtering?
(S12_filt_glom_g = tax_glom(S12_P_A_filt, taxrank = "Genus"))
(S16_filt_glom_g = tax_glom(S16_P_A_filt, taxrank = "Genus"))
(S18_filt_glom_g = tax_glom(S18_P_A_filt, taxrank = "Genus"))
(COI_filt_glom_g = tax_glom(COI_P_A_filt, taxrank = "Genus"))

# How many species are present after filtering?
(S12_filt_glom_s = tax_glom(S12_P_A_filt, taxrank = "Species"))
(S16_filt_glom_s = tax_glom(S16_P_A_filt, taxrank = "Species"))
(S18_filt_glom_s = tax_glom(S18_P_A_filt, taxrank = "Species"))
(COI_filt_glom_s = tax_glom(COI_P_A_filt, taxrank = "Species"))



table(sample_data(S12_physeq_pruned)$Size.Fraction)
table(sample_data(S12_physeq_pruned)$Station)
table(sample_data(S12_physeq_pruned)$Year)
table(sample_data(S12_physeq_pruned)$Season)
table(sample_data(S12_physeq_pruned)$Depth)

table(sample_data(S12_physeq_pruned)$Size.Fraction)
table(sample_data(S12_physeq_pruned)$Station)
table(sample_data(S12_physeq_pruned)$Year)
table(sample_data(S12_physeq_pruned)$Season)
table(sample_data(S12_physeq_pruned)$Depth)


# Step 1: Extract sample IDs
sample_ids12S <- as.character(sample_names(S12_physeq_data))
sample_ids16S <- as.character(sample_names(S16_physeq_data))
sample_ids18S <- as.character(sample_names(S18_physeq_data))
sample_idsCOI <- as.character(sample_names(COI_physeq_data))

# Step 2: Find common sample IDs
common_samples <- Reduce(intersect, list(sample_ids12S, sample_ids16S, sample_ids18S, sample_idsCOI))

# Step 3: Prune samples in each phyloseq object
S12_physeq_pruned <- prune_samples(common_samples, S12_deep_data)
S16_physeq2_pruned <- prune_samples(common_samples, S16_deep_data)
S18_physeq1_pruned <- prune_samples(common_samples, S18_deep_data)
COI_physeq2_pruned <- prune_samples(common_samples, COI_deep_data)

# Step 4: Merge the pruned phyloseq objectsCOI_physeq2_pruned
merged_physeq <- merge_phyloseq(S12_physeq_pruned, S16_physeq2_pruned, S18_physeq1_pruned,COI_physeq2_pruned)
