library(tidyverse)

# Step 1: Define phyloseq objects and process unique taxa
phyloseq_list <- list(`12S` = S12_taxfilt_data, 
                      `16S` = S16_taxfilt_data, 
                      `18S` = S18_taxfilt_data, 
                      `COI` = COI_taxfilt_data)

# Function to count unique taxa at each taxonomic level
count_unique_taxa <- function(physeq, dataset_name) {
  tax_table(physeq) %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "Taxonomy_Level", values_to = "Taxon") %>%
    group_by(Taxonomy_Level) %>%
    summarise(Unique_Taxa_Count = n_distinct(Taxon), .groups = "drop") %>%
    mutate(Dataset = dataset_name)
}

# Process unique taxa for all datasets
unique_taxa_counts <- lapply(names(phyloseq_list), function(dataset) {
  count_unique_taxa(phyloseq_list[[dataset]], dataset)
}) %>%
  bind_rows() %>%
  filter(Taxonomy_Level %in% c("Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  mutate(
    Taxonomy_Level = factor(Taxonomy_Level, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")),
    Dataset = factor(Dataset, levels = c("12S", "16S", "18S", "COI"))
  )



# Step 2: Compute shared taxa
S12_tax <- tax_table(S12_taxfilt_data) %>%
  as.data.frame() %>%                   # Convert taxonomy table to a data frame
  select(Phylum, Class, Order, Family, Genus, Species) %>% # Select the desired columns first
  summarise(across(everything(), ~ list(unique(.)))) %>% 
  mutate(dataset = "12S")

S16_tax <- tax_table(S16_taxfilt_data) %>%
  as.data.frame() %>%                   # Convert taxonomy table to a data frame
  select(Phylum, Class, Order, Family, Genus, Species) %>% # Select the desired columns first
  summarise(across(everything(), ~ list(unique(.)))) %>% 
  mutate(dataset = "16S")

S18_tax <- tax_table(S18_taxfilt_data) %>%
  as.data.frame() %>%                   # Convert taxonomy table to a data frame
  select(Phylum, Class, Order, Family, Genus, Species) %>% # Select the desired columns first
  summarise(across(everything(), ~ list(unique(.)))) %>% 
  mutate(dataset = "18S")

COI_tax <- tax_table(COI_taxfilt_data) %>%
  as.data.frame() %>%                   # Convert taxonomy table to a data frame
  select(Phylum, Class, Order, Family, Genus, Species) %>% # Select the desired columns first
  summarise(across(everything(), ~ list(unique(.)))) %>% 
  mutate(dataset = "COI")

all_unique_tax <- bind_rows(S12_tax, S16_tax, S18_tax, COI_tax)


# Function to extract shared taxa and their sources for a given taxonomy level
find_shared_taxa <- function(data, level) {
  # Extract lists of unique taxa for the given level
  taxa_lists <- setNames(data[[level]], data$dataset)
  
  # Combine lists into a stackable format
  all_taxa_with_sources <- stack(taxa_lists)
  
  # Count occurrences of each taxon
  taxa_counts <- table(all_taxa_with_sources$values)
  
  # Filter for taxa appearing in more than one list
  shared_taxa <- names(taxa_counts[taxa_counts > 1])
  
  # Determine which lists each shared taxon comes from
  shared_taxa_sources <- lapply(shared_taxa, function(taxon) {
    unique(all_taxa_with_sources$ind[all_taxa_with_sources$values == taxon])
  })
  
  # Convert to a data frame for inspection
  data.frame(
    Taxonomic_Level = level,
    Shared_Taxon = shared_taxa,
    Source_Lists = I(shared_taxa_sources) # Use I() to keep lists unflattened
  )
}

# Apply to all taxonomy levels
taxonomy_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

shared_taxa_results <- bind_rows(
  lapply(taxonomy_levels, function(level) find_shared_taxa(all_unique_tax, level))
)

for_save <- shared_taxa_results %>%
  unnest(Source_Lists) %>% 
  filter(complete.cases(.))

write_csv(for_save, file = "outputs/taxa_detected_multiple_datasets.csv")


# Step 3: Prepare shared taxa counts for barplot
shared_taxa_counts <- shared_taxa_results %>%
  group_by(Taxonomic_Level) %>%   # Group by the Taxonomic_Level column
  summarise(
    Shared = n(),      # Count the number of entries per level
    .groups = "drop"              # Drop grouping for a clean output
  )



# Step 4: Adjust original dataset counts by subtracting shared taxa
# Unnest the Source_Lists column to count occurrences
shared_taxa_counts_per_dataset <- shared_taxa_results %>%
  unnest(Source_Lists) %>% # Unnest the Source_Lists to expand it into rows
  group_by(Taxonomic_Level, Source_Lists) %>% # Group by Taxonomic_Level and dataset
  summarise(
    Count = n(), # Count occurrences of each dataset in each Taxonomic_Level
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Source_Lists, values_from = Count, values_fill = 0) %>%  # Spread datasets into columns
  select(-`NA`)

unique_taxa_counts <- pivot_wider(unique_taxa_counts, names_from = Dataset, values_from = Unique_Taxa_Count, values_fill = 0)

adjusted_taxa_counts <- unique_taxa_counts[c("12S", "16S", "18S", "COI")]-shared_taxa_counts_per_dataset[c("12S", "16S", "18S", "COI")]
row.names(adjusted_taxa_counts) <- c("Class", "Family", "Genus", "Order", "Phylum", "Species")
adjusted_taxa_counts <- adjusted_taxa_counts %>% rownames_to_column(var = "Taxonomic_Level")

# Step 5: Combine adjusted counts with shared taxa counts
final_taxa_counts <- right_join(adjusted_taxa_counts, shared_taxa_counts, by = "Taxonomic_Level") %>%
  pivot_longer(cols = -Taxonomic_Level, names_to = "Dataset", values_to = "Count") %>%  # Make it plotable
  mutate(
    Taxonomic_Level = factor(Taxonomic_Level,levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")),
    Dataset = factor(Dataset,levels = c("Shared", "12S", "16S", "18S", "COI"))
  )



# Step 6: Plot the stacked bar graph
(plot <- ggplot(final_taxa_counts, aes(x = Taxonomic_Level, y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("gray50", "#4A5BE0", "#9c179e", "#ed7953", "#f0f921")) +
  labs(x = "Taxonomy Level", y = "Unique Taxa Count", fill = "Dataset of origin") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  ))

ggsave(
  filename = "outputs/tax_barplots/Unique_taxa_counts.png",
  plot = plot, width = 12, height = 4, dpi = 300
)
