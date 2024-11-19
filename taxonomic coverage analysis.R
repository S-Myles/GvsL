## Taxonomy coverage analysis
library(tidyverse)
library(phyloseq)


######
# 12S
######
# Import data
S12_tax <- tax_table(S12_physeq_data) %>% 
  as.data.frame()

# Calculate assigned vs unassigned proportions
S12_proportion_table <- S12_tax %>% 
  summarise(across(everything(),~ tibble(prop_non_na = sum(!is.na(.)) / n(),
                                         prop_na = sum(is.na(.)) / n()
                                         ))) %>%
  pivot_longer(cols = everything(), names_to = "taxonomic_level", values_to = "Proportion")
S12_proportion_table <- bind_cols(S12_proportion_table$taxonomic_level, S12_proportion_table$Proportion)
colnames(S12_proportion_table) <- c("Taxonomic_level", "Proportion_assigned", "Proportion_unknown")

# Make it plotable
S12_proportion_table_long <- S12_proportion_table %>%
  pivot_longer(
    cols = c(Proportion_assigned, Proportion_unknown),
    names_to = "proportion_type",
    values_to = "value"
  ) %>% 
  mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>% 
  mutate(proportion_type = factor(proportion_type, levels = c("Proportion_unknown", "Proportion_assigned")))

# Plot stacked proportions (up to 1) for each level
(plot <- ggplot(S12_proportion_table_long, aes(x = Taxonomic_level, y = value, fill = proportion_type)) +
  geom_bar(stat = "identity") +
#  labs(
#    x = "Taxonomic Level",
#    y = "Proportion"
#  ) +
  scale_fill_manual(
    values = c("Proportion_unknown" = "black", "Proportion_assigned" = "#4A5BE0"),    # Gene marker color
    name = "Proportion Type",
    labels = c("Unkown", "Assigned")
  ) +
  theme_minimal() +
  theme(
    #plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        #axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 14, face = "bold", vjust = 2),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none"
        ))

#Save
ggsave(
  filename = "outputs/Fig1/S12-raw_tax_coverage.png",
  plot = plot, width = 10, height = 2, dpi = 300
)
##########################################################################################################


##########
#16S
##########
S16_tax <- tax_table(S16_physeq_data) %>% 
  as.data.frame()

S16_proportion_table <- S16_tax %>% 
  summarise(across(everything(),~ tibble(prop_non_na = sum(!is.na(.)) / n(),
                                         prop_na = sum(is.na(.)) / n()
  ))) %>%
  pivot_longer(cols = everything(), names_to = "taxonomic_level", values_to = "Proportion")
S16_proportion_table <- bind_cols(S16_proportion_table$taxonomic_level, S16_proportion_table$Proportion)
colnames(S16_proportion_table) <- c("Taxonomic_level", "Proportion_assigned", "Proportion_unknown")


S16_proportion_table_long <- S16_proportion_table %>%
  pivot_longer(
    cols = c(Proportion_assigned, Proportion_unknown),
    names_to = "proportion_type",
    values_to = "value"
  ) %>% 
  mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>% 
  mutate(proportion_type = factor(proportion_type, levels = c("Proportion_unknown", "Proportion_assigned")))

(plot <- ggplot(S16_proportion_table_long, aes(x = Taxonomic_level, y = value, fill = proportion_type)) +
    geom_bar(stat = "identity") +
    #  labs(
    #    x = "Taxonomic Level",
    #    y = "Proportion"
    #  ) +
    scale_fill_manual(
      values = c("Proportion_unknown" = "black", "Proportion_assigned" = "#9c179e"),
      name = "Proportion Type",
      labels = c("Unkown", "Assigned")
    ) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      #axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(size = 14, face = "bold", vjust = 2),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.position = "none"
    ))

ggsave(
  filename = "outputs/Fig1/S16-raw_tax_coverage.png",
  plot = plot, width = 10, height = 2, dpi = 300
)
####################################################################################################




##########
#18S
##########
S18_tax <- tax_table(S18_physeq_data) %>% 
  as.data.frame()

S18_proportion_table <- S18_tax %>% 
  summarise(across(everything(),~ tibble(prop_non_na = sum(!is.na(.)) / n(),
                                         prop_na = sum(is.na(.)) / n()
  ))) %>%
  pivot_longer(cols = everything(), names_to = "taxonomic_level", values_to = "Proportion")
S18_proportion_table <- bind_cols(S18_proportion_table$taxonomic_level, S18_proportion_table$Proportion)
colnames(S18_proportion_table) <- c("Taxonomic_level", "Proportion_assigned", "Proportion_unknown")


S18_proportion_table_long <- S18_proportion_table %>%
  pivot_longer(
    cols = c(Proportion_assigned, Proportion_unknown),
    names_to = "proportion_type",
    values_to = "value"
  ) %>% 
  mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>% 
  mutate(proportion_type = factor(proportion_type, levels = c("Proportion_unknown", "Proportion_assigned")))

(plot <- ggplot(S18_proportion_table_long, aes(x = Taxonomic_level, y = value, fill = proportion_type)) +
    geom_bar(stat = "identity") +
    #  labs(
    #    x = "Taxonomic Level",
    #    y = "Proportion"
    #  ) +
    scale_fill_manual(
      values = c("Proportion_unknown" = "black", "Proportion_assigned" = "#ed7953"),
      name = "Proportion Type",
      labels = c("Unkown", "Assigned")
    ) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      #axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(size = 14, face = "bold", vjust = 2),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.position = "none"
    ))

ggsave(
  filename = "outputs/Fig1/S18-raw_tax_coverage.png",
  plot = plot, width = 10, height = 2, dpi = 300
)
####################################################################################################





##########
#18S
##########
COI_tax <- tax_table(COI_physeq_data) %>% 
  as.data.frame()

COI_proportion_table <- COI_tax %>% 
  summarise(across(everything(),~ tibble(prop_non_na = sum(!is.na(.)) / n(),
                                         prop_na = sum(is.na(.)) / n()
  ))) %>%
  pivot_longer(cols = everything(), names_to = "taxonomic_level", values_to = "Proportion")
COI_proportion_table <- bind_cols(COI_proportion_table$taxonomic_level, COI_proportion_table$Proportion)
colnames(COI_proportion_table) <- c("Taxonomic_level", "Proportion_assigned", "Proportion_unknown")


COI_proportion_table_long <- COI_proportion_table %>%
  pivot_longer(
    cols = c(Proportion_assigned, Proportion_unknown),
    names_to = "proportion_type",
    values_to = "value"
  ) %>% 
  mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>% 
  mutate(proportion_type = factor(proportion_type, levels = c("Proportion_unknown", "Proportion_assigned")))

(plot <- ggplot(COI_proportion_table_long, aes(x = Taxonomic_level, y = value, fill = proportion_type)) +
    geom_bar(stat = "identity") +
    #  labs(
    #    x = "Taxonomic Level",
    #    y = "Proportion"
    #  ) +
    scale_fill_manual(
      values = c("Proportion_unknown" = "black", "Proportion_assigned" = "#f0f921"),
      name = "Proportion Type",
      labels = c("Unkown", "Assigned")
    ) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      #axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(size = 14, face = "bold", vjust = 2),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.position = "none"
    ))

ggsave(
  filename = "outputs/Fig1/COI-raw_tax_coverage.png",
  plot = plot, width = 10, height = 2, dpi = 300
)
####################################################################################################
