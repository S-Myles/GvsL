library(ape)
########
# 12S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S12_physeq_data, 
                         phy_tree(read.tree("data/12S_aligned_masked_tree_rooted.nwk")))

# plot tree with qiime classified data
(p <- plot_tree(physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/12S-ASV_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(physeq, taxrank = "Species")

# plot tree
(p <- plot_tree(physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/12S-Species_tree.png", 
       plot = p, width = 10, height = 14)

#########
# Now plot same tree with Blastn results
# Blast tree needs new taxonomy
# Taxonomy reference table
S12_taxonomy <- read_tsv("./data/12S_blastn_taxized.txt") %>%
  column_to_rownames('ASV_ID') %>%
  as.matrix()


S12_blast_physeq <- merge_phyloseq(otu_table(S12_physeq_data), 
                                   sample_data(S12_physeq_data), 
                                   tax_table(S12_taxonomy),
                                   phy_tree(read.tree("data/12S_aligned_masked_tree_rooted.nwk")))

# plot tree with qiime classified data
(p <- plot_tree(S12_blast_physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/12S-ASV_BLAST_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(S12_blast_physeq, taxrank = "Species")

# plot tree
(p <- plot_tree(S12_blast_physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction") +
    #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/12S-Species_BLAST_tree.png", 
       plot = p, width = 10, height = 14)







########
# 16S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S16_physeq_data, 
                         phy_tree(read.tree("data/16S_filt_aligned_masked_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                #label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/16S-ASV_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(physeq, taxrank = "Phylum")

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Phylum", 
                size = "abundance", 
                shape = "Size.Fraction")+
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/16S-Phylum_tree.png", 
       plot = p, width = 10, height = 14)





########
# 18S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S18_physeq_data, 
                         phy_tree(read.tree("data/18S_filt_aligned_masked_tree_rooted.nwk")))
# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/18S-ASV_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(physeq, taxrank = "Phylum")

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Phylum", 
                size = "abundance", 
                shape = "Size.Fraction") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/18S-Phylum_tree.png", 
       plot = p, width = 10, height = 14)




########
# COI
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(COI_physeq_data, 
                         phy_tree(read.tree("data/COI_filt_aligned_masked_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/COI-ASV_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(physeq, taxrank = "Phylum")

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Phylum", 
                size = "abundance", 
                shape = "Size.Fraction") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/COI-Phylum_tree.png", 
       plot = p, width = 10, height = 14)


#########
# Now plot same tree with Blastn results
# Blast tree needs new taxonomy
# Taxonomy reference table
COI_taxonomy <- read_tsv("./data/COI_blastn_taxized.txt") %>%
  column_to_rownames('Feature ID') %>%
  as.matrix()


COI_blast_physeq <- merge_phyloseq(otu_table(COI_physeq_data), 
                                   sample_data(COI_physeq_data), 
                                   tax_table(COI_taxonomy),
                                   phy_tree(read.tree("data/COI_filt_aligned_masked_tree_rooted.nwk")))

# plot tree with qiime classified data
(p <- plot_tree(COI_blast_physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/COI-ASV_BLAST_tree.pdf", 
       plot = p, width = 22, height = 49)

# Optionally aggregate dataset
physeq <- tax_glom(COI_blast_physeq, taxrank = "Phylum")


# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Phylum", 
                size = "abundance", 
                shape = "Size.Fraction") +
    #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/COI-Phylum_BLAST_tree.png", 
       plot = p, width = 10, height = 14)
