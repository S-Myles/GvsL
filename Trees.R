########
# 12S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S12_physeq_data, 
                         phy_tree(read.tree("data/12S_aligned_masked_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Season",
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
                shape = "Season") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/12S-Species_tree.png", 
       plot = p, width = 6, height = 12)






########
# 16S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S16_physeq_data, 
                         phy_tree(read.tree("data/16S_aligned_masked_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                #label.tips = "Species", 
                size = "abundance", 
                shape = "Season",
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
                shape = "Season")+
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/16S-Phylum_tree.png", 
       plot = p, width = 16, height = 22)





########
# 18S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S18_physeq_data, 
                         phy_tree(read.tree("data/18S_aligned_masked_tree_rooted.nwk")))
# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Season",
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
                shape = "Season") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/18S-Phylum_tree.png", 
       plot = p, width = 10, height = 20)




########
# COI
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(COI_physeq_data, 
                         phy_tree(read.tree("data/COI_aligned_masked_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Season",
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
                shape = "Season") +
                #justify = "left") +
    theme(legend.position = "right"))

ggsave("outputs/Trees/COI-Phylum_tree.png", 
       plot = p, width = 6, height = 12)
