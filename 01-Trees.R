library(ape)
########
# 12S
########
# -------------- Trees
########

S12_tree <- read.tree("data/12S_tree_rooted.nwk")
plot(S12_tree)

# Find the MRCA of the two outgroup taxa (sharks)
outgroup_mrca <- getMRCA(S12_tree, tip = c("05fb9f9e7b4e1d7597209de9a6d444fe", 
                                           "c440dd4feac74c8e0c6a7946c1261bbe"))

# Root the tree at this MRCA
S12_tree_rerooted <- root(S12_tree, node = outgroup_mrca)
plot(S12_tree_rerooted)



# Add tree to dataset
physeq <- merge_phyloseq(S12_taxfilt_data, 
                         phy_tree(S12_tree_rerooted))

# plot tree with qiime classified data
(p <- plot_tree(physeq, 
                color = "Order", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right",     
          panel.border = element_blank())
  )

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
    theme(legend.position = "right",     
          panel.border = element_blank())
)

ggsave("outputs/Trees/12S-Species_tree.png", 
       plot = p, width = 10, height = 14)



########
# 16S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S16_taxfilt_data, 
                         phy_tree(read.tree("data/16S_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                #label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right",     
          panel.border = element_blank())
)

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
    theme(legend.position = "right",     
          panel.border = element_blank())
)

ggsave("outputs/Trees/16S-Phylum_tree.png", 
       plot = p, width = 10, height = 14)



########
# 18S
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(S18_taxfilt_data, 
                         phy_tree(read.tree("data/18S_tree_rooted.nwk")))
# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right",     
          panel.border = element_blank())
)

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
    theme(legend.position = "right",     
          panel.border = element_blank())
)

ggsave("outputs/Trees/18S-Phylum_tree.png", 
       plot = p, width = 10, height = 14)



########
# COI
########
# -------------- Trees
########
# Add tree to dataset
physeq <- merge_phyloseq(COI_taxfilt_data, 
                         phy_tree(read.tree("data/COI_tree_rooted.nwk")))

# plot tree
(p <- plot_tree(physeq, 
                color = "Phylum", 
                label.tips = "Species", 
                size = "abundance", 
                shape = "Size.Fraction",
                justify = "left") +
    theme(legend.position = "right",     
          panel.border = element_blank())
)

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
    theme(legend.position = "right",     
          panel.border = element_blank())
)

ggsave("outputs/Trees/COI-Phylum_tree.png", 
       plot = p, width = 10, height = 14)
