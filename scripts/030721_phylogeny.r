# Install packages
pacman::p_load("ggtree", "stringr", "ggrepel", "RColorBrewer", "scales")

# Read in the consensus tree
tr <- read.tree("data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa.contree")

# Fix  formatting of tip labels
tr$tip.label <- gsub("NRP_Polyketide", "NRPPolyketide", tr$tip.label)
tr$tip.label <- gsub("_Saccharide", "_Polyketide", tr$tip.label) # also a polyketide, removing saccharide column
tr$tip.label

# Ggtree
gtr <- ggtree(tr, layout = "rectangular", 
              color = "gray50", size = 0.5)

pal <- c("#BBD374", # NRPS
         "#F0D468", # Other
         "#CF8BBA", # T1 PKS
         "#EC8F6B") #RiPPs

metadat <- data.frame(label = gtr$data$label[gtr$data$isTip]) %>%
  dplyr::mutate(bgc_type = paste0(word(label, sep = "_", -1))) %>%
  dplyr::mutate(bgc_color1 = case_when(grepl("Polyketide", label) ~ "#CF8BBA",
                                      grepl("Other", label) ~ "#F0D468",
                                      grepl("NRP", label) ~ "#BBD374",
                                      grepl("study", label) ~ "#EC8F6B",
                                      grepl("polytheonamide|aeronamide", label) ~ NA_character_)) %>%
  dplyr::mutate(bgc_color2 = case_when(grepl("NRPPolyketide", label) ~ "#BBD374")) %>%
  dplyr::mutate(label_trim = gsub("_Polyketide|_NRPPolyketide|_Other|_NRP", "", label)) %>%
  dplyr::mutate(label_trim = gsub("_1", ".1", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("_", " ", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("WP ", "WP_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("BCRE ", "BCRE_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("STRAU ", "STRAU_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("AMYAL ", "AMYAL_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("Sare ", "Sare_", label_trim)) %>%
    dplyr::mutate(label_trim = gsub("A 201A", "A-201A", label_trim)) %>%
    dplyr::mutate(label_trim = gsub("A.102395", "A-102395", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("l 2 Amino 4 methoxy trans 3 butenoic acid",
                                  "l-2-Amino-4-methoxy-trans-3-butenoic acid", label_trim)) %>%
  dplyr::mutate(O_vs_N = case_when(grepl("study", label) ~ "N",
                                   grepl("polytheonamide|aeronamide", label) ~ NA_character_,
                TRUE ~ "O")) %>%
  dplyr::mutate(O_vs_N_color = case_when(grepl("AMYAL_RS47350|AAK83186|BAJ52700|PA2302", label) ~ bgc_color1,
                grepl("study", label) ~ bgc_color1,
                grepl("polytheonamide|aeronamide", label) ~ NA_character_,
                TRUE ~ "gray90")) %>%
  dplyr::mutate(tiplab_color = case_when(grepl("study", label) ~  "#EC8F6B",
                                         grepl("polytheonamide|aeronamide", label) ~ "gray60",
                                         TRUE ~ "black")) 
gdf <- gtr %<+% metadat

pdf("output/20210704_40_MT_tree.pdf", width = 12, height = 7)
gdf +
  geom_nodelab(aes(label = label), size = 3) + 
  geom_tippoint(x = 15.8, color = gdf$data$bgc_color1[gdf$data$isTip], size = 6) +
  geom_tippoint(x = 16.2, color = ifelse(!is.na(gdf$data$bgc_color2[gdf$data$isTip]), "#BBD374", "white"), size = 6) +
  geom_tiplab(aes(label = label_trim), color = gdf$data$tiplab_color[gdf$data$isTip],
              fontface = ifelse(grepl("study", gdf$data$label_trim[gdf$data$isTip]), "bold.italic", "italic"), align=T) + 
  geom_label2(x = 15.4, aes(subset = isTip), fill = gdf$data$O_vs_N_color[gdf$data$isTip], #alpha = 0.5,
             color = "black", label = gdf$data$O_vs_N[gdf$data$isTip])  +
  xlim_tree(16)
dev.off()


# Plot legend
pal <- c("#BBD374", # NRPS
         "#F0D468", # Other #
         "#CF8BBA", # T1 PKS
         "#EC8F6B") # RiPP
txt <- c("NRPS", "Other",
         "PKS", "RiPP")
# Legend
pdf(paste0("output/20210704_legend.pdf"),width=5,height=6)
plot.new()
legend("topleft", legend=txt, fill=pal,
       border=FALSE, bty="n", title = "BGC type")
dev.off()

