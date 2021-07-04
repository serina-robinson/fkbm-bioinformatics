# Install packages
pacman::p_load("ggtree", "stringr", "ggrepel", "RColorBrewer", "scales")

# Read in the consensus tree
tr <- read.tree("data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa.contree")
tr$tip.label[grep("study", tr$tip.label)]
# Fix some tip labels that are having problems
# tr$tip.label[tr$tip.label == "BGC0000769_AAN05766_1_NA_NA"] <- "BGC0000769_AAN05766_1_Mycobacterium_avium_subsp_hominis_A5_Saccharide"
tr$tip.label <- gsub("NRP_Polyketide", "NRPPolyketide", tr$tip.label)
tr$tip.label <- gsub("_Saccharide", "_Polyketide", tr$tip.label) # also polyketide, removing saccharide 
tr$tip.label
# Ggtree
gtr <- ggtree(tr, layout = "rectangular", 
              color = "gray50", size = 0.5)
gtr$data$label

pal <- c("#BBD374", # NRPS
         "#F0D468", # Other
         "gray80", #
        # "#909EC4",  #T2/T3 PKS 
         "#CF8BBA", # T1 PKS
         "#EC8F6B") #RiPPs
       #  "#8FBDA7") # Terpene

metadat <- data.frame(label = gtr$data$label[gtr$data$isTip]) %>%
 # dplyr::mutate(acc = paste0(word(label, sep = "_", 1), "_", word(label, sep = "_", 2))) %>%
  dplyr::mutate(bgc_type = paste0(word(label, sep = "_", -1))) %>%
  dplyr::mutate(bgc_color1 = case_when(grepl("Polyketide", label) ~ "#CF8BBA",
                                      grepl("Other", label) ~ "#F0D468",
                                      grepl("NRP", label) ~ "#BBD374",
                                      grepl("study|polytheonamide|aeronamide", label) ~ "#EC8F6B")) %>%
  dplyr::mutate(bgc_color2 = case_when(grepl("NRPPolyketide", label) ~ "#BBD374")) %>%
  dplyr::mutate(label_trim = gsub("_Polyketide|_NRPPolyketide|_Other|_NRP", "", label)) %>%
  dplyr::mutate(label_trim = gsub("_1", ".1", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("_", " ", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("WP ", "WP_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("BCRE ", "BCRE_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("STRAU ", "STRAU_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("AMYAL ", "AMYAL_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("Sare ", "Sare_", label_trim)) %>%
  dplyr::mutate(label_trim = gsub("l 2 Amino 4 methoxy trans 3 butenoic acid", "l-2-Amino-4-methoxy-trans-3-butenoic acid", label_trim)) 
metadat$bgc_color2

gdf <- gtr %<+% metadat

pal <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
pal2 <- brewer.pal(8,"Set1")[c(1:3, 5)]
#  "gray40",
pal3 <- c("#E41A1C", "gray40", "#377EB8", pal[c(2, 1)], "navyblue", pal[c(5, 4, 6)], "#8DA0CB", "gray40", "goldenrod4")


pdf("output/40_MT_tree_v4.pdf", width = 12, height = 7)
gdf +
  geom_nodelab(aes(label = label), size = 3) + 
  #size = ifelse(as.numeric(gdf$data$label[!gdf$data$isTip]) > 75, 8, 0.0001),
               #  color = ifelse(as.numeric(gdf$data$label[!gdf$data$isTip]) > 75, "gray50", NA)) +
  geom_tippoint(x = 15, color = gdf$data$bgc_color1[gdf$data$isTip], size = 6) +
  geom_tippoint(x = 15.3, color = ifelse(!is.na(gdf$data$bgc_color2[gdf$data$isTip]), "#BBD374", "white"), size = 6) +
  geom_tiplab(aes(label = label_trim), color = ifelse(grepl("study", gdf$data$label_trim[gdf$data$isTip]), "#EC8F6B", "black"),
              #size = ifelse(grepl("study", gdf$data$label_trim[gdf$data$isTip]), 8, 6),
              fontface = ifelse(grepl("study", gdf$data$label_trim[gdf$data$isTip]), "bold.italic", "italic"), align=T) + 
  xlim_tree(15)
dev.off()

# Make a legend




# Changes
# remove all positions not with 50% or more gaps
# HHBlitz HMM vs HMM

