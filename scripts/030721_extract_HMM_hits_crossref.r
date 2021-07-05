# Load packages
pacman::p_load("tidyverse", "Biostrings", 
               "DECIPHER", "readxl")

# Read in HMM hits
hits <- read_delim("data/HMMsearch_NCcutoff_40MTs.txt", delim = " ", col_names = F)
hits$X9 <- trimws(hits$X9)
accs <- word(hits$X9, sep = "\\|", -1) # 
accs

# Read in the MIBiG sequences
mibig <- readAAStringSet("data/mibig_prot_seqs_2.0.fasta")
inds <- grep(paste0(accs, collapse = "|"), names(mibig))
met_hits <- mibig[inds]
names(met_hits)
names(met_hits)

# Also cross-reference with organism and sequence type
metadat <- read_excel("data/bigfam_results.xlsx") %>%
  dplyr::select(-genome) 
metadat$BGC <- word(metadat$BGC, sep = "\\.gbk", 1)
metadat$taxon <- gsub("\\(Species\\)\\?", "", metadat$taxon)
metadat$taxon <- gsub("\\(Organism\\)\\?", "", metadat$taxon)

# Make unique tree identifiers
names(met_hits)
metdf <- data.frame(bgc = word(names(met_hits), sep = "\\|", 1),
                    nams = names(met_hits),
                    aa = met_hits) %>%
  left_join(., metadat, by = c("bgc" = "BGC")) %>%
  dplyr::mutate(acc = word(nams, sep = "\\|", 5))

# Read in manually-curated literature information about clusters
manual <- read_excel("data/39_MIBiG_Methyltransf_21_lit_dat_manually_curated.xlsx") %>%
  dplyr::filter(include_in_tree == "yes") %>%
  dplyr::filter(!duplicated(acc))

metdf[!metdf$bgc %in% manual$bgc,]
metdf$bgc %in% manual$bgc

merg <- manual %>%
  inner_join(., metdf, keep = F, by = "bgc", suffix = c("", ".y")) %>%
  select_at(
    vars(-ends_with(".y"))
  ) %>%
  dplyr::group_by(acc) %>%
  dplyr::slice(1) %>%
  ungroup()
merg

write_csv(merg, "data/29_MIBiG_hits_manually_curated_to_include.csv")
# Manually identified and replaced the missing taxon
# https://mibig.secondarymetabolites.org/repository/BGC0000769/index.html#r1c1
# BGC0000769
# Mycobacterium avium subsp. hominissuis A5
# Saccharide
# metdf$taxon[is.na(metdf$taxon)] <- "Mycobacterium avium subsp. hominissuis A5"           
# metdf$class[is.na(metdf$class)] <- "Saccharide: unknown" 
# metdf$dataset[is.na(metdf$dataset)] <- "mibig"
# write_csv(metdf, "data/39_MIBiG_Methyltransf_21_lit_dat.csv")
# Make a new AA set with updated names
aaset <- AAStringSet(merg$aa)

# names(aaset) <- paste(merg$bgc, merg$acc, merg$taxon, merg$natural_product, merg$class, sep = "_")
names(aaset) <- paste(merg$acc, merg$taxon, merg$natural_product, merg$class, sep = "_")
names(aaset)

names(aaset) <- gsub(" |-|,|\\.", "_", names(aaset))
names(aaset) <- gsub(":unknown", "", names(aaset))
names(aaset) <- gsub("__|___", "_", names(aaset))
names(aaset) <- gsub("NRP_Polyketide", "NRPPolyketide", names(aaset))

# Now add in the Eremiobacterota sequence
erem <- readAAStringSet("data/Eudoremicrobiaceae_EreM.fasta.txt")
names(erem)
poyE <- readAAStringSet("data/PoyE_AerE_renamed.fasta")

aacomb <- AAStringSet(c(poyE, erem, aaset))
length(aacomb)
names(aacomb)
writeXStringSet(aacomb, "data/32_MT_homologs_with_outgroup.fasta")

# Align using muscle
# muscle -in 32_MT_homologs_with_outgroup.fasta -out 32_MT_homologs_with_outgroup_aligned.afa

# Trim gaps using trimal
#./trimal -in ~/Documents/github/fkbm-bioinformatics/data/32_MT_homologs_with_outgroup_aligned.afa 
# -gt 0.5 -out ~/Documents/github/fkbm-bioinformatics/data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa

hmmal <- readAAStringSet("data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa")
hmmal <- AAStringSet(toupper(hmmal))
BrowseSeqs(hmmal, "32_MTs_homologs_PoyE_outgroup_aligned_50gaps_removed.html")

mthomologs <- readAAStringSet("data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa")
names(mthomologs) <- gsub(" 347 bp", "", names(mthomologs))
BrowseSeqs(mthomologs)
writeXStringSet(mthomologs, "data/32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle.afa")
names(mthomologs)[1]

# Tree properly using Iqtree
# iqtree -s HMMaligned_trimmed_40_MTs.fasta -m MFP 
# Akaike Information Criterion:           VT+F+R6
# Corrected Akaike Information Criterion: VT+F+R5
# Bayesian Information Criterion:         VT+F+R5
# Best-fit model: VT+F+R5 chosen according to BIC

# Now do bootstrapping using best-fit model
# iqtree -s HMMaligned_trimmed_40_MTs_for_bootstrap.fasta -m VT+F+R5 -B 1000 -T 2

# iqtree -s 32_MT_homologs_PoyE_outgroup_aligned_50gaps_removed_muscle_nams_fixed.afa 
# -o AFS60641.1_Candidatus_Entotheonella_factor_polytheonamide_PoyE_polytheonamide -m VT+F+R5 -B 5000 -T 2

# Analysis results written to: 
#   IQ-TREE report:                HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.iqtree
# Maximum-likelihood tree:       HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.treefile
# Likelihood distances:          HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.mldist
# 
# Ultrafast bootstrap approximation results written to:
#   Split support values:          HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.splits.nex
# Consensus tree:                HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.contree
# Screen log file:               HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.log


