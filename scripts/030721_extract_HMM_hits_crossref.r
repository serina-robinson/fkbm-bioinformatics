# Load packages
pacman::p_load("tidyverse", "Biostrings", 
               "DECIPHER", "readxl")

# Read in HMM hits
hits <- read_delim("data/MiBIG_Methyltransf_2.0_hits.txt", delim = " ", col_names = F)
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
metdf$acc
# Manually identified and replaced the missing taxon
# https://mibig.secondarymetabolites.org/repository/BGC0000769/index.html#r1c1
# BGC0000769
# Mycobacterium avium subsp. hominissuis A5
# Saccharide
metdf$taxon[is.na(metdf$taxon)] <- "Mycobacterium avium subsp. hominissuis A5"           
metdf$class[is.na(metdf$class)] <- "Saccharide"              


# Make a new AA set with updated names
aaset <- AAStringSet(metdf$aa)

names(aaset) <- paste(metdf$bgc, metdf$acc, metdf$taxon, metdf$class, sep = "_")
names(aaset) <- gsub(" |-|,|\\.", "_", names(aaset))
names(aaset) <- gsub(":unknown", "", names(aaset))
names(aaset) <- gsub("__|___", "_", names(aaset))
names(aaset)

# Now add in the Eremiobacterota sequence
erem <- readAAStringSet("data/Eudoremicrobiaceae_EreM.fasta.txt")

aacomb <- AAStringSet(c(erem, aaset))
writeXStringSet(aacomb, "data/40_methyltransferase_homologs.fasta")

# Tried alignment using DECIPHER
# width(aacomb)
# dec <- AlignSeqs(aacomb)
# subs <- AAStringSet(substr(dec, start = 320, stop = 860))
# width(subs)
# BrowseSeqs(subs, "data/40_MTs_aligned_trimmed.html") # need to find a local alignment or extract
# BrowseSeqs(dec, "data/40_MTs_aligned.html")

# Instead decided to align using HMMAlign
hmmal <- readAAStringSet("data/HMMaligned_sto_to_fasta_40_MTs.fasta")
hmmal <- AAStringSet(toupper(hmmal))
BrowseSeqs(hmmal, "40_MTs_hmmaligned.html")

# Subset the alignment to only include EreM domain
subhmm <- AAStringSet(substr(hmmal, start = 2450, stop = 2870))
BrowseSeqs(subhmm, "data/HMMaligned_trimmed_40_MTs.html")
writeXStringSet(subhmm, "data/HMMaligned_trimmed_40_MTs_for_bootstrap.fasta")

# Tree properly using Iqtree
# iqtree -s HMMaligned_trimmed_40_MTs.fasta -m MFP 
# Akaike Information Criterion:           VT+F+R6
# Corrected Akaike Information Criterion: VT+F+R5
# Bayesian Information Criterion:         VT+F+R5
# Best-fit model: VT+F+R5 chosen according to BIC

# Now do bootstrapping using best-fit model
# iqtree -s HMMaligned_trimmed_40_MTs_for_bootstrap.fasta -m VT+F+R5 -B 1000 -T 2

# Analysis results written to: 
#   IQ-TREE report:                HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.iqtree
# Maximum-likelihood tree:       HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.treefile
# Likelihood distances:          HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.mldist
# 
# Ultrafast bootstrap approximation results written to:
#   Split support values:          HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.splits.nex
# Consensus tree:                HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.contree
# Screen log file:               HMMaligned_trimmed_40_MTs_for_bootstrap.fasta.log

# Read in the tree 
