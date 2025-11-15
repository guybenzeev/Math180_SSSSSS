# ==============================================
# run_analysis.R
# Complete workflow:
#   1. Detect pol frame + PR/RT motifs
#   2. Build codon table from alignment
#   3. Map PR + RT resistance AA positions (HXB2 numbering)
#        -> pol AA indices
#        -> codon_index in codon_table
#   4. Add resistance labels and write CSV
# ==============================================

library(dplyr)
library(Biostrings)

# Your existing helpers
source("R/preprocess_alignment.R")        # defines build_codon_table()
source("R/utils/detect_pol_proteins.R")   # prints PR/RT motifs, returns best_frame + positions
source("R/utils/pol_mapping.R")           # map_aa_to_codon_index(), etc.

fasta_path <- "Data/hiv-db.fasta"

# -----------------------
# 1. Detect pol frame + motifs
# -----------------------
det_res <- detect_pol_proteins(fasta_path)
# det_res$best_frame   : 1–3
# det_res$PR$positions : position of PR motif PQITLWQRPLV in pol AA indexing
# det_res$RT$positions : position of RT motif YMDD in pol AA indexing

if (is.null(det_res$RT) || length(det_res$RT$positions) == 0) {
  stop("Could not find RT YMDD motif in reference – cannot anchor RT numbering.")
}
if (is.null(det_res$PR) || length(det_res$PR$positions) == 0) {
  stop("Could not find PR motif PQITLWQRPLV in reference – cannot anchor PR numbering.")
}

# --- RT anchor (YMDD = M184 in RT) ---
ymdd_pol_aa <- det_res$RT$positions[1]  # AA index in pol of YMDD motif
ymdd_rt_pos <- 184L                     # YMDD M184 in RT numbering
rt_offset    <- ymdd_pol_aa - ymdd_rt_pos

# --- PR anchor (PQITLWQRPLV starts at PR position 1 in protease) ---
pr_motif_pol_aa <- det_res$PR$positions[1]  # AA index in pol of the P in PQITLWQRPLV
pr_motif_pr_pos <- 1L                       # PR N-terminus is position 1
pr_offset       <- pr_motif_pol_aa - pr_motif_pr_pos

cat("Anchors:\n")
cat("  PR motif at pol AA", pr_motif_pol_aa, "≙ PR position", pr_motif_pr_pos, "\n")
cat("  RT YMDD at pol AA", ymdd_pol_aa, "≙ RT position", ymdd_rt_pos, "\n\n")

# -----------------------
# 2. Build codon table from alignment
# -----------------------
codon_table <- build_codon_table(fasta_path)

# Map pol AA index -> codon_index (your indexing)
aa_to_codon <- map_aa_to_codon_index(
  fasta_path = fasta_path,
  best_frame = det_res$best_frame
)

# -----------------------
# 3. Define PR + RT resistance positions (HXB2 protein numbering)
# -----------------------

# RT resistance positions (HXB2 RT numbering) – your original list:
rt_resistance_positions <- c(41, 65, 67, 69, 70, 74, 100, 101, 103, 106, 115, 138, 151, 181, 184, 188, 190, 210, 215, 219, 230)

# Example PR resistance positions (HXB2 PR numbering).
# You can tweak this list according to the exact set you want (HIVDB/IAS).
pr_resistance_positions <- c(
  30, 32, 33, 46, 47, 48, 50, 54,
  76, 82, 84, 88, 90
)

# --- Convert RT positions -> pol AA indices ---
rt_pol_aa_positions <- rt_offset + rt_resistance_positions

# --- Convert PR positions -> pol AA indices ---
pr_pol_aa_positions <- pr_offset + pr_resistance_positions

# Clip to the range covered by aa_to_codon
valid_rt <- rt_pol_aa_positions >= 1 & rt_pol_aa_positions <= length(aa_to_codon)
valid_pr <- pr_pol_aa_positions >= 1 & pr_pol_aa_positions <= length(aa_to_codon)

rt_pol_aa_positions <- rt_pol_aa_positions[valid_rt]
pr_pol_aa_positions <- pr_pol_aa_positions[valid_pr]

rt_resistance_positions <- rt_resistance_positions[valid_rt]
pr_resistance_positions <- pr_resistance_positions[valid_pr]

# --- Map pol AA -> codon_index ---
rt_codons <- aa_to_codon[rt_pol_aa_positions]
pr_codons <- aa_to_codon[pr_pol_aa_positions]

# Drop NA and deduplicate
rt_codons <- sort(unique(rt_codons[!is.na(rt_codons)]))
pr_codons <- sort(unique(pr_codons[!is.na(pr_codons)]))

# Combined resistance codons
resistance_codons <- sort(unique(c(rt_codons, pr_codons)))

cat("RT resistance mapping (RT numbering -> pol AA -> codon_index):\n")
print(
  tibble(
    RT_pos    = rt_resistance_positions,
    pol_aa    = rt_pol_aa_positions,
    codon_idx = aa_to_codon[rt_pol_aa_positions]
  )
)

cat("\nPR resistance mapping (PR numbering -> pol AA -> codon_index):\n")
print(
  tibble(
    PR_pos    = pr_resistance_positions,
    pol_aa    = pr_pol_aa_positions,
    codon_idx = aa_to_codon[pr_pol_aa_positions]
  )
)

# -----------------------
# 4. Add resistance labels and save
# -----------------------
codon_table <- codon_table %>%
  mutate(
    is_RT_resistance_site = as.integer(codon_index %in% rt_codons),
    is_PR_resistance_site = as.integer(codon_index %in% pr_codons),
    resistance_site       = as.integer(codon_index %in% resistance_codons)
  )

write.csv(codon_table, "Results/codon_features.csv", row.names = FALSE)
cat("\nWrote codon features with PR + RT resistance annotations to Results/codon_features.csv\n")
