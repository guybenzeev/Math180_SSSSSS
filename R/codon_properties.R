# ==============================================
# codon_properties.R
# Defines codon-level properties:
#   - Amino acid chemistry categories
#   - Codon degeneracy
#   - Synonymous/nonsynonymous mutability
# ==============================================

library(Biostrings)
library(dplyr)
library(tibble)

# -----------------------------
# 1. Amino acid chemistry map
# -----------------------------
aa_chemistry_map <- c(
  A = "hydrophobic", V = "hydrophobic", L = "hydrophobic",
  I = "hydrophobic", M = "hydrophobic",
  F = "aromatic", Y = "aromatic", W = "aromatic",
  S = "polar", T = "polar", N = "polar", Q = "polar",
  K = "positive", R = "positive", H = "positive",
  D = "negative", E = "negative",
  G = "special", P = "special", C = "special",
  "*" = "stop"
)

# -----------------------------
# 2. Codon degeneracy + mutability
# -----------------------------
make_codon_properties <- function() {
  # NO data() call here – just use the object
  gc <- Biostrings::GENETIC_CODE        # named vector: codon -> AA

  all_codons <- names(gc)
  aa_vec <- unname(gc)                  # amino acids in same order
  names(aa_vec) <- all_codons           # keep codon names if helpful

  # degeneracy: how many codons map to each amino acid
  deg_table <- table(aa_vec)
  codon_degeneracy <- as.integer(deg_table[aa_vec])

  nts <- c("A", "C", "G", "T")

  # List single-nt neighbors of a codon
  single_nt_neighbors <- function(codon) {
    codon_chars <- strsplit(codon, "")[[1]]
    neighbors <- character()
    for (pos in 1:3) {
      for (nt in nts) {
        if (nt != codon_chars[pos]) {
          mut <- codon_chars
          mut[pos] <- nt
          neighbors <- c(neighbors, paste0(mut, collapse = ""))
        }
      }
    }
    neighbors
  }

  syn_frac <- numeric(length(all_codons))
  nonsyn_frac <- numeric(length(all_codons))

  for (i in seq_along(all_codons)) {
    codon <- all_codons[i]
    ref_aa <- aa_vec[i]

    neigh <- single_nt_neighbors(codon)
    neigh_aa <- gc[neigh]

    syn   <- sum(neigh_aa == ref_aa)
    nonsyn <- sum(neigh_aa != ref_aa)
    total <- syn + nonsyn

    syn_frac[i]    <- syn / total
    nonsyn_frac[i] <- nonsyn / total
  }

  tibble(
    codon        = all_codons,
    aa           = aa_vec,
    degeneracy   = codon_degeneracy,
    syn_frac     = syn_frac,
    nonsyn_frac  = nonsyn_frac
  )
}

# Precompute global table
codon_props <- make_codon_properties()
