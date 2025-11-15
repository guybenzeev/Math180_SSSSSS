# ==============================================
# preprocess_alignment.R
# Reads aligned FASTA file and constructs:
#   - codon-level sequence features
#   - evolutionary features (entropy, conservation, dNdS proxy)
# Outputs: tibble codon_table
# ==============================================

library(Biostrings)
library(dplyr)
library(tidyr)
library(tibble)

source("R/codon_properties.R")   # brings in codon_props + aa_chemistry_map

# --------------------------------------
# Helper: translate codons (with gaps)
# --------------------------------------
codon_to_aa <- function(codon) {
  if (grepl("-", codon) || nchar(codon) != 3) return(NA_character_)
  if (!codon %in% names(Biostrings::GENETIC_CODE
      )) return(NA_character_)
  unname(Biostrings::GENETIC_CODE
         [[codon]])
}

# --------------------------------------
# Main function: FASTA → feature table
# --------------------------------------
build_codon_table <- function(fasta_path) {
  dna <- readDNAStringSet(fasta_path)
  seq_names <- names(dna)
  ref_seq <- as.character(dna[[1]])  # first sequence = HXB2
  aln_len <- nchar(ref_seq)

  if (aln_len %% 3 != 0)
    stop("Alignment length is not divisible by 3")

  n_codons <- aln_len / 3

  # Convert alignment to matrix
  seq_mat <- do.call(rbind, lapply(as.character(dna),
                                   function(s) strsplit(s, "")[[1]]))
  rownames(seq_mat) <- seq_names

  # Storage
  out <- tibble(
    codon_index = seq_len(n_codons),
    ref_codon = NA_character_,
    ref_aa = NA_character_,
    aa_chemistry = NA_character_,
    codon_degeneracy = NA_real_,
    syn_frac = NA_real_,
    nonsyn_frac = NA_real_,
    shannon_entropy = NA_real_,
    conservation = NA_real_,
    dNdS_proxy = NA_real_
  )

  for (k in seq_len(n_codons)) {
    cols <- ((k - 1) * 3 + 1):(k * 3)
    ref_codon <- paste0(seq_mat[1, cols], collapse = "")
    out$ref_codon[k] <- ref_codon

    valid_ref <- !grepl("-", ref_codon) && ref_codon %in% codon_props$codon

    if (valid_ref) {
      ref_aa <- codon_to_aa(ref_codon)
      out$ref_aa[k] <- ref_aa

      cp <- codon_props[codon_props$codon == ref_codon, ]
      out$codon_degeneracy[k] <- cp$degeneracy
      out$syn_frac[k] <- cp$syn_frac
      out$nonsyn_frac[k] <- cp$nonsyn_frac
      out$aa_chemistry[k] <- aa_chemistry_map[ref_aa]
    }

    codons_all <- character()
    aas_all <- character()

    for (j in seq_len(nrow(seq_mat))) {
      codon_j <- paste0(seq_mat[j, cols], collapse = "")
      if (grepl("-", codon_j)) next
      if (!(codon_j %in% names(Biostrings::GENETIC_CODE
            ))) next
      aa_j <- codon_to_aa(codon_j)
      if (is.na(aa_j)) next

      codons_all <- c(codons_all, codon_j)
      aas_all <- c(aas_all, aa_j)
    }

    if (length(aas_all) == 0 || !valid_ref) next

    tab <- table(aas_all)
    p <- tab / sum(tab)
    out$shannon_entropy[k] <- -sum(p * log2(p))

    out$conservation[k] <- mean(aas_all == out$ref_aa[k])

    syn <- sum(aas_all == ref_aa)
    nonsyn <- sum(aas_all != ref_aa)
    out$dNdS_proxy[k] <- (nonsyn + 0.5) / (syn + 0.5)
  }
  out
}
