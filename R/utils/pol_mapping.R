# ==============================================
# pol_mapping.R
# - build_codon_table(): alignment -> codon feature table
# - map_aa_to_codon_index(): AA index (in HXB2 pol) -> codon_index
# - map_motif_hits_to_codons(): PR/RT/IN motif AA positions -> codon_index
# ==============================================

library(Biostrings)
library(dplyr)
library(tidyr)
library(tibble)

source("R/codon_properties.R")   # codon_props + aa_chemistry_map
source("R/utils/detect_pol_proteins.R")  # your working detect_pol_proteins()

# --------------------------------------
# Helper: translate codons (with gaps)
# --------------------------------------
codon_to_aa <- function(codon) {
  if (grepl("-", codon) || nchar(codon) != 3) return(NA_character_)
  if (!codon %in% names(Biostrings::GENETIC_CODE)) return(NA_character_)
  unname(Biostrings::GENETIC_CODE[[codon]])
}

# --------------------------------------
# 1. Main: FASTA → codon-level feature table
# --------------------------------------
build_codon_table <- function(fasta_path) {
  dna <- readDNAStringSet(fasta_path)
  seq_names <- names(dna)
  
  # first sequence = HXB2 (aligned, may contain '-')
  ref_seq <- as.character(dna[[1]])
  aln_len <- nchar(ref_seq)

  if (aln_len %% 3 != 0)
    stop("Alignment length is not divisible by 3")

  n_codons <- aln_len / 3

  # alignment matrix: rows = sequences, cols = alignment positions
  seq_mat <- do.call(
    rbind,
    lapply(as.character(dna), function(s) strsplit(s, "")[[1]])
  )
  rownames(seq_mat) <- seq_names

  out <- tibble(
    codon_index       = seq_len(n_codons),
    ref_codon         = NA_character_,
    ref_aa            = NA_character_,
    aa_chemistry      = NA_character_,
    codon_degeneracy  = NA_real_,
    syn_frac          = NA_real_,
    nonsyn_frac       = NA_real_,
    shannon_entropy   = NA_real_,
    conservation      = NA_real_,
    dNdS_proxy        = NA_real_
  )

  for (k in seq_len(n_codons)) {
    cols <- ((k - 1) * 3 + 1):(k * 3)
    ref_codon <- paste0(seq_mat[1, cols], collapse = "")
    out$ref_codon[k] <- ref_codon

    valid_ref <- !grepl("-", ref_codon) &&
      ref_codon %in% codon_props$codon

    if (valid_ref) {
      ref_aa <- codon_to_aa(ref_codon)
      out$ref_aa[k] <- ref_aa

      cp <- codon_props[codon_props$codon == ref_codon, ]
      out$codon_degeneracy[k] <- cp$degeneracy
      out$syn_frac[k]         <- cp$syn_frac
      out$nonsyn_frac[k]      <- cp$nonsyn_frac
      out$aa_chemistry[k]     <- aa_chemistry_map[ref_aa]
    }

    codons_all <- character()
    aas_all    <- character()

    for (j in seq_len(nrow(seq_mat))) {
      codon_j <- paste0(seq_mat[j, cols], collapse = "")
      if (grepl("-", codon_j)) next
      if (!(codon_j %in% names(Biostrings::GENETIC_CODE))) next

      aa_j <- codon_to_aa(codon_j)
      if (is.na(aa_j)) next

      codons_all <- c(codons_all, codon_j)
      aas_all    <- c(aas_all, aa_j)
    }

    if (length(aas_all) == 0 || !valid_ref) next

    tab <- table(aas_all)
    p <- tab / sum(tab)
    out$shannon_entropy[k] <- -sum(p * log2(p))

    out$conservation[k] <- mean(aas_all == out$ref_aa[k])

    syn    <- sum(aas_all == ref_aa)
    nonsyn <- sum(aas_all != ref_aa)
    out$dNdS_proxy[k] <- (nonsyn + 0.5) / (syn + 0.5)
  }

  out
}

# --------------------------------------
# 2. Build mapping: ungapped nt index -> alignment column
# --------------------------------------
.build_ungapped_nt_to_align_col <- function(ref_seq_gappy) {
  ref_chars <- strsplit(ref_seq_gappy, "")[[1]]
  aln_len   <- length(ref_chars)
  
  # ungapped index for each alignment column (NA if gap)
  ungap_idx_for_col <- integer(aln_len)
  ungap_idx_for_col[] <- NA_integer_
  
  u <- 0L
  for (c in seq_len(aln_len)) {
    if (ref_chars[c] != "-") {
      u <- u + 1L
      ungap_idx_for_col[c] <- u
    }
  }
  
  # invert: ungapped index -> alignment column
  max_u <- max(ungap_idx_for_col, na.rm = TRUE)
  col_for_u <- integer(max_u)
  col_for_u[] <- NA_integer_
  
  for (c in seq_len(aln_len)) {
    u <- ungap_idx_for_col[c]
    if (!is.na(u)) {
      col_for_u[u] <- c
    }
  }
  
  list(
    ungap_idx_for_col = ungap_idx_for_col,
    col_for_u         = col_for_u
  )
}

# --------------------------------------
# 3. AA index (in HXB2 pol frame) -> codon_index
# --------------------------------------
map_aa_to_codon_index <- function(fasta_path, best_frame) {
  dna <- readDNAStringSet(fasta_path)
  ref_seq_gappy <- as.character(dna[[1]])  # HXB2 row in alignment
  aln_len <- nchar(ref_seq_gappy)
  
  # strip gaps (same as in detect_pol_proteins)
  ref_seq_nogap <- gsub("-", "", ref_seq_gappy)
  u_max <- nchar(ref_seq_nogap)
  
  # frame: 1,2,3 -> offset 0,1,2
  offset <- best_frame - 1L
  start_u <- 1L + offset
  
  # translate that ungapped sequence in this frame
  subseq_chr <- substr(ref_seq_nogap, start_u, u_max)
  aa <- translate(DNAString(subseq_chr))
  aa_str <- as.character(aa)
  L_aa <- nchar(aa_str)
  
  # mapping from ungapped nt index -> alignment column
  nt_map <- .build_ungapped_nt_to_align_col(ref_seq_gappy)
  col_for_u <- nt_map$col_for_u
  
  # mapping from alignment column -> codon_index
  # must match build_codon_table(): codon 1 = cols 1:3, codon 2 = 4:6, etc.
  col2codon <- ceiling(seq_len(aln_len) / 3)
  
  # AA index -> codon_index
  aa_to_codon <- integer(L_aa)
  aa_to_codon[] <- NA_integer_
  
  for (i in seq_len(L_aa)) {
    # ungapped nt indices for this AA in this frame
    u1 <- start_u + 3L * (i - 1L)
    u2 <- u1 + 1L
    u3 <- u1 + 2L
    
    if (u3 > u_max) break  # ran off end of ungapped sequence
    
    c1 <- col_for_u[u1]
    c2 <- col_for_u[u2]
    c3 <- col_for_u[u3]
    
    if (any(is.na(c(c1, c2, c3)))) {
      aa_to_codon[i] <- NA_integer_
      next
    }
    
    codon_idxs <- unique(col2codon[c(c1, c2, c3)])
    if (length(codon_idxs) == 1L) {
      aa_to_codon[i] <- codon_idxs
    } else {
      aa_to_codon[i] <- codon_idxs[1L]  # defensive
    }
  }
  
  aa_to_codon
}

# --------------------------------------
# 4. Convenience: map motif hits (PR/RT/IN) -> codon_index
# --------------------------------------
map_motif_hits_to_codons <- function(fasta_path, detect_res) {
  aa_to_codon <- map_aa_to_codon_index(
    fasta_path = fasta_path,
    best_frame = detect_res$best_frame
  )
  
  tibs <- list()
  
  if (!is.null(detect_res$PR)) {
    tibs[["PR"]] <- tibble(
      region      = "PR",
      aa_pos      = detect_res$PR$positions,
      codon_index = aa_to_codon[detect_res$PR$positions]
    )
  }
  
  if (!is.null(detect_res$RT)) {
    tibs[["RT"]] <- tibble(
      region      = "RT",
      aa_pos      = detect_res$RT$positions,
      codon_index = aa_to_codon[detect_res$RT$positions]
    )
  }
  
  if (!is.null(detect_res$IN)) {
    tibs[["IN"]] <- tibble(
      region      = "IN",
      aa_pos      = detect_res$IN$positions,
      codon_index = aa_to_codon[detect_res$IN$positions]
    )
  }
  
  if (length(tibs) == 0L) {
    return(tibble(region = character(), aa_pos = integer(), codon_index = integer()))
  }
  
  bind_rows(tibs)
}
