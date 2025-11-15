# ==============================================
# map_pol_to_codons.R
# Utilities to map AA positions (in HXB2 pol)
# → nucleotide positions (ungapped)
# → alignment columns
# → codon_index in build_codon_table()
# ==============================================

library(Biostrings)
library(dplyr)
library(tibble)

# ---- 1. Build mapping: ungapped nt index -> alignment column ----
.build_ungapped_nt_to_align_col <- function(ref_seq_gappy) {
  ref_chars <- strsplit(ref_seq_gappy, "")[[1]]
  aln_len <- length(ref_chars)
  
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
  
  # now invert: for each ungapped index u, which column c?
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
    col_for_u = col_for_u
  )
}

# ---- 2. Compute AA→codon_index map for the best frame ----
map_aa_to_codon_index <- function(fasta_path, best_frame) {
  dna <- readDNAStringSet(fasta_path)
  ref_seq_gappy <- as.character(dna[[1]])  # HXB2 in the alignment
  aln_len <- nchar(ref_seq_gappy)
  
  # strip gaps (same as detect_pol_proteins)
  ref_seq_nogap <- gsub("-", "", ref_seq_gappy)
  u_max <- nchar(ref_seq_nogap)
  
  # start position in ungapped coords for this frame
  # best_frame is 1,2,3; offset is 0,1,2
  offset <- best_frame - 1L
  start_u <- 1L + offset
  
  # translate to know how many AA we actually have in this frame
  subseq_chr <- substr(ref_seq_nogap, start_u, u_max)
  aa <- translate(DNAString(subseq_chr))
  aa_str <- as.character(aa)
  L_aa <- nchar(aa_str)
  
  # mapping from ungapped nt index → alignment column
  map_nt <- .build_ungapped_nt_to_align_col(ref_seq_gappy)
  col_for_u <- map_nt$col_for_u
  
  # mapping from alignment column → codon_index
  # must match build_codon_table(): codon 1 = cols 1:3, codon 2 = 4:6, etc
  col2codon <- ceiling(seq_len(aln_len) / 3)
  
  # AA index → codon_index (NA where mapping broken)
  aa_to_codon <- integer(L_aa)
  aa_to_codon[] <- NA_integer_
  
  for (i in seq_len(L_aa)) {
    # ungapped nt indices for this AA in this frame
    u1 <- start_u + 3L * (i - 1L)
    u2 <- u1 + 1L
    u3 <- u1 + 2L
    
    # safety: if beyond available nts, stop
    if (u3 > u_max) break
    
    # map to alignment columns
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
      # should not really happen, but be defensive
      aa_to_codon[i] <- codon_idxs[1L]
    }
  }
  
  aa_to_codon
}

# ---- 3. Convenience: map motif hits to codon indices ----
map_motif_hits_to_codons <- function(fasta_path, detect_res) {
  aa_to_codon <- map_aa_to_codon_index(
    fasta_path = fasta_path,
    best_frame = detect_res$best_frame
  )
  
  tibs <- list()
  
  if (!is.null(detect_res$PR)) {
    tibs[["PR"]] <- tibble(
      region = "PR",
      aa_pos = detect_res$PR$positions,
      codon_index = aa_to_codon[detect_res$PR$positions]
    )
  }
  
  if (!is.null(detect_res$RT)) {
    tibs[["RT"]] <- tibble(
      region = "RT",
      aa_pos = detect_res$RT$positions,
      codon_index = aa_to_codon[detect_res$RT$positions]
    )
  }
  
  if (!is.null(detect_res$IN)) {
    tibs[["IN"]] <- tibble(
      region = "IN",
      aa_pos = detect_res$IN$positions,
      codon_index = aa_to_codon[detect_res$IN$positions]
    )
  }
  
  if (length(tibs) == 0L) {
    return(tibble(region = character(), aa_pos = integer(), codon_index = integer()))
  }
  
  bind_rows(tibs)
}
