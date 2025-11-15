# ==============================================
# detect_pol_proteins.R  (no width(), robust)
# ==============================================

library(Biostrings)

# -----------------------------
# Helper: choose best frame
# -----------------------------
.choose_best_frame <- function(ref_seq_nogap) {
  # Work with plain character internally to avoid S4 generics issues
  if (!is.character(ref_seq_nogap)) {
    ref_seq_nogap <- as.character(ref_seq_nogap)
  }
  
  frame_info <- vector("list", 3)
  
  for (offset in 0:2) {
    if (nchar(ref_seq_nogap) - offset < 3) {
      frame_info[[offset + 1]] <- list(
        offset = offset,
        aa = AAString(""),
        n_stop = Inf,
        longest_orf = 0
      )
      next
    }
    
    subseq_chr <- substr(ref_seq_nogap, 1 + offset, nchar(ref_seq_nogap))
    aa <- translate(DNAString(subseq_chr))
    aa_str <- as.character(aa)
    aa_chars <- strsplit(aa_str, "")[[1]]
    
    if (length(aa_chars) == 0) {
      frame_info[[offset + 1]] <- list(
        offset = offset,
        aa = aa,
        n_stop = Inf,
        longest_orf = 0
      )
      next
    }
    
    # Count stop codons
    n_stop <- sum(aa_chars == "*")
    
    # Longest non-stop run
    longest_orf <- 0
    current_run <- 0
    for (x in aa_chars) {
      if (x == "*") {
        if (current_run > longest_orf) longest_orf <- current_run
        current_run <- 0
      } else {
        current_run <- current_run + 1
      }
    }
    if (current_run > longest_orf) longest_orf <- current_run
    
    frame_info[[offset + 1]] <- list(
      offset = offset,
      aa = aa,
      n_stop = n_stop,
      longest_orf = longest_orf
    )
  }
  
  # Pick best frame: fewest stops, then longest ORF
  stops <- sapply(frame_info, function(x) x$n_stop)
  best_candidates <- which(stops == min(stops))
  if (length(best_candidates) > 1) {
    orfs <- sapply(frame_info[best_candidates], function(x) x$longest_orf)
    best_idx <- best_candidates[which.max(orfs)]
  } else {
    best_idx <- best_candidates
  }
  
  frame_info[[best_idx]]
}

# -----------------------------
# Helper: motif finder
# -----------------------------
.find_motif <- function(aa_seq, motif, max_mismatch = 0L) {
  # aa_seq: AAString
  # motif: character, e.g. "PQITLWQRPLV"
  hits <- matchPattern(motif, aa_seq, max.mismatch = max_mismatch)
  if (length(hits) == 0) return(NULL)
  starts <- start(hits)
  list(positions = as.integer(starts))
}

# -----------------------------
# Main function
# -----------------------------
detect_pol_proteins <- function(fasta_path,
                                ref_pattern = "HXB2",
                                show_aa_n = 150) {
  
  dna <- readDNAStringSet(fasta_path)
  seq_names <- names(dna)
  
  ref_idx <- which(grepl(ref_pattern, seq_names, ignore.case = TRUE))[1]
  if (is.na(ref_idx))
    stop("Could not find reference matching: ", ref_pattern)
  
  ref_name <- seq_names[ref_idx]
  cat("Using reference:", ref_name, "\n\n")
  
  # Get aligned reference and strip gaps
  ref_seq <- as.character(dna[[ref_idx]])
  ref_seq_nogap <- gsub("-", "", ref_seq)
  
  if (nchar(ref_seq_nogap) < 3)
    stop("Reference sequence (gap-stripped) too short.")
  
  # Choose best frame
  best <- .choose_best_frame(ref_seq_nogap)
  aa_best <- best$aa
  aa_best_str <- as.character(aa_best)
  aa_best_chars <- strsplit(aa_best_str, "")[[1]]
  
  cat("Best frame (1–3):", best$offset + 1, "\n")
  cat("Stop codons in frame:", best$n_stop, "\n")
  cat("Longest ORF:", best$longest_orf, "aa\n\n")
  
  # Show first N amino acids
  n_show <- min(show_aa_n, length(aa_best_chars))
  cat("First", n_show, "AAs (best frame):\n")
  if (n_show > 0) {
    cat(paste0(aa_best_chars[seq_len(n_show)], collapse = ""), "\n\n")
  } else {
    cat("[no amino acids translated]\n\n")
  }
  
  aaA <- AAString(aa_best_str)
  
cat("Protein motif detection:\n")

# --- PR (Protease) ---
pr_motif <- "PQITLWQRPLV"   # canonical PR motif
pr_hit <- .find_motif(aaA, pr_motif, max_mismatch = 1L)
if (is.null(pr_hit)) {
  cat("  PR (Protease; motif", pr_motif, "): NOT FOUND\n")
} else {
  cat("  PR (Protease; motif", pr_motif, 
      ", ≤1 mismatch): FOUND at positions:",
      paste(pr_hit$positions, collapse = ", "), "\n")
}

# --- RT (Reverse Transcriptase YMDD) ---
rt_motif <- "YMDD"
rt_hit <- .find_motif(aaA, rt_motif, max_mismatch = 0L)
if (is.null(rt_hit)) {
  cat("  RT (Reverse transcriptase; motif", rt_motif, "): NOT FOUND\n")
} else {
  cat("  RT (Reverse transcriptase; motif", rt_motif, 
      "): FOUND at positions:",
      paste(rt_hit$positions, collapse = ", "), "\n")
}

# --- IN (Integrase HHCC) ---
in_motif <- "HHCC"
in_hit <- .find_motif(aaA, in_motif, max_mismatch = 1L)
if (is.null(in_hit)) {
  cat("  IN (Integrase; motif", in_motif, ", ≤1 mismatch): NOT FOUND\n")
} else {
  cat("  IN (Integrase; motif", in_motif, 
      ", ≤1 mismatch): FOUND at positions:",
      paste(in_hit$positions, collapse = ", "), "\n")
}

  invisible(list(
    ref_name = ref_name,
    best_frame = best$offset + 1,
    aa = aa_best_str,
    PR = pr_hit,
    RT = rt_hit,
    IN = in_hit
  ))
}
