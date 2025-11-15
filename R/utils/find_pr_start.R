## ----- SETTINGS -----
fasta_path <- "Data/hiv-db.fasta" 
pr_motif   <- "PQITL"               # HXB2 protease N-terminus

## ----- PACKAGES -----
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)

## ----- LOAD FASTA -----
seqs <- readDNAStringSet(fasta_path)

# pick HXB2-like reference by name
ref_idx <- grep("HXB2|LAI|IIIB|K03455", names(seqs), ignore.case = TRUE)
if (length(ref_idx) == 0) {
  stop("Couldn't find HXB2 reference in FASTA headers (no 'HXB2', 'LAI', 'IIIB', 'K03455').")
}
ref_idx <- ref_idx[1]
cat("Using reference:", names(seqs)[ref_idx], "\n")

ref_seq <- seqs[[ref_idx]]

## ----- CONVERT TO CHARACTER AND SPLIT INTO CODONS -----
# turn DNAString into a plain character
ref_char <- as.character(ref_seq)
ref_char <- toupper(ref_char)

L <- nchar(ref_char)

# we only take complete codons (multiples of 3)
last_full_base <- L - (L %% 3)
if (L %% 3 != 0) {
  warning("Reference length not divisible by 3; truncating trailing ",
          L %% 3, " base(s).")
}

starts <- seq(1, last_full_base, by = 3)
ends   <- starts + 2

codons <- substring(ref_char, starts, ends)

## ----- TRANSLATE CODONS TO AMINO ACIDS -----
translate_codon <- function(c) {
  # treat any codon with a gap as NA
  if (grepl("-", c)) return(NA_character_)
  as.character(translate(DNAString(c)))
}

aa_vec <- vapply(codons, translate_codon, character(1))

# keep only non-NA AAs to search motif robustly
aa_valid <- aa_vec[!is.na(aa_vec)]
aa_string <- paste0(aa_valid, collapse = "")

## ----- FIND PQITL MOTIF -----
pos_in_aa <- regexpr(pr_motif, aa_string, fixed = TRUE)
if (pos_in_aa[1] == -1) {
  stop("Motif '", pr_motif, "' not found in translated AA sequence.")
}

idx_in_valid <- as.integer(pos_in_aa[1])  # index in aa_valid

# map back to codon index (1-based among all codons)
codon_index_valid <- which(!is.na(aa_vec))
pr_start_idx <- codon_index_valid[idx_in_valid]

cat("\nProtease PR1 (P of PQITL) found at codon_index =", pr_start_idx, "\n")

## ----- SANITY CHECK -----
window_len <- 10
end_idx <- min(length(aa_valid), idx_in_valid + window_len - 1)
aa_window <- aa_valid[idx_in_valid:end_idx]

cat("AA window from PR1:",
    paste(aa_window, collapse = ""), "\n")
