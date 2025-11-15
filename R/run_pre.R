source("R/utils/detect_pol_proteins.R")
source("R/preprocess_alignment.R")    # <-- this defines build_codon_table()
source("R/utils/pol_mapping.R")       # <-- the new mapping helpers

fasta_path <- "Data/hiv-db.fasta"

# 1. Detect motifs and best frame
det_res <- detect_pol_proteins(fasta_path)

# 2. Build codon table from the same alignment
codon_table <- build_codon_table(fasta_path)

# 3. Map motif AA positions (PR, RT, IN) to codon_index
motif_codons <- map_motif_hits_to_codons(fasta_path, det_res)
motif_codons