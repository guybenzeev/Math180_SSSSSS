# ==============================================
# run_analysis.R
# Complete workflow:
#   1. Build codon table
#   2. Add resistance labels
#   3. Fit logistic regression model
# ==============================================

library(dplyr)
source("R/preprocess_alignment.R")

# -----------------------
# 1. Build codon table
# -----------------------
codon_table <- build_codon_table("Data/hiv-db.fasta")

# -----------------------
# 2. Resistance label list
# -----------------------

# Currently a list placeholder values, replace with actual resistance codon indices
resistance_codons <- c(41, 65, 70, 74, 103, 106, 151, 184, 215)

codon_table <- codon_table %>%
  mutate(resistance_site = as.integer(codon_index %in% resistance_codons))

# -----------------------
# 3. Logistic Regression
# -----------------------

# Save outputs
write.csv(codon_table, "Results/codon_features.csv", row.names = FALSE)

