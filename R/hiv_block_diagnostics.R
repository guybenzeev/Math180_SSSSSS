############################################################
# hiv_pol_block_diagnostics.R
#
# Goal: investigate why blocks 4 and 5 perform worse
#       than blocks 2 and 3 in the reduced logit model.
#
# Steps:
#   - Rebuild data_model with same predictors as reduced model
#   - Assign block_id (100-codon windows)
#   - For blocks 2–5:
#       * inspect class counts
#       * compare predictor distributions by resistance_site
#       * compute univariate AUC per predictor per block
#       * (optional) plot boxplots by block + resistance_site
############################################################

library(dplyr)
library(pROC)
library(ggplot2)

set.seed(180)

############################################################
## 1. Load and preprocess data ----
############################################################

codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

cat("Columns in codon_features.csv:\n")
print(names(codons))

# Ensure codon_index is numeric
if (!is.numeric(codons$codon_index)) {
  codons$codon_index <- as.numeric(codons$codon_index)
}

# Keep only complete, non-stop codons (same as before)
is_complete_codon <- !is.na(codons$ref_aa) &
  !grepl("-", codons$ref_codon) &
  nchar(codons$ref_codon) == 3 &
  codons$ref_aa != "*"

data <- codons[is_complete_codon, ]

cat("Number of codons after filtering incomplete/stop codons:",
    nrow(data), "\n")

# Outcome as numeric 0/1
if (!"resistance_site" %in% names(data)) {
  stop("Column 'resistance_site' not found in codon_features.csv.")
}
data$resistance_site <- as.numeric(data$resistance_site)

# Factor for aa_chemistry if present
if ("aa_chemistry" %in% names(data)) {
  data$aa_chemistry <- factor(data$aa_chemistry)
}

############################################################
## 2. Reduced predictor set + block IDs ----
############################################################

predictor_vars <- c(
  "codon_degeneracy",
  "syn_frac",
  "shannon_entropy",
  "conservation",
  "aa_chemistry"
)

missing_predictors <- setdiff(predictor_vars, names(data))
if (length(missing_predictors) > 0) {
  stop("Missing predictor columns: ",
       paste(missing_predictors, collapse = ", "))
}

model_cols <- c("codon_index", "resistance_site", predictor_vars)
data_model <- data[, model_cols]
data_model <- data_model[complete.cases(data_model), ]

cat("Number of codons after complete-case filtering:", nrow(data_model), "\n")
cat("Overall resistance_site distribution:\n")
print(table(data_model$resistance_site))

# Define 100-codon blocks (same as before)
min_idx <- min(data_model$codon_index, na.rm = TRUE)
max_idx <- max(data_model$codon_index, na.rm = TRUE)

block_starts <- seq(from = min_idx, to = max_idx, by = 100)
blocks <- data.frame(
  block_id = seq_along(block_starts),
  start    = block_starts,
  end      = pmin(block_starts + 99, max_idx)
)

cat("\nDefined blocks (100-codon windows):\n")
print(blocks)

# Attach block_id to each codon
data_model <- data_model %>%
  mutate(
    block_id = floor((codon_index - min_idx) / 100) + 1
  )

############################################################
## 3. Class counts per block ----
############################################################

cat("\n=== Class counts (resistance_site) per block ===\n")
class_counts <- data_model %>%
  group_by(block_id) %>%
  summarise(
    n_total = n(),
    n_resistance = sum(resistance_site == 1),
    n_nonres    = sum(resistance_site == 0),
    .groups = "drop"
  )

print(class_counts)

# Focus on blocks 2–5 (where we know there are positives)
focus_blocks <- 2:5
cat("\n=== Class counts for blocks 2–5 ===\n")
print(class_counts %>% filter(block_id %in% focus_blocks))

############################################################
## 4. Predictor distributions by block & resistance_site ----
############################################################

cat("\n=== Summary stats by block & resistance_site (blocks 2–5) ===\n")

numeric_predictors <- c("codon_degeneracy", "syn_frac",
                        "shannon_entropy", "conservation")

summary_stats <- data_model %>%
  filter(block_id %in% focus_blocks) %>%
  mutate(resistance_site = factor(resistance_site, levels = c(0, 1))) %>%
  group_by(block_id, resistance_site) %>%
  summarise(
    n = n(),
    across(
      all_of(numeric_predictors),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

print(summary_stats)

############################################################
## 5. Univariate AUC per predictor per block ----
############################################################

cat("\n=== Univariate ROC AUC per predictor per block (blocks 2–5) ===\n")

auc_rows <- list()

for (b in focus_blocks) {
  df_b <- data_model %>% filter(block_id == b)
  
  # Only compute if we have both classes
  if (length(unique(df_b$resistance_site)) < 2) {
    next
  }
  
  for (pred in numeric_predictors) {
    x <- df_b[[pred]]
    y <- df_b$resistance_site
    
    if (all(is.na(x))) next
    
    # pROC expects predictor where larger values -> more likely positive;
    # but if direction is reversed, AUC will be < 0.5 (still informative!)
    roc_obj <- try(pROC::roc(response = y, predictor = x, quiet = TRUE), silent = TRUE)
    if (inherits(roc_obj, "try-error")) next
    
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    auc_rows[[length(auc_rows) + 1]] <- data.frame(
      block_id = b,
      predictor = pred,
      AUC = auc_val
    )
  }
}

if (length(auc_rows) > 0) {
  auc_table <- do.call(rbind, auc_rows)
  auc_table <- auc_table %>%
    arrange(block_id, desc(AUC))
  print(auc_table)
} else {
  cat("No AUCs computed (no blocks had both classes?).\n")
}

############################################################
## 6. (Optional) Boxplots of predictors by resistance & block ----
############################################################

# This part is for visualization; comment out if you don't want plots.

cat("\nGenerating boxplots for numeric predictors by resistance_site, faceted by block_id (2–5)...\n")

for (pred in numeric_predictors) {
  p <- ggplot(
    data_model %>% filter(block_id %in% focus_blocks),
    aes(x = factor(resistance_site), y = .data[[pred]])
  ) +
    geom_boxplot() +
    facet_wrap(~ block_id, scales = "free_y") +
    labs(
      title = paste("Predictor", pred, "by resistance_site, blocks 2–5"),
      x = "Resistance site (0 = no, 1 = yes)",
      y = pred
    )
  
  print(p)
}


############################################################
# End of script
############################################################
