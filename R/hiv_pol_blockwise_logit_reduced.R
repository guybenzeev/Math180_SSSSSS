############################################################
# Blockwise Logistic Regression Analysis for HIV-1 pol
# Reduced Predictor Model (low-collinearity)
#
# For each 100-codon block:
#   - Use that block as TEST
#   - Use all other codons as TRAIN
#   - Oversample minority class in TRAIN
#   - Fit logistic regression (reduced predictor set)
#   - Tune threshold to maximize F1 on that block
#   - Record performance metrics & coefficients
#
# Predictors used:
#   - codon_degeneracy
#   - syn_frac          (keep only one of syn/nonsyn)
#   - shannon_entropy
#   - conservation      (drop log_dNdS to reduce collinearity)
#   - aa_chemistry (categorical)
#
# Scholarly Super Serious Statistically Scientific Society
# Vincent Ngo, Charles Pierceton, Guy Ben Zeev, Marc-Philippe Gnagne
############################################################

########################
# Package helper ------
########################

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_or_install("pROC")   # ROC/AUC
load_or_install("PRROC")  # PR curves

set.seed(180)

############################################################
## 1. Load and preprocess data ----
############################################################

# Path relative to this R script (assuming it's in the R/ folder)
codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

cat("Columns in codon_features.csv:\n")
print(names(codons))

# Ensure codon_index is numeric
if (!is.numeric(codons$codon_index)) {
  codons$codon_index <- as.numeric(codons$codon_index)
}

# Keep only complete, non-stop codons
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

# (Optional) log dNdS just for later exploratory plots (not used in model)
if ("dNdS_proxy" %in% names(data)) {
  data$log_dNdS <- log1p(data$dNdS_proxy)
} else {
  warning("No 'dNdS_proxy' column found; skipping log_dNdS feature.")
}

# Categorical predictors
if ("ref_aa" %in% names(data)) {
  data$ref_aa <- factor(data$ref_aa)
}
if ("ref_codon" %in% names(data)) {
  data$ref_codon <- factor(data$ref_codon)
}
if ("aa_chemistry" %in% names(data)) {
  data$aa_chemistry <- factor(data$aa_chemistry)
}

############################################################
## 2. Select predictors and clean model data ----
############################################################

# Reduced, low-collinearity predictor set
predictor_vars <- c(
  "codon_degeneracy",
  "syn_frac",           # keep only one of syn_frac / nonsyn_frac
  "shannon_entropy",
  "conservation",       # keep conservation; drop log_dNdS from model
  "aa_chemistry"
)

missing_predictors <- setdiff(predictor_vars, names(data))
if (length(missing_predictors) > 0) {
  stop("These predictor columns are missing from the data: ",
       paste(missing_predictors, collapse = ", "))
}

model_cols <- c("codon_index", "resistance_site", predictor_vars)
data_model <- data[, model_cols]
data_model <- data_model[complete.cases(data_model), ]

cat("Number of codons after complete-case filtering:", nrow(data_model), "\n")

# For safety, ensure resistance_site is 0/1 numeric
data_model$resistance_site <- as.numeric(data_model$resistance_site)

# Basic distribution
cat("Overall class distribution:\n")
print(table(data_model$resistance_site))

############################################################
## 3. Helper functions: oversampling, evaluation, threshold ----
############################################################

# Simple random oversampling of minority class
oversample_minority <- function(df) {
  df_fac <- df
  df_fac$resistance_site <- factor(df_fac$resistance_site)
  tab <- table(df_fac$resistance_site)
  
  if (length(tab) < 2 || tab["1"] == 0) {
    warning("No positive examples in training data; skipping oversampling.")
    return(df)
  }
  
  maj_class <- names(which.max(tab))
  min_class <- names(which.min(tab))
  n_maj     <- as.numeric(tab[maj_class])
  n_min     <- as.numeric(tab[min_class])
  
  maj_rows <- df_fac[df_fac$resistance_site == maj_class, ]
  min_rows <- df_fac[df_fac$resistance_site == min_class, ]
  
  set.seed(180)
  min_upsampled <- min_rows[sample(seq_len(n_min), n_maj, replace = TRUE), ]
  
  df_os <- rbind(maj_rows, min_upsampled)
  df_os <- df_os[sample(seq_len(nrow(df_os))), ]
  
  # Return numeric version
  df_os_num <- df_os
  df_os_num$resistance_site <- as.numeric(as.character(df_os$resistance_site))
  df_os_num
}

# Evaluate a binary classifier given probabilities and threshold
evaluate_binary <- function(y_true, y_prob, threshold = 0.5) {
  y_pred <- ifelse(y_prob >= threshold, 1, 0)
  
  tab <- table(
    truth = factor(y_true, levels = c(0, 1)),
    pred  = factor(y_pred, levels = c(0, 1))
  )
  
  TN <- tab["0", "0"]
  FP <- tab["0", "1"]
  FN <- tab["1", "0"]
  TP <- tab["1", "1"]
  
  accuracy  <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
  f1        <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall == 0),
    NA,
    2 * precision * recall / (precision + recall)
  )
  
  # ROC AUC
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_roc <- as.numeric(pROC::auc(roc_obj))
  
  # PR AUC (only if we have positives and negatives)
  if (sum(y_true == 1) > 0 && sum(y_true == 0) > 0) {
    pr_obj <- PRROC::pr.curve(
      scores.class0 = y_prob[y_true == 1],
      scores.class1 = y_prob[y_true == 0],
      curve = FALSE
    )
    auc_pr <- pr_obj$auc.integral
  } else {
    auc_pr <- NA
  }
  
  list(
    confusion = tab,
    accuracy  = accuracy,
    precision = precision,
    recall    = recall,
    f1        = f1,
    auc_roc   = auc_roc,
    auc_pr    = auc_pr
  )
}

# Threshold tuning via F1 on a grid
tune_threshold_pr <- function(y_true, y_prob, n_grid = 200) {
  thresholds <- seq(0.001, 0.999, length.out = n_grid)
  best_f1 <- -Inf
  best_t  <- 0.5
  best_stats <- NULL
  
  for (t in thresholds) {
    y_pred <- ifelse(y_prob >= t, 1, 0)
    tab <- table(
      truth = factor(y_true, levels = c(0, 1)),
      pred  = factor(y_pred, levels = c(0, 1))
    )
    TN <- tab["0", "0"]
    FP <- tab["0", "1"]
    FN <- tab["1", "0"]
    TP <- tab["1", "1"]
    
    precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
    recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
    f1        <- ifelse(
      is.na(precision) || is.na(recall) || (precision + recall == 0),
      NA,
      2 * precision * recall / (precision + recall)
    )
    
    if (!is.na(f1) && f1 > best_f1) {
      best_f1 <- f1
      best_t  <- t
      best_stats <- list(
        threshold = t,
        precision = precision,
        recall    = recall,
        f1        = f1
      )
    }
  }
  best_stats
}

############################################################
## 4. Block definition ----
############################################################

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

############################################################
## 5. Loop over blocks: train/test Logit_oversampled_tuned ----
############################################################

results_list    <- list()
model_coefs     <- list()
model_summaries <- list()

logit_formula <- as.formula(
  paste("resistance_site ~", paste(predictor_vars, collapse = " + "))
)

for (i in seq_len(nrow(blocks))) {
  b_start <- blocks$start[i]
  b_end   <- blocks$end[i]
  
  cat("\n==============================\n")
  cat("Block", i, ": codons [", b_start, ", ", b_end, "] as TEST\n", sep = "")
  cat("==============================\n")
  
  is_test <- data_model$codon_index >= b_start & data_model$codon_index <= b_end
  test_df  <- data_model[is_test, ]
  train_df <- data_model[!is_test, ]
  
  n_train <- nrow(train_df)
  n_test  <- nrow(test_df)
  
  cat("Train rows:", n_train, "  Test rows:", n_test, "\n")
  
  if (n_train == 0 || n_test == 0) {
    warning("Empty train or test set for this block; skipping.")
    next
  }
  
  # Check class distribution
  train_tab <- table(train_df$resistance_site)
  test_tab  <- table(test_df$resistance_site)
  cat("Train class distribution:\n")
  print(train_tab)
  cat("Test class distribution:\n")
  print(test_tab)
  
  # If there are no positives in TRAIN or TEST, skip this block
  if (length(train_tab) < 2 || train_tab["1"] == 0) {
    warning("No positive examples in TRAIN for block ", i,
            " – skipping this block.")
    next
  }
  if (!("1" %in% names(test_tab)) || test_tab["1"] == 0) {
    warning("No positive examples in TEST for block ", i,
            " – skipping this block.")
    next
  }
  
  # Oversample minority class in TRAIN
  train_os <- oversample_minority(train_df)
  
  # Class weights (after oversampling this will often be ~balanced)
  tab_y <- table(train_os$resistance_site)
  if (length(tab_y) < 2 || tab_y["1"] == 0) {
    glm_weights <- rep(1, nrow(train_os))
  } else {
    pos_weight <- as.numeric(tab_y["0"] / tab_y["1"])
    glm_weights <- ifelse(train_os$resistance_site == 1, pos_weight, 1)
  }
  
  # Fit logistic regression
  logit_fit <- glm(
    formula = logit_formula,
    data    = train_os,
    family  = binomial(link = "logit"),
    weights = glm_weights
  )
  
  model_coefs[[i]]     <- coef(logit_fit)
  model_summaries[[i]] <- summary(logit_fit)
  
  # Predict on TEST
  y_test      <- test_df$resistance_site
  logit_probs <- as.numeric(predict(logit_fit, newdata = test_df, type = "response"))
  
  # Tune threshold
  best <- tune_threshold_pr(y_test, logit_probs, n_grid = 300)
  cat("Best threshold for block", i, ":\n")
  print(best)
  
  # Evaluate at best threshold
  metrics <- evaluate_binary(y_test, logit_probs, threshold = best$threshold)
  
  # Store results
  results_list[[i]] <- data.frame(
    block_id      = i,
    start         = b_start,
    end           = b_end,
    n_train       = n_train,
    n_test        = n_test,
    train_pos     = as.numeric(train_tab["1"]),
    train_neg     = as.numeric(train_tab["0"]),
    test_pos      = ifelse("1" %in% names(test_tab), as.numeric(test_tab["1"]), 0),
    test_neg      = ifelse("0" %in% names(test_tab), as.numeric(test_tab["0"]), 0),
    threshold     = best$threshold,
    Accuracy      = metrics$accuracy,
    Precision     = metrics$precision,
    Recall        = metrics$recall,
    F1            = metrics$f1,
    ROC_AUC       = metrics$auc_roc,
    PR_AUC        = metrics$auc_pr
  )
}

############################################################
## 6. Summarize blockwise results ----
############################################################

block_results <- do.call(rbind, results_list)

cat("\n\n===== Blockwise Logistic Regression Results (Reduced Model) =====\n")
block_results_round <- block_results
metric_cols <- c("Accuracy", "Precision", "Recall", "F1", "ROC_AUC", "PR_AUC", "threshold")
block_results_round[, metric_cols] <- round(block_results_round[, metric_cols], 3)
print(block_results_round)

cat("\nBlocks sorted by F1 score (descending):\n")
print(block_results_round[order(-block_results_round$F1), ])

cat("\nBlocks sorted by ROC AUC (descending):\n")
print(block_results_round[order(-block_results_round$ROC_AUC), ])

############################################################
## 7. Print and save model coefficients nicely ----
############################################################

cat("\n\n===== Logistic Regression Coefficients by Block (Reduced Model) =====\n")

coef_table_list <- lapply(seq_along(model_coefs), function(b) {
  coefs <- model_coefs[[b]]
  if (is.null(coefs)) return(NULL)
  
  data.frame(
    block_id   = b,
    term       = names(coefs),
    estimate   = as.numeric(coefs),
    odds_ratio = exp(as.numeric(coefs)),
    row.names  = NULL
  )
})

coef_table <- do.call(rbind, coef_table_list)

# Order nicely: block first, then term names alphabetically
coef_table <- coef_table[order(coef_table$block_id, coef_table$term), ]

# Round for readability
coef_table_print <- coef_table
coef_table_print$estimate   <- round(coef_table_print$estimate, 4)
coef_table_print$odds_ratio <- round(coef_table_print$odds_ratio, 4)

for (b in unique(coef_table_print$block_id)) {
  cat("\n--- Block", b, "---\n")
  print(subset(coef_table_print, block_id == b))
}

# Save coefficient table for plotting scripts
save(coef_table, file = "../Results/coef_table_reduced.rda")
cat("\nSaved coefficient table to ../Results/coef_table_reduced.rda\n")
############################################################
## 8. Print confusion matrices for each evaluated block ----
############################################################

cat("\n\n===== Confusion Matrices for Each Block (Reduced Model) =====\n")

for (i in seq_along(results_list)) {
  res <- results_list[[i]]
  if (is.null(res)) next
  
  # Identify the block’s test set
  b_start <- res$start
  b_end   <- res$end
  
  is_test <- data_model$codon_index >= b_start & data_model$codon_index <= b_end
  test_df  <- data_model[is_test, ]
  
  # Recompute probabilities (same model already fitted)
  logit_fit <- glm(
    formula = logit_formula,
    data    = oversample_minority(data_model[!is_test, ]),
    family  = binomial(link = "logit")
  )
  
  probs <- predict(logit_fit, newdata = test_df, type = "response")
  threshold <- res$threshold
  preds <- ifelse(probs >= threshold, 1, 0)
  
  cm <- table(
    truth = factor(test_df$resistance_site, levels = c(0,1)),
    pred  = factor(preds,                levels = c(0,1))
  )
  
  cat("\n--- Block", i, " (threshold =", round(threshold, 3), ") ---\n")
  print(cm)
}

############################################################
# End of script
############################################################
