############################################################
# Neural Network Analysis for HIV-1 pol codon features
#
# - Load codon_features.csv
# - Filter incomplete/stop codons
# - Choose train/test split (even/odd or block)
# - Oversample minority class in TRAIN
# - Fit 1-hidden-layer neural network (nnet)
# - Tune threshold (F1) on TEST
# - Report performance
#
# Scholarly Super Serious Statistically Scientific Society
# Vincent Ngo, Charles Pierceton, Guy Ben Zeev, Marc-Philippe Gnagne
############################################################

########################
# 0. Config -----------
########################

CONFIG <- list(
  split_type      = "block",   # "even_odd" or "block"
  block_range     = c(200, 300),  # used if split_type == "block"
  use_oversampling = TRUE,
  tune_threshold   = TRUE,
  hidden_units     = 5,           # size of hidden layer
  decay            = 1e-4,        # L2 penalty
  maxit            = 500          # training iterations
)

########################
# 1. Package helper ----
########################

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_or_install("nnet")
load_or_install("pROC")
load_or_install("PRROC")

set.seed(180)

############################################################
## 2. Load and preprocess data ----
############################################################

# Assuming this script is in the R/ folder:
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

# log dNdS
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
## 3. Select predictors and build model data ----
############################################################

predictor_vars <- c(
  "codon_degeneracy",
  "syn_frac",
  "nonsyn_frac",
  "aa_chemistry",
  "shannon_entropy",
  "conservation",
  if ("log_dNdS" %in% names(data)) "log_dNdS" else "dNdS_proxy"
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

# Ensure resistance_site is 0/1 numeric
data_model$resistance_site <- as.numeric(data_model$resistance_site)

cat("Overall class distribution:\n")
print(table(data_model$resistance_site))

############################################################
## 4. Train/test split ----
############################################################

if (CONFIG$split_type == "even_odd") {
  cat("\nUsing EVEN-ODD split: even codon_index = TRAIN, odd = TEST.\n")
  is_train <- (data_model$codon_index %% 2 == 0)
} else if (CONFIG$split_type == "block") {
  br <- CONFIG$block_range
  cat("\nUsing BLOCK split: codons in [", br[1], ",", br[2],
      "] = TEST, others = TRAIN.\n", sep = "")
  is_train <- !(data_model$codon_index >= br[1] & data_model$codon_index <= br[2])
} else {
  stop("Unknown CONFIG$split_type: ", CONFIG$split_type)
}

train_df <- data_model[is_train, ]
test_df  <- data_model[!is_train, ]

cat("Training rows:", nrow(train_df), "\n")
cat("Test rows    :", nrow(test_df), "\n")

cat("\nClass distribution in TRAIN:\n")
print(table(train_df$resistance_site))
cat("Class distribution in TEST:\n")
print(table(test_df$resistance_site))

if (length(table(train_df$resistance_site)) < 2 ||
    !("1" %in% names(table(train_df$resistance_site))) ||
    table(train_df$resistance_site)["1"] == 0) {
  stop("No positive examples in TRAIN – cannot train a NN classifier.")
}

if (length(table(test_df$resistance_site)) < 2 ||
    !("1" %in% names(table(test_df$resistance_site))) ||
    table(test_df$resistance_site)["1"] == 0) {
  warning("No positives in TEST – recall/F1/PR AUC will be undefined or NA.")
}

############################################################
## 5. Helper functions: oversampling, eval, threshold ----
############################################################

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
  
  df_os_num <- df_os
  df_os_num$resistance_site <- as.numeric(as.character(df_os$resistance_site))
  df_os_num
}

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
  
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_roc <- as.numeric(pROC::auc(roc_obj))
  
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

tune_threshold_pr <- function(y_true, y_prob, n_grid = 200) {
  thresholds <- seq(0.001, 0.999, length.out = n_grid)
  best_f1 <- -Inf
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
## 6. Oversample TRAIN (optional) ----
############################################################

if (CONFIG$use_oversampling) {
  cat("\nApplying simple random oversampling of minority class in TRAIN...\n")
  cat("Before oversampling:\n")
  print(table(train_df$resistance_site))
  train_df <- oversample_minority(train_df)
  cat("After oversampling:\n")
  print(table(train_df$resistance_site))
}

############################################################
## 7. Build NN design matrices (model.matrix + scaling) ----
############################################################

# Use model.matrix to handle factors -> numeric design matrix
X_train <- model.matrix(
  ~ . - resistance_site,
  data = train_df[, c("resistance_site", predictor_vars)]
)
X_test <- model.matrix(
  ~ . - resistance_site,
  data = test_df[, c("resistance_site", predictor_vars)]
)

# Align columns (just in case)
common_cols <- intersect(colnames(X_train), colnames(X_test))
X_train <- X_train[, common_cols, drop = FALSE]
X_test  <- X_test[, common_cols, drop = FALSE]

# Standardize (mean 0, sd 1) using TRAIN statistics
train_means <- apply(X_train, 2, mean)
train_sds   <- apply(X_train, 2, sd)
train_sds[train_sds == 0] <- 1  # avoid division by zero

X_train_scaled <- scale(X_train, center = train_means, scale = train_sds)
X_test_scaled  <- scale(X_test,  center = train_means, scale = train_sds)

y_train <- factor(train_df$resistance_site)
y_test  <- test_df$resistance_site

############################################################
## 8. Fit neural network (nnet) ----
############################################################

cat("\nFitting neural network with nnet...\n")
cat("Hidden units:", CONFIG$hidden_units,
    "  decay:", CONFIG$decay,
    "  maxit:", CONFIG$maxit, "\n")

train_nn_df <- data.frame(y = y_train, X_train_scaled)

nn_fit <- nnet::nnet(
  y ~ .,
  data   = train_nn_df,
  size   = CONFIG$hidden_units,
  decay  = CONFIG$decay,
  maxit  = CONFIG$maxit,
  trace  = FALSE
)

############################################################
## 9. Predict & evaluate on TEST ----
############################################################

test_nn_df <- data.frame(X_test_scaled)
nn_probs_raw <- predict(nn_fit, newdata = test_nn_df, type = "raw")

# For binary classification, nnet returns a vector or 1-column matrix
nn_probs <- as.numeric(nn_probs_raw)

if (CONFIG$tune_threshold) {
  best <- tune_threshold_pr(y_test, nn_probs, n_grid = 300)
  cat("\nBest threshold for NN (max F1 on TEST):\n")
  print(best)
  threshold <- best$threshold
} else {
  threshold <- 0.5
  best <- list(threshold = threshold, precision = NA, recall = NA, f1 = NA)
}

metrics <- evaluate_binary(y_test, nn_probs, threshold = threshold)

cat("\n==============================\n")
cat("Model: Neural Network (nnet)\n")
cat("==============================\n")
cat("Threshold used:", threshold, "\n")
cat("Confusion matrix (truth x pred):\n")
print(metrics$confusion)

cat("\nAccuracy :", round(metrics$accuracy, 3), "\n")
cat("Precision:", round(metrics$precision, 3), "\n")
cat("Recall   :", round(metrics$recall, 3), "\n")
cat("F1 score :", round(metrics$f1, 3), "\n")
cat("ROC AUC  :", round(metrics$auc_roc, 3), "\n")
cat("PR  AUC  :", round(metrics$auc_pr, 3), "\n")

############################################################
# End of script
############################################################
