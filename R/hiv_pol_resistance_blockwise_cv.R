############################################################
# HIV-1 pol codon-level resistance site analysis
# BLOCKWISE CROSS-VALIDATION (100-codon blocks)
#
# Models:
#  - Logistic regression (glm) + oversampling + optional weights + tuned threshold
#  - Regularized logistic regression (glmnet)
#  - Decision tree (rpart)
#  - Random forest (randomForest) balanced sampling
############################################################

########################
# CONFIG ------
########################
config <- list(
  use_oversampling   = TRUE,
  use_weights_logit  = TRUE,
  tune_threshold     = TRUE,

  use_glm            = TRUE,
  use_glmnet         = TRUE,
  use_tree           = TRUE,
  use_rf             = TRUE,

  block_size         = 100,     # size of each held-out block
  min_test_pos       = 1,       # skip blocks with < this many positives in TEST
  min_train_pos      = 1,       # skip blocks with < this many positives in TRAIN
  seed               = 180
)

cat("CONFIG:\n"); print(config)
set.seed(config$seed)

########################
# Package helper ----
########################
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

load_or_install("pROC")
load_or_install("glmnet")
load_or_install("rpart")
load_or_install("randomForest")
load_or_install("PRROC")

############################################################
## 1. Load data ----
############################################################
codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)
cat("Columns in codon_features.csv:\n"); print(names(codons))

if (!is.numeric(codons$codon_index)) codons$codon_index <- as.numeric(codons$codon_index)

# Keep only complete codons, not stop codons
is_complete_codon <- !is.na(codons$ref_aa) &
  !grepl("-", codons$ref_codon) &
  nchar(codons$ref_codon) == 3 &
  codons$ref_aa != "*"

data <- codons[is_complete_codon, ]
cat("Number of codons after filtering incomplete/stop codons:", nrow(data), "\n")

if (!"resistance_site" %in% names(data)) stop("Column 'resistance_site' not found.")
data$resistance_site <- as.numeric(data$resistance_site)

# Reduced predictor set (matches your reduced model)
if ("aa_chemistry" %in% names(data)) data$aa_chemistry <- factor(data$aa_chemistry)

predictor_vars <- c(
  "codon_degeneracy",
  "syn_frac",
  "shannon_entropy",
  "conservation",
  "aa_chemistry"
)

missing_predictors <- setdiff(predictor_vars, names(data))
if (length(missing_predictors) > 0) {
  stop("Missing predictors: ", paste(missing_predictors, collapse = ", "))
}

model_cols <- c("codon_index", "resistance_site", predictor_vars)
data_model <- data[, model_cols]
data_model <- data_model[complete.cases(data_model), ]
cat("Rows after complete-case filtering:", nrow(data_model), "\n")
cat("Overall class distribution:\n"); print(table(data_model$resistance_site))

############################################################
## 2. Helpers ----
############################################################
oversample_minority <- function(df) {
  df_fac <- df
  df_fac$resistance_site <- factor(df_fac$resistance_site)
  tab <- table(df_fac$resistance_site)

  if (length(tab) < 2 || tab["1"] == 0) return(df)

  maj_class <- names(which.max(tab))
  min_class <- names(which.min(tab))
  n_maj <- as.numeric(tab[maj_class])
  n_min <- as.numeric(tab[min_class])

  maj_rows <- df_fac[df_fac$resistance_site == maj_class, ]
  min_rows <- df_fac[df_fac$resistance_site == min_class, ]

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

  TN <- tab["0","0"]; FP <- tab["0","1"]; FN <- tab["1","0"]; TP <- tab["1","1"]
  accuracy  <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
  f1        <- ifelse(is.na(precision) || is.na(recall) || (precision + recall == 0),
                      NA, 2 * precision * recall / (precision + recall))

  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_roc <- as.numeric(pROC::auc(roc_obj))

  # PR AUC if possible
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

  list(accuracy=accuracy, precision=precision, recall=recall, f1=f1, auc_roc=auc_roc, auc_pr=auc_pr, confusion=tab)
}

tune_threshold_pr <- function(y_true, y_prob, n_grid = 300) {
  thresholds <- seq(0.001, 0.999, length.out = n_grid)
  best_f1 <- -Inf; best_t <- 0.5; best_stats <- NULL

  for (t in thresholds) {
    y_pred <- ifelse(y_prob >= t, 1, 0)
    tab <- table(
      truth = factor(y_true, levels = c(0, 1)),
      pred  = factor(y_pred, levels = c(0, 1))
    )
    TN <- tab["0","0"]; FP <- tab["0","1"]; FN <- tab["1","0"]; TP <- tab["1","1"]
    precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
    recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
    f1        <- ifelse(is.na(precision) || is.na(recall) || (precision + recall == 0),
                        NA, 2 * precision * recall / (precision + recall))

    if (!is.na(f1) && f1 > best_f1) {
      best_f1 <- f1; best_t <- t
      best_stats <- list(threshold=t, precision=precision, recall=recall, f1=f1)
    }
  }
  best_stats
}

############################################################
## 3. Define blocks ----
############################################################
min_idx <- min(data_model$codon_index, na.rm = TRUE)
max_idx <- max(data_model$codon_index, na.rm = TRUE)

block_starts <- seq(from = min_idx, to = max_idx, by = config$block_size)
blocks <- data.frame(
  block_id = seq_along(block_starts),
  start    = block_starts,
  end      = pmin(block_starts + config$block_size - 1, max_idx)
)

cat("\nDefined blocks:\n"); print(blocks)

############################################################
## 4. CV loop over blocks ----
############################################################
logit_formula <- as.formula(paste("resistance_site ~", paste(predictor_vars, collapse = " + ")))

results <- list()

for (i in seq_len(nrow(blocks))) {
  b_start <- blocks$start[i]
  b_end   <- blocks$end[i]

  is_test <- data_model$codon_index >= b_start & data_model$codon_index <= b_end
  test_df  <- data_model[is_test, ]
  train_df <- data_model[!is_test, ]

  train_tab <- table(train_df$resistance_site)
  test_tab  <- table(test_df$resistance_site)

  cat("\n==============================\n")
  cat("Block", i, "TEST codons [", b_start, ",", b_end, "]\n", sep = "")
  cat("==============================\n")
  cat("Train class dist:\n"); print(train_tab)
  cat("Test class dist:\n"); print(test_tab)

  # Skip blocks without enough positives
  train_pos <- ifelse("1" %in% names(train_tab), as.numeric(train_tab["1"]), 0)
  test_pos  <- ifelse("1" %in% names(test_tab),  as.numeric(test_tab["1"]),  0)
  if (train_pos < config$min_train_pos || test_pos < config$min_test_pos) {
    cat("Skipping block", i, "due to insufficient positives (train_pos=", train_pos, ", test_pos=", test_pos, ")\n", sep="")
    next
  }

  # Oversample train for models that benefit from it
  train_os <- if (config$use_oversampling) oversample_minority(train_df) else train_df

  # Factor versions for tree/RF
  train_fac <- train_os
  test_fac  <- test_df
  train_fac$resistance_site <- factor(train_fac$resistance_site)
  test_fac$resistance_site  <- factor(test_fac$resistance_site)

  y_test <- test_df$resistance_site

  # Storage for block metrics
  block_row <- data.frame(
    block_id=i, start=b_start, end=b_end,
    train_pos=train_pos, test_pos=test_pos,
    stringsAsFactors = FALSE
  )

  ########################
  # (A) Logistic regression (glm)
  ########################
  if (config$use_glm) {
    if (config$use_weights_logit) {
      tab_y <- table(train_os$resistance_site)
      if (length(tab_y) < 2 || tab_y["1"] == 0) {
        glm_weights <- rep(1, nrow(train_os))
      } else {
        pos_weight <- as.numeric(tab_y["0"] / tab_y["1"])
        glm_weights <- ifelse(train_os$resistance_site == 1, pos_weight, 1)
      }
    } else {
      glm_weights <- rep(1, nrow(train_os))
    }

    glm_fit <- glm(
      formula = logit_formula,
      data    = train_os,
      family  = binomial(link = "logit"),
      weights = glm_weights
    )

    glm_probs <- as.numeric(predict(glm_fit, newdata = test_df, type = "response"))

    if (config$tune_threshold) {
      best <- tune_threshold_pr(y_test, glm_probs, n_grid = 300)
      th <- ifelse(is.null(best), 0.5, best$threshold)
    } else {
      th <- 0.5
    }

    m <- evaluate_binary(y_test, glm_probs, threshold = th)

    block_row$glm_threshold <- th
    block_row$glm_acc <- m$accuracy
    block_row$glm_prec <- m$precision
    block_row$glm_rec <- m$recall
    block_row$glm_f1 <- m$f1
    block_row$glm_rocauc <- m$auc_roc
    block_row$glm_prauc <- m$auc_pr
  }

  ########################
  # (B) Regularized logistic regression (glmnet)
  ########################
  if (config$use_glmnet) {
    X_train <- model.matrix(~ . - codon_index - resistance_site, data = train_os)
    y_train <- train_os$resistance_site
    X_test  <- model.matrix(~ . - codon_index - resistance_site, data = test_df)

    cv_fit <- cv.glmnet(x = X_train, y = y_train, family = "binomial", alpha = 0.5)
    fit <- glmnet(x = X_train, y = y_train, family = "binomial", alpha = 0.5, lambda = cv_fit$lambda.min)
    probs <- as.numeric(predict(fit, newx = X_test, type = "response"))

    m <- evaluate_binary(y_test, probs, threshold = 0.5)

    block_row$glmnet_acc <- m$accuracy
    block_row$glmnet_prec <- m$precision
    block_row$glmnet_rec <- m$recall
    block_row$glmnet_f1 <- m$f1
    block_row$glmnet_rocauc <- m$auc_roc
    block_row$glmnet_prauc <- m$auc_pr
    block_row$glmnet_lambda <- cv_fit$lambda.min
  }

  ########################
  # (C) Decision tree (rpart)
  ########################
  if (config$use_tree) {
    tree_fit <- rpart(
      formula = logit_formula,
      data    = train_fac,
      method  = "class",
      control = rpart.control(cp = 0.01, minsplit = 20, maxdepth = 5)
    )

    # probability of class "1"
    tree_probs <- as.numeric(predict(tree_fit, newdata = test_fac, type = "prob")[, "1"])
    m <- evaluate_binary(y_test, tree_probs, threshold = 0.5)

    block_row$tree_acc <- m$accuracy
    block_row$tree_prec <- m$precision
    block_row$tree_rec <- m$recall
    block_row$tree_f1 <- m$f1
    block_row$tree_rocauc <- m$auc_roc
    block_row$tree_prauc <- m$auc_pr
  }

  ########################
  # (D) Random forest (balanced)
  ########################
  if (config$use_rf) {
    tab_rf <- table(train_fac$resistance_site)
    min_class_size <- min(tab_rf)

    rf_fit <- randomForest(
      formula    = logit_formula,
      data       = train_fac,
      ntree      = 700,
      mtry       = max(1, floor(sqrt(length(predictor_vars)))),
      importance = TRUE,
      strata     = train_fac$resistance_site,
      sampsize   = rep(min_class_size, length(tab_rf))
    )

    rf_probs <- as.numeric(predict(rf_fit, newdata = test_fac, type = "prob")[, "1"])
    m <- evaluate_binary(y_test, rf_probs, threshold = 0.5)

    block_row$rf_acc <- m$accuracy
    block_row$rf_prec <- m$precision
    block_row$rf_rec <- m$recall
    block_row$rf_f1 <- m$f1
    block_row$rf_rocauc <- m$auc_roc
    block_row$rf_prauc <- m$auc_pr
  }

  results[[length(results) + 1]] <- block_row
}

############################################################
## 5. Summarize results ----
############################################################
if (length(results) == 0) stop("No blocks evaluated (likely due to min_pos filters).")

block_results <- do.call(rbind, results)

cat("\n\n===== BLOCKWISE RESULTS =====\n")
print(block_results)

# Optional: simple averages (over evaluated blocks only)
metric_cols <- names(block_results)[grepl("_(acc|prec|rec|f1|rocauc|prauc)$", names(block_results))]
avg <- sapply(block_results[, metric_cols, drop = FALSE], function(x) mean(x, na.rm = TRUE))

cat("\n===== AVERAGE METRICS (evaluated blocks only) =====\n")
print(round(avg, 3))

# Save for plotting
write.csv(block_results, "../Results/blockwise_model_results.csv", row.names = FALSE)
cat("\nWrote ../Results/blockwise_model_results.csv\n")
