############################################################
# HIV-1 pol codon-level resistance site analysis
# with simple oversampling, class weights,
# threshold tuning, and PR curve for logistic regression
#
# Scholarly Super Serious Statistically Scientific Society
# Vincent Ngo, Charles Pierceton, Guy Ben Zeev, Marc-Philippe Gnagne
############################################################

########################
# CONFIG SWITCHES ------
########################

config <- list(
  use_oversampling   = TRUE,   # simple random oversampling of minority class
  use_weights_logit  = TRUE,   # class weights in logistic regression
  tune_threshold     = TRUE,   # choose optimal threshold via F1 on PR grid
  use_tree           = TRUE,
  use_rf             = TRUE,
  
  # NEW: train/test split strategy
  # "parity" = even codon_index train, odd test (original method)
  # "block"  = hold out a continuous block of codon indices as test
  split_type   = "block",      # choose "parity" or "block"
  
  # if split_type == "block", this range (inclusive) is used as TEST set
  block_range  = c(200, 300)   # e.g. codons 200–300 as test
)

cat("CONFIG:\n")
print(config)

########################
# 0. Package helper ----
########################

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_or_install("pROC")        # ROC/AUC
load_or_install("glmnet")      # regularized logistic regression
load_or_install("rpart")       # decision tree
load_or_install("rpart.plot")  # optional tree plots
load_or_install("randomForest")# RF / bagging-style trees
load_or_install("PRROC")       # precision–recall curves

set.seed(180)

############################################################
## 1. Load data ----
############################################################

# Path relative to R/ folder
codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

cat("Columns in codon_features.csv:\n")
print(names(codons))

############################################################
## 2. Basic cleaning & filtering ----
############################################################

if (!is.numeric(codons$codon_index)) {
  codons$codon_index <- as.numeric(codons$codon_index)
}

# Keep only complete codons, not stop codons
is_complete_codon <- !is.na(codons$ref_aa) &
  !grepl("-", codons$ref_codon) &
  nchar(codons$ref_codon) == 3 &
  codons$ref_aa != "*"

data <- codons[is_complete_codon, ]

cat("Number of codons after filtering incomplete/stop codons:",
    nrow(data), "\n")

if (!"resistance_site" %in% names(data)) {
  stop("Column 'resistance_site' not found in codon_features.csv.")
}

# Outcome as numeric 0/1
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
## 3. Train/test split: parity vs block ----
############################################################

if (config$split_type == "parity") {
  cat("\nUsing PARITY split: even codon_index = train, odd = test.\n")
  
  is_train <- (data$codon_index %% 2) == 0
  
} else if (config$split_type == "block") {
  br <- config$block_range
  if (length(br) != 2) {
    stop("config$block_range must be a numeric vector of length 2, e.g. c(200,300)")
  }
  block_start <- min(br)
  block_end   <- max(br)
  
  cat("\nUsing BLOCK split: test codons in [",
      block_start, ",", block_end, "], others = train.\n", sep = "")
  
  is_test_block <- data$codon_index >= block_start & data$codon_index <= block_end
  
  # ensure at least one positive in train and test, otherwise warn
  # (we'll still proceed, but this helps you catch bad ranges)
  tmp_train <- data[!is_test_block, ]
  tmp_test  <- data[is_test_block, ]
  if (!any(tmp_train$resistance_site == 1)) {
    warning("No resistance sites in TRAIN for this block split. Consider a different block_range.")
  }
  if (!any(tmp_test$resistance_site == 1)) {
    warning("No resistance sites in TEST for this block split. Consider a different block_range.")
  }
  
  is_train <- !is_test_block
  
} else {
  stop("config$split_type must be either 'parity' or 'block'.")
}

train <- data[is_train, ]
test  <- data[!is_train, ]

cat("Training codons:", nrow(train), "\n")
cat("Test codons:", nrow(test), "\n")

############################################################
## 4. Define predictors and modeling datasets ----
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

train_df <- train[, model_cols]
test_df  <- test[,  model_cols]

train_df <- train_df[complete.cases(train_df), ]
test_df  <- test_df[complete.cases(test_df), ]

cat("Training rows after complete-case filtering:", nrow(train_df), "\n")
cat("Test rows after complete-case filtering:", nrow(test_df), "\n")

# Numeric version
train_num <- train_df
test_num  <- test_df

# Factor version
train_fac <- train_df
test_fac  <- test_df
train_fac$resistance_site <- factor(train_fac$resistance_site)
test_fac$resistance_site  <- factor(test_fac$resistance_site)

cat("\nClass distribution in TRAIN (original):\n")
print(table(train_fac$resistance_site))
cat("Class distribution in TEST:\n")
print(table(test_fac$resistance_site))

############################################################
## 5. Simple random oversampling of minority class ----
############################################################

if (config$use_oversampling) {
  cat("\nApplying simple random oversampling of minority class...\n")
  
  tab <- table(train_fac$resistance_site)
  cat("Class distribution before oversampling:\n")
  print(tab)
  
  if (length(tab) < 2 || tab["1"] == 0) {
    warning("No positive examples in training data; cannot oversample.")
    train_fac_os <- train_fac
    train_num_os <- train_num
  } else {
    maj_class <- names(which.max(tab))
    min_class <- names(which.min(tab))
    n_maj     <- as.numeric(tab[maj_class])
    n_min     <- as.numeric(tab[min_class])
    
    maj_rows <- train_fac[train_fac$resistance_site == maj_class, ]
    min_rows <- train_fac[train_fac$resistance_site == min_class, ]
    
    # Upsample minority class to match majority size
    set.seed(180)
    min_upsampled <- min_rows[sample(seq_len(n_min), n_maj, replace = TRUE), ]
    
    train_fac_os <- rbind(maj_rows, min_upsampled)
    train_fac_os <- train_fac_os[sample(seq_len(nrow(train_fac_os))), ]
    
    cat("Class distribution after oversampling:\n")
    print(table(train_fac_os$resistance_site))
    
    # Numeric version
    train_num_os <- train_fac_os
    train_num_os$resistance_site <- as.numeric(
      as.character(train_fac_os$resistance_site)
    )
  }
} else {
  cat("\nSkipping oversampling. Using original training data.\n")
  train_fac_os <- train_fac
  train_num_os <- train_num
}

############################################################
## 6. Evaluation helpers ----
############################################################

evaluate_binary <- function(y_true, y_prob, threshold = 0.5, model_name = "") {
  # y_true numeric 0/1; y_prob numeric in [0,1]
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
  f1        <- ifelse(is.na(precision) || is.na(recall) || (precision + recall == 0),
                      NA,
                      2 * precision * recall / (precision + recall))
  
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_roc <- as.numeric(pROC::auc(roc_obj))
  
  cat("\n==============================\n")
  cat("Model:", model_name, "\n")
  cat("==============================\n")
  cat("Threshold used:", threshold, "\n")
  cat("Confusion matrix (truth x pred):\n")
  print(tab)
  cat("\nAccuracy :", round(accuracy, 3), "\n")
  cat("Precision:", round(precision, 3), "\n")
  cat("Recall   :", round(recall, 3), "\n")
  cat("F1 score :", round(f1, 3), "\n")
  cat("ROC AUC  :", round(auc_roc, 3), "\n\n")
  
  invisible(list(
    confusion = tab,
    accuracy  = accuracy,
    precision = precision,
    recall    = recall,
    f1        = f1,
    auc_roc   = auc_roc
  ))
}

# Threshold tuning via PR / F1
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
    f1        <- ifelse(is.na(precision) || is.na(recall) || (precision + recall == 0),
                        NA,
                        2 * precision * recall / (precision + recall))
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
## 7. Logistic regression (glm) with weights ----
############################################################

logit_formula <- as.formula(
  paste("resistance_site ~", paste(predictor_vars, collapse = " + "))
)

cat("\nFitting logistic regression (glm) with oversampling and optional class weights...\n")

if (config$use_weights_logit) {
  tab_y <- table(train_num_os$resistance_site)
  if (length(tab_y) < 2 || tab_y["1"] == 0) {
    warning("No positive examples in training data after oversampling; weights disabled.")
    glm_weights <- rep(1, nrow(train_num_os))
  } else {
    pos_weight <- as.numeric(tab_y["0"] / tab_y["1"])
    glm_weights <- ifelse(train_num_os$resistance_site == 1, pos_weight, 1)
    cat("Using class weights for glm. Positive class weight:", pos_weight, "\n")
  }
} else {
  glm_weights <- rep(1, nrow(train_num_os))
  cat("Not using class weights for glm.\n")
}

logit_fit <- glm(
  formula = logit_formula,
  data    = train_num_os,
  family  = binomial(link = "logit"),
  weights = glm_weights
)

# Predict probabilities on test set
logit_probs <- predict(logit_fit, newdata = test_num, type = "response")

# PR curve + AUPRC
y_test <- test_num$resistance_site
pr_obj <- PRROC::pr.curve(
  scores.class0 = logit_probs[y_test == 1],
  scores.class1 = logit_probs[y_test == 0],
  curve = TRUE
)
auc_pr <- pr_obj$auc.integral
cat("Logistic regression AUPRC (area under PR curve):",
    round(auc_pr, 3), "\n")

# Tune threshold based on F1
if (config$tune_threshold) {
  best <- tune_threshold_pr(y_test, logit_probs, n_grid = 300)
  cat("Best threshold (max F1) for logistic regression:\n")
  print(best)
  best_threshold <- best$threshold
} else {
  best_threshold <- 0.5
}

metrics_logit <- evaluate_binary(
  y_true    = y_test,
  y_prob    = logit_probs,
  threshold = best_threshold,
  model_name = "Logistic Regression (oversampling + weights + tuned threshold)"
)

############################################################
## 8. Regularized logistic regression (glmnet) ----
############################################################

cat("\nFitting regularized logistic regression (glmnet)...\n")

X_train <- model.matrix(
  ~ . - codon_index - resistance_site,
  data = train_num_os
)
y_train <- train_num_os$resistance_site

X_test <- model.matrix(
  ~ . - codon_index - resistance_site,
  data = test_num
)

cv_fit <- cv.glmnet(
  x      = X_train,
  y      = y_train,
  family = "binomial",
  alpha  = 0.5
)

cat("Best lambda (min):", cv_fit$lambda.min, "\n")

glmnet_fit <- glmnet(
  x      = X_train,
  y      = y_train,
  family = "binomial",
  alpha  = 0.5,
  lambda = cv_fit$lambda.min
)

glmnet_probs <- as.numeric(predict(glmnet_fit, newx = X_test, type = "response"))

metrics_glmnet <- evaluate_binary(
  y_true     = y_test,
  y_prob     = glmnet_probs,
  threshold  = 0.5,
  model_name = "Regularized Logistic (glmnet, alpha = 0.5, threshold=0.5)"
)

############################################################
## 9. Decision tree (optional) ----
############################################################

if (config$use_tree) {
  cat("\nFitting decision tree (rpart)...\n")
  
  tree_fit <- rpart(
    formula = logit_formula,
    data    = train_fac_os,
    method  = "class",
    control = rpart.control(
      cp       = 0.01,
      minsplit = 20,
      maxdepth = 5
    )
  )
  
  # rpart.plot(tree_fit) # optional
  
  tree_probs <- predict(tree_fit, newdata = test_fac, type = "prob")[, "1"]
  
  metrics_tree <- evaluate_binary(
    y_true     = y_test,
    y_prob     = tree_probs,
    threshold  = 0.5,
    model_name = "Decision Tree (rpart)"
  )
} else {
  metrics_tree <- list(
    accuracy = NA, precision = NA, recall = NA,
    f1 = NA, auc_roc = NA
  )
}

############################################################
## 10. Random forest (optional) ----
############################################################

if (config$use_rf) {
  cat("\nFitting random forest (balanced sampling)...\n")
  
  tab_rf <- table(train_fac_os$resistance_site)
  min_class_size <- min(tab_rf)
  cat("Balanced RF with sampsize per class =", min_class_size, "\n")
  
  rf_fit <- randomForest(
    formula   = logit_formula,
    data      = train_fac_os,
    ntree     = 700,
    mtry      = floor(sqrt(length(predictor_vars))),
    importance = TRUE,
    strata    = train_fac_os$resistance_site,
    sampsize  = rep(min_class_size, length(tab_rf))
  )
  
  rf_probs <- as.numeric(predict(rf_fit, newdata = test_fac, type = "prob")[, "1"])
  
  metrics_rf <- evaluate_binary(
    y_true     = y_test,
    y_prob     = rf_probs,
    threshold  = 0.5,
    model_name = "Random Forest (balanced)"
  )
  
  cat("\nRandom forest variable importance:\n")
  print(importance(rf_fit))
} else {
  metrics_rf <- list(
    accuracy = NA, precision = NA, recall = NA,
    f1 = NA, auc_roc = NA
  )
}

############################################################
## 11. Compare models ----
############################################################

compare_df <- data.frame(
  Model     = c("Logit_oversampled_tuned", "Logit_glmnet", "Tree", "RandomForest"),
  Accuracy  = c(metrics_logit$accuracy,
                metrics_glmnet$accuracy,
                metrics_tree$accuracy,
                metrics_rf$accuracy),
  ROC_AUC   = c(metrics_logit$auc_roc,
                metrics_glmnet$auc_roc,
                metrics_tree$auc_roc,
                metrics_rf$auc_roc),
  Precision = c(metrics_logit$precision,
                metrics_glmnet$precision,
                metrics_tree$precision,
                metrics_rf$precision),
  Recall    = c(metrics_logit$recall,
                metrics_glmnet$recall,
                metrics_tree$recall,
                metrics_rf$recall),
  F1        = c(metrics_logit$f1,
                metrics_glmnet$f1,
                metrics_tree$f1,
                metrics_rf$f1)
)

cat("\n\n===== Model performance comparison (test set) =====\n")
compare_df_round <- compare_df
compare_df_round[ , -1] <- round(compare_df_round[ , -1], 3)
print(compare_df_round)

############################################################
# End of script
############################################################
