############################################################
# HIV-1 pol codon-level resistance site analysis with Bagging
# Scholarly Super Serious Statistically Scientific Society
# Vincent Ngo, Charles Pierceton, Guy Ben Zeev, Marc-Philippe Gnagne
############################################################

########################
# CONFIG SWITCHES ------
########################

config <- list(
  use_weights_logit  = TRUE,  # weighted logistic regression
  use_oversampling   = TRUE,  # simple random oversampling of minority class
  use_balanced_bagging = TRUE # bagging with balanced bootstrap sampling
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

load_or_install("pROC")        # AUC
load_or_install("glmnet")      # regularized logistic regression
load_or_install("rpart")       # decision tree
load_or_install("rpart.plot")  # tree plotting (optional)
load_or_install("randomForest")# for bagged trees

set.seed(180)  # reproducibility

############################################################
## 1. Load data ----
############################################################

# Adjust path if needed (relative to this script)
codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

cat("Columns in codon_features.csv:\n")
print(names(codons))

############################################################
## 2. Basic cleaning & filtering ----
############################################################

if (!is.numeric(codons$codon_index)) {
  codons$codon_index <- as.numeric(codons$codon_index)
}

# Keep only complete codons with defined amino acid and not stop codons
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

# Coerce outcome to numeric 0/1
data$resistance_site <- as.numeric(data$resistance_site)

# log-transform dNdS_proxy if present
if ("dNdS_proxy" %in% names(data)) {
  data$log_dNdS <- log1p(data$dNdS_proxy)
} else {
  warning("No 'dNdS_proxy' column found; skipping log_dNdS feature.")
}

# Categorical variables
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
## 3. Train/test split by codon_index parity ----
############################################################

is_train <- (data$codon_index %% 2) == 0

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

# Numeric version (for glm / glmnet)
train_num <- train_df
test_num  <- test_df

# Factor version (for tree / bagging / oversampling)
train_fac <- train_df
test_fac  <- test_df
train_fac$resistance_site <- factor(train_fac$resistance_site)
test_fac$resistance_site  <- factor(test_fac$resistance_site)

############################################################
## 5. Optional oversampling (for imbalance) ----
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
    n_maj <- as.numeric(tab[maj_class])
    n_min <- as.numeric(tab[min_class])
    
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
## 6. Evaluation helper ----
############################################################

evaluate_binary <- function(y_true, y_prob, threshold = 0.5, model_name = "") {
  # y_true should be numeric 0/1
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
  auc     <- as.numeric(pROC::auc(roc_obj))
  
  cat("\n==============================\n")
  cat("Model:", model_name, "\n")
  cat("==============================\n")
  cat("Confusion matrix (truth x pred):\n")
  print(tab)
  cat("\nAccuracy :", round(accuracy, 3), "\n")
  cat("Precision:", round(precision, 3), "\n")
  cat("Recall   :", round(recall, 3), "\n")
  cat("F1 score :", round(f1, 3), "\n")
  cat("AUC      :", round(auc, 3), "\n\n")
  
  invisible(list(
    confusion = tab,
    accuracy  = accuracy,
    precision = precision,
    recall    = recall,
    f1        = f1,
    auc       = auc
  ))
}

############################################################
## 7. Logistic regression (glm, weighted) ----
############################################################

logit_formula <- as.formula(
  paste("resistance_site ~", paste(predictor_vars, collapse = " + "))
)

cat("\nFitting logistic regression (glm)...\n")

if (config$use_weights_logit) {
  tab_y <- table(train_num_os$resistance_site)
  if (length(tab_y) < 2 || tab_y["1"] == 0) {
    warning("No positive examples in training data after oversampling/filters; weights disabled.")
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

summary(logit_fit)

logit_probs <- predict(logit_fit, newdata = test_num, type = "response")

metrics_logit <- evaluate_binary(
  y_true    = test_num$resistance_site,
  y_prob    = logit_probs,
  threshold = 0.5,
  model_name = if (config$use_weights_logit)
    "Logistic Regression (weighted glm)"
  else
    "Logistic Regression (unweighted glm)"
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
y_test <- test_num$resistance_site

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
  model_name = "Regularized Logistic (glmnet, alpha = 0.5)"
)

############################################################
## 9. Decision tree (rpart) ----
############################################################

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

# Optional plot:
# rpart.plot(tree_fit, main = "Decision Tree for Resistance Site")

tree_probs <- predict(tree_fit, newdata = test_fac, type = "prob")[, "1"]

metrics_tree <- evaluate_binary(
  y_true     = test_num$resistance_site,
  y_prob     = tree_probs,
  threshold  = 0.5,
  model_name = "Decision Tree (rpart)"
)

############################################################
## 10. Bagged trees (randomForest with mtry = p) ----
############################################################

cat("\nFitting bagged trees (randomForest with mtry = number of predictors)...\n")

bag_formula <- logit_formula
p <- length(predictor_vars)

if (config$use_balanced_bagging) {
  tab_bag <- table(train_fac_os$resistance_site)
  min_class_size <- min(tab_bag)
  cat("Balanced bagging with stratified sampling; sampsize per class =",
      min_class_size, "\n")
  
  bag_fit <- randomForest(
    formula   = bag_formula,
    data      = train_fac_os,
    ntree     = 700,
    mtry      = p,  # all predictors at each split = bagging
    importance = TRUE,
    strata    = train_fac_os$resistance_site,
    sampsize  = rep(min_class_size, length(tab_bag))
  )
} else {
  cat("Standard bagging (no explicit class balancing).\n")
  bag_fit <- randomForest(
    formula   = bag_formula,
    data      = train_fac_os,
    ntree     = 500,
    mtry      = p,
    importance = TRUE
  )
}

bag_probs <- as.numeric(predict(bag_fit, newdata = test_fac, type = "prob")[, "1"])

metrics_bag <- evaluate_binary(
  y_true     = test_num$resistance_site,
  y_prob     = bag_probs,
  threshold  = 0.5,
  model_name = if (config$use_balanced_bagging)
    "Bagging (balanced)"
  else
    "Bagging (unbalanced)"
)

cat("\nBagging variable importance:\n")
print(importance(bag_fit))
# Optional:
# varImpPlot(bag_fit)

############################################################
## 11. Compare models ----
############################################################

compare_df <- data.frame(
  Model     = c("Logit_glm", "Logit_glmnet", "Tree", "Bagging"),
  Accuracy  = c(metrics_logit$accuracy,
                metrics_glmnet$accuracy,
                metrics_tree$accuracy,
                metrics_bag$accuracy),
  AUC       = c(metrics_logit$auc,
                metrics_glmnet$auc,
                metrics_tree$auc,
                metrics_bag$auc),
  Precision = c(metrics_logit$precision,
                metrics_glmnet$precision,
                metrics_tree$precision,
                metrics_bag$precision),
  Recall    = c(metrics_logit$recall,
                metrics_glmnet$recall,
                metrics_tree$recall,
                metrics_bag$recall),
  F1        = c(metrics_logit$f1,
                metrics_glmnet$f1,
                metrics_tree$f1,
                metrics_bag$f1)
)

cat("\n\n===== Model performance comparison (test set) =====\n")
compare_df_round <- compare_df
compare_df_round[ , -1] <- round(compare_df_round[ , -1], 3)
print(compare_df_round)

############################################################
# End of script
############################################################
