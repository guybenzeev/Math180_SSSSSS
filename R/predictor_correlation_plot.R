############################################################
# Correlation Heatmap of Predictors for HIV-1 pol Features
#
# Produces a correlation chart of all continuous predictors
# used in the blockwise logistic regression project.
#
# Scholarly Super Serious Statistically Scientific Society
############################################################

# Install/load required packages
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_or_install("corrplot")
load_or_install("ggplot2")

set.seed(180)

############################################################
## 1. Load and preprocess data ----
############################################################

codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

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

# Compute log_dNdS if available
if ("dNdS_proxy" %in% names(data)) {
  data$log_dNdS <- log1p(data$dNdS_proxy)
}

############################################################
## 2. Select predictor variables ----
############################################################

predictor_vars <- c(
  "codon_degeneracy",
  "syn_frac",
  "nonsyn_frac",
  "shannon_entropy",
  "conservation",
  if ("log_dNdS" %in% names(data)) "log_dNdS" else "dNdS_proxy"
)

# Subset to predictors only + drop rows with missing values
predictor_df <- data[, predictor_vars]
predictor_df <- predictor_df[complete.cases(predictor_df), ]

############################################################
## 3. Compute correlation matrix ----
############################################################

cor_matrix <- cor(predictor_df, use = "pairwise.complete.obs")

############################################################
## 4. Plot correlation heatmap ----
############################################################

# Print to RStudio viewer
corrplot(cor_matrix,
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 0.9,
         number.cex = 0.7,
         addCoef.col = "black",
         mar = c(1,1,1,1),
         title = "Correlation of Predictors for HIV-1 pol Codon Features")

# Save to PDF
pdf("../Results/predictor_correlation_heatmap.pdf", width = 7, height = 7)
corrplot(cor_matrix,
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 0.9,
         number.cex = 0.7,
         addCoef.col = "black",
         mar = c(1,1,1,1),
         title = "Correlation of Predictors for HIV-1 pol Codon Features")
dev.off()

cat("\nSaved correlation heatmap to ../Results/predictor_correlation_heatmap.pdf\n")
