############################################################
# Multi-panel HIV pol Boxplots (3×2 grid + squished bottom)
############################################################

load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

load_or_install("ggplot2")
load_or_install("dplyr")
load_or_install("patchwork")

set.seed(180)

########################
# Load + Clean Data
########################

codons <- read.csv("../Results/codon_features.csv", stringsAsFactors = FALSE)

if (!is.numeric(codons$codon_index)) {
  codons$codon_index <- as.numeric(codons$codon_index)
}

is_complete_codon <- !is.na(codons$ref_aa) &
  !grepl("-", codons$ref_codon) &
  nchar(codons$ref_codon) == 3 &
  codons$ref_aa != "*"

data <- codons[is_complete_codon, ]
data$resistance_site <- as.numeric(data$resistance_site)

if ("dNdS_proxy" %in% names(data)) {
  data$log_dNdS <- log1p(data$dNdS_proxy)
}

predictor_vars <- c(
  "codon_degeneracy", "syn_frac", "nonsyn_frac",
  "aa_chemistry", "shannon_entropy", "conservation",
  if ("log_dNdS" %in% names(data)) "log_dNdS" else "dNdS_proxy"
)

model_cols <- c("codon_index", "resistance_site", predictor_vars)
data_model <- data[, model_cols, drop = FALSE]
data_model <- data_model[complete.cases(data_model), ]

data_model$resistance_site_factor <- factor(
  data_model$resistance_site,
  c(0, 1),
  c("Non-resistance site", "Resistance site")
)

########################
# Top: 6 continuous-variable boxplots (3 × 2 grid)
########################

continuous_vars <- c(
  "shannon_entropy", "conservation", "syn_frac",
  "nonsyn_frac", if ("log_dNdS" %in% names(data_model)) "log_dNdS" else "dNdS_proxy",
  "codon_degeneracy"
)

plots <- lapply(continuous_vars, function(var) {
  ggplot(data_model, aes(x = resistance_site_factor, y = .data[[var]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.3) +
    labs(
      x = "Resistance site",
      y = var,
      title = paste("Boxplot of", var)
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title   = element_text(size = 10),
      axis.text.x  = element_text(size = 8, angle = 20, hjust = 1),
      axis.text.y  = element_text(size = 8)
    )
})

# Create 3×2 layout
combined_top <- wrap_plots(plots, ncol = 3)

########################
# Bottom: AA chemistry (squished)
########################

p_chemistry <- ggplot(data_model, aes(
  x = aa_chemistry,
  y = shannon_entropy,
  fill = resistance_site_factor
)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Amino acid chemistry group",
    y = "Shannon entropy",
    title = "Entropy by AA chemistry group",
    fill = "Resistance site"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 9),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

########################
# Combine (top large, bottom squished)
########################

combined_all <- combined_top / p_chemistry +
  plot_layout(heights = c(2, 1))   # top gets 4× more space

########################
# Save to PDF
########################

pdf("../Results/boxplots_by_resistance_site.pdf",
    width = 12, height = 9)
print(combined_all)
dev.off()
