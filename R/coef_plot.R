############################################################
# Coefficient Visualization for Blockwise Logistic Models
# Creates a faceted plot of log(odds ratios) for each block
############################################################

library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------------------------
# 1. Load coefficient table from the reduced model script
#    (Assumes you saved it as `coef_table` like earlier)
# ----------------------------------------------------------

# If sourcing separately, read the object saved from your script.
# Otherwise, comment this out.
load("../Results/coef_table_reduced.rda")

# Remove intercept rows
coef_clean <- coef_table %>%
  filter(term != "(Intercept)") %>%
  mutate(
    log_odds = log(odds_ratio),
    block_id = factor(block_id),
    term = factor(term),
    effect_dir = ifelse(log_odds > 0, "Positive", "Negative")
  )

# ----------------------------------------------------------
# 2. Nice human-readable predictor labels
# ----------------------------------------------------------

pretty_names <- c(
  "syn_frac" = "Synonymous Fraction",
  "conservation" = "Conservation",
  "shannon_entropy" = "Shannon Entropy",
  "codon_degeneracy" = "Codon Degeneracy",
  "aa_chemistryhydrophobic" = "Hydrophobic",
  "aa_chemistryaromatic" = "Aromatic",
  "aa_chemistrypolar" = "Polar",
  "aa_chemistrypositive" = "Positive Charge",
  "aa_chemistrynegative" = "Negative Charge",
  "aa_chemistryspecial" = "Special (G,P,C)"
)

coef_clean$term_label <- pretty_names[as.character(coef_clean$term)]

# ----------------------------------------------------------
# 3. Plot: log(odds ratio) per predictor per block
# ----------------------------------------------------------

p <- ggplot(coef_clean, aes(x = term_label, y = log_odds, fill = effect_dir)) +
  geom_col() +
  facet_wrap(~ block_id, ncol = 2) +
  scale_fill_manual(
    values = c(Positive = "#EF4444", Negative = "#3B82F6"),
    name   = "Effect Direction"
  ) +
  labs(
    title = "Log-Odds Coefficients by Predictors Across Blocks",
    x = "Predictor",
    y = "Log(odds ratio)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# ----------------------------------------------------------
# 4. Display in RStudio & save to PDF
# ----------------------------------------------------------

print(p)

ggsave("../Results/blockwise_logodds_plot.pdf", plot = p, width = 12, height = 8)

cat("\nSaved coefficient plot to ../Results/blockwise_logodds_plot.pdf\n")
