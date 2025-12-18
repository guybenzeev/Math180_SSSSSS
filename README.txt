HIV pol Codon-Level Resistance Analysis

This repository contains an end-to-end R analysis pipeline for identifying codon-level sequence and evolutionary features associated with drug resistance sites in the HIV-1 pol gene. The workflow includes sequence preprocessing, feature extraction, block-wise logistic regression, cross-validation, and visualization.

Repository Structure and File Descriptions

Data/
hiv-db.fasta
FASTA file containing aligned HIV-1 pol gene sequences used as the primary input for the analysis.

R/
utils/
Helper functions used by multiple scripts.

boxplots_resistance.R
Generates boxplots comparing codon-level predictors between resistance and non-resistance sites.

codon_properties.R
Defines and computes codon-level features, including amino acid chemical class, codon degeneracy, and synonymous/nonsynonymous mutability. Produces a codon-level feature table.

coef_plot.R
Creates coefficient (log-odds) plots for the fitted logistic regression models.

hiv_block_diagnostics.R
Runs diagnostic checks and summaries for block-wise logistic regression models.

hiv_pol_blockwise_logit_reduced.R
Fits block-wise logistic regression models predicting resistance sites using a reduced set of codon-level predictors.

hiv_pol_resistance_blockwise_cv.R
Performs block-wise cross-validation to evaluate model performance and robustness across regions of the pol gene.

predictor_correlation_plot.R
Computes and visualizes correlations among predictors as a heatmap.

preprocess_alignment.R
Preprocesses the sequence alignment, maps codons to reference positions, and prepares data structures for feature extraction.

run_pre.R
Wrapper script that runs all preprocessing and feature extraction steps and saves intermediate outputs.

Build_Data_Table.R.R
Main driver script that collapsses the hiv-db.fasta into codon_features.csv.

Results/
blockwise_logodds_plot.pdf
Visualization of logistic regression coefficients across blocks.

blockwise_model_results.csv
Table of model coefficients and performance metrics from block-wise analyses.

boxplots_by_resistance_site.pdf
Boxplots comparing predictors between resistance and non-resistance codons.

boxplots_by_resistance_site_one_page.pdf
Condensed version of the resistance boxplots on a single page.

codon_features.csv
Final codon-level feature table used for modeling.

predictor_correlation_heatmap.pdf
Heatmap showing correlations among predictors.