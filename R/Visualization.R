#' Statistical Visualization Functions for cfDNA Detection Analysis
#'
#' Functions for visualizing detection rates across different biological
#' and clinical variables using linear mixed-effects models to account
#' for patient-level clustering.
#'
#' @author Manasvita Vashisth
#' @organization Fred Hutchinson Cancer Center

# Dependencies ----------------------------------------------------------------
library(tidyverse)
library(lme4)          # For linear mixed-effects models
library(scales)        # For percentage formatting
library(viridis)       # For color palettes
library(hrbrthemes)    # For theme_ipsum()

# Main Plotting Functions -----------------------------------------------------

#' Plot Private SNV Detection Rate by Metastatic Organ Site
#'
#' Creates a boxplot showing the detection rate of private SNVs in cfDNA 
#' stratified by metastatic organ site. Uses a linear mixed-effects model
#' (LMER) to test for statistically significant differences between organ
#' sites while accounting for patient-level clustering (multiple samples
#' from the same patient).
#'
#' The statistical model used is:
#'   PrivateSNVs_detected ~ Organ + (1|patient)
#' 
#' Where:
#'   - Fixed effect: Organ (tests if organ site affects detection rate)
#'   - Random effect: (1|patient) (accounts for correlation between samples
#'     from the same patient)
#'
#' @param track Data frame. Detection rate results from run_power_analysis().
#'   Must include columns: Sample, PrivateSNVs_detected, SharedSNVs_detected,
#'   FounderSNVs_detected.
#' @param metadata Data frame. Sample metadata with clinical/biological annotations.
#'   Must include columns: Sample (matching track$Sample), Organ, patient.
#'   Additional columns may include: Phenotype, collection_date, etc.
#' 
#' @return Numeric. P-value from the LMER model testing whether organ site
#'   has a statistically significant effect on private SNV detection rate.
#'   This is the p-value for the fixed effect of the second coefficient
#'   (first organ site vs. reference organ site).
#' 
#' @details
#' The function performs the following steps:
#'   1. Merges detection rates with metadata
#'   2. Fits a linear mixed-effects model with organ as fixed effect
#'   3. Creates a boxplot with overlaid data points
#'   4. Annotates each point with its detection percentage
#'   5. Returns the p-value for statistical significance testing
#'
#' @note The p-value returned is for the second coefficient (first contrast),
#'   which compares the first non-reference organ to the reference organ.
#'   For comprehensive analysis, examine the full model summary.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load detection results
#' track <- read.table("results/detection_rates.tsv", header = TRUE)
#' 
#' # Load metadata
#' metadata <- read.table("data/sample_metadata.tsv", header = TRUE)
#' 
#' # Generate plot and get p-value
#' p_value <- plot_by_organ_site(track, metadata)
#' cat(sprintf("P-value for organ effect: %.4f\n", p_value))
#' }
plot_by_organ_site <- function(track, metadata) {
  
  # Data Preparation ----------------------------------------------------------
  
  # Merge detection rate data with sample metadata
  # This combines the power analysis results with clinical annotations
  # left_join ensures we keep all samples from track, even if some lack metadata
  data <- track %>%
    left_join(metadata, by = 'Sample')
  
  # Check for missing metadata
  n_missing <- sum(is.na(data$Organ))
  if (n_missing > 0) {
    warning(sprintf("%d samples missing organ annotation", n_missing))
  }
  
  # Statistical Analysis ------------------------------------------------------
  
  # Fit linear mixed-effects model (LMER)
  # Model specification:
  #   Response variable: PrivateSNVs_detected (proportion, 0-1)
  #   Fixed effect: Organ (tests if detection rate differs by organ)
  #   Random effect: (1|patient) (random intercept per patient)
  #
  # The random effect accounts for the fact that multiple samples from the
  # same patient are not independent - they share patient-specific factors
  # (e.g., overall tumor burden, genetics, treatment history)
  #
  # This is crucial because ignoring patient clustering would violate the
  # independence assumption and could lead to inflated Type I error rates
  model <- lmer(
    PrivateSNVs_detected ~ Organ + (1|patient),
    data = data
  )
  
  # Extract p-value for the organ effect
  # coef(summary(model)) returns a matrix with:
  #   Rows: model coefficients (intercept, then each organ level)
  #   Columns: Estimate, Std.Error, df, t value, Pr(>|t|)
  # 
  # [, "Pr(>|t|)"] extracts the p-value column
  # [2] gets the second row (first organ contrast vs. reference)
  #
  # Note: This gives the p-value for ONE specific contrast.
  # For a global test across all organs, consider using anova() or
  # a likelihood ratio test comparing the full model to a null model
  p_value <- coef(summary(model))[, "Pr(>|t|)"][2]
  
  # Print model summary for user inspection
  cat("\n=== Linear Mixed-Effects Model Summary ===\n")
  print(summary(model))
  cat("\nP-value for first organ contrast:", p_value, "\n\n")
  
  # Visualization -------------------------------------------------------------
  
  # Create boxplot with individual data points overlaid
  p <- ggplot(data = data, aes(x = Organ, y = PrivateSNVs_detected, fill = Organ)) +
    
    # Add boxplots showing distribution for each organ
    # Box shows: median (center line), IQR (box edges), whiskers (1.5×IQR)
    geom_boxplot() +
    
    # Color scheme: viridis palette (colorblind-friendly, perceptually uniform)
    # discrete = TRUE: use distinct colors for categorical data
    # alpha = 0.6: semi-transparent to see overlapping points
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    
    # Overlay individual data points (jittered to avoid overlap)
    # This shows the actual sample distribution and sample size
    # size = 0.4: small points to avoid crowding
    # alpha = 0.9: slightly transparent
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    
    # Add text labels showing percentage for each data point
    # sprintf formats the proportion as a percentage with 1 decimal place
    # hjust = -0.2: offset text slightly to the right of points
    geom_text(
      aes(label = sprintf("%.1f%%", PrivateSNVs_detected * 100)),
      hjust = -0.2,
      fontface = "bold",
      size = 3  # Text size
    ) +
    
    # Format y-axis as percentages
    # limits = c(0, 1): constrain to valid proportion range
    scale_y_continuous(
      labels = percent,
      limits = c(0, 1)
    ) +
    
    # Apply clean, professional theme
    # theme_ipsum() provides a minimalist, publication-ready appearance
    theme_ipsum() +
    
    # Additional theme customizations
    theme(
      legend.position = "none",        # Remove legend (redundant with x-axis)
      plot.title = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate organ names if needed
    ) +
    
    # Add descriptive labels
    ggtitle("Detection Rate of Private SNVs in cfDNA by Metastatic Organ Site") +
    xlab("Metastatic Organ Site") +
    ylab("Detection Rate in cfDNA")
  
  # Display the plot
  print(p)
  
  # Return p-value for downstream statistical reporting
  return(p_value)
}


#' Compare Detection Rates Across Multiple Variables
#'
#' Wrapper function to generate plots and statistics for multiple grouping
#' variables at once (e.g., organ, phenotype, treatment status).
#'
#' @param track Data frame. Detection rate results
#' @param metadata Data frame. Sample metadata with multiple annotation columns
#' @param variables Character vector. Column names in metadata to test
#' @return Named list of p-values for each variable
#' @export
#' @examples
#' \dontrun{
#' results <- compare_detection_rates(
#'   track = detection_results,
#'   metadata = sample_info,
#'   variables = c("Organ", "Phenotype","CCP score")
#' )
#' }
compare_detection_rates <- function(track, metadata, variables) {
  
  # Store p-values for each variable
  p_values <- list()
  
  for (var in variables) {
    
    cat(sprintf("\n\n=== Analyzing variable: %s ===\n", var))
    
    # Check if variable exists in metadata
    if (!var %in% names(metadata)) {
      warning(sprintf("Variable '%s' not found in metadata. Skipping.", var))
      next
    }
    
    # Prepare data with current variable
    data <- track %>%
      left_join(metadata, by = 'Sample')
    
    # Build model formula dynamically
    formula <- as.formula(paste("PrivateSNVs_detected ~", var, "+ (1|patient)"))
    
    # Fit model
    model <- lmer(formula, data = data)
    
    # Extract p-value
    p_values[[var]] <- coef(summary(model))[, "Pr(>|t|)"][2]
    
    # Create plot
    p <- ggplot(data = data, aes_string(x = var, y = "PrivateSNVs_detected", fill = var)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
      geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
      scale_y_continuous(labels = percent, limits = c(0, 1)) +
      theme_ipsum() +
      theme(legend.position = "none") +
      ggtitle(paste("Detection Rate by", var)) +
      xlab(var) +
      ylab("Detection Rate in cfDNA")
    
    print(p)
  }
  
  # Print summary of results
  cat("\n\n=== Summary of P-values ===\n")
  p_df <- data.frame(
    Variable = names(p_values),
    P_value = unlist(p_values),
    Significant = unlist(p_values) < 0.05
  )
  print(p_df)
  
  return(p_values)
}
#' Plot Relationship Between cfDNA Tumor Fraction and SNV Detection
#'
#' Creates a scatter plot with linear regression showing how cfDNA tumor 
#' fraction (tumor purity in plasma) correlates with the detection rate of
#' private SNVs. This visualization helps determine if samples with higher
#' tumor content in cfDNA show better variant detection rates.
#'
#' @param track Data frame. Detection rate results from run_power_analysis().
#'   Must include columns: Sample (or sample), PrivateSNVs_detected.
#' @param sample_file_path Character. Path to sample manifest file (TSV format).
#'   Required columns: sample, cfdna, patient, tumor_purity, cfdna_tf.
#'   
#' @return ggplot object. Scatter plot with linear regression line showing
#'   relationship between cfDNA tumor fraction and detection rate.
#'
#' @details
#' This analysis addresses the question: "Do samples with higher circulating
#' tumor DNA fractions have better variant detection rates?"
#' 
#' Expected biological relationship:
#'   - Higher cfDNA tumor fraction → More tumor DNA molecules in plasma
#'   - More tumor DNA → Higher probability of sequencing tumor variants
#'   - Result: Positive correlation between cfdna_tf and detection rate
#'
#' Clinical implications:
#'   - Validates that detection rates are biology-driven (not technical artifacts)
#'   - Identifies samples with low tumor fraction that may need deeper sequencing
#'   - Informs minimum tumor fraction requirements for reliable detection
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' p <- plot_reduction_in_SNVs(
#'   track = detection_results,
#'   sample_file_path = "data/sample_manifest.tsv"
#' )
#' print(p)
#' 
#' # Save plot
#' ggsave("cfdna_tf_vs_detection.pdf", p, width = 8, height = 6)
#' 
#' # Extract regression statistics
#' model <- lm(PrivateSNVs_detected ~ cfdna_tf, data = p$data)
#' summary(model)
#' }
plot_reduction_in_SNVs <- function(track, sample_file_path) {
  
  # Load Sample Manifest ------------------------------------------------------
  
  # Read sample manifest containing tumor/cfDNA pairing and purity estimates
  # This file links each cfDNA sample to its tumor purity and tumor fraction
  # 
  # Expected columns:
  #   - sample: tumor sample identifier
  #   - cfdna: cfDNA sample identifier  
  #   - patient: patient identifier (for grouping)
  #   - tumor_purity: estimated purity of tumor biopsy (0-1)
  #   - cfdna_tf: circulating tumor DNA fraction in plasma (0-1)
  #
  # cfdna_tf (cfDNA tumor fraction) represents:
  #   - The proportion of cell-free DNA that originates from tumor cells
  #   - Typically much lower than tumor_purity (e.g., 0.01-0.10 vs 0.5-0.9)
  #   - Key determinant of variant detection sensitivity
  #   - Can be estimated from shallow WGS, copy number, or variant allele frequencies
  tumor <- as.data.frame(fread(
    sample_file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c(".", "NA")
  ))
  
  cat(sprintf("Loaded manifest for %d samples\n", nrow(tumor)))
  
  # Merge Detection Rates with Tumor Fraction Data ---------------------------
  
  # Join the detection rate results (track) with the tumor fraction data (tumor)
  # This combines:
  #   - Detection rates (outcome variable: PrivateSNVs_detected)
  #   - cfDNA tumor fraction (predictor variable: cfdna_tf)
  #
  # Using left_join ensures we keep all samples from track,

data=left_join(track,tumor,by='sample')
ggplot(data, aes(x=cfdna_tf, y=PrivateSNVs_detected)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()
}
