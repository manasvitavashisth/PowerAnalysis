#' Power Analysis for cfDNA Variant Validation
#'
#' This module performs binomial power analysis to determine if cfDNA sequencing
#' depth is sufficient to detect tumor-informed variants. It processes Mutect2
#' force calling output and PyClone clonality data to calculate detection rates
#' for founder, shared, and private SNVs.
#'
#' @author Manasvita Vashisth
#' @organization Fred Hutchinson Cancer Center

# Dependencies ----------------------------------------------------------------
library(tidyverse)
library(data.table)

# Constants -------------------------------------------------------------------

# Default minimum number of alternate reads required for variant detection
DEFAULT_MIN_ALT_READS <- 3

# Default probability cutoff for considering a site to have sufficient power
# (0.8 = 80% probability of detecting min_alt_reads or more)
DEFAULT_PROB_CUTOFF <- 0.8

# Minimum read depth threshold for including a site in analysis
MIN_READ_DEPTH <- 3


# Main Functions --------------------------------------------------------------

#' Calculate Proportion of SNVs Detected in cfDNA
#'
#' Determines what fraction of tumor SNVs (that have sufficient sequencing power)
#' are actually detected in cfDNA for a specific SNV category (founder/shared/private).
#' 
#' The detection rate is calculated as:
#'   (Number of SNVs with >0 alt reads in cfDNA) / (Total SNVs with sufficient power)
#'
#' @param tumor Data frame. PyClone output containing SNVs with sufficient power
#'   to be detected. Must include columns: cluster_id, mutation_id.
#' @param cfdna Data frame. cfDNA variants with at least >1 alternate allele 
#'   detected at tumor-informed sites. Must include column: varID.
#' @param cluster Data frame. Maps cluster IDs to SNV categories. 
#'   Must include columns: cluster_id, Cluster (category type).
#' @param type Character. SNV category to calculate proportion for 
#'   (e.g., "Founder", "Shared", "Private").
#' 
#' @return Numeric. Proportion of SNVs detected (0 to 1).
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' proportion_founder <- calc_proportion(
#'   tumor = pyclone_data,
#'   cfdna = detected_variants,
#'   cluster = cluster_assignments,
#'   type = "Founder"
#' )
#' }
calc_proportion <- function(tumor, cfdna, cluster, type) {
  
  # Get all cluster IDs that match the requested SNV type (e.g., all "Founder" clusters)
  cluster_type <- cluster$cluster_id[cluster$Cluster == type]
  
  # Filter tumor data to only include SNVs from the specified cluster type
  tumor_type <- tumor[tumor$cluster_id %in% cluster_type, ]
  
  # Filter cfDNA data to only include variants that match the tumor SNVs of this type
  # This gives us the subset of tumor SNVs that were actually detected in cfDNA
  cfdna_type <- cfdna[cfdna$varID %in% tumor_type$mutation_id, ]
  
  # Remove duplicate mutations from tumor data (some may appear multiple times
  # if they're in multiple samples/clusters)
  tumor_type <- tumor_type[!duplicated(tumor_type$mutation_id), ]
  
  # Calculate proportion: detected SNVs / total SNVs with sufficient power
  proportion <- nrow(cfdna_type) / nrow(tumor_type)
  
  return(proportion)
}


#' Run Power Analysis on cfDNA Samples
#'
#' Main function to perform comprehensive power analysis across multiple cfDNA samples.
#' For each sample, this function:
#'   1. Loads Mutect2 force calling results
#'   2. Loads PyClone clonality information
#'   3. Calculates which tumor SNVs have sufficient sequencing depth for detection
#'   4. Applies binomial power analysis to identify adequately powered sites
#'   5. Determines detection rates for founder, shared, and private SNVs
#'
#' The binomial power calculation uses:
#'   P(X ≥ min_alt_reads) = 1 - pbinom(min_alt_reads - 1, depth, p)
#'   where p = tumor_VAF × (cfDNA_purity / tumor_purity)
#'
#' @param mutect_force_calling_path Character. Path to directory containing 
#'   Mutect2 force calling output folders (one per cfDNA sample). Each folder
#'   should contain 'mutations_unfiltered.hg38_multianno.txt'.
#' @param output_file_path Character. Full path for output TSV file containing
#'   detection rates for all samples.
#' @param sample_file_path Character. Path to sample manifest file (TSV format).
#'   Required columns: sample (tumor sample ID), cfdna (cfDNA sample ID), 
#'   patient (patient ID), tumor_purity (0-1), cfdna_tf (cfDNA tumor fraction, 0-1).
#' @param snv_list Character. Path to directory containing PyClone output files.
#'   Expected structure: {patient_id}/pyclone/{patient_id}_pyclone_output.tsv
#'   PyClone output should include: mutation_id, sample_id, cluster_id, vaf.
#' @param min_alt_reads Integer. Minimum number of alternate reads required
#'   to call a variant as detected (default: 3).
#' @param prob_cutoff Numeric. Minimum probability threshold for considering
#'   a site to have sufficient power (default: 0.8 = 80%).
#' @param snv_type Character. Path to file mapping cluster IDs to SNV categories
#'   (Founder/Shared/Private). Required columns: sample, cluster_id, Cluster.
#' 
#' @return Writes results to output_file_path. Returns invisible data frame with
#'   detection rates for each sample and SNV category.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' run_power_analysis(
#'   mutect_force_calling_path = "results/mutect_force_calling/",
#'   output_file_path = "output/detection_rates.tsv",
#'   sample_file_path = "data/sample_manifest.tsv",
#'   snv_list = "results/pyclone/",
#'   min_alt_reads = 3,
#'   prob_cutoff = 0.8,
#'   snv_type = "data/cluster_categories.tsv"
#' )
#' }
run_power_analysis <- function(mutect_force_calling_path, 
                               output_file_path,
                               sample_file_path,
                               snv_list,
                               min_alt_reads = DEFAULT_MIN_ALT_READS,
                               prob_cutoff = DEFAULT_PROB_CUTOFF,
                               snv_type) {
  
  # Load Input Data -----------------------------------------------------------
  
  cat("Loading input data...\n")
  
  # Read sample manifest containing tumor/cfDNA pairing and purity estimates
  # Expected columns: sample, cfdna, patient, tumor_purity, cfdna_tf
  tumor <- as.data.frame(fread(
    sample_file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c(".", "NA")
  ))
  
  # Read cluster/clonal identity mapping
  # Expected columns: sample, cluster_id, Cluster (Founder/Shared/Private)
  clusters <- as.data.frame(fread(
    snv_type,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c(".", "NA")
  ))
  
  # Validate Input Data -------------------------------------------------------
  
  # Ensure all required columns are present
  validate_input_data(tumor, clusters)
  
  # Extract sample identifiers
  sample <- tumor$sample
  
  # Initialize Results Table --------------------------------------------------
  
  # Create tracking table to store detection rates for each sample
  # Columns: Sample ID, and proportion detected for each SNV category
  track <- as.data.frame(matrix(
    data = NA,
    nrow = length(sample),
    ncol = 4
  ))
  colnames(track) <- c(
    'Sample',
    'Reduction_in_SNVs',
    'PrivateSNVs_detected',
    'SharedSNVs_detected',
    'FounderSNVs_detected'
  )
  track$Sample <- sample
  
  cat("\nProcessing samples...\n")
  
  # Process Each Sample -------------------------------------------------------
  
  # Iterate through each tumor sample and its corresponding cfDNA sample
  for (i in seq_along(sample)) {
     
    # Load Mutect2 Force Calling Results --------------------------------------
    
    # Check if file exists
    if (!file.exists(mutect_force_calling_path)) {
      warning(sprintf("Mutect file not found for sample %s: %s", 
                     tumor$cfdna[i], mutect_file))
      next
    }
    
    # Read Mutect2 force calling output
    # This contains read depth and allele counts at all tumor-informed sites
    mutect <- as.data.frame(fread(
      mutect_force_calling_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      na.strings = c(".", "NA")
    ))
    
    # Parse Read Depth from Mutect2 Output ------------------------------------
    
    # Extract total read depth from the FORMAT field (Otherinfo13)
    # FORMAT field structure: GT:AD:AF:DP (e.g., "0/1:20,5:0.2:25")
    # We extract the 4th field (DP = total depth)
    mutect$mutect_force_call_depth <- sapply(
      strsplit(mutect$Otherinfo13, ':'),
      function(x) ifelse(length(x) >= 4, x[4], NA)
    )
    
    # Parse Alternate Allele Depth --------------------------------------------
    
    # Extract alternate allele depth from the AD field (Allelic Depth)
    # AD field structure: "ref_depth,alt_depth" (e.g., "20,5")
    # We extract the second value (alt_depth)
    mutect$alt_depth <- sapply(
      strsplit(mutect$Otherinfo13, ':'),
      function(x) ifelse(length(x) >= 4, x[2], NA)
    )
    
    # Extract just the alternate depth (after the comma)
    # sub() replaces everything up to and including the comma with nothing
    mutect$alt_depth <- as.numeric(sub(".*,\\s*", "", mutect$alt_depth))
    
    # Create variant ID for matching with PyClone data
    # Format: chr_position_ref_alt (e.g., "chr1_12345_A_T")
    mutect$varID <- paste(mutect$Chr, mutect$Start, mutect$Ref, mutect$Alt, sep = '_')
    
    # Load PyClone Clonality Data ---------------------------------------------
    
    # Get cluster assignments for this specific sample
    s_clusters <- clusters[clusters$sample == sample[i], ]
    
    # Construct path to PyClone output for this patient
    pyclone_file <- file.path(
      snv_list,
      tumor$patient[i],
      'pyclone',
      paste0(tumor$patient[i], '_pyclone_output.tsv')
    )
    
    # Check if file exists
    if (!file.exists(pyclone_file)) {
      warning(sprintf("PyClone file not found for patient %s: %s", 
                     tumor$patient[i], pyclone_file))
      next
    }
    
    # Read PyClone output containing clonality and VAF information
    pyclone <- as.data.frame(fread(
      pyclone_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      na.strings = c(".", "NA")
    ))
    cat(sprintf("  Loaded %d PyClone variants\n", nrow(pyclone)))
    
    # Merge PyClone with Mutect Data ------------------------------------------
    
    # Join PyClone clonality data with Mutect force calling depth data
    # This gives us VAF + clonality info + cfDNA read depth for each variant
    pyclone <- left_join(pyclone, mutect, by = 'varID')
    
    # Convert depth to numeric (in case it was read as character)
    pyclone$mutect_force_call_depth <- as.numeric(pyclone$mutect_force_call_depth)
    
    # Filter by Minimum Read Depth --------------------------------------------
    
    # Only keep sites with sufficient read depth (>3 reads by default)
    # Sites with very low depth don't provide reliable information
    pyclone_reduced <- pyclone[pyclone$mutect_force_call_depth > MIN_READ_DEPTH, ]
    
    # Remove sites with missing depth information
    pyclone_reduced <- pyclone_reduced[!is.na(pyclone$mutect_force_call_depth), ]
    
    # Calculate Expected Allele Frequency in cfDNA ----------------------------
    
    # Calculate the "dilution factor" for tumor variants in cfDNA
    # This accounts for the lower tumor content in cfDNA vs. tumor tissue
    # Formula: expected_cfDNA_VAF = tumor_VAF × (cfDNA_purity / tumor_purity)
    fraction <- tumor$cfdna_tf[i] / tumor$tumor_purity[i]
    
    
    # Apply Binomial Power Analysis -------------------------------------------
    
    # Calculate detection power for each variant using binomial distribution
    # Power = P(observing ≥ min_alt_reads | depth, expected_VAF)
    #       = 1 - P(observing < min_alt_reads)
    #       = 1 - pbinom(min_alt_reads - 1, depth, prob)
    # 
    # Where prob = tumor_VAF × (cfDNA_purity / tumor_purity)
    #
    # Only keep variants where power ≥ prob_cutoff (default: 0.8)
    # These are sites where we have sufficient depth to reliably detect
    # the variant if it's truly present in cfDNA
    
    pyclone_reduced <- pyclone_reduced[
      (1 - pbinom(
        (min_alt_reads - 1),              # x - 1 (for P(X ≥ x))
        pyclone_reduced$mutect_force_call_depth,  # n (number of trials = read depth)
        pyclone_reduced$vaf * fraction            # p (expected probability of alt allele)
      )) >= prob_cutoff,                   # Only keep if power ≥ threshold
    ]
     
    # Calculate proportion of SNV sites reduced due to insufficient read depth
    track[i,'Reduction_in_SNVs']=(nrow(pyclone)-nrow(pyclone_reduced))/nrow(pyclone)
    # Identify variants that were actually detected in cfDNA
    # "Detected" = has at least 1 alternate read (alt_depth > 0)
    ctdna_detected <- pyclone_reduced[pyclone_reduced$alt_depth > 0, ]
    
     
    # Calculate Detection Rates by SNV Category -------------------------------
    
    # Calculate what proportion of each SNV category was detected
    # This answers: "Of all the founder/shared/private SNVs we had power to 
    # detect, what fraction did we actually detect in cfDNA?"
    
    # Private SNVs (found in only one metastatic site)
    track[i, 'PrivateSNVs_detected'] <- calc_proportion(
      tumor = pyclone_reduced,
      cfdna = ctdna_detected,
      cluster = s_clusters,
      type = 'Private'
    )
     
    # Shared SNVs (found in multiple but not all metastatic sites)
    track[i, 'SharedSNVs_detected'] <- calc_proportion(
      tumor = pyclone_reduced,
      cfdna = ctdna_detected,
      cluster = s_clusters,
      type = 'Shared'
    )
    
    # Founder SNVs (found in all metastatic sites - likely early/trunk mutations)
    track[i, 'FounderSNVs_detected'] <- calc_proportion(
      tumor = pyclone_reduced,
      cfdna = ctdna_detected,
      cluster = s_clusters,
      type = 'Founder'
    )
  }
  
  # Save Results --------------------------------------------------------------
  
  cat("\n=== Writing results ===\n")
  
  # Write detection rates to output file
  write.table(
    track,
    file = output_file_path,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
  
  cat(sprintf("Results saved to: %s\n", output_file_path))
  
  # Print Summary Statistics --------------------------------------------------
  
  cat("\n=== Summary Statistics ===\n")
  
  # Calculate mean detection rates across all samples
  mean_private <- mean(track$PrivateSNVs_detected, na.rm = TRUE)
  mean_shared <- mean(track$SharedSNVs_detected, na.rm = TRUE)
  mean_founder <- mean(track$FounderSNVs_detected, na.rm = TRUE)
  
  cat(sprintf("Mean detection rates across %d samples:\n", 
              sum(!is.na(track$PrivateSNVs_detected))))
  cat(sprintf("  Private SNVs:  %.1f%%\n", mean_private * 100))
  cat(sprintf("  Shared SNVs:   %.1f%%\n", mean_shared * 100))
  cat(sprintf("  Founder SNVs:  %.1f%%\n", mean_founder * 100))
  
  # Return results invisibly
  invisible(track)
}


# Helper Functions ------------------------------------------------------------

#' Validate Input Data
#'
#' Checks that input data frames contain all required columns.
#' Stops execution with informative error if validation fails.
#'
#' @param tumor Data frame. Sample manifest with purity estimates
#' @param clusters Data frame. Cluster to SNV category mapping
#' @keywords internal
validate_input_data <- function(tumor, clusters) {
  
  # Define required columns for sample manifest
  # - sample: tumor sample identifier
  # - cfdna: corresponding cfDNA sample identifier
  # - patient: patient identifier (links to PyClone data)
  # - cfdna_tf: cfDNA tumor fraction (estimated proportion of tumor DNA in cfDNA)
  # - tumor_purity: tumor purity estimate (proportion of tumor cells in tumor sample)
  required_tumor_cols <- c("sample", "cfdna", "patient", "cfdna_tf", "tumor_purity")
  
  # Define required columns for cluster assignments
  # - sample: sample identifier (matches tumor$sample)
  # - cluster_id: PyClone cluster identifier
  # - Cluster: SNV category (Founder/Shared/Private)
  required_clusters_cols <- c("sample", "cluster_id", "Cluster")
  
  # Check tumor manifest has all required columns
  if (!all(required_tumor_cols %in% names(tumor))) {
    missing_cols <- setdiff(required_tumor_cols, names(tumor))
    stop(
      "Sample manifest (tumor data) missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nExpected columns: ",
      paste(required_tumor_cols, collapse = ", ")
    )
  }
  
  # Check cluster data has all required columns
  if (!all(required_clusters_cols %in% names(clusters))) {
    missing_cols <- setdiff(required_clusters_cols, names(clusters))
    stop(
      "Cluster data missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nExpected columns: ",
      paste(required_clusters_cols, collapse = ", ")
    )
  }
  
  # Validate purity values are in valid range (0, 1]
  if (any(tumor$tumor_purity <= 0 | tumor$tumor_purity > 1, na.rm = TRUE)) {
    warning("Some tumor_purity values are outside valid range (0, 1]")
  }
  
  if (any(tumor$cfdna_tf < 0 | tumor$cfdna_tf > 1, na.rm = TRUE)) {
    warning("Some cfdna_tf values are outside valid range [0, 1]")
  }
  
  # Check for valid SNV category types
  valid_categories <- c("Founder", "Shared", "Private")
  invalid_categories <- setdiff(unique(clusters$Cluster), valid_categories)
  if (length(invalid_categories) > 0) {
    warning(
      "Found unexpected SNV categories: ",
      paste(invalid_categories, collapse = ", "),
      "\nExpected: ",
      paste(valid_categories, collapse = ", ")
    )
  }
}


# Utility Functions -----------------------------------------------------------

#' Calculate Power for Single Variant
#'
#' Helper function to calculate detection power for a single variant
#' using the binomial distribution.
#'
#' @param depth Integer. Read depth at the locus
#' @param vaf Numeric. Variant allele frequency in tumor (0-1)
#' @param cfDNA_purity Numeric. Tumor fraction in cfDNA (0-1)
#' @param tumor_purity Numeric. Tumor purity in tissue sample (0-1)
#' @param min_alt_reads Integer. Minimum alternate reads for detection
#' @return Numeric. Probability of detecting min_alt_reads or more (0-1)
#' @export
#' @examples
#' \dontrun{
#' power <- calculate_single_variant_power(
#'   depth = 100,
#'   vaf = 0.5,
#'   cfDNA_purity = 0.05,
#'   tumor_purity = 0.7,
#'   min_alt_reads = 3
#' )
#' }
calculate_single_variant_power <- function(depth,
                                          vaf,
                                          cfDNA_purity,
                                          tumor_purity,
                                          min_alt_reads = DEFAULT_MIN_ALT_READS) {
  
  # Calculate expected probability in cfDNA
  expected_prob <- vaf * (cfDNA_purity / tumor_purity)
  
  # Ensure probability is in valid range
  expected_prob <- pmax(0, pmin(1, expected_prob))
  
  # Calculate power using binomial distribution
  # P(X ≥ k) = 1 - P(X < k) = 1 - pbinom(k-1, n, p)
  power <- 1 - pbinom(min_alt_reads - 1, depth, expected_prob)
  
  return(power)
}


#' Generate Summary Report
#'
#' Creates a formatted summary report of detection rates.
#'
#' @param results_file Character. Path to results TSV file
#' @param output_file Character. Path for summary report (optional)
#' @export
generate_summary_report <- function(results_file, output_file = NULL) {
  
  # Read results
  results <- fread(results_file)
  
  # Calculate statistics
  summary_stats <- list(
    n_samples = nrow(results),
    mean_private = mean(results$PrivateSNVs_detected, na.rm = TRUE),
    mean_shared = mean(results$SharedSNVs_detected, na.rm = TRUE),
    mean_founder = mean(results$FounderSNVs_detected, na.rm = TRUE),
    sd_private = sd(results$PrivateSNVs_detected, na.rm = TRUE),
    sd_shared = sd(results$SharedSNVs_detected, na.rm = TRUE),
    sd_founder = sd(results$FounderSNVs_detected, na.rm = TRUE)
  )
  
  # Create report text
  report <- sprintf(
    "Detection Rate Summary Report
=============================

Number of samples: %d

Mean Detection Rates:
  Private SNVs:  %.1f%% (SD: %.1f%%)
  Shared SNVs:   %.1f%% (SD: %.1f%%)
  Founder SNVs:  %.1f%% (SD: %.1f%%)

Statistical Comparison:
  Founder vs Private: %.1f%% higher detection
  Founder vs Shared:  %.1f%% higher detection
  Shared vs Private:  %.1f%% higher detection
",
    summary_stats$n_samples,
    summary_stats$mean_private * 100, summary_stats$sd_private * 100,
    summary_stats$mean_shared * 100, summary_stats$sd_shared * 100,
    summary_stats$mean_founder * 100, summary_stats$sd_founder * 100,
    (summary_stats$mean_founder - summary_stats$mean_private) * 100,
    (summary_stats$mean_founder - summary_stats$mean_shared) * 100,
    (summary_stats$mean_shared - summary_stats$mean_private) * 100
  )
  
  # Print to console
  cat(report)
  
  # Save to file if requested
  if (!is.null(output_file)) {
    writeLines(report, output_file)
    cat(sprintf("\nReport saved to: %s\n", output_file))
  }
  
  invisible(summary_stats)
}


