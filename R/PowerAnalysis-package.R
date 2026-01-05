#' PowerAnalysis: cfDNA Variant Validation Using Power Analysis
#'
#' @description
#' PowerAnalysis provides a comprehensive toolkit for validating tumor-informed
#' variants in cell-free DNA (cfDNA) samples. The package uses binomial power
#' analysis to determine minimum sequencing depth requirements and distinguish
#' true biological negatives from under-sampled sites.
#'
#' @section Main Functions:
#' 
#' **Data Preparation:**
#' * [prepare_bed_file()] - Generate BED/VCF files for force calling
#' 
#' **Power Analysis:**
#' * [run_power_analysis()] - Main analysis pipeline
#' * [calc_proportion()] - Calculate detection proportions
#' 
#' **Visualization:**
#' * [plot_by_organ_site()] - Detection by metastatic organ
#' * [plot_reduction_in_SNVs()] - Tumor fraction vs detection
#' * [compare_detection_rates()] - Comprehensive summary plots
#' 
#' @section Workflow:
#' 
#' 1. Prepare BED files from tumor VCFs
#' 2. Force call cfDNA samples at tumor-informed sites (GATK)
#' 3. Run power analysis to identify adequately powered sites
#' 4. Calculate detection rates by SNV category
#' 5. Generate publication-quality visualizations
#' 
#' @section Citation:
#' 
#' To cite PowerAnalysis in publications, use:
#' 
#'   Vashisth M (2025). PowerAnalysis: Tumor-Informed cfDNA Variant 
#'   Validation Using Binomial Power Analysis. R package version 1.0.0.
#'   https://github.com/manasvitavashisth/PowerAnalysis
#'
#' @docType package
#' @name PowerAnalysis-package
#' @aliases PowerAnalysis
#' @keywords internal
"_PACKAGE"

# Import commonly used functions to avoid namespace issues
#' @import ggplot2
#' @import dplyr
#' @importFrom data.table fread
NULL
