#' Preparing BED files at tumor informed SNVs
#'
#' This module prepares BED files as input for mutect force calling or Collect Allelic Counts
#'
#' @author Manasvita Vashisth
#' @organization Fred Hutchinson Cancer Center

# Dependencies ----------------------------------------------------------------
library(tidyverse)
library(data.table)

# Constants -------------------------------------------------------------------

# VCF header information
VCF_HEADER <- '##fileformat=VCFv4.3\n##fileDate=20250510\n##source=Ensembl\n##reference=GRCh38'

# Main Functions --------------------------------------------------------------

#' Prepare BED file
#'
#' @param tumor_file_path. Path to Concatenated Tumor SNV file for all patients
#' @param output_file_path. Path for output BED files
#' @param patient_file_path. Path to list of patients and corresponding cfDNA samples
#' @param validate_input Logical. Perform input validation (default: TRUE).
#' @export

prepare_bed_file <- function(tumor_file_path,output_file_path,patient_file_path,validate_input = TRUE) {
  # Validate inputs -----------------------------------------------------------
  if (validate_input) {
    validate_file_paths(tumor_file_path, patient_file_path)
  }
  
  # Load input data -----------------------------------------------------------
  cat("Loading input files...\n")

  # Read tumor SNV data
  tumor=as.data.frame(fread(tumor_file_path,header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  # Read patient list
  patient=as.data.frame(fread(patient_file_path,header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  
  # Ensure output directory exists
  if (!dir.exists(output_file_path)) {
    dir.create(output_file_path, recursive = TRUE)
    cat(sprintf("  Created output directory: %s\n", output_file_path))
  }
  for(i in 1:nrow(patient))
  {
    bed=tumor[tumor$patient==patient$Patient[i],]
    bed=bed[!duplicated(bed$varID),c('Chr','Start','sample_name','Ref','Alt')]
    bed=na.omit(bed)
    bed=bed[order(bed$Chr,bed$Start),]
    colnames(bed)=c('#CHROM','POS','ID','REF','ALT')
    bed$QUAL=rep('.',nrow(bed))
    bed$FILTER=rep('.',nrow(bed))
    bed$INFO=rep('.',nrow(bed))
    writeLines(VCF_HEADER,paste0(output_file_path,patient$sample[i],'.vcf'))
    write.table(bed,file=paste0(output_file_path,patient$sample[i],'.vcf'),row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t",append = TRUE)
    
  }

}

# Helper Functions ------------------------------------------------------------

#' Validate File Paths
#'
#' @param tumor_file_path Character. Tumor file path
#' @param patient_file_path. Path to list of patients and corresponding cfDNA samples
#' @keywords internal
validate_file_paths <- function(tumor_file_path, patient_file_path) {
  
  # Check input files exist
  if (!file.exists(tumor_path)) {
    stop("Tumor file not found: ", tumor_file_path)
  }
  if (!file.exists(patient_path)) {
    stop("Patient file not found: ", patient_file_path)
  }
}
