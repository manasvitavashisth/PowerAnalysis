#' Power Analysis for cfDNA Variant Validation
#'
#' This module performs binomial power analysis to determine if cfDNA sequencing
#' depth is sufficient to detect tumor-informed variants.
#'
#' @author Manasvita Vashisth
#' @organization Fred Hutchinson Cancer Center

# Dependencies ----------------------------------------------------------------
library(tidyverse)
library(data.table)

# Constants -------------------------------------------------------------------
DEFAULT_MIN_ALT_READS <- 3
DEFAULT_PROB_CUTOFF <- 0.8

# Main Functions --------------------------------------------------------------

#' Calculate proportion of SNVs detected in cfDNA
#'
#' @param tumor Data Frame with tumor SNVs for required category
#' @param cfdna Data Frame with atleast >1 alt allele detected in tumor informed sites
#' @return Data frame with power analysis results
#' @export
#' @examples
#' calc_expected_prob(0.5, 0.7, 0.05)
calc_expected_prob <- function(tumor, cfdna) {
  if (tumor_purity == 0) {
    warning("Tumor purity is zero. Returning NA.")
    return(NA)
  }
  
  expected_prob <- vaf_tumor * (cfdna_purity / tumor_purity)
  
  # Ensure probability is within valid range
  expected_prob <- pmin(pmax(expected_prob, 0), 1)
  
  return(expected_prob)
}

#' Calculate Expected Probability of Observing Mutant Allele
#'
#' @param mutect_force_calling_path. Path to folder with output of Mutect force Calling
#' @param output_file_path. Path for output file
#' @param sample_file_path. Path to list of tumor samples, patient ids, corresponding cfDNA samples, tumor purity and cfDNA tumor fraction (See data folder for example format)
#' @param snv_list. Path to pyclone output or list of SNVs in the patient and cluster/clonal identity
#' @param snv_type. Path to file with cluster/clonal identity (Founder, Shared, Private)
#' @export

run_power_analysis <- function(mutect_force_calling_path, 
                               output_file_path,
                               sample_file_path,
                               snv_list,
                               min_alt_reads = DEFAULT_MIN_ALT_READS,
                               prob_cutoff = DEFAULT_PROB_CUTOFF,
                               snv_type,) {
  
  
  
  # Read sample data file
  tumor=as.data.frame(fread(sample_file_path,header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  # Read cluster/clonal identity data file
  clusters=as.data.frame(fread(snv_type,header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))

  # Validate inputs
  validate_input_data(tumor, clusters)
  
  # Store list of each tumor sample in a new variable
  sample=tumor$sample
  
 for(i in 1:length(sample)
{
  mutect=as.data.frame(fread(paste0(mutect_force_calling_path,tumor$cfdna[i],'/mutations_unfiltered.hg38_multianno.txt'),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  mutect$mutect_force_call_depth=sapply(strsplit(mutect$Otherinfo13, ':'), function(x) ifelse(length(x) >= 4, x[4], NA))
  mutect$alt_depth=sapply(strsplit(mutect$Otherinfo13, ':'), function(x) ifelse(length(x) >= 4, x[2], NA))
  mutect$alt_depth=as.numeric(sub(".*,\\s*", "", mutect$alt_depth))
  mutect$varID=paste(mutect$Chr,mutect$Start,mutect$Ref,mutect$Alt,sep='_')
  s_clusters=clusters[clusters$sample==sample[i],]
  pyclone=as.data.frame(fread(paste0(snv_list,tumor$patient[i],'/pyclone/',tumor$patient[i],'_pyclone_output.tsv'),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  pyclone=left_join(pyclone,m,by='varID')
  pyclone$mutect_force_call_depth=as.numeric(pyclone$mutect_force_call_depth)
  pyclone=pyclone[comb$mutect_force_call_depth>3,]
  pyclone=pyclone[!is.na(pyclone$mutect_force_call_depth),]
  fraction=tumor$cfdna_tf[i]/tumor$tumor_purity[i]
  pyclone=pyclone[(1-pbinom((min_alt_reads-1),pyclone$mutect_force_call_depth,pyclone$vaf*fraction))>=prob_cutoff,]
  samples=unique(pyclone$sample_id)
  ctDNAp=pyclone[pyclone$alt_depth>0,]
  
  shared=s_clusters$cluster[s_clusters$Cluster=='Shared']
  pyclone_shared=pyclone[pyclone$cluster_id %in% shared,]
  ct_s=ctDNAp[ctDNAp$varID %in% pyclone_shared$mutation_id,]
  pyclone_shared=pyclone_shared[unique(pyclone_shared$mutation_id),]
  tf$overlap_shared[i]=nrow(ct_s)/nrow(pyclone_shared)
  
  private=tumor[tumor$sample_name==sample[i],]
  track[i,2]=nrow(private)
  private_overlap=ctDNA[ctDNA$varID %in% private$varID & ctDNA$patient==private$patient[1],]
  ctDNA1=ctDNA[ctDNA$patient==sample_patient[i],]
  track[i,3]=nrow(private_overlap)
  fraction=ctdna_tf[ctdna_tf$patient==sample_patient[i],2]/tumor_tf[tumor_tf$Sample==sample[i],2]
  #private=private[private$vaf*fraction[1,1]>=0.1,]
  comb=left_join(private,m,by='varID')
  comb$mutect_force_call_depth=as.numeric(comb$mutect_force_call_depth)
  comb1=comb[comb$mutect_force_call_depth>3,]
  comb1=comb1[!is.na(comb1$mutect_force_call_depth),]
  comb1=comb1[(1-pbinom(2,comb1$mutect_force_call_depth,comb1$vaf*fraction[1,1]))>=0.8,]
  #comb3=comb1[(1-binom.test(x=2,n=comb1$depth.y,p=comb1$vaf*fraction[1,1],alternative = 'greater'))>=0.8,]
  track[i,4]=nrow(comb1)
  comb2=comb1[comb1$alt_depth>0 & comb1$Alt.x==comb1$Alt.y,]
  track[i,5]=nrow(comb2)
  comb3=comb1[comb1$varID %in% ctDNA1$varID,]
  track[i,6]=nrow(comb3)
}
track[,7]=track[,3]/track[,2]
track[,8]=track[,5]/track[,4]
track[,9]=track[,6]/track[,4]
colnames(track)=c('sample','TumorSNVs','UnadjustedOverlap','ObservableSNVs','Overlap_mutect_force_call_alt_depth','Overlap_ctDNA_calling','Prop_Unadjusted','Prop_mutect_force_Call','Prop_ctDNA')

write.table(track,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PatientPrivateMutations/ctDNA_Tumor_Overlap.txt',row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

# Helper Functions ------------------------------------------------------------

#' Validate Input Data
#' @keywords internal
validate_input_data <- function(tumor, clusters) {
  required_tumor_cols <- c("sample","cfdna","patient","cfdna_tf","tumor_purity")
  required_clusters_cols <- c("sample","cluster_id","cluster_type")
  
  if (!all(required_tumor_cols %in% names(tumor))) {
    stop("Tumor data missing required columns: ", 
         paste(setdiff(required_tumor_cols, names(tumor)), collapse = ", "))
  }
  
  if (!all(required_clusters_cols %in% names(clusters))) {
    stop("cfDNA data missing required columns: ", 
         paste(setdiff(required_clusters_cols, names(clusters)), collapse = ", "))
  }
}
