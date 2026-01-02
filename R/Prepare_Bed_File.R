setwd("/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis")
ctdna_tf=read_excel('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/TumorFraction.xlsx')
ctDNA=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_ctDNA.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
ctdna_tf=ctdna_tf[ctdna_tf$Tumor_Fraction>=0.2,]
patient=str_extract(ctdna_tf$Sample, "[^_]+")
tumor=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_tumor_wo_perSample_PowerAnalysis.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
patient=unique(ctDNA$patient)
sample=unique(ctDNA$sample_name)
patient=patient[-c(7,8)]
sample=sample[-c(7,8)]
for(i in 1:length(patient))
{
  bed=tumor[tumor$patient==patient[i],]
  bed=bed[!duplicated(bed$varID),c('Chr','Start','sample_name','Ref','Alt')]
  bed=na.omit(bed)
  bed=bed[order(bed$Chr,bed$Start),]
  colnames(bed)=c('#CHROM','POS','ID','REF','ALT')
  bed$QUAL=rep('.',nrow(bed))
  bed$FILTER=rep('.',nrow(bed))
  bed$INFO=rep('.',nrow(bed))
  writeLines('##fileformat=VCFv4.3\n##fileDate=20240821\n##source=Ensembl\n##reference=GRCh38',paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',ctdna_tf$Sample[i],'.vcf'))
  write.table(bed,file=paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',ctdna_tf$Sample[i],'.vcf'),row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t",append = TRUE)

}

for(i in 1:length(patient))
{
  bed=tumor[tumor$patient==patient[i],]
  bed=bed[!duplicated(bed$varID),c('Chr','Start','sample_name','Ref','Alt')]
  bed=na.omit(bed)
  bed=bed[order(bed$Chr,bed$Start),]
  colnames(bed)=c('#CHROM','POS','ID','REF','ALT')
  bed$QUAL=rep('.',nrow(bed))
  bed$FILTER=rep('.',nrow(bed))
  bed$INFO=rep('.',nrow(bed))
  writeLines('##fileformat=VCFv4.3\n##fileDate=20250510\n##source=Ensembl\n##reference=GRCh38',paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/bed_files/',sample[i],'.vcf'))
  write.table(bed,file=paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/bed_files/',sample[i],'.vcf'),row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t",append = TRUE)
  
}

a=paste0('gatk IndexFeatureFile -I /fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',ctdna_tf$Sample,'.vcf')
write.table(a,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/index.sh',row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

bam=read_excel('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/TopOff_BAM.xlsx')
a=paste0('    ',bam$Sample,':\n        path: ',bam$BAM,'\n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',bam$Sample,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/config/samples.yaml',sep='\n')

a=paste0('gatk IndexFeatureFile -I /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/bed_files/',sample,'.vcf')
write.table(a,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/bed_files/index.sh',row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

a=paste0('    ',sample,':\n        path: /fh/working/ha_g/projects/ProstateTAN/ctDNA_BamsSL/',sample,'/',sample,'.bam \n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/bed_files/',sample,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/config/samples.yaml',sep='\n')


bam=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/Pyclone/all_info.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
bam=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/Pyclone/all_info_merged.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))

bam=bam[bam$sample_id %in% ctdna_tf$Sample,]
a=paste0('    ',ctdna_tf$Sample,':\n        path: ',bam$bam,'\n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',ctdna_tf$Sample,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/config/samples.yaml',sep='\n')

a=paste0('    ',ctdna_tf$Sample,':\n        path: ',bam$bam,'\n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',ctdna_tf$Sample,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/config/samples.yaml',sep='\n')

a=paste0('    ',bam$sample_id,':\n        path: ',bam$bam,'\n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/bed_files/',bam$sample_id,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/PowerAnalysis/config/samples.yaml',sep='\n')


