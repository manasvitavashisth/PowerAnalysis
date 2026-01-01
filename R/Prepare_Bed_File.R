tumor=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_tumor.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
sample_number=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/Number_of_samples.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
sample_number=sample_number[sample_number$number_of_samples>1,]
sample=sample_number$sample
bam=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/Pyclone_Power_Analysis/all_info.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))


for(i in 1:length(sample))
{
  bed=tumor[tumor$patient==sample_number$patient[i],]
  bed=bed[!duplicated(bed$varID),c('Chr','Start','sample_name','Ref','Alt')]
  bed=na.omit(bed)
  bed=bed[order(bed$Chr,bed$Start),]
  other_sample=sample_number$sample[sample_number$patient==sample_number$patient[i] & sample_number$sample!=sample_number$sample[i]]
  colnames(bed)=c('#CHROM','POS','ID','REF','ALT')
  bed$QUAL=rep('.',nrow(bed))
  bed$FILTER=rep('.',nrow(bed))
  bed$INFO=rep('.',nrow(bed))
  for(j in 1:length(other_sample))
  {
    writeLines('##fileformat=VCFv4.3\n##fileDate=20240821\n##source=Ensembl\n##reference=GRCh38',paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',other_sample[j],'.vcf'))
    write.table(bed,file=paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',other_sample[j],'.vcf'),row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t",append = TRUE)
    }
}
a=paste0('bgzip /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',sample,'.vcf')
write.table(a,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/zip.sh',row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

a=paste0('gatk IndexFeatureFile -I /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',sample,'.vcf')
write.table(a,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/index.sh',row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#Made a bed file with all unique SNVs in the patient
for(i in 1:length(sample))
{
  bed=tumor[tumor$patient==sample_number$patient[i],]
  bed=bed[!duplicated(bed$varID),c('Chr','Start')]
  other_sample=sample_number$sample[sample_number$patient==sample_number$patient[i] & sample_number$sample!=sample_number$sample[i]]
  colnames(bed)=c('chromosome','start')
  bed$end=bed$start
  bed$start=bed$start-1
  for(j in 1:length(other_sample))
  {
    write.table(bed,file=paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',other_sample[j],'.bed'),row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
}
bam=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/Pyclone/all_info.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
bam$titan=paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/optimalTITAN/',file)
bam=bam[bam$sample_id %in% sample,]
a=paste0('    ',sample,':\n        path: ',bam$bam,'\n        intervals_file: /fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/bed_files/',sample,'.vcf')
writeLines(a,'/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis/config/samples.yaml',sep='\n')
