#a=list.files('/fh/working/ha_g/projects/ProstateTAN/ctDNA_BamsSL/')
#write.table(paste0(a,': /fh/working/ha_g/projects/ProstateTAN/ctDNA_BamsSL/',a,'/',a,'.bam'),file='/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/samples_bampath.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)

#Compare Force calling mutect to samples
library(lmerTest)
library(tidyverse)
library(vcfR)
library(UpSetR)
library(ggplot2)
library(dplyr)
library(limma)
library(VennDiagram)
library(RColorBrewer)
library(ggVennDiagram)
library(maftools)
library(data.table)
library(stringr)
library(hrbrthemes)
library(viridis)
library(geomtextpath)

tumor=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_tumor_wo_perSample_PowerAnalysis.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
number=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/Number_of_samples_Per_Patient.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
#tumor=as.data.frame(fread('/Volumes/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_tumor_wo_perSample_PowerAnalysis.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
#number=as.data.frame(fread('/Volumes/projects/ProstateTAN/analysis_mutational/Number_of_samples_Per_Patient.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
ctDNA=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_ctDNA/SNVs/MasterSNVFile_ctDNA.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))

number=number[number$Number_of_Samples>1,]
#sample_ctDNA=unique(ctDNA$sample_name)
#sample_ctDNA=sample_ctDNA[order(sample_ctDNA)]
#patient_ctDNA=str_extract(sample_ctDNA, "[^_]+")
#sample_ctDNA=sample_ctDNA[patient_ctDNA %in% number$Patient]
ctdna_tf=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/Sample_Ploidy_Purity_ctDNA.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
ctdna_tf$patient=str_extract(ctdna_tf$Sample, "[^_]+")
ctdna_tf=ctdna_tf[order(ctdna_tf$Sample),]
ctdna_tf=ctdna_tf[ctdna_tf$Purity>=0.1,]
sample_ctDNA=ctdna_tf$Sample
patient_ctDNA=str_extract(sample_ctDNA, "[^_]+")
sample_ctDNA=sample_ctDNA[patient_ctDNA %in% number$Patient]
ctdna_tf=ctdna_tf[ctdna_tf$patient %in% number$Patient,]
patient_ctDNA=str_extract(sample_ctDNA, "[^_]+")

tumor_tf=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/rerun2025_Sample_Ploidy_Purity.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
tumor_tf=tumor_tf[tumor_tf$Purity>=0.1,]
tumor=tumor[tumor$patient %in% patient_ctDNA,]
sample=unique(tumor$sample_name)
sample=sample[sample %in% tumor_tf$Sample]
sample=sample[order(sample)]
sample_patient=str_extract(sample, "[^_]+")
track=as.data.frame(matrix(data=NA,nrow=length(sample),ncol=14))
colnames(track)[1]='Sample'
track$Sample=sample

patient=unique(sample_patient)
for(i in 1:length(patient))
{
  pt=tumor[tumor$patient==patient[i],]
  pt_sample=unique(pt$sample_name)
  s1=tumor[tumor$sample_name==pt_sample[1],]
  s2
}

for(i in 1:length(sample))
{
  m=as.data.frame(fread(paste0('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/results_force_call/',ctdna_tf[ctdna_tf$patient==sample_patient[i],'Sample'],'/mutations_unfiltered.hg38_multianno.txt'),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  m$mutect_force_call_depth=sapply(strsplit(m$Otherinfo13, ':'), function(x) ifelse(length(x) >= 4, x[4], NA))
  m$alt_depth=sapply(strsplit(m$Otherinfo13, ':'), function(x) ifelse(length(x) >= 4, x[2], NA))
  m$alt_depth=as.numeric(sub(".*,\\s*", "", m$alt_depth))
  m$varID=paste(m$Chr,m$Start,m$Ref,m$Alt,sep='_')
  private=tumor[tumor$sample_name==sample[i],]
  track[i,2]=nrow(private)
  ctDNA1=ctDNA[ctDNA$patient==private$patient[1],]
  private_overlap=private[private$varID %in% ctDNA1$varID,]
  
  track[i,3]=nrow(private_overlap)
  fraction=ctdna_tf[ctdna_tf$patient==sample_patient[i],'Purity']/tumor_tf[tumor_tf$Sample==sample[i],2]
  #private=private[private$vaf*fraction[1,1]>=0.1,]
  comb=left_join(private,m,by='varID')
  comb$mutect_force_call_depth=as.numeric(comb$mutect_force_call_depth)
  comb_unadj=comb[comb$alt_depth>0 & comb$Alt.x==comb$Alt.y,]
  comb1=comb[comb$mutect_force_call_depth>3,]
  comb1=comb1[!is.na(comb1$mutect_force_call_depth),]
  comb1=comb1[(1-pbinom(2,comb1$mutect_force_call_depth,comb1$vaf*fraction))>=0.8,]
  #comb3=comb1[(1-binom.test(x=2,n=comb1$depth.y,p=comb1$vaf*fraction[1,1],alternative = 'greater'))>=0.8,]
  track[i,4]=nrow(comb1)
  comb2=comb1[comb1$alt_depth>0 & comb1$Alt.x==comb1$Alt.y,]
  track[i,5]=nrow(comb2)
  comb3=comb1[comb1$varID %in% ctDNA1$varID,]
  track[i,6]=nrow(comb3)
  track[i,7]=nrow(comb_unadj)
 pt=tumor[(tumor$patient==sample_patient[i]) & (tumor$sample_name!=sample[i]),]
 private_snv=private[private$varID %notin% pt$varID,]
 private_snv1=comb1[comb1$varID %in% private_snv$varID,]
 track[i,8]=nrow(private_snv1)
 comb4=comb2[comb2$varID %in% private_snv1$varID,]
 track[i,9]=nrow(comb4)
 
}
track[,10]=track[,3]/track[,2]
track[,11]=track[,5]/track[,4]
track[,12]=track[,6]/track[,4]
track[,13]=track[,7]/track[,2]
track[,14]=track[,9]/track[,8]

colnames(track)=c('sample','TumorSNVs','UnadjustedOverlap','ObservableSNVs','Overlap_mutect_force_call_alt_depth','Overlap_ctDNA_calling','Overlap_mutect_force_call_unadj','PrivateSNVs','PrivateSNVOverlap','Prop_Unadjusted','Prop_mutect_force_Call','Prop_ctDNA','Prop_mutect_unadj','Prop_private')
track$organ=str_extract(track$sample, "(?<=_)[^_]*(?=_)")
track$patient=str_extract(track$sample, "[^_]+")
track1=left_join(track,ctdna_tf,by='patient')
tumor_tf$sample=tumor_tf$Sample
track1=left_join(track1,tumor_tf,by='sample')
track1$difference=track1$TumorSNVs-track1$ObservableSNVs

ggplot(track1, aes(x=Purity.x, y=Prop_Unadjusted)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE) +
  theme_ipsum()

model=lm(Purity.x~Prop_private, data=track1)
summary(model)
model=lm(Purity.x~Prop_mutect_force_Call, data=track1)
summary(model)
model=lm(Purity.x~Prop_Unadjusted, data=track1)
summary(model)
model=lm(Purity.y~Prop_Unadjusted, data=track1)
summary(model)
model=lm(Purity.y~Prop_private, data=track1)
summary(model)

write.table(track,file='/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/ctDNA_Tumor_Overlap.txt',row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
track=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational/PowerAnalysis_new/ctDNA_Tumor_Overlap.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
tf=as.data.frame(fread('/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/Sample_Ploidy_Purity_ctDNA.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
colnames(tf)[1]='sample'
track2=left_join(track, tf, by='sample')
ggplot(data=track2, aes(x=Purity, y=Prop_mutect_force_Call)) +
  geom_point() +
  theme_ipsum()

boxplot(Prop_private ~ organ, data = track)
boxplot(Prop_mutect_force_Call ~ organ, data = track)

phenotype=as.data.table(fread('/fh/fast/ha_g/projects/ProstateTAN/IHC/phenotype_new.txt',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
phenotype$sample=phenotype$Sample
track2=left_join(track1, phenotype, by='sample')
track2=track2[!is.na(track2$phenotype_IHC_included),]
ggplot(data=track2, aes(x=phenotype_IHC_included, y=Prop_Unadjusted, fill=phenotype_IHC_included)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Representation of ctDNA based on phenotype") +
  xlab("Phenotype")
model <- lmer(Prop_Unadjusted ~ phenotype_IHC_included + (1|patient), data = track2[track2$phenotype_IHC_included %in% c('ARneg_NEpos','ARpos_NEneg'),])
coef(summary(model))[, "Pr(>|t|)"][2]

track2=track2[!is.na(track2$CCP.31),]
ggplot(track2, aes(x=CCP.31, y=Prop_private)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE) +
  theme_ipsum()

a=lm(CCP.31~Prop_mutect_force_Call,data=track2)
summary(a)
df_long <- track[,c(1,10,11)] %>%
  pivot_longer(cols = c('Prop_Unadjusted','Prop_mutect_force_Call'), names_to = "variable", values_to = "value")
ggplot(df_long, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Value", fill = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(c(0,1))

ggplot(data=track, aes(x=organ, y=Prop_private, fill=organ)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Representation of ctDNA based on met site organ") +
  xlab("Organ")

ggplot(data=track, aes(x=organ, y=Prop_mutect_force_Call, fill=organ)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Representation of ctDNA based on met site organ") +
  xlab("Organ")

track_pt=track[track$patient=='18-115',]
ggplot(track_pt, aes(x=sample, y=Prop_mutect_force_Call,fill=sample)) + 
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set1") +
  theme_ipsum()

barplot(track$Prop_mutect_force_Call[1:4])

track2=track[track$organ %in% c('LN','LIVR'),]

model <- lmer(Prop_private ~ organ + (1|patient), data = track2)
coef(summary(model))[, "Pr(>|t|)"][2]
