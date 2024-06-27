#Write a BED file with union of all the SNVs in a patient
for(i in 1:length(patient))
{
  snv_bed=tumor[tumor$patient==patient[i],]
  snv_bed=snv_bed[!duplicated(snv_bed$varID),]
  snv_bed=snv_bed[,c('Chr', 'Start', 'End')]
  snv_bed$Start=snv_bed$Start-1
  write.table(snv_bed,file=paste0('Path to File') ,row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

#Write a .sh file quering gatk CollectAllelicCounts for all the patient SNV in the ctDNA file
ctDNA_tf=as.data.frame(fread('Path to BAM files for all samples',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
a=paste0('gatk CollectAllelicCounts -I ctDNA_BAM_file -R hg38_reference -L Bed_file_with_private_SNVs_from_the_patient -O Output_file_path')
file_conn <- file('gatk.sh', open = "w")
# Write the strings to separate lines in the .sh file
writeLines(a, file_conn)
# Close the file connection
close(file_conn)
#Run this file using sh Filename.sh on the terminal

# Function to read output file of CollectAllelicCounts and ignore rows starting with '@', then read the rest into a data frame
readFileIgnoreAtRows <- function(filePath) {
  # Step 1: Read the entire file into a character vector
  fileContent <- readLines(filePath)
  
  # Initialize an empty vector to hold the non-@ lines
  nonAtLines <- c()
  
  # Step 2: Filter out lines starting with '@'
  for(line in fileContent) {
    if (!grepl("^@", line)) { # Using grepl to check if the line does NOT start with '@'
      nonAtLines <- c(nonAtLines, line)
    }
  }
  # Convert the list of lines into a single string, ensuring tabs are preserved
  dataString <- paste(nonAtLines, collapse = "\n")
  
  # Step 3: Convert the string into a data frame, handling tab separation
  # Splitting the string into lines based on newline characters
  dataLines <- strsplit(dataString, "\n")[[1]]
  
  # Converting the list of trimmed lines into a matrix, using tab as separator
  dataMatrix <- do.call(rbind, lapply(dataLines, function(line) strsplit(line, "\\t")[[1]]))
  
  # Converting the matrix to a data frame
  dataFrame <- as.data.frame(dataMatrix, stringsAsFactors = FALSE)
  colnames(dataFrame)=dataFrame[1,]
  dataFrame=dataFrame[-1,]
  return(dataFrame)
}

#Example Plotting functions, please feel free to build your own depending on your data framework

unfiltered <- readFileIgnoreAtRows("ctDNAFile")
private=as.data.frame(fread('Private Mutations File',header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))

track=as.data.frame(matrix(data=NA,nrow=length(sample),ncol=1))
colnames(track)='Sample'
track$Sample=sample

for(i in 1:length(sample))
{
  unfiltered <- readFileIgnoreAtRows(paste0("ctDNA file"))
  private=as.data.frame(fread(paste0('private mtuation file'),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  track[i,2]=nrow(private)
  private_overlap=ctDNA[ctDNA$varID %in% private$varID & ctDNA$patient==private$patient[1],]
  track[i,3]=nrow(private_overlap)
  fraction=ctdna_tf[ctdna_tf$patient==sample_patient[i],2]/tumor_tf[tumor_tf$Sample==sample[i],3]
  private=private[private$vaf*fraction[1,1]>=0.1,]
  colnames(unfiltered)=c('Chr','Start','depth','alt_depth','Ref','Alt')
  unfiltered$Start=as.integer(unfiltered$Start)
  comb=inner_join(private,unfiltered,by=c('Chr'='Chr','Start'='Start'))
  comb$depth.y=as.integer(comb$depth.y)
  comb1=comb[comb$depth.y>3,]
  comb1=comb1[!is.na(comb1$depth.y),]
  comb1=comb1[(1-pbinom(2,comb1$depth.y,comb1$vaf*fraction[1,1]))>=0.8,]
  track[i,4]=nrow(comb1)
  comb2=comb1[comb1$alt_depth>0 & comb1$Alt.x==comb1$Alt.y,]
  track[i,5]=nrow(comb2)
}
track[,6]=track[,3]/track[,2]
track[,7]=track[,5]/track[,4]
colnames(track)=c('sample','PrivateMutations','UnadjustedOverlap','ObservableMutations','AdjustedOverlap','Unadjusted%','Adjusted%')

df_long <- track[,c(1,6,7)] %>% 
  pivot_longer(cols = c('Unadjusted%', 'Adjusted%'), names_to = "Type", values_to = "Value")

# Create the scatter plot
ggplot(df_long, aes(x = sample, y = Value, color = Type)) +
  geom_point(size = 3) +  # Adjusts the size of the points
  labs(title = "Scatter Plot of Unadjusted and Adjusted Percentages by Sample",
       x = "Sample",
       y = "Percentage",
       color = "Type") +
  theme_minimal() +  # Optional: Applies a minimal theme
  scale_color_manual(values = c("Unadjusted%" = "blue", "Adjusted%" = "red"))  # Customizes colors
