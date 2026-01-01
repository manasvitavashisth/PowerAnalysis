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
