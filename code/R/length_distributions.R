dir.create("_Length_distribution", path = "B_Cutadapt/_stats/_Length_distribution")

trimmed_files <- list.files(here("B_Cutadapt", "_data"), full.names = TRUE) # Create list of trimmed filenames

start_time <- Sys.time() # start timer to see how long the following for loop will take

# Loop over all trimmed files. Create length distribution pdf with name of the sample (cutoff file path and suffix)
for (file in trimmed_files) {
  Length_distribution(
    sequFile = file, 
    out = paste("B_Cutadapt/_stats/_Length_distribution/", 
                substr(file, start = nchar("/Users/cedarmackaness/gc/B_Cutadapt/_data/ "), 
                       nchar(file) - nchar("_PE_cut.fastq")), ".pdf", sep=""), fastq=TRUE)
}

end_time <- Sys.time() # end timer once for loop is completed

end_time - start_time