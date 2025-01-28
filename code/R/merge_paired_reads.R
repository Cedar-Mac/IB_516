# Load packages
packages <- c('JAMP', 'tidyverse', 'vegan', 'here')
lapply(packages, library, character.only = T)
source("code/R/JAMP_folder_rename.R")


merge_reads <- function(){
# Make sure that all reads are paired, sudo check that no files were left behind
first_reads <- list.files(
  path = here("data", "test_data"), 
  pattern = ".*R1.fastq", 
  full.names = TRUE)
second_reads <- list.files(
  path = here("data", "test_data"), 
  pattern = ".*R2.fastq", 
  full.names = TRUE)
unmerged_file_names <- data.frame(first_reads, second_reads) 

# Remove the _R1.fastq suffix and _R2.fastq suffix and verify names match
unmerged_file_names$first_reads <- substring(
  text = unmerged_file_names$first_reads, 
  first = 1, 
  last = nchar(unmerged_file_names$first_reads)-9)

unmerged_file_names$second_reads <- substring(
  text = unmerged_file_names$second_reads, 
  first = 1, 
  last = nchar(unmerged_file_names$second_reads)-9)

unmerged_file_names$ERROR <- ifelse(
  test = unmerged_file_names$first_reads != unmerged_file_names$second_reads, 
  yes = "ERROR", 
  no = "fine")

# Use lapply to iterate over unmerged_file_names df.
unmerged_file_names %>%
  lapply(Merge_PE(file1 = first_reads, file2 = second_reads, exe = "vsearch"))

# Make sure all reads merged successfully (there should be half the number of resulting files).
if (length(list.files(here("A_merge_PE", "_data"))) == 
    length(list.files(here("input_files", "_data"))) / 2) {
      print("Success!") 
  }

}


merge_reads()
