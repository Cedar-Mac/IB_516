# denoise files filtered at EE = 0.5
ee0_5_files <- "maxee_haplotypes_maxee0_5_2023-08-25/_data"
files_for_denoise_maxee0_5 <- list.files(ee0_5_files, full.names = TRUE)

for (i in c("1", "3", "5", "7", "9", "11", "13", "15")) {
  Denoise(files = files_for_denoise_maxee0_5,
          # remove sequences with less than 10 occurrences within a site
          minsize = 10,
          # remove sequences with less than 0.01% relative abundance for each site
          minrelsize = 0.001,
          # set alpha parameter (see for loop)
          unoise_alpha = i)
  
  JAMP_folder_rename(
    newname = paste("denoised_haplotypes_maxee0_5", "_alpha_", i, sep=""))
}

# denoise files filtered at EE = 1.0
ee1_0_files <- "maxee_haplotypes_maxee1_0_2023-08-25/_data"
files_for_denoise_maxee1_0 <- list.files(ee1_0_files, full.names = TRUE)

for (i in c("1", "3", "5", "7", "9", "11", "13", "15")) {
  Denoise(files = files_for_denoise_maxee1_0,
          # remove sequences with less than 10 occurrences within a site
          minsize = 10,
          # remove sequences with less than 0.1% relative abundance for each site
          minrelsize = 0.001,
          # set alpha parameter (see for loop)
          unoise_alpha = i)
  
  JAMP_folder_rename(
    newname = paste("denoised_haplotypes_maxee1_0", "_alpha_", i, sep=""))
}