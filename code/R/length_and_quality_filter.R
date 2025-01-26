#Filter trimmed reads to only match 142 bp
Cutadapt_files <- list.files("B_Cutadapt/_data", full.names = TRUE)
Minmax(files = Cutadapt_files,
       min=(142),
       max=(142))
JAMP_folder_rename(newname = "minmax_haplotypes")

# filter by quality scores (ee)
min_max_path = paste(recent_newname, "/_data", sep = "")
haplotype_minmax_files <- list.files(min_max_path, full.names = TRUE)

Max_ee(files = haplotype_minmax_files,
       max_ee=0.5)
JAMP_folder_rename(newname = "maxee_haplotypes_maxee0_5")

Max_ee(files = haplotype_minmax_files,
       max_ee=1.0)
JAMP_folder_rename(newname = "maxee_haplotypes_maxee1_0")