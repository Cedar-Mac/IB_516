# Import primer sequences (just COI primer sequence, not including universal tail)
primers_file <- "input_files/fwhF2_EPTDr2n_primer_sequences.csv"
primers <- read.csv(primers_file)

# Identify name of forward and reverse primer
f.primer <- "fwhF2_"
r.primer <- "EPTDr2n_"

# Return TRUE or FALSE if each primer is Forward or Reverse
is_f <- grepl(f.primer, primers$Name)
is_r <- grepl(r.primer, primers$Name)

# Create character string with forward and reverse primers
fwhF2 <- primers$Sequence[primers$Name=="fwhF2"]
EPTDr2n <- primers$Sequence[primers$Name=="EPTDr2n"]

primer_sequences <- cbind(fwhF2, EPTDr2n)

# Trim primers 
merged_files <- list.files(here("A_merge_PE", "_data"), full.names = TRUE)

Cutadapt(files = merged_files, 
         forward = fwhF2,
         reverse = EPTDr2n,
         bothsides=T)

# By using "bothsides=T", forward or reverse primers are detected on both ends. This is not nessesary for fusion primers.