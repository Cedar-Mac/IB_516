all_taxa_sequences <- read_csv("all_taxa_sequences.csv", 
                               col_names = TRUE,
                               show_col_types = FALSE)

assign_taxon_names <- function(file, include_streams=FALSE) {
  # read in a haplo_table file
  file_data <- read_csv(file,
                   col_names = TRUE, 
                   show_col_types = FALSE)
  
  # join the haplo_table with the taxon information. Use the sequence columns to merge to avoid ambiguous OTU numbers.
  # select important columns.
  if (include_streams == FALSE) {
    complete_data <- left_join(file_data, all_taxa_sequences, by = join_by(sequences == seq)) %>% 
      select(Phylum, Class, Order, Family, Genus, Species, haplotype, sequences) %>% 
      filter(if_all(c(Phylum, Class, Order), ~ !is.na(.)))
    # can choose to include site information (all the streams)
  } else {
    complete_data <- left_join(file_data, all_taxa_sequences, by = join_by(sequences == seq)) %>% 
      select(-c("sort", "haplotype", "OTU.x", "OTU.y")) %>% 
      filter(if_all(c(Phylum, Class, Order), ~ !is.na(.)))
  }

  
  return(complete_data)
}