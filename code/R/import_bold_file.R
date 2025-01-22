import_bold_file <- function(file=""){
  bold_file <- readxl::read_xlsx(file, 
                                 col_names = TRUE) 
  
  # Only take the best match (first row) for each OTU.
  best_matches_only <- bold_file %>%
    filter_at(vars(contains("You ")), any_vars(!is.na(.))) 
  
  # Remove caret from OTU names in taxa list.
  best_matches_only$`You searched for` <- substr(best_matches_only$`You searched for`, 2, 9)
  
  # Rename OTU column
  colnames(best_matches_only)[1] <- "OTU"
  
  return(best_matches_only)
}

