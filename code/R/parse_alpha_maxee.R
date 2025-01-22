# A function that will add two columns to a dataframe.
# The first column is the maxee filtering parameter,
# the second column is the alpha value.
parse_alpha_maxee <- function(data=data.frame()){
  if ("method" %in% colnames(data)){
    # create new column with maxee values
    data$maxee <- substring(text = data$method, 
                            first = nchar("maxee") + 1, 
                            last = nchar("maxeeX_X"))
    data$maxee <- as.factor(data$maxee)
    
    # create new column with alpha values
    data$alpha <- substring(text = data$method, 
                            first = nchar("maxeeX_X_alpha_") + 1, 
                            last = nchar("maxeeX_X_alpha_XX")) %>%
      parse_number()
    return(data)
    
  } else{ print("No method column! :(") }

}

