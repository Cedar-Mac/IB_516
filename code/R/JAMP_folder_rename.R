# Function to rename JAMP folders

JAMP_folder_rename <- function(newname=""){
  dirinfo <- file.info(list.dirs(recursive = FALSE)) %>% filter(str_detect(row.names(.), "./.git", negate=T)) #get directories in WD that is not .git/
  recent_dir <- row.names(dirinfo)[which.max(dirinfo$ctime)]
  recent_newname <- paste(newname, Sys.Date(), sep="_")
  file.rename(recent_dir, recent_newname) #rename file with a more descriptive name
  assign("recent_newname", recent_newname, envir = globalenv())
}

