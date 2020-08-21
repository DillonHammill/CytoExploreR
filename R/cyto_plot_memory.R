## CYTO_PLOT_MEMORY ------------------------------------------------------------

#' Recall a list of saved cyto_plot arguments
#' @noRd
.cyto_plot_args_recall <- function(){
  
  temp_dir <- tempdir()
  temp_files <- list.files(temp_dir)
  file_match <- grepl("cyto_plot_memory", temp_files) & 
    grepl(".rds", temp_files)
  
  # FILE EXISTS
  if(any(file_match)){
    file <- temp_files[file_match]
    file <- file[1]
    return(readRDS(paste0(temp_dir,
                          .Platform$file.sep,
                          file)))
    # FILES DOES NOT EXIST
  }else{
    return(NULL)
  }
  
}

#' Save a list of cyto_plot arguments
#' @noRd
.cyto_plot_args_save <- function(args){
  
  # TEMP DIR
  temp_dir <- tempdir()
  
  # REMOVE ANY EXISTING FILE
  .cyto_plot_args_remove()
  
  # SAVE DEPARSED ARGUMENTS 
  saveRDS(args,
          paste0(temp_dir,
                 .Platform$file.sep,
                 "cyto_plot_memory.rds"),
          compress =  "xz")
  
}

#' Delete cyto_plot arguments saved to RDS tempfile
#' @noRd
.cyto_plot_args_remove <- function(){
  
  temp_dir <- tempdir()
  temp_files <- list.files(temp_dir)
  file_match <- grepl("cyto_plot_memory", temp_files) & 
    grepl(".rds", temp_files)
  
  # FILE EXISTS
  if(any(file_match)){
    file.remove(paste0(temp_dir,
                       .Platform$file.sep,
                       temp_files[file_match]))
  }
  
}
