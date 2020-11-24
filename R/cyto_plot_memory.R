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

#' Inherit saved cyto_plot arguments
#' @noRd
.cyto_plot_args_inherit <- function(args, 
                                    memory) {
  
  # MEMORY
  if(missing(memory)) {
    memory <- .cyto_plot_args_recall()
  }
  
  # NO MEMORY
  if(is.null(memory)) {
    return(args)
  # MEMORY
  } else {
    # MEMORY INDEX
    if(is.null(names(memory))) {
      memory_ind <- 1
    } else {
      memory_ind <- match(NA, names(memory))[1]
    }
    # UPDATE MEMORY NAMES - DONE
    names(memory)[memory_ind] <- "DONE" # sets other names to NA
    # RESET CYTO_PLOT_MEMORY INDICATOR - ALL DONE
    if(!any(is.na(names(memory)))) {
      names(memory) <- NULL
    }
    # LIST OF CYTO_PLOT ARGUMENT LISTS
    args[names(memory[[memory_ind]])] <- memory[[memory_ind]]
    # UPDATE MEMORY
    .cyto_plot_args_save(memory)
    # UPDATED ARGUMENTS
    return(args)
  }
}
