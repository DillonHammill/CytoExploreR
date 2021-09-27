## CYTO_INFINITY_SCREEN --------------------------------------------------------

#' Use infinityFlow to integrate multiplexed parallel cytometry datasets
#'
#' \code{cyto_infinity_screen()} is a CytoExploreR wrapper around the
#' \code{infinity_flow()} function in the \code{infinityFlow} package that uses
#' machine learning to integrate multiplexed parallel cytometry datasets. Users
#' should note that \code{cyto_infinity_screen()} takes on the path to the FCS
#' files rather than a \code{cytoset} or \code{GatingSet}. Furthermore,
#' \code{cyto_infinity_screen()} does not return a \code{cytoset} or
#' \code{GatingSet}, instead the new FCS files written to \code{save_as} should
#' be read in afterwards using \code{cyto_setup} as per usual.
#'
#' @param path points to the location where the FCS files (one per well) are
#'   stored, replaces the \code{path_to_fcs} argument in \code{infinity_flow()}.
#' @param backbone_selection_file name of a CSV file containing columns
#'   \code{"name"}, \code{"desc"} and \code{"type"}, where type indicates
#'   whether the named parameter is a part of the \code{"backbone"}, or instead
#'   an \code{"exploratory"} parameter. Parameters can be optionally excluded by
#'   setting \code{type} to \code{"discard"}. Users need not create this file by
#'   hand as CytoExploreR will allow the user to interactively edit this
#'   information when required.
#' @param annotation Named character vector. Elements should be the targets of
#'   the exploratory antibodies, names should be the name of the FCS file where
#'   that exploratory antibody was measured. Users need not supply this
#'   information manually as CytoExploreR will allow users to interactively
#'   enter this information.
#' @param isotype Named character vector. Elements should be the isotype used in
#'   each of the wells (e.g. IgG2) and the corresponding isotype should be
#'   present in \emph{annotation} (e.g. Isotype_IgG2, with this capitalization
#'   exactly). Autofluorescence measurements should be listed here as "Blank".
#'   Users need not supply this information manually as CytoExploreR will allow
#'   users to interactively enter this information.
#' @param save_as name of the directory to which the output FCS files should be
#'   written, replaces the \code{path_to_output} argument in
#'   \code{infinity_flow()}.
#' @param ... additional arguments passed to \code{infinity_flow()}, refer to
#'   the \code{infinityFlow} package documentation for more details.
#'
#' @return raw and background-corrected imputed expression data for every
#'   Infinity antibody, new FCS files written to \code{save_as} along with plots
#'   of UMAP projections of the backbone antibody panel.
#'
#' @importFrom DataEditR data_edit
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @references Becht E, Headley M, Newell E, Gottardo R (2020). Infinity Flow:
#'   comprehensive single-cell protein profiling via massively parallel flow
#'   cytometry and machine learning. J Immunol 204:159.2
#'
#' @seealso \code{\link{cyto_setup}}
#'
#' @export
cyto_infinity_screen <- function(path = NULL,
                                 backbone_selection_file = NULL,
                                 annotation = NA,
                                 isotype = NA,
                                 save_as = NULL,
                                 ...) {
  
  # REQUIRE INFINITYFLOW
  cyto_require(
    "infinityFlow",
    source = "BioC",
    repo = "ebecht/infinityFlow",
    ref = paste0(
      "Becht E, Headley M, Newell E, Gottardo R (2020). Infinity Flow: ",
      "comprehensive single-cell protein profiling via massively parallel ",
      "flow cytometry and machine learning. J Immunol 204:159.2"
    )
  )
  
  # PATH 
  if(is.null(path)) {
    stop(
      "Supply the name of the directory containing the FCS files to 'path'!"
    )
  }
  
  # FCS FILENAMES
  fcs_files <- list.files(
    path, 
    pattern = ".fcs",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # READ IN A FILE - EXTRACT CHANNELS/MARKERS
  cs <- cyto_load(
    path,
    select = basename(fcs_files[1])
  )
  
  # CHANNELS
  channels <- cyto_channels(cs)
  
  # PREPARE INFINITY TARGETS & ISOTYPES
  if(any(is.na(isotype)) | any(is.na(annotation))) {
    # CHECK IF TARGET FILE EXISTS
    files <- list.files(".", pattern = ".csv")
    ind <- grep(
      "InfinityFlow-Targets.csv$",
      files,
      ignore.case = TRUE
    )
    # TARGET FILE EXISTS
    if(length(ind) != 0) {
      # READ FIRST TARGET FILE
      targets <- read_from_csv(
        files[ind[1]]
      )
    # INTERACTIVELY CREATE & EDIT TARGET FILE
    } else {
      # MESSAGE
      message(
        paste0(
          "Interactively edit the InfinityFlow antibody targets and isotypes",
          "..."
        )
      )
      # CREATE TARGETS
      targets <- data.frame(
        "name" = basename(fcs_files),
        "Infinity_target" = rep(
          annotation, 
          length(fcs_files)
        ),
        "Infinity_isotype" = rep(
          isotype, 
          length.out = length(fcs_files)
        ),
        stringsAsFactors = FALSE
      )
      # EDIT TARGETS
      targets <- data_edit(
        targets,
        logo = CytoExploreR_logo(),
        title = "InfinityFlow Targets & Isotypes Editor",
        col_names = FALSE,
        row.names = FALSE,
        quiet = TRUE,
        hide = TRUE,
        viewer = "pane"
      )
      # WRITE TARGETS TO FILE
      write_to_csv(
        targets,
        paste0(
          format(
            Sys.Date(),
            "%d%m%y"
          ),
          "-InfinityFlow-Targets.csv"
        )
      )
    }
    # UPDATE ISOTYPE & ANNOTATION
    isotype <- targets[, "Infinity_isotype"]
    names(isotype) <- targets[, "name"]
    annotation <- targets[, "Infinity_target"]
    names(annotation) <- targets[, "name"]
  }
  
  # BACKBONE SELECTION FILE - CREATE INTERACTIVELY
  if(is.null(backbone_selection_file)) {
    # CHECK BACKBONE FILE EXISTS
    files <- list.files(".", pattern = ".csv")
    ind <- grep(
      "InfinityFlow-Backbone.csv$",
      files,
      ignore.case = TRUE
    )
    # BACKBONE FILE EXISTS 
    if(length(ind) != 0) {
      backbone_selection_file <- files[ind[1]]
    # INTERACTIVELY CREATE & EDIT BACKBONE FILE
    } else {
      # MESSAGE
      message(
        paste0(
          "Interactively edit the InfinityFlow backbone selection file",
          "..."
        )
      )
      # CREATE BACKBONE SELECTION FILE
      backbone <- data.frame(
        "name" = channels,
        "desc" = cyto_markers_extract(cs, channels  = channels),
        "type" = "backbone",
        stringsAsFactors = FALSE
      )
      # INTERACTIVELY EDIT BACKBONE SELECTION
      backbone <- data_edit(
        backbone,
        logo = CytoExploreR_logo(),
        title = "InfinityFlow Backbone Selection Editor",
        col_names = FALSE,
        col_readonly = c("name", "desc"),
        quiet = TRUE,
        hide = TRUE,
        viewer = "pane",
        col_options = list("type" = c("backbone", "exploratory", "discard"))
      )
      # BACKBONE SELECTION TO FILE
      backbone_selection_file = paste0(
        format(
          Sys.Date(),
          "%d%m%y"
        ),
        "-InfinityFlow-Backbone.csv"
      )
      # WRITE NEW BACKBONE SELECTION FILE
      write_to_csv(
        backbone,
        backbone_selection_file
      )
    }
  }

  # SAVE_AS - APPEND INFINITYFLOW
  if(is.null(save_as)) {
    save_as <- paste0(
      basename(
        path
      ),
      "-InfinityFlow"
    )
  }
  
  # COMBINE ARGUMENTS
  args <- list(
    path_to_fcs = path,
    path_to_output = save_as,
    backbone_selection_file = backbone_selection_file,
    annotation = annotation,
    isotype = isotype,
    ...
  )
  
  # RUN INFINITYFLOW
  res <- do_call(
    "infinityFlow::infinity_flow",
    args
  )
  
  # RETURN DATA
  invisible(res)
  
}