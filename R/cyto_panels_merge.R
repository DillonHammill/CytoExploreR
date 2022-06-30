## CYTO_PANELS_MERGE -----------------------------------------------------------

#' Combine multiple overlapping antibody panels
#'
#' @param x an object of class \code{\link[flowWorkspace:cytoset]{cytoset}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population to extract from each sample prior
#'   to combining the antibody panels, set to the "root" node to use all events
#'   in each sample by default.
#' @param select vector or list of sample selection criteria passed to
#'   \code{cyto_select()} to extract the samples that should have their panels
#'   merged.
#' @param channels names of the channels or markers over which the panel
#'   integration is to be performed, set to all channels by default.
#' @param merge_by a vector of experiment variables to dictate how samples
#'   should be split into groups prior to merging their panels.
#' @param events indicates the proportion \code{events < 1} or the number of
#'   events \code{events > 1} to keep in the merged samples per group, set to
#'   the average number of events within each group by default. \code{events}
#'   cannot exceed the maximum number of events within each group. Setting
#'   \code{events = 1} will keep all events within each merged group.
#' @param save_as name of a CSV to which the edited panel details should be
#'   saved, defaults to \code{Panel-Details.csv} prefixed with the date unless
#'   \code{file} is specified or a file is found when searching the current
#'   directory for panel details. Setting this argument to NA will prevent
#'   writing of edited panel details to file. Custom file names passed to
#'   \code{save_as} should end in \code{"Panel-Details.csv"} so that
#'   CytoExploreR can automatically find and read in these details when
#'   required.
#' @param file name of csv file containing the panel information associated with
#'   each sample in the supplied data. This file must contain a column called
#'   \code{"name"} which contains the names of the samples.
#' @param pool logical to indicate whether the integrated data should be
#'   combined into a single cytoset if the groups share the same channels, set
#'   to TRUE by default.
#' @param k number of neighbourhoods to consider when performing kNN regression,
#'   set to 10 by default.
#' @param ... additional arguments passed to \code{FNN::knn.reg()} to offer
#'   control over the algorithm to use when constructing the kNN trees.
#'   
#' @return a list of cytosets containing the merged samples for each group in
#'   \code{merge_by}.
#'
#' @importFrom BiocGenerics nrow
#'
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_merge_by}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' \dontrun{
#' # read in data - add patient_id variable to details 
#' gs <- cyto_setup(
#'   "FCS_Files"
#' )
#' 
#' # merge panels by patient
#' cs_list <- cyto_panels_merge(
#'   gs,
#'   parent = "root",
#'   channels = c(
#'   "FSC-A", 
#'   "SSC-A",
#'   "PE-A",
#'   "PE-Cy7-A",
#'   "FITC-A",
#'   "APC-A",
#'   "APC-Cy7-A"
#'   ),
#'   merge_by = "patient_id",
#'   events = 50000,
#'   pool = TRUE
#' )
#' 
#' # add merged fdata to new GatingSet
#' gs <- GatingSet(cs_list[[1]])
#' }
#'
#' @export
cyto_panels_merge <- function(x,
                              parent = "root",
                              select = NULL,
                              channels = NULL,
                              merge_by = NULL,
                              events = NULL,
                              save_as = NULL,
                              file = NULL,
                              pool = TRUE,
                              k = 10,
                              ...) {
  
  # TODO: TIDY UP - TRANSFER TRANSFORMERS COMPENSATION GATES
  # TODO: DEFAULT VALUE FOR K COMPUTED BY EVENTS
  # TODO: ADD INFORMATIVE MESSAGES
  
  # LOAD REQUIRED PACKAGES
  cyto_require(
    "FNN",
    source = "CRAN"
  )
  
  # MERGE_BY
  if(is.null(merge_by)) {
    stop(
      "Supplying names of the grouping variables to merge_by!"
    )
  }
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(
      x,
      select
    )
  }
  
  # DEFAULT CHANNELS
  if(is.null(channels)) {
    channels <- cyto_channels(
      x[[1]],
      exclude = c(
        "Time",
        "Event-?ID",
        "Sample-?ID",
        "Original"
      )
    )
  # SUPPLIED CHANNELS
  } else {
    channels <- cyto_channels_extract(
      x[[1]],
      channels = channels
    )
  }
  
  # MARKERS
  markers <- cyto_markers(x)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # TODO: LOOSEN MATCHING TO DETAILS.CSV?
  # LOCATE PANEL DETAILS
  pd_new <- cyto_file_search(
    "Panel-Details.*\\.csv$",
    colnames = "name",
    rownames = rownames(pd),
    ignore.case = TRUE,
    data.table = FALSE,
    type = "panel details",
    files = file
  )
  
  # INHERIT PANEL DETAILS
  if(length(pd_new) > 0) {
    pd <- pd_new[[1]]
    file <- names(pd_new)
  }
  
  # MARKES -> CHANNELS
  if(any(markers %in% colnames(pd))) {
    colnames(pd)[colnames(pd) %in% markers] <- cyto_channels_extract(
      x,
      channels = colnames(pd)[colnames(pd) %in% markers]
    )
  }
  
  # APPEND MISSING CHANNELS
  m <- channels[!channels %in% colnames(pd)]
  if(length(m) > 0) {
    pd <- cbind(
      pd,
      matrix(
        NA,
        ncol = length(m),
        nrow = nrow(pd),
        dimnames = list(
          rownames(pd),
          m
        )
      )
    )
  }
  
  # ASSIGN BACKBONE MARKERS
  if(length(markers) > 0) {
    lapply(
      channels,
      function(z) {
        if(.all_na(pd[,z])) {
          # FSC|SSC
          if(.grepl("FSC|SSC", z, ignore.case = TRUE)) {
            pd[, z] <<- z
          } else if(z %in% names(markers)) {
            pd[, z] <<- markers[match(z, names(markers))]
          }
        }
      }
    )
  }
  
  # SAVE_AS
  if(is.null(save_as)) {
    save_as <- file
    if(is.null(save_as)) {
      save_as <- cyto_file_name(
        paste0(
          format(
            Sys.Date(), 
            "%d%m%y"
          ), 
          "-Panel-Details.csv"
        )
      )
    }
  }
  
  # EDIT PANEL DETAILS
  if(interactive() & cyto_option("CytoExploreR_interactive")) {
    # ROWNAMES
    cyto_names <- rownames(pd)
    rownames(pd) <- NULL
    # EDIT EXPERIMENT DETAILS
    pd <- data_edit(
      pd,
      logo = CytoExploreR_logo(),
      title = "Panel Details Editor",
      row_edit = FALSE,
      # col_readonly = "name",
      quiet = TRUE,
      hide = TRUE,
      viewer = "pane"
    )
    # REPLACE ROW NAMES - REQUIRED TO UPDATE CYTO_DETAILS
    rownames(pd) <- cyto_names
  }
  
  # SET ROWNAMES
  if(is.null(rownames(pd))) {
    rownames(pd) <- pd[, "name"]
  }
  
  # UPDATE DETAILS - PD MAY CONTAIN EXTRA INFORMATION
  cyto_details(x) <- pd[match_ind(
    rownames(
      cyto_details(x)
    ), 
    rownames(pd)
  ), , drop = FALSE]
  
  # SAVE UPDATED DETAILS - CANNOT SAVE ABOVE AS ROWNAMES REMOVED
  if(!.all_na(save_as)) {
    write_to_csv(
      pd,
      file = save_as,
      row.names = TRUE
    )
  }
  
  # UPDATE DETAILS - PD MAY CONTAIN EXTRA INFORMATION
  cyto_details(x) <- pd[match_ind(
    rownames(
      cyto_details(x)
    ), 
    rownames(pd)
  ), , drop = FALSE]
  
  # SPLIT DATA INTO GROUPS
  x <- cyto_group_by(
    x,
    group_by = merge_by
  )
  
  # COMPUTE MEAN EVENTS PER GROUP
  if(is.null(events)) {
    events <- LAPPLY(
      x,
      function(z) {
        ceiling(
          mean(
            unlist(
              BiocGenerics::nrow(z)
            ),
            na.rm = TRUE
          )
        )
      }
    )
  # EVENTS PER GROUP
  } else {
    events <- rep(
      events, 
      length.out = length(x)
    )
  }
  
  # TODO: ADD PROGRESS BARS
  # PANEL INTEGRATION
  x <- structure(
    lapply(
      seq_along(x),
      function(z) {
        # TODO: ADD MESSAGE RE EVENTS
        # GROUP CYTOSET | GATINGSET
        y <- x[[z]]
        # GROUP PANEL DETAILS
        pd_group <- cyto_details(y)
        # BACKBONE CHANNELS
        backbone_channels <- LAPPLY(
          channels,
          function(v) {
            # EMPTY
            if(.all_na(pd_group[, v, drop = TRUE]) | 
               all(pd_group[, v, drop = TRUE] %in% "NA")) {
              pd_group <<- pd_group[, -match(v, colnames(pd_group))]
              return(NULL)
            }
            # BACKBONE
            if(length(unique(pd_group[, v, drop = TRUE])) == 1) {
              return(v)
            }
            return(NULL)
          }
        )
        # CHECK BACKBONE CHANNELS
        if(length(backbone_channels) == 0) {
          stop(
            paste0(
              names(x)[z], 
              " does not contain any overlapping backbone markers!"
            )
          )
        # KNN REQUIRE > 2D
        } else if(length(backbone_channels) < 3) {
          warning(
            paste0(
              "Panel integration works best when there are at least 3 ",
              "overlapping backbone parameters."
            )
          )
        }
        # BACKBONE MARKERS
        backbone_markers <- as.character(
          pd_group[1, backbone_channels, drop = FALSE]
        )
        # EXTRACT BACKBONE DATA
        backbone_data <- cyto_data_extract(
          y,
          parent = parent,
          channels = backbone_channels,
          format = "cytoset",
          copy = FALSE,
          events = events[z],
          coerce = TRUE
        )[[1]]
        # QUANTILE NORMALISATION
        backbone_data <- cyto_quantile_norm(
          backbone_data,
          channels = backbone_channels
        )
        # COMPUTE ORIGINAL CHANNEL LIMITS
        backbone_limits <- apply(
          cyto_apply(
            backbone_data,
            FUN = "cyto_stat_range",
            channels = backbone_channels,
            input = "matrix",
            copy = FALSE
          ),
          2,
          "range"
        )
        # RANGE NORMALISATION
        backbone_data <- cyto_apply(
          backbone_data,
          FUN = function(cf) {
            cyto_exprs(cf) <- cyto_stat_rescale(
              cyto_exprs(
                cf,
                channels = backbone_channels,
                drop = FALSE
              ),
              scale = c(0, 1),
              limits = backbone_limits
            )
            return(cf)
          },
          channels = backbone_channels,
          input = "cytoframe",
          copy = FALSE,
          simplify = TRUE
        )
        # OPEN CHANNELS
        open_channels <- colnames(pd_group)[
          colnames(pd_group) %in% channels &
            !colnames(pd_group) %in% backbone_channels
        ]
        # IMPUTED MISSING CHANNELS PER SAMPLE
        imputed_data <- structure(
          lapply(
            seq_along(y),
            function(v) {
              # SAMPLE PANEL DETAILS
              pd_sample <- cyto_details(y[v])
              # LOCATE CHANNELS TO IMPUTE
              impute_channels <- open_channels[
                !as.character(
                    pd_sample[1, open_channels, drop = FALSE]
                  ) %in% c(NA, "NA", "")
              ]
              # BYPASS SAMPLES WITHOUT ADDITIONAL MARKERS
              if(length(impute_channels) == 0) {
                return(NULL)
              }
              # IMPUTE MARKERS
              impute_markers <- as.character(
                pd_sample[1, impute_channels, drop = FALSE]
              )
              # BUILD KNN MODEL FOR EACH MARKER
              do.call(
                "cbind",
                structure(
                  lapply(
                    impute_channels,
                    function(w) {
                      # EXTRACT SAMPLE TRAINING DATA
                      train_data <- cyto_data_extract(
                        y,
                        select = v,
                        format = "matrix",
                        channels = c(backbone_channels, w),
                        events = 1, # TODO: SEPARATE SAMPLING HERE TO TRAIN KNN?
                        copy = FALSE
                      )[[1]][[1]]
                      # BUILD KNN MODEL & IMPUTE MARKER
                      cyto_func_call(
                        "FNN::knn.reg",
                        args = list(
                          train = train_data[, -match(w, colnames(train_data))],
                          test = cyto_data_extract(
                            backbone_data,
                            channels = backbone_channels,
                            format = "matrix",
                            events = 1,
                            copy = FALSE
                          )[[1]][[1]],
                          y = train_data[, match(w, colnames(train_data))],
                          k = k,
                          ...
                        )
                      )$pred
                    }
                  ),
                  names = impute_markers
                )
              )
            }
          ),
          names = cyto_names(y)
        )
        # REMOVE IMPUTED DATA FROM BACKBONE SAMPLES
        imputed_data[unlist(lapply(imputed_data, "is.null"))] <- NULL
        # COMBINE IMPUTED DATA
        imputed_data <- do.call("cbind", imputed_data)
        # TODO: HANDLE DUPLICATED IMPUTED PARAMETERS
        # BACKBONE DATA -> ORIGINAL SCALE
        backbone_data <- cyto_apply(
          backbone_data,
          FUN = function(cf) {
            cyto_exprs(cf) <- cyto_stat_rescale(
              cyto_exprs(
                cf,
                channels = backbone_channels,
                drop = FALSE
              ),
              scale = backbone_limits,
              limits = c(0, 1)
            )
            return(cf)
          },
          channels = backbone_channels,
          input = "cytoframe",
          copy = FALSE,
          simplify = TRUE
        )
        # REPLACE BACKBONE CHANNEL NAMES WITH MARKERNAMES
        cyto_channels(backbone_data) <- backbone_markers
        # COMBINE BACKBONE & IMPUTED DATA
        cyto_cbind(
          backbone_data,
          imputed_data
        )
      }
    ),
    names = names(x)
  )

  # POOL SAMPLES - MATCHING PANELS
  if(pool) {
    if(all(LAPPLY(x, function(z){
      setequal(cyto_channels(x[[1]]),
               cyto_channels(z))
    }))) {
      x <- list(
        cytoset(
          structure(
            lapply(x, `[[`, 1)
          ),
          names = names(x)
        )
      )
    }
  }
  
  # RETURN LIST OF CYTOSETS - INTEGRATED DATA
  return(x)
  
}