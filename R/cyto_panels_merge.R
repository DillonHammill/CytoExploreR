## CYTO_PANELS_MERGE -----------------------------------------------------------

#' Use machine learning to unite multiple panels on the same samples
#'
#' @param x object of class \code{\link[flwoWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} or a list of these
#'   objects each consisting of a different panel.
#' @param select named list containing experimental variables to be used to
#'   select samples using \code{\link{cyto_select}}.
#' @param exclude vector of channels or markers that should be excluded during
#'   modelling, \code{Time}, \code{Event-ID} and \code{Sample-ID} are
#'   automatically excluded.
#' @param merge_by the name of an experiment variable by which panels should be
#'   united (e.g. "ID").
#' @param ... additional arguments passed to \code{xgboost::xgboost}.
#'
#' @return a new \code{cytoset} or \code{Gatingset} containing data combined
#'   from multiple panels.
#'
#' @importFrom purrr transpose
#'
#' @seealso \code{\link{cyto_select}}
#' @seealso \code{\link{cyto_infinity_screen}}
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
cyto_panels_merge <- function(x,
                              parent = "root",
                              select = NULL,
                              exclude = NULL,
                              merge_by = NULL,
                              trans = NA,
                              inverse = FALSE,
                              save_as = NULL,
                              ...) {
  
  # LOAD XGBOOST ---------------------------------------------------------------
  
  # REQUIRE XGBOOST PACKAGE
  cyto_require(
    "xgboost"
  )
  
  # LIST OF LENGTH SAMPLES - LIST OF LENGTH PANELS
  # [[1]] SAMPLE 1
  #   [[1]] Sample_1_Panel_A.fcs - CYTOSET
  #   [[2]] Sample_1_Panel_B.fcs - CYTOSET
  #   [[3]] Sample_1_Panel_C.fcs - CYTOSET
  
  # DATA IS NOT EXTRACTED HERE...
  
  # FORMAT DATA ----------------------------------------------------------------
  
  # DATA - LIST OF CYTOSETS OR GATINGSETS PER PANEL
  if(cyto_class(x, "list", TRUE)) {
    # LIST OF CYTOSETS PER SAMPLE
    cs_list <- structure(
      lapply(
        x,
        function(z) {
          # TRANSFORMERS
          if(.all_na(trans)) {
            trans <<- cyto_transformers_extract(z)
          }
          # GROUP
          cyto_group_by(
            z,
            select = select,
            group_by = merge_by
          )
        }
      ),
      names = names(x)
    )
    # TRANSPOSE
    cs_list <- transpose(cs_list)
  # DATA - CYTOSET OR GATINGSET
  } else {
    # GROUPS
    cs_list <- cyto_group_by(
      x,
      select = select,
      group_by = merge_by
    )
    # SPLIT GROUPS
    cs_list <- structure(
      lapply(
        cs_list,
        function(z) {
          # TRANSFORMERS
          if(.all_na(trans)) {
            trans <<- cyto_transformers_extract(z)
          }
          # SPLIT GROUP BY NAME
          cyto_group_by(
            z,
            group_by = "name"
          )
        }
      ),
      names = names(cs_list)
    )
  }
  
  # CHANNELS TO EXCLUDE
  exclude = c(
    "Time",
    "Event-ID",
    "Sample-ID",
    exclude
  )
  
  # ANTIBODY PANELS
  panels <- structure(
    lapply(
      cs_list[[1]],
      cyto_markers,
      exclude = exclude,
      append = TRUE
    )
  )
  
  # NORMALISATION & PREDICTIONS ------------------------------------------------
  
  # LOOP THROUGH SAMPLES
  res <- structure(
    lapply(
      cs_list,
      function(z) {
        # Z IS A NAMED LIST OF (SINGLE) CYTOSETS
        
      }
    )
  )
  
  # WRITE NEW MERGED FCS FILES -------------------------------------------------
  
  # SAVE NEW FCS FILES
  if(!is.null(save_as)) {
    cyto_save(
      cs,
      trans = trans,
      inverse = inverse,
      save_as = save_as
    )
  }
  
}