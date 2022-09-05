## CYTO_SOM --------------------------------------------------------------------

#' Train a self-organising map (SOM) on cytometry data
#'
#' \code{cyto_som()} provides an intuitive way to map data stored in a
#' \code{cytoset}, \code{GatingHierarchy} or \code{GatingSet} into a compressed
#' self-organising map (SOM). The SOM is trained on the combined data from the
#' population specified by \code{parent} in the required \code{channels} using
#' the requested proportion or number of \code{events}. The codebook vectors
#' generated from the trained SOM are then passed to \code{cyto_map()} to obtain
#' dimension-reduced co-ordinates using the method specified by \code{map}. The
#' codes and dimension-reduced co-ordinates from the the trained SOM are stored
#' in a new cytoframe for each sample with the number of events assigned to each
#' node being stored in a new channel called \code{"SOM_counts"}. Data can also
#' be mapped to an existing SOM trained by \code{cyto_som()} if a SOM
#' \code{cytoset} or \code{GatingSet} is passed to \code{som}.
#'
#' @param x object of class \code{\link[flowWorkspace:cytoset]{cytoset}},
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param select passed to \code{cyto_select()} to control which samples are
#'   passed to the SOM, set to all samples by default.
#' @param parent name of the parent population upon which the SOM should be
#'   trained, set to the \code{"root"} node by default.
#' @param channels name(s) of the channel(s) or marker(s) to use when training
#'   the SOM, no default set so channels must be supplied manually.
#' @param grid vector of length 2 indicating the dimensions of the SOM grid, set
#'   to \code{c(14,14)} by default.
#' @param rlen number of time the data should be presented to the network, set
#'   to 10 by default.
#' @param alpha learning rate, default is to decline linearly from 0.05 to 0.01
#'   over \code{rlen} updates.
#' @param radius radius of the neighbourhood for SOM training, set to start with
#'   a value that covers two thirds of all unit-to-unit distances.
#' @param events indicates the number or proportion of events to use when
#'   training the SOM, set to 1 by default to use all available events.
#' @param map passed to the \code{type} argument of cyto_map() to perform
#'   dimensionality reduction on the SOM code vectors, set to \code{"t-SNE"} by
#'   default. Alternatively, a vector of channel names indicating the
#'   dimension-reduced channels in a supplied \code{som} or the supplied data
#'   that should be inherited by the new SOM. Refer to \code{\link{cyto_map}} or
#'   details.
#' @param trans object of class \code{transformerList} containing the
#'   transformers applied to the supplied data, only required for objects of
#'   class \code{cytoset}.
#' @param som a \code{cytoset}, \code{GatingHierarchy} or \code{GatingSet}
#'   containing a pre-computed SOM upon which the supplied data should be
#'   mapped. If the parameters defined by \code{map} are located, the existing
#'   dimension-reduced co-ordinates will be transferred to the new SOM. The data
#'   will be mapped to the existing SOM by computing the distance to the nearest
#'   node using the parameters defined by \code{channels}.
#' @param plot logical indicating whether the constructed SOM should be plotted
#'   in the channels produced by \code{map}, set to TRUE by default.
#' @param gatingTemplate name of a gatingTemplate csv file to set as the active
#'   \code{gatingTemplate} to store downstream gates on the constructed SOM, set
#'   to NULL by default to bypass \code{gatingTemplate} assignment.
#' @param ... additional arguments passed to \code{cyto_map()} for
#'   dimensionality reduction.
#'
#' @return a cytoset or GatingSet containing the code vectors for every sample
#'   with associated counts and dimension-reduced co-ordinates.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace cytoset GatingSet
#'
#' @seealso \code{\link{cyto_map}}
#'
#' @examples
#' \dontrun{
#' library(CytoExploreRData)
#'
#' # Activation GatingSet
#' gs <- GatingSet(Activation)
#' gs <- cyto_compensate(gs)
#' gs <- cyto_transform(gs)
#' gs <- cyto_gatingTemplate_apply(gs, Activation_gatingTemplate)
#'
#' # Build SOM
#' gs_som <- cyto_som(
#'   gs,
#'   select = 1:32,
#'   parent = "Live Cells",
#'   channels = c(
#'     "FSC-A",
#'     "SSC-A",
#'     cyto_markers(gs)
#'   ),
#'   grid = c(8,8)
#' )
#'
#' # Visualise SOMs
#' cyto_plot_som(
#'   gs,
#'   parent = "Live Cells",
#'   channels = c("t-SNE-1", "t-SNE-2")
#' )
#'
#' # Visualise SOM maps
#' cyto_plot_som_map(
#'   gs,
#'   select = c(1, 32),
#'   parent = "Live Cells",
#'   channels = c("t-SNE-1", "t-SNE-2"),
#'   point_col = c(
#'     "FSC-A",
#'     "SSC-A",
#'     cyto_markers(gs)
#'   )
#' )
#' }
#'
#' @export
cyto_som <- function(x,
                     select = NULL,
                     parent = "root",
                     channels = NULL,
                     grid = c(14,14),
                     rlen = 10,
                     alpha = c(0.05, 0.01),
                     events = 1,
                     map = "t-SNE",
                     trans = NA,
                     som = NULL,
                     plot = TRUE,
                     gatingTemplate = NULL,
                     ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # CYTOSET | GATINGSET ONLY
  if(!cyto_class(x, c("flowSet", "GatingSet"))) {
    stop(
      "cyto_som() only accepts cytoset, GatingHierachy or GatingSet objects!"
    )
  }
  
  # CHANNELS
  if(is.null(channels)) {
    stop(
      paste0(
        "Supply the name(s) of the channel(s) or marker(s) to 'channels' ",
        "to construct a SOM."
      )
    )
  }
  
  # SELECT
  if(!is.null(select)) {
    x <- cyto_select(
      x,
      select
    )
  }
  
  # TRANSFORMERS
  if(.all_na(trans)) {
    trans <- cyto_transformers_extract(x)
  }
  
  # CHANNELS MISSING
  if(is.null(channels)) {
    stop(
      paste0(
        "Supply the name(s) of the channel(s) or marker(s) to 'channels' ",
        "to train a SOM."
      )
    )
  # CONVERT CHANNELS
  } else {
    channels  <- cyto_channels_extract(
      x,
      channels = channels
    )
  }
  
  # MARKER ASSIGNMENTS
  markers <- cyto_markers(x)
  
  # EXPERIMENT DETAILS
  pd <- cyto_details(x)
  
  # GRID
  grid <- rep(grid, length.out = 2)
  
  # REQUIRE KOHONEN
  cyto_require(
    "kohonen",
    source = "CRAN",
    ref = paste0(
      "Wehrens, R., & Buydens, L. M. C. (2007). Self- and Super-organizing ",
      "Maps in R: The kohonen Package. Journal of Statistical Software, ",
      "21(5), 1–19."
    )
  )
  
  # SOM TRAINING ---------------------------------------------------------------
  
  # SOM SUPPLIED
  if(!is.null(som)) {
    # SOM CODES SUPPLIED
    if(cyto_class(som, c("flowSet", "GatingSet"))) {
      # SOM TRANSFORMERS
      som_trans <- cyto_transformers_extract(
        som
      )
      # SOM -> LINEAR SCALE
      codes <- cyto_data_extract(
        som,
        select = 1,
        parent = "root",
        channels = NULL,
        format = "matrix",
        trans = if(.all_na(som_trans)) {
          trans
        } else {
          som_trans
        },
        inverse = TRUE,
        copy = FALSE
      )[[1]][[1]]
    # SOM CODES SUPPLIED
    } else {
      stop(
        paste0(
          "'som' must be a cytoset or GatingSet containing the codes for a ",
          "SOM trained using cyto_som()!"
        )
      )
    }
    # VALID SOM
    if(all(channels %in% colnames(codes)) & nrow(codes) == prod(grid)) {
      # DIMENSION REDUCTION SUPPLIED
      if(length(map) > 1) {
        map <- cyto_channels_extract(
          som,
          channels = map
        )
        if(!all(map %in% colnames(codes))) {
          stop(
            paste0(
              "The names of the dimension-reduced paremeters supplied to 'map'",
              " don't match any parameters in 'som'!"
            )
          )
        }
        coords <- codes[, map]
        codes <- codes[, channels]
      # NO DIMENSIONALITY REDUCTION
      } else {
        codes <- codes[, channels]
        coords <- NULL
      }
      # RECONSTRUCT SOM FROM CODES
      som <- structure(
        list(
          codes = list(codes),
          data = list(codes),
          whatmap = 1,
          user.weights = 1,
          dist.fcts = "sumofsquares",
          maxNA.fraction = 0
        ),
        class = "kohonen"
      )
    # INVALID SOM
    } else {
      warning(
        paste0(
          "The supplied SOM does not match the required paramters or grid size",
          "! Training a new SOM..."
        )
      )
      som <- NULL
    }
  }
  
  # SOM REQUIRED
  if(is.null(som)) {
    # DATA TO TRAIN NEW SOM -> LINEAR SCALE
    x_train <- cyto_data_extract(
      x,
      parent = parent,
      channels = if(length(map) > 1) {
        c(channels, map)
      } else {
        channels
      },
      trans = trans,
      inverse = TRUE,
      format = "matrix",
      events = events,
      coerce = TRUE,
      barcode = FALSE,
      markers = FALSE,
      copy = TRUE
    )[[1]][[1]]
    # SOM GRID
    grid <- cyto_func_call(
      "kohonen::somgrid",
      list(
        xdim = grid[1],
        ydim = grid[2]
      )
    )
    # TRAIN SOM
    som <- cyto_func_call(
      "kohonen::som",
      list(
        x_train[, channels, drop = FALSE],
        grid = grid,
        rlen = rlen,
        alpha = alpha,
        keep.data = FALSE
      )
    )
    # MATCH ABOVE
    som$data <- list(som$codes[[1]])
    # COMPATIBILITY WITH SUPPLIED SOM
    som_trans <- trans
    # SOM CODE VECTORS
    codes <- som$codes[[1]]
    # DIMENSIONALITY REDUCTION
    if(length(map) > 1) {
      # MAP -> CHANNELS
      map <- cyto_channels_extract(
        x,
        channels = map
      )
      # MAP TO DATA EMBEDDING
      if(all(map %in% colnames(x_train))) {
        coords <- do.call(
          "cbind",
          structure(
            lapply(
              map,
              function(z) {
                cyto_impute(
                  train = x_train[, channels, drop = FALSE],
                  test = codes,
                  channels = channels,
                  labels = x_train[, z, drop = TRUE],
                  k = 5,
                  scale = NULL
                )
              }
            ),
            names = map
          )
        )
        colnames(coords) <- map
      } else {
        coords <- NULL
      }
    } else {
      coords <- NULL
    }
  }
  
  # CODES CYTOSET --------------------------------------------------------------
  
  # NUMBER OF CODES
  N <- nrow(codes)
  
  # CONVERT CODES TO CYTOSET ON LINEAR SCALE
  codes <- cytoset(
    structure(
      list(
        as(
          codes,
          "cytoframe"
        )
      ),
      names = "SOM"
    )
  )
  
  # DIMENSIONALITY REDUCTION ---------------------------------------------------
  
  # CYTO_MAP() ON CODES CYTOSET
  if(is.null(coords)) {
    codes <- cyto_map(
      codes,
      channels = channels,
      type = map,
      scale = FALSE,   # NO SCLAING REQUIRED
      events = 1,
      merge_by = "all",
      trans = trans,
      inverse = FALSE,   # CODES ON LINEAR SCALE
      plot = FALSE,
      ...
    )
    map <- cyto_channels(
      codes,
      exclude = c(
        channels,
        "Sample-ID",
        "Event-ID"
      )
    )
  # APPEND CO-ORDINATES TO CODES CYTOSET
  } else {
    codes <- cyto_cbind(
      codes,
      coords
    )
  }
  
  # PLOT -----------------------------------------------------------------------
  
  # PLOT SOM 
  if(plot) {
    cyto_plot(
      codes,
      channels = map,
      point_shape = 21,
      point_size = 2.5
    )
  }
  
  # SOM MAPPING ----------------------------------------------------------------
  
  # MAP DATA TO SOM
  message(
    "Mapping data to the trained SOM..."
  )
  
  # CREATE LINEAR CYTOSET
  x_som <- cytoset(
    structure(
      lapply(
        seq_along(x),
        function(z) {
          # SAMPLE PROGRESS
          message(
            cyto_names(x)[z]
          )
          # EXTRACT DATA -> LINEAR SCALE
          exprs <- cyto_data_extract(
            x,
            select = z,
            parent = parent,
            channels = channels,
            trans = trans,
            inverse = TRUE,
            format = "matrix",
            events = 1,
            coerce = FALSE,
            barcode = FALSE,
            markers = FALSE,
            copy = TRUE
          )[[1]][[1]]
          # MAP DATA TO SOM
          res <- cyto_func_call(
            "predict",
            list(
              som,
              newdata = data.matrix(
                exprs
              )
            )
          )$unit.classif
          # COUNTS
          cnts <- structure(
            rep(0, N),
            names = 1:N
          )
          tbl <- table(res)
          cnts[as.numeric(names(tbl))] <- tbl
          # PREPARE DATA - LINEAR CODES + DIM REDUCTION + COUNTS + PROPORTIONS
          as(
            cbind(
              cyto_exprs(
                codes[[1]],
                channels = c(channels, map),
                drop = FALSE
              ),
              cols = matrix(
                c(
                  unname(cnts),
                  unname(cnts)/sum(cnts)
                ),
                byrow = FALSE,
                ncol = 2,
                dimnames = list(
                  NULL,
                  c("SOM_counts", "SOM_freq")
                )
              )
            ),
            "cytoframe"
          )
        }
      ),
      names = cyto_names(x)
    )
  )
  cyto_details(x_som) <- pd[, 1:ncol(pd), drop = FALSE]
  
  # BARCODING FOR CYTO_PLOT()
  x_som <- cyto_barcode(
    x_som,
    type = "events"
  )
  
  # TRANSFER MARKER ASSSIGNMENTS
  if(length(markers) > 0) {
    markers <- markers[
      names(markers) %in% channels
    ]
    if(length(markers) > 0) {
      cyto_markers(x_som) <- markers
    }
  }
  
  # CODES ON TRANSFORMED SCALE + TRANSFORMERS ATTACHED
  if(!.all_na(trans)) {
    x_som <- cyto_transform(
      if(cyto_class(x, "GatingSet")) {
        GatingSet(x_som)
      }else {
        x_som
      },
      trans = som_trans,  # APPLY DATA TRANSFORMERS
      channels = channels,
      inverse = FALSE,
      plot = FALSE,
      quiet = TRUE
    )
  # LINEAR GATINGSET
  } else {
    if(cyto_class(x, "GatingSet")) {
      x_som <- GatingSet(x_som)
    }
  }
  
  # GATINGTEMPLATE
  if(!is.null(gatingTemplate)) {
    cyto_gatingTemplate_active(
      gatingTemplate
    )
  }
  
  # RETURN CYTOSET | GATINGSET
  return(x_som)

}

## CYTO_GATE_SOM ---------------------------------------------------------------

#' Map data to a gated self-organising map (SOM)
#'
#' \code{cyto_som()} provides a way to compress the original data into a more
#' manageable format for downstream processing. \code{cyto_gate_som()} transfers
#' the gates from a gated SOM \code{GatingSet} produced by \code{cyto_som()}
#' back onto the original data stored in a separate \code{GatingSet}.
#'
#' @param x an object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param parent name of the parent population which should be gated using the
#'   SOM, ideally the same population used to train the SOM in
#'   \code{cyto_som()}.
#' @param som object of class \code{GatingSet} containing the codebook vectors
#'   for the trained SOM.
#' @param alias name(s) of the populations gated on the SOM that should be
#'   transferred back onto the original data, set to the \code{descendants} of
#'   the \code{"root"} node by default.
#' @param channels name(s) of the channel(s) or marker(s) that should be used to
#'   map the data to the SOM, ideally the same set of channels used to train the
#'   SOM in \code{cyto_som()}.
#' @param hidden logical to indicate whether hidden nodes be included as
#'   descendants in \code{alias} when \code{alias} is not manually supplied, set
#'   to TRUE by default.
#' @param ... additional arguments passed to \code{cyto_som()} to train a new
#'   SOM when a \code{som} is not supplied manually.
#'
#' @return a GatingSet containing the original data with gates inherited from
#'   the gated SOM.
#'
#' @importFrom flowWorkspace gs_pop_add
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_som}}
#'
#' @export
cyto_gate_som <- function(x,
                          parent = "root",
                          som = NULL,
                          alias = NULL,
                          channels = NULL,
                          hidden = TRUE,
                          ...) {
  
  # CHECKS ---------------------------------------------------------------------
  
  # GATINGSET ONLY
  if(!cyto_class(x, "GatingSet")) {
    stop(
      "cyto_gate_som() only supports GatingHierarchy or GatingSet objects!"
    )
  }
  
  # SOM MISSING
  if(is.null(som)) {
    # TRAIN A NEW SOM?
    opt <- cyto_enquire(
      paste(
        "'som' must be a trained SOM returned by cyto_som(). Do you want to ",
        "train a new SOM? (Y|N)"
      ),
      options = c("T", "Y", "TRUE")
    )
    # TRAIN NEW SOM
    if(opt) {
      som <- cyto_som(
        x,
        parent = parent,
        ...
      )
    # SOM REQUIRED
    } else {
      return(NULL)
    }
  # SOM CHECKS
  } else {
    # GATINGSET REQUIRED
    if(!cyto_class(som, "GatingSet") & length(som) != length(x)) {
      stop(
        "'som' must be a GatingSet object created by cyto_som()!"
      )
    }
  }
  
  # CHANNELS MISSING
  if(is.null(channels)) {
    stop(
      paste0(
        "Supply the name(s) of the channels or markers used to train the SOM ",
        "to 'channels'."
      )
    )
  # CONVERT CHANNELS
  } else {
    channels <- cyto_channels_extract(
      x,
      channels = channels
    )
  }
  
  # GATINGTEMPLATE WARNING
  warning(
    paste0(
      "'cyto_gate_som()' does not transfer gates into a gatingTemplate for ",
      "'x'. It is recommended to export your SOM GatingSet using cyto_save()."
    )
  )  
  
  # KOHONEN
  cyto_require(
    "kohonen",
    source = "CRAN",
    ref = paste0(
      "Wehrens, R., & Buydens, L. M. C. (2007). Self- and Super-organizing ",
      "Maps in R: The kohonen Package. Journal of Statistical Software, ",
      "21(5), 1–19."
    )
  )
  
  # SOM GATING -----------------------------------------------------------------

  # CREATE GATES
  gate <- structure(
    lapply(
      seq_along(x),
      function(z) {
        # ALIAS
        if(is.null(alias)) {
          alias <- cyto_nodes_kin(
            x, 
            nodes = "root",
            type = "descendants",
            hidden = hidden
          )
        }
        # SOM CODES -> LINEAR SCALE
        som_codes <- cyto_data_extract(
          som,
          select = z,
          parent = "root",
          channels = channels,
          format = "matrix",
          markers = FALSE,
          inverse = TRUE,
          copy = TRUE
        )[[1]][[1]]
        # UNGATED SOM
        if(length(alias) == 0) {
          # DEFAULT SOM LABELS
          som_labels <- paste0(
            "SOM-",
            1:nrow(som_codes)
          )
          # GATED SOM
        } else {
          # PREPARE ALIAS
          if(any(alias %in% "root")) {
            ind <- which(alias %in% "root")
            alias[ind] <- paste0(
              "SOM-",
              ind
            )
          }
          # EXTRACT LABELS
          som_labels <- cyto_gate_indices(
            som,
            select = z,
            parent = "root",
            nodes = alias,
            labels = TRUE
          )[[1]]
        }
        # RECONSTRUCT SOM
        SOM <- structure(
          list(
            codes = list(som_codes),
            data = list(som_codes),
            whatmap = 1,
            user.weights = 1,
            dist.fcts = "sumofsquares",
            maxNA.fraction = 0
          ),
          class = "kohonen"
        )
        # EXTRACT DATA -> LINEAR SCALE
        exprs <- cyto_data_extract(
          x,
          select = z,
          parent = parent,
          channels = channels,
          format = "matrix",
          markers = FALSE,
          inverse = TRUE,
          copy = TRUE
        )[[1]][[1]]
        # MAP DATA TO SOM
        som_ind <- cyto_func_call(
          "predict",
          list(
            SOM,
            newdata = data.matrix(
              exprs
            )
          )
        )$unit.classif
        # PREPARE GATE
        res <- as(
          factor(
            som_labels[som_ind],
            levels = unique(
              som_labels
            )
          ),
          "filterResult"
        )
        res@filterId <- "cyto_gate_som"
        res@frameId <- cyto_names(x)[z]
        return(res)
      }
    ),
    names= cyto_names(x)
  )

  # ADD GATES TO GATINGSET
  gs_pop_add(
    x,
    parent = parent,
    gate = gate
  )
  
  # # RECOMPUTE
  # suppressMessages(
  #   recompute(gs)
  # )
  
  # RETURN SOM GATED GATINGSET
  return(x)
  
}
