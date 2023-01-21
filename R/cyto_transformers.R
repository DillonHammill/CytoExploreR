## CYTO_TRANSFORMERS -----------------------------------------------------------

# CytoExploreR provides wrapper functions for data transformations as defined
# within the cytoverse suite of cytometry packages. cyto_transformers_define()
# makes it easy to customise and visualise data transformations on the fly
# without changing the underlying data. Additional helper functions
# cyto_transformer_combine and cyto_transformer_extract are supplied to make it
# easier to work with transformerList objects.

# The old cyto_transformer_log(), cyto_transformer_arcsinh(),
# cyto_transformer_biex() and cyto_transformer_logicle() APIs have now been
# merged into a single consistent API called cyto_transformers_define(). Now
# users can specify transformations on a per channel basis using the 'type'
# argument. Maintaining support for the old cyto_transformer functions is
# cumbersome and as such these APIs are now defunct.

#' @noRd
#' @export
cyto_transformer_log <- function(...){
  .Defunct("cyto_transformers_define")
}

#' @noRd
#' @export
cyto_transformer_biex <- function(...){
  .Defunct("cyto_transformers_define")
}

#' @noRd
#' @export
cyto_transformer_arcsinh <- function(...){
  .Defunct("cyto_transformers_define")
}

#' @noRd
#' @export
cyto_transformer_logicle <- function(...){
  .Defunct("cyto_transformers_define")
}

# For consistency, cyto_transformer_combine and cyto_transformer_extract have
# been renamed to instead use the cyto_transformers_ prefix. As such, the old
# cyto_transformer_ APIs are now defunct as well.

#' @noRd
#' @export
cyto_transformer_combine <- function(...){
  .Defunct("cyto_transformers_combine")
}

#' @noRd
#' @export
cyto_transformer_extract <- function(...){
  .Defunct("cyto_transformers_extract")
}

## CYTO_TRANSFORMERS_DEFINE ----------------------------------------------------

#' Define transformers for channels
#'
#' \code{cyto_transformers_define()} is a convenient function that allows users
#' to optimise cytoverse transformer definitions for channels, as
#' transformations are not applied directly to the supplied data and data
#' transformations are automatically visualised. The type of transformation to
#' apply can be controlled through the \code{type} argument and links to
#' relevant documentation for additional arguments are supplied below. Separate
#' transformer definitions can be combined into the same transformerList using
#' \code{\link{cyto_transformers_combine}}. Once optimisied, transformers can be
#' applied to the data using \code{\link{cyto_transform}}.
#'
#' @param x an object of class \code{flowFrame}, \code{flowSet},
#'   \code{GatingHierarchy} or \code{GatingSet}.
#' @param parent name of the parent population to use for generating the
#'   transformation(s).
#' @param channels name(s) of channel(s)/marker(s) for which transformation
#'   functions must be generated. Set to all the fluorescent channels by default
#'   if no channels/markers are supplied.
#' @param type a vector indicating the type of transformations to apply to each
#'   channel, options include \code{"log"}, \code{"asinh"}, \code{"asinh_Gml2"},
#'   \code{"biexponential"} or \code{"logicle"}. To set different transformation
#'   types for each channel, supply a vector of transformation types named with
#'   each of the channels. Users can also define their own custom transformation
#'   functions by supplying a list of transformation and inverse transformation
#'   functions for each channel (e.g. list("PE-A" = list(transform = "asinh",
#'   inverse = "sinh"))). All channels will have the same settings for each
#'   transformation type, multiple calls to \code{cyto_transformers_define()}
#'   can be made to modify these settings per channel. Multiple transformerLists
#'   can be joined together using \code{cyto_transformers_combine()}.
#' @param select list of selection criteria passed to \code{cyto_select} to
#'   select a subset of samples for \code{cyto_plot}.
#' @param plot logical indicating whether the results of transformations should
#'   be plotted using \code{cyto_plot}.
#' @param events number or frequency of events to display in plots, set to 50000
#'   events by default. See \code{\link{cyto_plot}} for details.
#' @param progress logical indicating whether to display the progress bar(s)
#'   when computing transformer definitions, set to TRUE by default.
#' @param ... additional arguments passed to
#'   \code{\link[flowWorkspace:flowjo_log_trans]{flowjo_log_trans}},
#'   \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}},
#'   \code{\link[flowWorkspace:flowjo_biexp]{flowjo_biexp}} or
#'   \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}} and
#'   \code{\link{cyto_plot}}.
#'
#' @return a \code{transformerList} object containing the transformation
#'   definitions for each channel.
#'
#' @importFrom flowCore logicleTransform inverseLogicleTransform
#'   CytoExploreR_.estimateLogicle
#' @importFrom flowWorkspace flowjo_log_trans asinh_Gml2 flowjo_biexp flow_trans
#'   transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link[flowWorkspace:flowjo_log_trans]{flowjo_log_trans}}
#' @seealso \code{\link[flowWorkspace:asinh_Gml2]{asinh_Gml2}}
#' @seealso \code{\link[flowWorkspace:flowjo_biexp]{flowjo_biexp}}
#' @seealso \code{\link[flowCore:logicleTransform]{estimateLogicle}}
#' @seealso \code{\link[flowWorkspace:estimateLogicle]{estimateLogicle}}
#'
#' @export
cyto_transformers_define <- function(x,
                                     parent = "root",
                                     channels = NULL,
                                     type = "logicle",
                                     select = NULL,
                                     events = 50000,
                                     plot = TRUE,
                                     progress = TRUE,
                                     ...) {
  
  # NOTE: AXES_LIMITS CAUSE ARGUMENT CONFLICT WITH A FOR LOGICLE
  
  # RESET PROGRESS BAR ON EXIT
  on.exit({
    cyto_option("CytoExploreR_progress", NULL)
  })
  
  # PREPARE CHANNELS
  if(is.null(channels)) {
    channels <- cyto_fluor_channels(x)
  } else {
    channels <- cyto_channels_extract(x, channels)
  }
  
  # PREPARE ARGUMENTS
  args <- list(...)
  
  # TYPE - LIST OF CUSTOM TRANSFORMERS
  if(cyto_class(type, "list", TRUE)) {
    # LIST OF TRANSFORMERS - LIST(TRANSFORM = F(), INVERSE = G())
    if(!cyto_class(type[[1]], "list", TRUE)) {
      type <- structure(
        rep(
          list(type),
          length(channels)
        ),
        names = channels
      )
    }
    # CHECK TYPE FORMAT
    if(!all(LAPPLY(type, function(w){
      cyto_class(w, "list", TRUE) & length(w) == 2
    }))) {
      stop(
        paste0(
          "Custom transformers should be supplied to 'type' in this form: ",
          "\n",
          "list( \n",
          "  'PE-A' = list('transform' = f, 'inverse' = g) \n",
          ") \n",
          "where: \n",
          "'f' = name of the transform function \n",
          "'g' = name of the inverse transform function"
        )
      )
    }
  }
  
  # PREPARE TYPE
  if(is.null(names(type))) {
    type <- rep(type, length(channels))
    names(type) <- channels
  }
  
  # PROGRESS BAR
  if(progress) {
    pb <- cyto_progress(
      label = "cyto_transformers_define()",
      total = length(channels)
    )
  }
  
  # TRANSFORMATION DEFINITIONS
  transformer_list <- structure(
    lapply(
      channels, 
      function(z) {
        # TYPE VECTOR
        if(!cyto_class(type, "list", TRUE)) {
          # LOG TRANSFORM
          if(grepl("^log$", type[z], ignore.case = TRUE)) {
            trans <- cyto_func_execute("flowjo_log_trans", args)
          # ARCSINH TRANSFORM
          } else if(grepl("^a", type[z], ignore.case = TRUE)) {
            # ARCSINH GML2 FLOWWORKSPACE
            if(grepl("^arc", type[z], ignore.case = TRUE) |
               grepl("g", type[z], ignore.case = TRUE)) {
              # TRANSFORMERS
              args[["inverse"]] <- FALSE
              trans <- cyto_func_execute("asinh_Gml2", args)
              # INVERSE TRANSFORMERS
              args[["inverse"]] <- TRUE
              inv_trans <- cyto_func_execute("asinh_Gml2", args)
              # COMBINE TRANSFORMERS
              trans <- flow_trans(
                "arcsinh_Gml2",
                trans@.Data,
                inv_trans@.Data
              )
            # ARCSINH - FLOWVS ESTIMATE COFACTOR
            } else {
              # FLOWVS - CANNOT TURN OFF PLOTS & MESSAGES USE CAT
              # COFACTOR SUPPLIED MANUALLY
              if("cofactor" %in% names(args)) {
                cf <- args[["cofactor"]]
              # ESTIMATE COFACTOR USING FLOWVS
              } else {
                # FLOWVS
                cyto_require("flowVS",
                             source = "BioC",
                             repo = NULL,
                             version = NULL,
                             ref = paste0(
                               "Azad A, Rajwa B, Pothen A (2016). flowVS:",
                               " channel-specific variance stabilisation in",
                               " flow cytometry, BMC Bioinformatics 17(291)."))
                # ESTIMATE COFACTOR TO STABILISE VARAIANCE
                message(
                  paste0(
                    "Using flowVS to estimate cofactor for ", z, "..."
                  )
                )
                invisible(
                  capture.output(
                    cf <- cyto_func_call(
                      "flowVS::estParamFlowVS",
                      list(
                        cyto_data_extract(
                          x,
                          parent = parent,
                          select = select,
                          format = "cytoset",
                          channels = z
                        )[[1]],
                        z
                      )
                    )
                  )
                )
              }
              # TRANSFORMER DEFINITIONS
              asinh_trans <- function(x, cofactor = cf) {
                asinh(x/cofactor)
              }
              sinh_trans <- function(x, cofactor = cf) {
                sinh(x/cofactor)
              }
              # TRANSFORMERS
              trans <- flow_trans(
                "arcsinh",
                asinh_trans,
                sinh_trans
              )
            }
          # BIEXPONENTIAL TRANSFORM
          } else if(grepl("^biex", type[z], ignore.case = TRUE)) {
            # DEFAULT WIDTHBASIS
            if(is.null(args[["widthBasis"]])) {
              args[["widthBasis"]] <- -100
            }
            # MAXVALUE - HARD CODED
            if(is.null(args[["maxValue"]])) {
              args[["maxValue"]] <- 262144
            }
            # CHECK INSTRUMENT RANGE
            rng <- range(
              cyto_apply(
                x,
                parent = parent,
                select = select,
                channels = z,
                input = "cytoframe",
                FUN = "range",
                type = "instrument"
              ),
              na.rm = TRUE
            )
            # MAXVALUE ADJUST - INSTRUMENT RANGE
            if(max(rng) > args[["maxValue"]]) {
              args[["maxValue"]] <- max(rng)
            }
            # TRANSFORMERS
            args[["inverse"]] <- FALSE
            trans <- cyto_func_execute("flowjo_biexp", args)
            # INVERSE TRANSFORMERS
            args[["inverse"]] <- TRUE
            inv_trans <- cyto_func_execute("flowjo_biexp", args)
            trans <- flow_trans(
              "biexponential",
              trans@.Data,
              inv_trans@.Data
            )
          # LOGICLE TRANSFORM
          } else if(grepl("^logicle$", type[z], ignore.case = TRUE)) {
            
            # ESTIMATE LOGICLE PARAMETERS - SEE FLOWCORE:::.lgclTrans()
            # r = data range (user supply data/instrument) - default instrument
            # t = top of scale
            # a = additional negative decades
            # m = full length of transformed display in decades
            # q = quantile of negative data value (w adjustment)
            # w = linear decades
            
            # RANGE
            if(!"r" %in% names(args)) {
              # DON'T UPDATE GLOBALLY - CHANNEL SPECIFIC RANGES
              # TODO: DO WE NEED TO COMPUTE INSTRUMENT RANGE ON EVERY SAMPLE?
              args$r <- range(
                cyto_apply(
                  x,
                  parent = parent,
                  select = select,
                  channels = z,
                  input = "cytoframe",
                  FUN = "range",
                  type = "instrument"
                ),
                na.rm = TRUE
              )
            }
            # TOP OF SCALE
            if(!"t" %in% names(args)) {
              args$t <- max(args$r)
            }
            # COMPUTE M FROM RANGE
            if(!"m" %in% names(args)) {
              # FLOWCORE USES DEFAULT 4.5 & ADDS +1 FOR ESTIMATE
              args$m <- log10(args$t) - 1 # FLOWCORE USES + 1
              # ROUND UP TO NEAREST HALF DECADE
              d <- args$m - floor(args$m)
              if(d > 0) {
                if(d < 0.5) {
                  args$m <- args$m + (0.5 - d)
                } else {
                  args$m <- args$m + (1 - d)
                }
              }
            }
            # COMPUTE W - (R MOST NEGATIVE VALUES FOR DISPLAY)
            if(!"w" %in% names(args)) {
              args$w <- 0
              # EXTRACT DATA < ZERO
              d <- do.call(
                "c",
                cyto_apply(
                  x,
                  parent = parent,
                  select = select,
                  channels = z,
                  input = "column",
                  FUN = function(v) {
                    v[v < 0]
                  },
                  simplify = FALSE
                )
              )
              # COMPUTE MINIMUM VALUE
              if(length(d) > 0) {
                args$p <- .Machine$double.eps + quantile(
                  d,
                  0.05
                )
                args$w <- (args$m - log10(args$t/abs(args$p))) / 2
              }
            }
            # CHECK W > 0
            if(args$w < 0) {
              stop(
                paste0(
                  "'w' is negative - try increasing m above ",
                  args$m
                )
              )
            }
            # # ADDITIONAL NEGATIVE RANGE
            # if(!"a" %in% names(args)) {
            #   args$a <- args$m
            # }
            # LOGICLE TRANSFORM - DROP R AND P ARGUMENTS
            trans <- cyto_func_execute(
              "logicleTransform",
              c(list("logicle"), args)
            )
            # INVERSE LOGICLE TRANSFORM
            inv <- inverseLogicleTransform(
              trans
            )
            # LOGICLE TRANSFORMERS
            trans <- flow_trans(
              "logicle",
              trans@.Data,
              inv@.Data
            )

            # # ESTIMATELOGICLE CYTOEXPLORER WRAPPER DOESN'T EXPOSE ARGUMENTS
            # estimateLogicle_args <- c("t", "m", "a", "q")
            # trans <- do.call("CytoExploreR_.estimateLogicle",
            #                  c("x" = cyto_data_extract(
            #                    x,
            #                    coerce = TRUE)[[1]][[1]],
            #                    list("channels" = z),
            #                    args[names(args) %in% estimateLogicle_args]))
            # inv_trans <- inverseLogicleTransform(trans[[z]])
            # flow_trans(
            #   "logicle",
            #   trans[[z]]@.Data,
            #   inv_trans@.Data
            # )
            
          # UNSUPPORTED TRANSFORM
          } else {
            stop(
              paste(
                type[z], "is not a supported transformation type!"
              )
            )
          }
        # TYPE - CUSTOM TRANSFORMERS
        } else {
          # PREPARE TRANSFORM FUNCTION
          trans <- cyto_func_match(
            if(is.null(names(type[[z]]))) {
              type[[z]][[1]]
            } else {
              type[[z]][[grep("^t", names(type[[z]]), ignore.case = TRUE)]]
            }
          )
          # PREPARE INVERSE TRANSFORM FUNCTION
          inv <- cyto_func_match(
            if(is.null(names(type[[z]]))) {
              type[[z]][[2]]
            } else {
              type[[z]][[grep("^i", names(type[[z]]), ignore.case = TRUE)]]
            }
          )
          # PREPARE TRANSFORMERS
          trans <- flow_trans(
            "custom",
            trans, # TRANSFORM
            inv    # INVERSE TRANSFORM
          )
        }
        # INCREMENT PROGRESS BAR
        if(progress){
          cyto_progress(pb)
        }
        return(trans)
      }
    ), names = channels
  )
  
  # TRANSFORMERLIST
  transformer_list <- transformerList(
    channels,
    transformer_list
  )
  
  # PLOT TRANSFORMATIONS
  if(plot) {
    # EXTRACT DATA FOR TRANSFORMATION
    cs <- cyto_data_extract(
      x,
      parent = parent,
      select = select,
      channels = channels,
      format = "cytoset",
      events = events,
      coerce = TRUE,
      copy = TRUE
    )[[1]]
    # APPLY TRANSFORMERS
    cs <- cyto_transform(
      cs,
      trans = transformer_list,
      plot = FALSE,
      quiet = TRUE
    )
    # PLOT DATA TRANSFORMATIONS
    tryCatch(
      cyto_plot_profile(
        cs,
        channels = channels,
        axes_trans = transformer_list,
        axes_limits = "machine",
        events = 1 # SAMPLED ABOVE FOR FASTER TRANSFORMATIONS
      ), 
      error = function(e){
        if(.grepl("figure margins too large", e$message)) {
          message("Insufficient plotting space to display transformations!")
          if(!"RStudioGD" %in% names(dev.cur())) {
            dev.off()
          }
        } else {
          e
        }
      }
    )
  }
  
  # RETURN TRANSFORMERLIST
  return(transformer_list)
  
}

## CYTO_TRANSFORMER_EXTRACT ----------------------------------------------------

#' Extract transformers from a GatingHierarchy or GatingSet
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingHierarchy-class]{GatingHierarchy}}
#'   or \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @return transformerList or NA.
#'
#' @importFrom flowWorkspace gh_get_transformations
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples 
#' 
#' library(CytoExploreRData)
#' 
#' # Activation flowSet
#' fs <- Activation
#' 
#' # Activation GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Apply transformations
#' gs <- cyto_transform(gs, channels = c("CD4", "CD8"))
#'
#' # Extract transformations
#' trans <- cyto_transformers_extract(gs)
#'
#' @export
cyto_transformers_extract <- function(x){
  
  # LIST
  if(cyto_class(x, "list", TRUE)) {
    x <- x[[1]]
  }
  
  # GATINGSET
  if(cyto_class(x, "GatingSet", TRUE)){
    x <- x[[1]]
  }
  
  # GATINGHIERARCHY
  if(cyto_class(x, "GatingHierarchy", TRUE)){
    transformers <- tryCatch(
      gh_get_transformations(
        x, 
        only.function = FALSE
      ),
      error = function(e) {
        warning(e$message)
        return(NULL)
      }
    )
    if(length(transformers) != 0){
      transformers <- cyto_transformers_combine(transformers)
    } else {
      transformers <- NA
    }
  # FLOWFRAME/FLOWSET
  } else {
    transformers <- NA
  }
  
  # RETURN TRANSFORMERLIST
  return(transformers)
}

## CYTO_TRANSFORMERS_COMBINE ---------------------------------------------------

#' Combine cyto_transformer definitions
#'
#' \code{cyto_transformers_combine} makes it easy to combine transformation
#' definitions obtained from \code{\link{cyto_transformers_define}} prior to
#' applying these transformations to the data using
#' \code{\link{cyto_transform}}.
#'
#' @param ... objects of class
#'   \code{\link[flowWorkspace:transformerList]{transformerList}} to be combined
#'   into a single \code{transformerList} object.
#'
#' @importFrom flowWorkspace transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_transformers_define}}
#'
#' @export
cyto_transformers_combine <- function(...) {
  
  # Combine transformerList objects
  transformer_list <- c(...)
  
  # Convert to transformerList
  transformer_list <- transformerList(
    names(transformer_list),
    transformer_list
  )
  
  # Return transformations in a single transformerList for cyto_transform
  return(transformer_list)
  
}

## CYTO_TRANSFORMERS_DEPARSE ---------------------------------------------------

#' Convert transformerList into character string
#'
#' @param trans object of class transformerList containing transformer
#'   definitions.
#'
#' @return a character string defining the transformerList object.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
cyto_transformers_deparse <- function(trans) {
  
  # SEPARATOR REQUIRED FOR SPLITTING -> '***'
  paste0(
    deparse(
      trans
    ),
    collapse = "***"
  )
  
}

## CYTO_TRANSFORMERS_PARSE -----------------------------------------------------

#' Convert a character string of transformers into transformerList
#'
#' @param trans a charcter string defining a transformerList as created using
#'   \code{cyto_transformers_deparse()}.
#'   
#' @return a transformerList object containing the transformer definitions.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
cyto_transformers_parse <- function(trans) {
  
  eval(
    parse(
      text = strsplit(
        trans,
        "***",
        fixed = TRUE
      )[[1]]
    )
  )
  
}

## .CYTO_TRANSFORMERS_COMPLETE -------------------------------------------------

#' Get complete transformerList for compensation functions
#' @importFrom methods is
#' @noRd
.cyto_transformers_complete <- function(x, ...){
  UseMethod(".cyto_transformers_complete")
}

#' @noRd
.cyto_transformers_complete.flowFrame <- function(x,
                                                 axes_trans = NA) {
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # Transformations supplied
  if(!.all_na(axes_trans)){
    # Must be transformerList
    if(!is(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    # Some transformations are missing
    if(!all(channels %in% names(axes_trans))){
      # Missing channels
      trans_chans <- channels[channels %in% names(axes_trans)]
      linear_chans <- channels[!channels %in% names(axes_trans)]
      # Get remaining transformers
      trans <- cyto_transformers_define(
        x,
        channels = linear_chans,
        type = "biex",
        plot = FALSE
      )
      # Combine transformers
      axes_trans <- cyto_transformers_combine(
        axes_trans[trans_chans], 
        trans
      )
    }
    # No transformations supplied
  }else{
    # Get biex transformers
    axes_trans <- cyto_transformers_define(
      x, 
      channels = channels, 
      type = "biex",
      plot = FALSE
    )
    
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformers_complete.flowSet <- function(x, 
                                               axes_trans = NA){
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # Transformations supplied
  if(!.all_na(axes_trans)){
    # Must be transformerList
    if(!is(axes_trans, "transformerList")){
      stop("'axes_trans' must be an object of class transformerList.")
    }
    # Some stranformations are missing
    if(!all(channels %in% names(axes_trans))){
      # Missing channels
      trans_chans <- channels[channels %in% names(axes_trans)]
      linear_chans <- channels[!channels %in% names(axes_trans)]
      # Get remaining transformers
      trans <- cyto_transformers_define(
        x,
        channels = linear_chans,
        type = "biex",
        plot = FALSE
      )
      # Combine transformers
      axes_trans <- cyto_transformers_combine(
        axes_trans[trans_chans], 
        trans
      )
    }
    # No transformations supplied
  }else{
    # Get biex transformers
    axes_trans <- cyto_transformers_define(
      x, 
      channels = channels,
      type = "biex",
      plot = FALSE
    )
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformers_complete.GatingHierarchy <- function(x,
                                                       axes_trans = NA){
  
  # Fluorescent channels
  channels <- cyto_fluor_channels(x)
  
  # Extract transformations from GatingHierarchy
  gh_trans <- gh_get_transformations(x,
                                     only.function = FALSE)
  if(length(gh_trans) != 0){
    gh_trans <- cyto_transformers_combine(gh_trans)
  }
  
  # In case NULL axes_trans
  if(is.null(axes_trans)){
    axes_trans <- NA
  }
  
  # No transformations found in GatingHierarchy
  if(length(gh_trans) == 0){
    # No transformations supplied
    if(.all_na(axes_trans)){
      # Get biex transformers for all channels
      axes_trans <- cyto_transformers_define(
        x, 
        channels = channels, 
        type = "biex",
        plot = FALSE
      )
      # Transformations supplied
    }else if(!.all_na(axes_trans)){
      # All channels covered
      if(all(channels %in% names(axes_trans))){
        axes_trans <- axes_trans
        # Some channels covered
      }else if(any(channels %in% names(axes_trans))){
        # Channels not covered
        trans_chans <- channels[channels %in% names(axes_trans)]
        linear_chans <- channels[!channels %in% names(axes_trans)]
        # Get biex transformers for missing channels
        trans <- cyto_transformers_define(
          x,
          channels = linear_chans,
          type = "biex",
          plot = FALSE
        )
        # Combine transformers into transformerList
        axes_trans <- cyto_transformers_combine(
          axes_trans[trans_chans],
          trans
        )
        # No channels covered
      }else if(!any(channels %in% names(axes_trans))){
        # Get biex transformers for all channels
        axes_trans <- cyto_transformers_define(
          x,
          channels = channels,
          type = "biex",
          plot = FALSE
        )
      }
    }
    
    # Transformations found in GatingHierarchy
  }else{
    # All the transformations have been applied to GatingHierarchy
    if(all(channels %in% names(gh_trans))){
      # Regardless of axes_trans use complete gh_trans
      axes_trans <- cyto_transformers_combine(gh_trans[channels])
      # Some transformations have been applied to GatingHierachy
    }else if(any(channels %in% names(gh_trans))){
      # Channels not covered
      gh_trans_chans <- channels[channels %in% names(gh_trans)]
      linear_chans <- channels[!channels %in% names(gh_trans)]
      # No transformations supplied
      if(.all_na(axes_trans)){
        # Get biex transformers for missing channels
        trans <- cyto_transformers_define(
          x,
          channels = linear_chans,
          type = "biex",
          plot = FALSE
        )
        # Combine transformers into transformerList
        axes_trans <- cyto_transformers_combine(gh_trans[gh_trans_chans], trans)
        # Some transformations supplied
      }else if(!.all_na(axes_trans)){
        # All missing transformations covered
        if(all(linear_chans %in% names(axes_trans))){
          axes_trans <- cyto_transformers_combine(
            gh_trans[gh_trans_chans],
            axes_trans[linear_chans]
          )
          # Some missing transformations covered
        }else if(any(linear_chans %in% names(axes_trans))){
          # Channels not covered
          trans_chans <- linear_chans[linear_chans%in% names(axes_trans)]
          extra_chans <- linear_chans[!linear_chans%in% names(axes_trans)]
          # Get biex transformers for missing channels
          trans <- cyto_transformers_define(
            x,
            channels = extra_chans,
            type = "biex",
            plot = FALSE
          )
          # Combine transformers into transformerList
          axes_trans <- cyto_transformers_combine(
            gh_trans[gh_trans_chans],
            axes_trans[trans_chans],
            trans
          )
          # No missing transformations covered
        }else if(!any(linear_chans %in% names(axes_trans))){
          # Get biex transformers for all channels
          axes_trans <- cyto_transformers_define(
            x,
            channels = linear_chans,
            type = "biex",
            plot = FALSE
          )
        }
      }
      # Applied transformations are not in fluorescent channels
    }else if(!any(channels %in% names(gh_trans))){
      # No transformations supplied
      if(.all_na(axes_trans)){
        # Get biex transformers for all channels
        axes_trans <- cyto_transformers_define(
          x,
          channels = channels,
          type = "biex",
          plot = FALSE
        )
        # Some transformations supplied
      }else if(!.all_na(axes_trans)){
        # All channels covered by axes_trans
        if(all(channels %in% names(axes_trans))){
          axes_trans <- cyto_transformers_combine(
            axes_trans[channels]
          )
          # Some channels are covered by axes_trans
        }else if(any(channels %in% names(axes_trans))){
          # Channels not covered
          trans_chans <- channels[channels %in% names(axes_trans)]
          linear_chans <- channels[!channels %in% names(axes_trans)]
          # Get biex transformers for missing channels
          trans <- cyto_transformers_define(
            x,
            channels = linear_chans,
            type = "biex",
            plot = FALSE
          )
          # Combine transformers into transformerList
          axes_trans <- cyto_transformers_combine(
            axes_trans[trans_chans],
            trans
          )
          # No channels covered by axes_trans
        } else if(!any(channels %in% names(axes_trans))){
          # Get biex transformers for all channels
          axes_trans <- cyto_transformers_define(
            x,
            channels = channels,
            type = "biex",
            plot = FALSE
          )
        }
      }
    }
  }
  
  # Return complete transformerList
  return(axes_trans)
  
}

#' @noRd
.cyto_transformers_complete.GatingSet <- function(x,
                                                 axes_trans = NA){
  
  # Extract GatingHierarchy
  x <- x[[1]]
  
  # Make call to GatingHierarchy method
  axes_trans <- .cyto_transformers_complete(x,
                                           axes_trans = axes_trans)
  
  # Return complete transformerList
  return(axes_trans)
  
}
