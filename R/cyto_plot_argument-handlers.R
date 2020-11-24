## ARGUMENT HANDLERS -----------------------------------------------------------

#' New and improved argument splitter (repeat per plot now)
#' @param x list of named arguments to repeat.
#' @param plots numeric indicating the number of plots.
#' @param layers vector indicating the number of layers in each plot.
#'
#' @importFrom methods is
#' @importFrom utils type.convert
#'
#' @noRd
.cyto_plot_args_split <- function(x,
                                  plots,
                                  layers) {
  
  # Calls within a cyto_plot wrapper require manual entries for plots and
  # layers. The argument list must contain a gate and negate argument that must
  # be of length plots.
  
  # NUMBER OF PLOTS - NP -------------------------------------------------------
  
  # CALLED WITHIN CYTO_PLOT METHOD
  if (missing(plots)) {
    if (all(LAPPLY(x$x, "is", "flowFrame"))) {
      NP <- 1
      MTD <- "flowFrame"
    } else if (all(LAPPLY(x$x, function(z) {
      LAPPLY(z, function(y) {
        is(y, "flowFrame")
      })
    }))) {
      NP <- length(x$x)
      MTD <- "flowSet"
    }
    # CALLED WITHIN CYTO_PLOT WRAPPER
  } else {
    NP <- plots
    if (NP == 1) {
      MTD <- "flowFrame"
    } else {
      MTD <- "flowSet"
    }
  }
  
  # LAYERS PER PLOT - L --------------------------------------------------------
  
  # CALLED WITHIN CYTO_PLOT METHOD - VECTOR OF LENGTH NP
  if (missing(layers)) {
    if (MTD == "flowFrame") {
      L <- length(x$x)
    } else if (MTD == "flowSet") {
      L <- LAPPLY(x$x, "length")
    }
    # CALLED WITHIN CYTO_PLOT WRAPPER
  } else {
    if (length(layers) != NP) {
      message("Using the same number of layers per plot to split arguments.")
      layers <- rep(layers, length.out = NP)
    }
    L <- layers
  }
  
  # TOTAL LAYERS - TL ----------------------------------------------------------
  
  # VECTOR OF LENGTH NP
  TL <- sum(L)
  
  # GATE COUNT PER PLOT - GC ---------------------------------------------------
  
  # CALLED WITHIN CYTO_PLOT METHOD/WRAPPER
  if (MTD == "flowFrame") {
    if (all(LAPPLY(x$gate, function(z) {
      .all_na(z)
    }))) {
      GC <- 1
    } else {
      GC <- length(x$gate)
    }
  } else if (MTD == "flowSet") {
    if (.all_na(x$gate)) {
      x$gate <- split(rep(NA, length.out = NP), seq_len(NP))
    } else if (length(x$gate) != NP) {
      stop("Gates must be supplied per plot to split cyto_plot arguments.")
    }
    GC <- LAPPLY(x$gate, function(z) {
      if (.all_na(z)) {
        return(1)
      } else {
        return(length(z))
      }
    })
  }
  
  # GATED POPULATIONS PER LAYER - GP -------------------------------------------
  if (MTD == "flowFrame") {
    if (GC != 0) {
      GP <- LAPPLY(x$gate, function(z) {
        if (class(z) == "quadGate") {
          return(4)
        } else {
          return(1)
        }
      })
      GP <- sum(GP)
    } else {
      GP <- 1
    }
  } else if (MTD == "flowSet") {
    if (!all(GC == 0)) {
      GP <- lapply(x$gate, function(z) {
        LAPPLY(z, function(y) {
          if (class(y) == "quadGate") {
            return(4)
          } else {
            return(1)
          }
        })
      })
      GP <- unlist(LAPPLY(GP, "sum"))
      # NO GATES
    } else {
      GP <- rep(1, NP)
    }
  }
  
  # TOTAL GATED POPULATIONS PER LAYER - TGP ------------------------------------
  TGP <- sum(GP)
  
  # GATED/NEGATED POPULATIONS PER LAYER - GNP ----------------------------------
  
  # PREPARE NEGATE ARGUMENT
  x$negate <- rep(x$negate, length.out = NP)
  if(MTD == "flowSet"){
    x$negate <- split(x$negate, seq_len(NP))
  }
  if (MTD == "flowFrame") {
    if (GC != 0 &
        x$negate == TRUE &
        !"quadGate" %in% LAPPLY(x$gate, "is")) {
      GNP <- GP + 1
    } else {
      GNP <- GP
    }
  } else if (MTD == "flowSet") {
    GNP <- LAPPLY(seq_len(NP), function(z) {
      if (GC[z] != 0 &
          x$negate[z] == TRUE &
          !"quadGate" %in% LAPPLY(x$gate[[z]], "is")) {
        return(GP[z] + 1)
      } else {
        return(GP[z])
      }
    })
  }
  
  # TOTAL GATED/NEGATED POPULATIONS - TGNP -------------------------------------
  TGNP <- sum(GNP)
  
  # CYTO_PLOT ARGUMENTS --------------------------------------------------------
  
  # MISSING LISTED AS NULL -> CONVERTED TO EMPTY
  args <- formals("cyto_plot.flowFrame")
  lapply(names(args), function(z) {
    if (all(class(args[[z]]) == "name")) {
      args[[z]] <<- ""
    }
  })
  
  # FLOWSET ARGUMENTS (LABEL_TEXT_SIZE ETC)
  if(MTD == "flowSet"){
    fs_args <- formals("cyto_plot.flowSet")
    lapply(names(fs_args), function(z) {
      if (all(class(fs_args[[z]]) == "name")) {
        fs_args[[z]] <<- ""
      }
    })
    args <- fs_args[names(args)]
  }
  
  # REMOVE ARGUMENTS THAT SHOULD PREPARED ALREADY
  args <- args[-match(
    c(
      "x",
      "channels",
      "overlay",
      "gate", # arguments that require manual preparation
      "negate",
      "axes_trans",
      "axes_text",
      "..."
    ),
    names(args)
  )]
  
  # ARGUMENT NAMES
  arg_names <- names(args)
  
  # ARGUMENTS PER PLOT ---------------------------------------------------------
  
  # NEGATE HANDLED ABOVE
  plot_args <- c(
    "display",
    "axes_limits",
    "axes_limits_buffer",
    "margins",
    "popup",
    "xlab",
    "ylab",
    "xlim",
    "ylim",
    "hist_stat",
    "hist_stack",
    "hist_smooth",
    "hist_bins",
    "hist_bandwidth",
    "hist_cols",
    "point_cols",
    "point_col_scale",
    "point_fast",
    arg_names[grepl("title", arg_names)], # all title arguments
    arg_names[grepl("axes_text", arg_names)], # all axes_text arguments
    arg_names[grepl("axes_label", arg_names)], # all axes_label arguments
    arg_names[grepl("border", arg_names)], # all border arguments
    arg_names[grepl("grid", arg_names)], # all grid arguments
    "label",
    "label_position",
    "legend"
  )
  
  # UPDATE ARG_NAMES
  arg_names <- arg_names[-match(plot_args, arg_names)]
  
  # ARGUMENTS PER LAYER ------------------------------------------------------
  
  layer_args <- c(
    arg_names[grepl("contour", arg_names)], # contour arguments
    arg_names[grepl("hist_fill", arg_names)], # hist_fill arguments
    arg_names[grepl("hist_line", arg_names)], # hist_line arguments
    arg_names[grepl("legend_", arg_names)], # legend aes arguments
    arg_names[grepl("point_", arg_names)] # point args
  )
  
  # UPDATE ARG_NAMES
  arg_names <- arg_names[-match(layer_args, arg_names)]
  
  # ARGUMENTS PER POPULATION ---------------------------------------------------
  
  pop_args <- arg_names[grepl("label_", arg_names)]
  
  # UPDATE ARG_NAMES
  arg_names <- arg_names[-match(pop_args, arg_names)]
  
  # ARGUMENTS PER GATED POPULATION ---------------------------------------------
  
  gate_pop_args <- arg_names[grepl("gate_fill", arg_names)]
  
  # UPDATE ARG_NAMES
  arg_names <- arg_names[-match(gate_pop_args, arg_names)]
  
  # ARGUMENTS PER GATE ---------------------------------------------------------
  
  gate_args <- c(
    arg_names[grepl("gate", arg_names)]
  )
  
  # PREPARE ARGUMENTS ----------------------------------------------------------
  
  # ARGUMENT NAMES
  x_names <- names(x)
  
  # REPEAT AND FORMAT ARGUMENTS
  x <- lapply(seq_along(x), function(z) {
    # ARGUMENT
    arg <- names(x[z])
    # ARGUMENT EXISTS
    if (arg %in% names(args)) {
      # DEFAULT FOR ARGUMENT
      if (arg == "label_text" & !all(LAPPLY(x[[z]], ".empty"))) {
        arg_default <- NA
      } else {
        arg_default <- args[[arg]]
        if(is.call(arg_default)){
          arg_default <- eval(arg_default)
        }
      }
      # DEFAULT ARGUMENT LENGTH
      arg_default_length <- length(arg_default)
      # ARGUMENTS PER PLOT - EXPECTED ARGUMENT LENGTHS
      if (arg %in% plot_args) {
        # WATCH OUT COLOURS
        if(arg %in% c("point_cols",
                      "point_col_scale",
                      "hist_cols")){
          arg_default_length <- length(x[[arg]])
        }
        # ARGUMENT LENGTH PER PLOT
        arg_plot_length <- rep(arg_default_length, NP)
        # TOTAL ARGUMENT LENGTH
        arg_total_length <- sum(arg_plot_length)
        # ARGUMENTS PER LAYER - EXPECTED ARGUMENT LENGTHS  
      }else if(arg %in% layer_args){
        # ARGUMENT LENGTH PER PLOT
        arg_plot_length <- arg_default_length * L
        # TOTAL ARGUMENT LENGTH
        arg_total_length <- sum(arg_plot_length)
        # ARGUMENTS PER GATE - EXPECTED ARGUMENT LENGTHS
      }else if(arg %in% gate_args){
        # ARGUMENT LENGTH PER PLOT
        arg_plot_length <- arg_default_length * GC
        # TOTAL ARGUMENT LENGTH
        arg_total_length <- sum(arg_plot_length)
        # ARGUMENTS PER GATED POPULATION - EXPECTED ARGUMENT LENGTHS
      }else if(arg %in% gate_pop_args){
        # ARGUMENT LENGTH PER PLOT
        arg_plot_length <- arg_default_length * GP 
        # TOTAL ARGUMENT LENGTH
        arg_total_length <- sum(arg_plot_length)
        # ARGUMENTS PER POPULATION - EXPECTED ARGUMENT LENGTHS
      }else if(arg %in% pop_args){
        # ARGUMENT LENGTH PER PLOT
        arg_plot_length <- arg_default_length * GNP * L
        # TOTAL ARGUMENT LENGTH
        arg_total_length <- sum(arg_plot_length)
      }
      # FILL ARGUMENT WITH PLACE HOLDER & SPLIT
      arg_repeat <- rep(c(
        x[[arg]],
        rep("*-*", arg_total_length)
      ),
      length.out = arg_total_length
      )
      # SPLIT INTO PLOTS
      arg_split <- split(
        arg_repeat,
        rep(seq_len(NP),
            times = arg_plot_length
        )
      )
      # PLOT ARGUMENTS - COMPLETE/INCOMPLETE/MISSING
      arg_complete <- which(LAPPLY(arg_split, function(z) {
        ind <- z == "*-*"
        if(any(is.na(ind))){
          ind[is.na(ind)] <- FALSE
        }
        return(!any(ind))
      }))
      arg_incomplete <- which(LAPPLY(arg_split, function(z) {
        ind <- z == "*-*"
        if(any(is.na(ind))){
          ind[is.na(ind)] <- FALSE
        }
        return(any(ind) & !all(ind))
      }))
      arg_missing <- which(LAPPLY(arg_split, function(z) {
        ind <- z == "*-*"
        if(any(is.na(ind))){
          ind[is.na(ind)] <- FALSE
        }
        return(all(ind))
      }))
      # MISSING TREATED AS INCOMPLETE (DIFFERENT PLOTS)
      if(!length(unique(arg_plot_length)) == 1){
        arg_incomplete <- sort(c(arg_incomplete, arg_missing))
        arg_missing <- NULL
      }
      # FILL INCOMPLETE
      if(length(arg_incomplete) != 0){
        lapply(arg_incomplete, function(z) {
          # REPEAT ARGUMENT (CHARACTERISTIC - MULTI-LEVEL)
          if(length(arg_split[[z]]) > arg_default_length &
             any(LAPPLY(c("fill_alpha",
                          "col_alpha",
                          "line",
                          "shape",
                          "size",
                          "font",
                          "label_fill",
                          "gate_fill"), function(r){
                            grepl(r, arg, ignore.case = TRUE)
                          }))){
            arg_ind <- which(arg_split[[z]] == "*-*")
            arg_length <- length(arg_ind)
            arg_split[[z]][arg_ind] <<- rep(arg_split[[z]][-arg_ind],
                                            length.out = arg_length)
            arg_complete <<- c(arg_complete, z)
            # FILL ARGUMENT WITH DEFAULT
          }else{
            arg_ind <- which(arg_split[[z]] == "*-*")
            arg_length <- length(arg_ind)
            arg_replace <- rep(arg_default,
                               length.out = arg_length)
            arg_split[[z]][arg_ind] <<- arg_replace
            arg_complete <<- c(arg_complete, z)
          }
        })
      }
      # REPLACE EMPTY 
      if(length(arg_missing) != 0){
        arg_split[arg_missing] <- arg_split[
          rep(sort(arg_complete), length.out = length(arg_missing))
        ]
      }
      # RETURN ARGUMENTS SPLIT BY PLOT
      if(MTD == "flowFrame"){
        arg_split <- arg_split[[1]]
      }
      return(arg_split)
      # OTHER ARGUMENTS
    }else{
      # CHANNELS
      if(arg == "channels"){
        if(class(x[[z]]) == "list"){
          arg_split <- rep(x[[z]], length.out = NP)
          return(arg_split)
        }else{
          arg_split <- rep(list(x[[z]]), length.out = NP)
          if(MTD == "flowFrame"){
            arg_split <- arg_split[[1]]
          }
          return(arg_split)
        }
        # AXES_TRANS
      }else if(arg == "axes_trans"){
        # TRANSFORMERLIST
        if (is(x[[z]], "transformerList")) {
          arg_split <- rep(list(x[[z]]), length.out = NP)
          # LIST OF TRANSFORMERLISTS
        } else {
          arg_split <- rep(x[[z]], length.out = NP)
        }
        if(MTD == "flowFrame"){
          arg_split <- arg_split[[1]]
        }
        return(arg_split)
      # AXES_TEXT
      }else if(arg == "axes_text") {
        # LIST ALREADY PREPARED
        if(!is(x[[z]], "list")) {
          # REPEAT
          if(length(x[[z]]) < 2) {
            arg_split <- split(rep(c(x[[z]], TRUE), length.out = NP * 2), 
                               rep(seq_len(NP), each = 2))
          } else {
            arg_split <- split(rep(x[[z]], length.out = NP * 2), 
                               rep(seq_len(NP), each = 2))
          }
          if(MTD == "flowFrame"){
            arg_split <- arg_split[[1]] # don't extract lists
          }
        } else {
          arg_split <- x[[z]]
        }
        return(arg_split)
      }else{
        return(x[[z]])
      }
    }
  })
  names(x) <- x_names
  
  # NUMERIC CONVERSION
  x <- lapply(seq_along(x), function(z){
    arg <- names(x[z])
    if(arg %in% names(args)){
      if(MTD == "flowFrame"){
        w <- x[[arg]]
        if(is.vector(w)) {
          empty_ind <- which(LAPPLY(w, ".empty"))
          if(length(empty_ind) == 0) {
            return(type.convert(w, as.is = TRUE))
          } else if(length(empty_ind) == length(w)) {
            return(w)
          } else {
            w[-empty_ind] <- type.convert(w[-empty_ind], as.is = TRUE)
            return(w)
          }
        } else {
          return(w)
        }
      }else if(MTD == "flowSet"){
        lapply(x[[arg]], function(w){
          if(is.vector(w)) {
            empty_ind <- which(LAPPLY(w, ".empty"))
            if(length(empty_ind) == 0) {
              return(type.convert(w, as.is = TRUE))
            } else if(length(empty_ind) == length(w)) {
              return(w)
            } else {
              w[-empty_ind] <- type.convert(w[-empty_ind], as.is = TRUE)
              return(w)
            }
          } else {
            return(w)
          }
        })
      }
    } else {
      return(x[[z]])
    }
  })
  names(x) <- x_names
  
  # RETURN SPLIT ARGUMENTS
  return(x)
}

