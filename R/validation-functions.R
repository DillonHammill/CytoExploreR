#' Check Supplied Channels
#'
#' \code{cyto_channel_check} will check whether the supplied channels are valid
#' for the \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowset}} or
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} .
#' \code{cyto_channel_check} will also return the channels if the marker names
#' are supplied.
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param channels vector of channel names (e.g. c("PE-A","APC-A")) or marker
#'   names.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to TRUE by default to accept 1 or 2 channels only.
#'
#' @return vector of valid channel names.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
setGeneric(
  name = "cyto_channel_check",
  def = function(x, channels, plot) {
    standardGeneric("cyto_channel_check")
  }
)

#' Check Supplied Channels - flowFrame Method
#'
#' \code{cyto_channel_check} will check whether the supplied channels are valid
#' for the supplied \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#'
#' @param x an object of class
#'   \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param channels vector of channel names (e.g. c("PE-A","APC-A")) or marker
#'   names.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to TRUE by default to accept 1 or 2 channels.
#'
#' @return vector of valid channel names.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Add samples to flowSet
#' fs <- Activation
#' 
#' # Assign marker names
#' chnls <- c(
#'   "PE-A",
#'   "Alexa Fluor 488-A",
#'   "Alexa Fluor 700-A",
#'   "Alexa Fluor 647-A",
#'   "7-AAD-A"
#' )
#' markers <- c("Va2", "CD8", "CD4", "CD44", "CD69")
#' names(markers) <- chnls
#' markernames(fs) <- markers
#' 
#' cyto_channel_check(fs[[1]],
#'   channels = c("CD4", "CD8"),
#'   plot = TRUE
#' )
#' 
#' cyto_channel_check(fs[[1]],
#'   channels = c("CD4", "CD8", "CD44", "CD69"),
#'   plot = FALSE
#' )
#' @export
setMethod(cyto_channel_check,
  signature = "flowFrame",
  definition = function(x,
                          channels,
                          plot = TRUE) {

    # Incorrect channels length
    if (plot == TRUE) {
      if (!length(channels) %in% c(1, 2)) {
        stop("Invalid number of supplied channels.")
      }
    }

    chans <- BiocGenerics::colnames(x)
    fr.data <- pData(parameters(x))

    # Channel Indices supplied
    if (is.numeric(channels)) {
      channels <- chans[channels]
    }

    # Check if channels match colnames of flowFrame
    if (all(channels %in% chans)) {

      # Supplied channels are valid
    } else if (!all(channels %in% chans)) {
      lapply(channels, function(channel) {
        if (channel %in% chans) {


        } else if (channel %in% fr.data$desc) {
          channels[channels %in% channel] <<- as.character(
            fr.data$name[match(channel, fr.data$desc)]
          )
        } else if (!channel %in% chans & !channel %in% fr.data$desc) {
          stop(paste(channel, "is not a valid channel for this flowFrame."))
        }
      })
    }

    return(channels)
  }
)

#' Check Supplied Channels - flowSet Method
#'
#' \code{cyto_channel_check} will check whether the supplied channels are valid
#' for the supplied \code{\link[flowCore:flowSet-class]{flowSet}}.
#'
#' @param x an object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param channels vector of channel names (e.g. c("PE-A","APC-A")) or marker
#'   names.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to TRUE by default to accept 1 or 2 channels.
#'
#' @return vector of valid channel names.
#'
#' @importFrom flowWorkspace pData
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Add samples to flowSet
#' fs <- Activation
#' 
#' # Assign marker names
#' chnls <- c(
#'   "PE-A",
#'   "Alexa Fluor 488-A",
#'   "Alexa Fluor 700-A",
#'   "Alexa Fluor 647-A",
#'   "7-AAD-A"
#' )
#' markers <- c("Va2", "CD8", "CD4", "CD44", "CD69")
#' names(markers) <- chnls
#' markernames(fs) <- markers
#' 
#' cyto_channel_check(fs,
#'   channels = c("CD4", "CD8"),
#'   plot = TRUE
#' )
#' cyto_channel_check(fs,
#'   channels = c("CD4", "CD8", "CD44", "CD69"),
#'   plot = FALSE
#' )
#' @export
setMethod(cyto_channel_check,
  signature = "flowSet",
  definition = function(x,
                          channels,
                          plot = TRUE) {

    # Incorrect channels length
    if (plot == TRUE) {
      if (!length(channels) %in% c(1, 2)) {
        stop("Invalid number of supplied channels.")
      }
    }

    # Assign x to fs
    fs <- x

    chans <- BiocGenerics::colnames(fs[[1]])
    fr.data <- pData(parameters(fs[[1]]))

    # Channel Indices supplied
    if (is.numeric(channels)) {
      channels <- chans[channels]
    }

    # Check if channels match colnames of flowFrame
    if (all(channels %in% chans)) {

      # Supplied channels are valid
    } else if (!all(channels %in% chans)) {
      lapply(channels, function(channel) {
        if (channel %in% chans) {


        } else if (channel %in% fr.data$desc) {
          channels[channels %in% channel] <<- as.character(
            fr.data$name[match(channel, fr.data$desc)]
          )
        } else if (!channel %in% chans & !channel %in% fr.data$desc) {
          stop(paste(channel, "is not a valid channel for this flowFrame."))
        }
      })
    }

    return(channels)
  }
)

#' Check Supplied Channels - GatingSet Method
#'
#' \code{cyto_channel_check} will check whether the supplied channels are valid
#' for the supplied \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#'
#' @param x an object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel names (e.g. c("PE-A","APC-A")) or marker
#'   names.
#' @param plot logical indicating whether the channels will be used to construct
#'   a plot, set to TRUE by default to accept 1 or 2 channels.
#'
#' @return vector of valid channel names.
#'
#' @importFrom flowWorkspace pData getData
#' @importFrom flowCore parameters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Add samples to flowSet
#' fs <- Activation
#' 
#' # Assign marker names
#' chnls <- c(
#'   "PE-A",
#'   "Alexa Fluor 488-A",
#'   "Alexa Fluor 700-A",
#'   "Alexa Fluor 647-A",
#'   "7-AAD-A"
#' )
#' markers <- c("Va2", "CD8", "CD4", "CD44", "CD69")
#' names(markers) <- chnls
#' markernames(fs) <- markers
#' 
#' # fs to GatingSet
#' gs <- GatingSet(fs)
#' 
#' cyto_channel_check(gs,
#'   channels = c("CD4", "CD8"),
#'   plot = TRUE
#' )
#' 
#' cyto_channel_check(gs,
#'   channels = c("CD4", "CD8", "CD44", "CD69"),
#'   plot = FALSE
#' )
#' @export
setMethod(cyto_channel_check,
  signature = "GatingSet",
  definition = function(x,
                          channels,
                          plot = TRUE) {

    # Incorrect channels length
    if (plot == TRUE) {
      if (!length(channels) %in% c(1, 2)) {
        stop("Invalid number of supplied channels.")
      }
    }

    # Assign x to gs
    gs <- x
    fs <- getData(gs, "root")

    chans <- BiocGenerics::colnames(fs[[1]])
    fr.data <- pData(parameters(fs[[1]]))

    # Channel Indices supplied
    if (is.numeric(channels)) {
      channels <- chans[channels]
    }

    # Check if channels match colnames of flowFrame
    if (all(channels %in% chans)) {

      # Supplied channels are valid
    } else if (!all(channels %in% chans)) {
      lapply(channels, function(channel) {
        if (channel %in% chans) {


        } else if (channel %in% fr.data$desc) {
          channels[channels %in% channel] <<- as.character(
            fr.data$name[match(channel, fr.data$desc)]
          )
        } else if (!channel %in% chans & !channel %in% fr.data$desc) {
          stop(paste(channel, "is not a valid channel for this flowFrame."))
        }
      })
    }

    return(channels)
  }
)

#' Check Gate Type(s) Supplied to gate_draw.
#'
#' @param type vector indicating the types of gates to construct using
#'   \code{gate_draw}.
#' @param alias names of the populations to be gated.
#'
#' @return Stop gating process if type is incorrect or returns \code{type} as
#'   full lower case name(s). If a single type is supplied for multiple
#'   populations, the same type will be used for all populations.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_gate_type_check(type = "r", alias = c("A", "B", "C"))
#' @noRd
.cyto_gate_type_check <- function(type, alias) {
  if (all(type %in% c("q", "Q", "quadrant", "Quadrant")) & length(alias) != 4) {
    stop("Supply the names of 4 poulations to alias for quadrant gates.")
  }

  gts <- c(
    "polygon",
    "Polygon",
    "p",
    "P",
    "rectangle",
    "Rectangle",
    "r",
    "R",
    "interval",
    "Interval",
    "i",
    "I",
    "threshold",
    "Threshold",
    "t",
    "T",
    "boundary",
    "Boundary",
    "b",
    "B",
    "ellipse",
    "Ellipse",
    "e",
    "E",
    "quadrant",
    "Quadrant",
    "q",
    "Q",
    "web",
    "Web",
    "w",
    "W"
  )

  if (!all(type %in% gts)) {
    if (length(type[type %in% gts == FALSE]) >= 2) {
      stop(
        paste(
          paste(type[type %in% gts == FALSE], collapse = " & "),
          "are not valid gate types for gate_draw!"
        )
      )
    } else {
      stop(paste(
        type[type %in% gts == FALSE],
        "is not a valid gate type for gate_draw!"
      ))
    }
  }

  type[type %in% c("polygon", "Polygon", "p", "P")] <- "polygon"

  type[type %in% c("rectangle", "Rectangle", "r", "R")] <- "rectangle"

  type[type %in% c("interval", "Interval", "i", "I")] <- "interval"

  type[type %in% c("threshold", "Threshold", "t", "T")] <- "threshold"

  type[type %in% c("boundary", "Boundary", "b", "B")] <- "boundary"

  type[type %in% c("ellipse", "Ellipse", "e", "E")] <- "ellipse"

  type[type %in% c("quadrant", "Quadrant", "q", "Q")] <- "quadrant"

  type[type %in% c("web", "Web", "w", "W")] <- "web"

  # Repeat type to equal length of alias
  if (length(type) != length(alias) &
    type[1] != "quadrant" &
    type[1] != "web") {
    type <- rep(type, length(alias))
  }

  return(type)
}

#' Check Alias Supplied to gate_draw
#'
#' @param alias vector indicating the names of the populations to be gated.
#' @param type vector indicating the type(s) of gate(s) to be constructed.
#'
#' @return Stops the gating process if alias is missing or of the incorrect
#'   length given the gate type.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .cyto_alias_check(alias = c("A", "B", "C", "D"), type = "q")
#' @noRd
.cyto_alias_check <- function(alias = NULL, type) {
  if (is.null(alias)) {
    stop("Supply names of populations to 'alias' for checking.")
  }

  if (type[1] == "quadrant" & length(alias) != 4) {
    stop("Supply 4 population names to 'alias' to construct quadrant gates.")
  }

  if (length(type) > 1) {
    if (length(alias) != length(type)) {
      stop("Length of alias must be the same length as type for multi-gates.")
    }
  }
}

#' Check Operating System & Open New Graphics Device
#'
#' \code{checkOSGD} is used internally by cyto_plot to open an OS-specific
#' interactive garphics device to facilitate gate drawing. Mac users will need
#' to install \href{https://www.xquartz.org/}{XQuartz} for this functionality.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @import grDevices
#'
#' @examples
#' \dontrun{
#' # Open platform-specific graphics device
#' .cyto_plot_window()
#' }
#' 
#' @noRd
.cyto_plot_window <- function() {
  if (.Platform$OS.type == "windows") {
    grDevices::windows()
  } else if (.Platform$OS.type == "unix") {
    if (Sys.info()["sysname"] == "Linux") {
      X11()
    } else if (Sys.info()["sysname"] == "Darwin") {
      quartz()
    }
  }
}

#' Check File Exists in Working Directory
#'
#' @param name filename including file extension to be checked.
#'
#' @return TRUE/FALSE if file exists in the current working directory.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' .file_wd_check("gatingTemplate.csv")
#' @noRd
.file_wd_check <- function(name) {
  if (length(which(list.files() == name)) != 0) {

    # File exists in working directory
    return(TRUE)
  } else if (length(which(list.files() == name)) == 0) {

    # File does not exist in working directory
    return(FALSE)
  }
}

#' Check gatingTemplate for Existing Entry
#'
#' @param parent name of the parent population.
#' @param alias name of the population of interest.
#' @param gatingTemplate csv file name of the gatingTemplate.
#'
#' @return stops the gating process if an entry already exists in the
#'   gatingTemplate for the supplied alias.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom utils read.csv
#'
#' @noRd
.cyto_gatingTemplate_check <- function(parent, alias, gatingTemplate) {
  if (inherits(gatingTemplate, "gatingTemplate")) {
    stop("'gatingTemplate' should the name of the gatingTemplate csv file.")
  } else {
    if (getOption("CytoRSuite_wd_check") == TRUE) {
      if (.file_wd_check(gatingTemplate)) {
        gt <- read.csv(gatingTemplate, header = TRUE)

        # Parent and alias entries match file
        if (any(gt$parent %in% parent & gt$alias %in% alias)) {
          message(
            paste(
              paste(gt$alias, collapse = " & "),
              "already exists in", gatingTemplate, "."
            )
          )
          stop("Supply another gatingTemplate or edit gate(s) using gate_edit.")
        }
      }
    }
  }
}

#' Check Transformation List Object
#'
#' @param trans object of class
#'   \code{\link[flowCore:transformList-class]{transformList}} or
#'   \code{\link[flowWorkspace]{transformerList}} generated by
#'   \code{\link[flowCore:logicleTransform]{estimateLogicle}} which was used to
#'   transform the fluorescent channels of the samples.
#' @param inverse logical indicating whether the returned transformList should
#'   contain the inverse transformations.
#'
#' @return NULL or object of class
#'   \code{\link[flowCore:transformList-class]{transformList}}
#'
#' @importFrom flowCore inverseLogicleTransform transformList
#' @importFrom methods is
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' library(CytoRSuiteData)
#' 
#' # Load in samples to flowSet
#' fs <- Activation
#' 
#' # Add fs to GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Get transformList using estimateLogicle
#' trans <- estimateLogicle(fs[[4]], cyto_fluor_channels(fs))
#' 
#' # Get transformList containing inverse transformarions
#' inv <- cyto_trans_check(trans, inverse = TRUE)
#' 
#' # Convert transformerList into transformList
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(gs))
#' 
#' # Convert transformerList into inverse transformList
#' inv <- cyto_trans_check(trans, inverse = FALSE)
#' @export
cyto_trans_check <- function(trans = NULL, inverse = FALSE) {
  if (is.null(trans)) {
    return(NULL)
  } else {
    if (!class(trans)[1] %in% c("transformList", "transformerList")) {
      stop("'trans' should be of class transformList or transformerList.")
    } else {
      if (is(trans, "transformList")) {
        if (inverse) {
          trans <- inverseLogicleTransform(trans)
          return(trans)
        } else {
          return(trans)
        }
      } else if (is(trans, "transformerList")) {
        if (inverse) {
          trans <- transformList(names(trans), lapply(trans, `[[`, "inverse"))
          return(trans)
        } else {
          trans <- transformList(names(trans), lapply(trans, `[[`, "transform"))
          return(trans)
        }
      }
    }
  }
}

#' Check Overlays Supplied to cyto_plot
#'
#' \code{.cyto_overlay_check} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}.
#'
#' @param x object of class \code{flowFrame} or \code{flowSet}.
#' @param overlay object to overlay.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
setGeneric(
  name = ".cyto_overlay_check",
  def = function(x, ...) {
    standardGeneric(".cyto_overlay_check")
  }
)

#' Check Overlays Supplied to cyto_plot
#'
#' \code{.cyto_overlay_check} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}. This flowFrame method will return a list of
#' flowFrames to overlay.
#'
#' @param x object of class \code{flowFrame}.
#' @param overlay object to overlay.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore sampleFilter fsApply Subset
#'
#' @noRd
setMethod(.cyto_overlay_check,
  signature = "flowFrame",
  definition = function(x,
                          overlay,
                          display = NULL) {

    # Assign x to fr
    fr <- x

    # Check overlay class
    if (class(overlay) == "flowFrame") {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        n <- BiocGenerics::nrow(overlay)
        overlay <- Subset(
          overlay,
          sampleFilter(size = display * n)
        )
      }
      overlay <- list(overlay)
    } else if (class(overlay) == "flowSet") {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- fsApply(overlay, function(x) {
          Subset(x, sampleFilter(size = display * BiocGenerics::nrow(x)))
        })
      }

      overlay <- lapply(seq(1, length(overlay), 1), function(x) overlay[[x]])
    } else if (all(unlist(lapply(overlay, class)) == "flowFrame")) {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- lapply(overlay, function(x) {
          Subset(x, sampleFilter(size = display * BiocGenerics::nrow(x)))
        })
      }

      overlay <- overlay
    } else if (all(unlist(lapply(overlay, class)) == "flowSet") &
      length(overlay) == 1) {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- lapply(overlay, function(x) {
          fsApply(x, function(y) {
            Subset(y, sampleFilter(size = display * BiocGenerics::nrow(y)))
          })
        })
      }

      overlay <- lapply(overlay, function(x) {
        lapply(seq(1, length(x), 1), function(y) x[[y]])
      })[[1]]
    } else {
      stop(paste("'overlay' must be a flowFrame, flowSet,",
         "list of flowFrames or a list containing a flowSet.",
         sep = " "))
    }

    # return is a list of flowFrames to overlay
    return(overlay)
  }
)

#' Check Overlays Supplied to cyto_plot
#'
#' \code{.cyto_overlay_check} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}. This flowSet method will return a list of flowFrame
#' lists to overlay.
#'
#' @param x object of class \code{flowSet}.
#' @param overlay object to overlay.
#' @param display  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowCore sampleFilter Subset fsApply
#'
#' @noRd
setMethod(.cyto_overlay_check,
  signature = "flowSet",
  definition = function(x,
                          overlay,
                          display = NULL) {

    # Assign x to fs
    fs <- x

    # Check overlay class
    if (class(overlay) == "flowFrame") {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        n <- BiocGenerics::nrow(overlay)
        overlay <- Subset(overlay, sampleFilter(size = display * n))
      }

      overlay <- lapply(rep(list(overlay), length(fs)), "list")
    } else if (class(overlay) == "flowSet") {
      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- fsApply(overlay, function(x) {
          Subset(x, sampleFilter(size = display * BiocGenerics::nrow(x)))
        })
      }

      overlay <- lapply(lapply(
        seq(1, length(overlay), 1),
        function(x) overlay[[x]]
      ), "list")
    } else if (all(unlist(lapply(overlay, class)) == "flowFrame")) {
      if (length(overlay) == 1) {
        overlay <- lapply(rep(list(overlay[[1]]), length(fs)), "list")
      } else {
        if (length(overlay) != length(fs)) {
          stop(paste("Supplied list of flowFrames must be of the", 
             "same length as the flowSet.", sep = " "))
        }
        overlay <- lapply(overlay, "list")
      }

      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- lapply(overlay, function(x) {
          Subset(x, sampleFilter(size = display * BiocGenerics::nrow(x)))
        })
      }
    } else if (all(unlist(lapply(overlay, class)) == "flowSet")) {
      if (!all(unlist(lapply(overlay, length)) == length(fs))) {
        stop(paste("Each flowSet in supplied list should be of the", 
           "same length as the supplied flowSet.", sep = " "))
      }

      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- lapply(overlay, function(x) {
          fsApply(x, function(y) {
            Subset(y, sampleFilter(size = display * BiocGenerics::nrow(y)))
          })
        })
      }

      # list of flowFrame lists
      overlay <- lapply(overlay, function(x) {
        lapply(seq(1, length(x), 1), function(y) x[[y]])
      })
      overlay <- lapply(seq_along(fs), function(x) {
        lapply(overlay, `[[`, x)
      })
    } else if (all(do.call("rbind", lapply(overlay, function(x) {
      lapply(x, class)
    })) == "flowFrame")) {
      if (length(overlay) != length(fs)) {
        stop(paste("'overlay' should be a list of flowFrames lists to overlay", 
           "on each flowFrame in the flowSet.", sep = " "))
      }

      if (!is.null(display) & getOption("CytoRSuite_overlay_display")) {
        overlay <- lapply(overlay, function(x) {
          lapply(x, function(y) {
            Subset(y, sampleFilter(size = display * BiocGenerics::nrow(y)))
          })
        })
      }
    } else {
      stop(paste("'overlay' must be a flowFrame, flowSet,",
        "list of flowFrames or a list of flowSets.", sep = " "))
    }

    # return is a list of flowFrame lists to overlay
    # 1 flowFrame list per flowFrame in fs
    return(overlay)
  }
)

#' Check Overlays Supplied to cyto_plot
#'
#' \code{.cyto_overlay_check} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}. This flowSet method will return a list of flowFrame
#' lists to overlay.
#'
#' @param x object of class \code{GatingHierarchy}.
#' @param overlay object to overlay.
#' @param display  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace getData
#'
#' @noRd
setMethod(.cyto_overlay_check,
  signature = "GatingHierarchy",
  definition = function(x,
                          overlay,
                          display = NULL) {

    # Assign x to gh
    gh <- x

    # Extract flowFrame
    fr <- getData(gh, "root")

    # Overlay should be character vector of population names
    if (inherits(overlay, "flowFrame") |
      inherits(overlay, "flowSet") |
      inherits(overlay, "list")) {

    } else if (inherits(overlay, "character")) {
      if (!all(overlay %in% basename(getNodes(gh)))) {
        stop("'overlay' does not exist in the GatingHierarchy.")
      } else {
        nms <- overlay
        overlay <- lapply(overlay, function(x) {
          getData(gh, x)
        })
        names(overlay) <- nms
      }
    }

    # .cyto_overlay_check to convert overlay to correct format
    .cyto_overlay_check(fr,
      overlay = overlay,
      display = display
    )
  }
)

#' Check Overlays Supplied to cyto_plot
#'
#' \code{.cyto_overlay_check} will check whether the supplied overlay is
#' supported and convert it into an appropriate format for use in
#' \code{\link{cyto_plot}}. This flowSet method will return a list of flowFrame
#' lists to overlay.
#'
#' @param x object of class \code{GatingSet}.
#' @param overlay object to overlay.
#' @param display  numeric indicating the number of events to plot, set to all
#'   events by default. Reducing the sample size can significantly increase
#'   plotting speed on less powerful machines.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @importFrom flowWorkspace getData
#'
#' @noRd
setMethod(.cyto_overlay_check,
  signature = "GatingSet",
  definition = function(x,
                          overlay,
                          display = NULL) {

    # Assign x to gh
    gs <- x

    # Extract flowFrame
    fs <- getData(gs, "root")

    # Overlay should be character vector of population names
    if (inherits(overlay, "flowFrame") |
      inherits(overlay, "flowSet") |
      inherits(overlay, "list")) {

    } else if (inherits(overlay, "character")) {
      if (!all(overlay %in% basename(getNodes(gs)))) {
        stop("overlay' does not exist in the GatingHierarchy.")
      } else {
        nms <- overlay
        overlay <- lapply(overlay, function(x) {
          getData(gs, x)
        })
        names(overlay) <- nms
      }
    }

    # .cyto_overlay_check to convert overlay to correct format
    .cyto_overlay_check(fs,
      overlay = overlay,
      display = display
    )
  }
)

#' Check Statistic for ComputeStats
#'
#' @param stat cyto_stats_compute statistic.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_stat_check <- function(stat) {
  if (!stat %in% c(
    "mean",
    "Mean",
    "median",
    "Median",
    "mode",
    "Mode",
    "count",
    "Count",
    "percent",
    "Percent",
    "freq",
    "Freq",
    "geo mean",
    "Geo mean",
    "Geo Mean",
    "CV",
    "cv"
  )) {
    stop("Supplied statistic not supported.")
  }

  if (stat %in% c("mean", "Mean")) {
    stat <- "mean"
  } else if (stat %in% c("median", "Median")) {
    stat <- "median"
  } else if (stat %in% c("mode", "Mode")) {
    stat <- "mode"
  } else if (stat %in% c("count", "Count")) {
    stat <- "count"
  } else if (stat %in% c("percent", "Percent", "freq", "Freq")) {
    stat <- "freq"
  } else if (stat %in% c("geo mean", "Geo mean", "Geo Mean")) {
    stat <- "geo mean"
  } else if (stat %in% c("cv", "CV")) {
    stat <- "CV"
  }

  return(stat)
}

#' Check gate object supplied to cyto_plot
#'
#' @param gate gate object(s) to add to plot.
#' @param smp number of samples
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @noRd
.cyto_gate_check <- function(gate, smp) {

  # A single gate supplied
  if (inherits(gate, "rectangleGate") |
    inherits(gate, "polygonGate") |
    inherits(gate, "ellipsoidGate")) {
    gate <- rep(list(list(gate)), smp)

    # A filters object supplied
  } else if (inherits(gate, "filters")) {
    gate <- rep(list(filters), smp)

    # A list of gate objects
  } else if (inherits(gate, "list")) {

    # List of individual gates
    if (all(unlist(lapply(gate, class)) %in%
      c("rectangleGate", "polygonGate", "ellipsoidGate"))) {

      # Assume 1 gate per sample
      if (length(gate) == smp) {
        gate <- lapply(gate, "list")

        # Assume list contains gates for each sample
      } else {
        gate <- rep(list(gate), smp)
      }

      # List of filters objects
    } else if (all(unlist(lapply(gate, class)) == "filters")) {

      # Assume 1 filters object per sample
      gate <- rep(gate, length.out = smp)
    }

    # filtersList
  } else if (inherits(gate, "filtersList")) {
    gate <- gate
  }

  return(gate)
}
