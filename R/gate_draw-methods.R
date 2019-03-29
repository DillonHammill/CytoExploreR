#' gate_draw
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' \code{gate_draw} is a convenient wrapper for the gating functions shipped
#' with \code{cytoRSuite} to facilitate analysis of flow cytometry by gate
#' drawing. Using \code{gate_draw} users can specify the type of gate(s) to be
#' constructed through the \code{type} argument and \code{gate_draw} will
#' automatically handle plotting the data and make calls to the relevant gating
#' function(s) to construct the gates around populations of interest.
#' \code{gate_draw} has methods for
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}} and
#' \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} objects, refer to
#' their respective help pages for more information. The flowFrame and flowSet
#' methods simply return the constructed gates as a list of
#' \code{\link[flowCore:filters-class]{filters}}, whilst the GatingSet method
#' automatically applies the constructed gates to the GatingSet and saves the
#' constructed gates in an \code{openCyto}
#' \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}for future use.
#' See \code{\link{gate_edit}} and \code{\link{gate_remove}} to manipulate
#' constructed gates and modify their entries in the gatingTemplate.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}},
#'   \code{\link[flowCore:flowSet-class]{flowSet}} or
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param ... additional method-specific arguments.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{gate_draw,flowFrame-method}}
#' @seealso \code{\link{gate_draw,flowSet-method}}
#' @seealso \code{\link{gate_draw,GatingSet-method}}
#'
#' @export
setGeneric(
  name = "gate_draw",
  def = function(x, ...) {
    standardGeneric("gate_draw")
  }
)

#' gate_draw flowFrame Method.
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types include \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param density_smooth smoothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowFrame-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} list containing the
#'   drawn gate objects.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,flowFrame-method}}
#' @seealso \code{\link{gate_draw,flowSet-method}}
#' @seealso \code{\link{gate_draw,GatingSet-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # draw gates using gate_draw - add contour lines & overlay control
#' gt <- gate_draw(fs[[4]],
#'   channels = c("FSC-A", "SSC-A"),
#'   alias = "Cells",
#'   type = "ellipse",
#'   contour_lines = 15,
#'   overlay = fs[[1]]
#' )
#' 
#' # gt is a filters object containing the contructed ellipsoidGate
#' gt[[1]]
#' }
#' 
#' @export
setMethod(gate_draw,
          signature = "flowFrame",
          definition = function(x,
                                alias = NULL,
                                channels = NULL,
                                type = "polygon",
                                display = 1 / length(x),
                                axis = "x",
                                density_smooth = 1.5,
                                label = TRUE,
                                plot = TRUE, ...) {
            
            # Turn off sampling of overlay
            options("CytoRSuite_overlay_display" = FALSE)
            
            # Assign x to fr
            fr <- x
            
            # Check type argument is valid
            type <- .cyto_gate_type_check(type = type, alias = alias)
            
            # Set default type for 1D gates to interval
            if (length(channels) == 1 & all(type %in% "polygon")) {
              type <- rep("interval", length(type))
            }
            
            # Check alias is supplied correctly
            .cyto_alias_check(alias = alias, type = type)
            
            # Check supplied channel(s) are valid
            channels <- cyto_channel_check(fr,
                                           channels = channels,
                                           plot = TRUE
            )
            
            # Make one call to drawPlot
            if (plot == TRUE) {
              if (getOption("CytoRSuite_interact") == FALSE) {
                if (length(channels) == 2) {
                  cyto_plot(fr,
                            channels = channels,
                            display = display,
                            popup = FALSE,
                            legend = FALSE, ...
                  )
                } else {
                  cyto_plot(fr,
                            channels = channels,
                            popup = FALSE,
                            legend = FALSE, ...
                  )
                }
              } else {
                if (length(channels) == 2) {
                  cyto_plot(fr,
                            channels = channels,
                            display = display,
                            popup = TRUE,
                            legend = FALSE, ...
                  )
                } else {
                  cyto_plot(fr,
                            channels = channels,
                            popup = TRUE,
                            legend = FALSE, ...
                  )
                }
              }
            }
            
            # Reset sampling of overlay
            options("CytoRSuite_overlay_display" = TRUE)
            
            # Construct gates save as filters object
            if (length(type) == 1 & type[1] == "quadrant") {
              gates <- gate_quadrant_draw(
                fr = fr,
                channels = channels,
                alias = alias,
                display = display,
                plot = FALSE,
                label = label, ...
              )
            } else if (length(type) == 1 & type[1] == "web") {
              gates <- gate_web_draw(
                fr = fr,
                channels = channels,
                alias = alias,
                display = display,
                plot = FALSE,
                label = label, ...
              )
            } else {
              gates <- mapply(function(type, alias) {
                if (type == "polygon") {
                  gate_polygon_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "rectangle") {
                  gate_rectangle_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "interval") {
                  gate_interval_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    axis = axis,
                    label = label, ...
                  )
                } else if (type == "threshold") {
                  gate_threshold_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "boundary") {
                  gate_boundary_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "ellipse") {
                  gate_ellipse_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    abel = label, ...
                  )
                }
              }, type, alias)
            }
            
            gates <- filters(gates)
            return(gates)
          }
)

#' gate_draw flowSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param select vector containing the indices of samples within gs to use for
#'   plotting. For large \code{flowSet} objects \code{select} is set to 20
#'   random \code{flowFrame} objects to improve processing speed.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param density_smooth smoothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,flowSet-method}}.
#'
#' @return a \code{\link[flowCore:filters-class]{filters}} list containing the
#'   drawn gate objects.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore filters
#' @importFrom methods as
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,flowSet-method}}
#' @seealso \code{\link{gate_draw,flowFrame-method}}
#' @seealso \code{\link{gate_draw,GatingSet-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' 
#' # draw gates using gate_draw - add contour lines & overlay control
#' gate_draw(fs,
#'   channels = c("FSC-A", "SSC-A"),
#'   alias = "Cells",
#'   type = "polygon",
#'   contour_lines = 15,
#'   overlay = fs[[1]]
#' )
#' 
#' # gt is a filters object containing the contructed polygonGate
#' gt[[1]]
#' }
#' 
#' @export
setMethod(gate_draw,
          signature = "flowSet",
          definition = function(x,
                                alias = NULL,
                                channels = NULL,
                                type = "polygon",
                                display = 1 / length(x),
                                select = NULL,
                                axis = "x",
                                density_smooth = 1.5,
                                label = TRUE,
                                plot = TRUE, ...) {
            
            # Turn off sampling of overlay
            options("CytoRSuite_overlay_display" = FALSE)
            
            # Assign x to fs
            fs <- x
            
            # Restrict to samples based on indices
            if (!is.null(select)) {
              if (class(select) != "numeric") {
                stop(
                  "'select' must contain the numeric indices of the samples to plot."
                )
              }
              
              # Extract samples using selectFrames
              fs <- fs[select]
              
            # Large flowSets restricted to 20 random flowFrames for speed
            }else if(is.null(select) & length(fs) > 20){
              
              fs <- fs[sample(seq(1,length(fs)), 20)]
              
            }
            fr <- as(fs, "flowFrame")
            
            # Check type argument is valid
            type <- .cyto_gate_type_check(type = type, alias = alias)
            
            # Set default type for 1D gates to interval
            if (length(channels) == 1 & all(type %in% "polygon")) {
              type <- rep("interval", length(type))
            }
            
            # Check alias is supplied correctly
            .cyto_alias_check(alias = alias, type = type)
            
            # Check supplied channel(s) are valid
            channels <- cyto_channel_check(fr,
                                           channels = channels,
                                           plot = TRUE
            )
            
            # Make one call to drawPlot
            if (plot == TRUE) {
              if (getOption("CytoRSuite_interact") == FALSE) {
                if (length(channels) == 2) {
                  cyto_plot(fs,
                            channels = channels,
                            display = display,
                            popup = FALSE,
                            legend = FALSE,
                            group_by = "all", ...
                  )
                } else {
                  cyto_plot(fs,
                            channels = channels,
                            popup = FALSE,
                            legend = FALSE,
                            group_by = "all", ...
                  )
                }
              } else {
                if (length(channels) == 2) {
                  cyto_plot(fs,
                            channels = channels,
                            display = display,
                            popup = TRUE,
                            legend = FALSE,
                            group_by = "all", ...
                  )
                } else {
                  cyto_plot(fs,
                            channels = channels,
                            popup = TRUE,
                            legend = FALSE,
                            group_by = "all", ...
                  )
                }
              }
            }
            
            # Reset sampling of overlay
            options("CytoRSuite_overlay_display" = TRUE)
            
            # Construct gates save as filters object
            if (length(type) == 1 & type[1] == "quadrant") {
              gates <- gate_quadrant_draw(
                fr = fr,
                channels = channels,
                alias = alias,
                display = display,
                plot = FALSE,
                label = label, ...
              )
            } else if (length(type) == 1 & type[1] == "web") {
              gates <- gate_web_draw(
                fr = fr,
                channels = channels,
                alias = alias,
                display = display,
                plot = FALSE,
                label = label, ...
              )
            } else {
              gates <- mapply(function(type, alias) {
                if (type == "polygon") {
                  gate_polygon_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "rectangle") {
                  gate_rectangle_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "interval") {
                  gate_interval_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    axis = axis,
                    label = label, ...
                  )
                } else if (type == "threshold") {
                  gate_threshold_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "boundary") {
                  gate_boundary_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                } else if (type == "ellipse") {
                  gate_ellipse_draw(
                    fr = fr,
                    channels = channels,
                    alias = alias,
                    plot = FALSE,
                    label = label, ...
                  )
                }
              }, type, alias)
            }
            
            gates <- filters(gates)
            return(gates)
          }
)

#' gate_draw GatingSet Method
#'
#' Manually draw gates around populations for analysis of flow cytometry data.
#'
#' @param x object of class
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}}.
#' @param channels vector of channel names to use for plotting, can be of length
#'   1 for 1-D density histogram or length 2 for 2-D scatter plot.
#' @param parent name of the \code{parent} population to extract for gating.
#' @param alias the name(s) of the populations to be gated. If multiple
#'   population names are supplied (e.g. \code{c("CD3,"CD4)}) multiple gates
#'   will be returned. \code{alias} is \code{NULL} by default which will halt
#'   the gating routine.
#' @param type vector of gate type names used to construct the gates. Multiple
#'   gate types are supported but should be accompanied with an \code{alias}
#'   argument of the same length (i.e. one \code{type} per \code{alias}).
#'   Supported gate types are \code{polygon, rectangle, ellipse, threshold,
#'   boundary, interval, quadrant and web} which can be abbreviated as upper or
#'   lower case first letters as well. Default \code{type} is \code{"interval"}
#'   for 1D gates and \code{"polygon"} for 2D gates.
#' @param gatingTemplate name of \code{gatingTemplate} csv file to be saved.
#' @param group_by vector of pData column names (e.g.
#'   c("Treatment","Concentration") indicating how the samples should be grouped
#'   prior to gating, set to the length of x by default to construct a single
#'   gate for all samples. If group_by is supplied a different gate will be
#'   constructed for each group.
#' @param display numeric [0,1] to control the percentage of events to be
#'   plotted. Specifying a value for \code{display} can substantial improve
#'   plotting speed for less powerful machines.
#' @param select vector containing the indices of samples within gs to use for
#'   plotting. For large \code{flowSet} objects \code{select} is set to 20
#'   random \code{flowFrame} objects to improve processing speed.
#' @param axis indicates whether the \code{"x"} or \code{"y"} axis should be
#'   gated for 2-D interval gates.
#' @param label logical indicating whether to include
#'   \code{\link{cyto_plot_label}} for the gated population(s), \code{TRUE} by
#'   default.
#' @param density_smooth smoothing factor passed to
#'   \code{\link[stats:density]{density}} for 1-D plots (defaults to 1.5).
#' @param plot logical indicating whether a plot should be drawn, set to
#'   \code{TRUE} by default.
#' @param ... additional arguments for \code{\link{cyto_plot,GatingSet-method}}.
#'
#' @return drawn gates are applied to the
#'   \code{\link[flowWorkspace:GatingSet-class]{GatingSet}} and saved to a
#'   \code{\link[openCyto:gatingTemplate-class]{gatingTemplate}}.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom flowWorkspace getData
#' @importFrom openCyto add_pop
#' @importFrom methods as
#' @importFrom utils read.csv write.csv
#' @importFrom flowCore filters
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @seealso \code{\link{cyto_plot,GatingSet-method}}
#' @seealso \code{\link{gate_draw,flowFrame-method}}
#' @seealso \code{\link{gate_draw,flowSet-method}}
#'
#' @examples
#' \dontrun{
#' library(CytoRSuiteData)
#' 
#' # Load in samples
#' fs <- Activation
#' gs <- GatingSet(fs)
#' 
#' # Apply compensation
#' gs <- compensate(gs, fs[[1]]@description$SPILL)
#' 
#' # Transform fluorescent channels
#' trans <- estimateLogicle(gs[[4]], cyto_fluor_channels(fs))
#' gs <- transform(gs, trans)
#' 
#' # Gate using gate_draw
#' gating(Activation_gatingTemplate, gs)
#' 
#' # draw gates using gate_draw - add contour lines & overlay control
#' gate_draw(gs,
#'   parent = "Dendritic Cells",
#'   channels = c("Alexa Fluor 488-A", "Alexa Fluor 700-A"),
#'   alias = c("CD8+ DC", "CD4+ DC"),
#'   gatingTemplate = "Example-gatingTemplate.csv",
#'   type = "rectangle",
#'   contour_lines = 15
#' )
#' 
#' # Constructed gate applied directly to GatingSet
#' getNodes(gs)
#' }
#' 
#' @export
setMethod(gate_draw,
          signature = "GatingSet",
          definition = function(x,
                                parent = "root",
                                alias = NULL,
                                channels = NULL,
                                type = "polygon",
                                gatingTemplate = NULL,
                                group_by = NULL,
                                display = NULL,
                                select = NULL,
                                axis = "x",
                                density_smooth = 1.5,
                                label = TRUE,
                                plot = TRUE, ...) {
            
            # Turn off sampling of overlay
            options("CytoRSuite_overlay_display" = FALSE)
            
            # Assign x to gs
            gs <- x
            smp <- length(gs)
            
            # Extract pData information
            pd <- pData(gs)
            
            # Check whether a gatingTemplate ready exists for this population
            if (!is.null(gatingTemplate)) {
              
              # Check whether gate already exists in gatingTemplate
              .cyto_gatingTemplate_check(parent, alias, gatingTemplate)
            }
            
            fs <- flowWorkspace::getData(x, parent)
            
            # grouping required
            if (is.null(group_by)) {
              group_by <- NA
              grps <- list(fs)
              pd$group_by <- rep(1, smp)
            } else if (all(is.na(group_by))) {
              group_by <- NA
              grps <- list(fs)
              pd$group_by <- rep(1, smp)
            } else if (all(is.character(group_by))) {
              if (!all(group_by %in% colnames(pd))) {
                stop("Names supplied to group_by do not exist in pData(x).")
              }
              
              pd$group_by <- do.call(paste, pd[, group_by, drop = FALSE])
              grps <- lapply(unique(pd$group_by), function(x) {
                fs[which(pd$group_by %in% x)]
              })
            }
            
            # Gate each group - list of filters
            fltrsLst <- lapply(grps, function(fs) {
              
              # Restrict to samples matching pData requirements
              if (!is.null(select)) {
                if (class(select) != "numeric") {
                  stop("'select' contain the numeric indices of the samples to plot.")
                }
                
                # Extract samples using selectFrames
                fs <- fs[select]
                
              # Large flowSets restricted to 20 random flowFrames for speed
              }else if(is.null(select) & length(fs) > 20){
                
                fs <- fs[sample(seq(1,length(fs)), 20)]
                
              }
              fr <- as(fs, "flowFrame")
              
              # Events to display
              if (is.null(display)) {
                display <- 1 / length(fs)
              }
              
              # Remove "Original" column introduced by coercion
              if (is.na(match("Original", BiocGenerics::colnames(fr))) == FALSE) {
                fr <- suppressWarnings(
                  fr[, -match("Original", BiocGenerics::colnames(fr))]
                )
              }
              
              # Check type argument is valid
              type <- .cyto_gate_type_check(type = type, alias = alias)
              
              # Set default type for 1D gates to interval
              if (length(channels) == 1 & all(type %in% "polygon")) {
                type <- rep("interval", length(type))
              }
              
              # Check alias is supplied correctly
              .cyto_alias_check(alias = alias, type = type)
              
              # Check supplied channel(s) are valid
              channels <- cyto_channel_check(fr,
                                             channels = channels,
                                             plot = TRUE
              )
              
              # title
              if (all(is.na(group_by))) {
                if (parent == "root") {
                  pnt <- "All Events"
                } else {
                  pnt <- parent
                }
                
                title <- paste("Combined Events", "\n", pnt)
              } else if (all(is.character(group_by))) {
                if (parent == "root") {
                  pnt <- "All Events"
                } else {
                  pnt <- parent
                }
                
                title <- paste(pd[pd$name %in% 
                                    sampleNames(fs), "group_by"][1], "\n", pnt)
              }
              
              # Make one call to cyto_plot
              if (plot == TRUE) {
                if (getOption("CytoRSuite_interact") == FALSE) {
                  if (length(channels) == 2) {
                    cyto_plot(
                      x = gs[sampleNames(fs)],
                      parent = parent,
                      channels = channels,
                      display = display,
                      popup = FALSE,
                      legend = FALSE,
                      group_by = "all",
                      title = title, ...
                    )
                  } else {
                    cyto_plot(
                      x = gs[sampleNames(fs)],
                      parent = parent,
                      channels = channels,
                      popup = FALSE,
                      legend = FALSE,
                      group_by = "all",
                      title = title, ...
                    )
                  }
                } else {
                  if (length(channels) == 2) {
                    cyto_plot(
                      x = gs[sampleNames(fs)],
                      parent = parent,
                      channels = channels,
                      display = display,
                      popup = TRUE,
                      legend = FALSE,
                      group_by = "all",
                      title = title, ...
                    )
                  } else {
                    cyto_plot(
                      x = gs[sampleNames(fs)],
                      parent = parent,
                      channels = channels,
                      popup = TRUE,
                      legend = FALSE,
                      group_by = "all",
                      title = title, ...
                    )
                  }
                }
              }
              
              # Reset sampling of overlay
              options("CytoRSuite_overlay_display" = TRUE)
              
              # Construct gates save as filters object
              if (length(type) == 1 & type[1] == "quadrant") {
                gates <- gate_quadrant_draw(
                  fr = fr,
                  channels = channels,
                  alias = alias,
                  display = display,
                  plot = FALSE,
                  label = label, ...
                )
              } else if (length(type) == 1 & type[1] == "web") {
                gates <- gate_web_draw(
                  fr = fr,
                  channels = channels,
                  alias = alias,
                  display = display,
                  plot = FALSE,
                  label = label, ...
                )
              } else {
                gates <- mapply(function(type, alias) {
                  if (type == "polygon") {
                    gate_polygon_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      label = label, ...
                    )
                  } else if (type == "rectangle") {
                    gate_rectangle_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      label = label, ...
                    )
                  } else if (type == "interval") {
                    gate_interval_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      axis = axis,
                      label = label, ...
                    )
                  } else if (type == "threshold") {
                    gate_threshold_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      label = label, ...
                    )
                  } else if (type == "boundary") {
                    gate_boundary_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      label = label, ...
                    )
                  } else if (type == "ellipse") {
                    gate_ellipse_draw(
                      fr = fr,
                      channels = channels,
                      alias = alias,
                      plot = FALSE,
                      label = label, ...
                    )
                  }
                }, type, alias)
              }
              
              gates <- filters(gates)
            })
            
            # Name gates with group_by info
            if (all(is.na(group_by))) {
              
              # group number
              names(fltrsLst) <- unique(pd$group_by)
            } else if (all(is.character(group_by))) {
              
              # merge columns
              names(fltrsLst) <- unique(do.call(paste, pd[, group_by, drop = FALSE]))
            }
            gates <- fltrsLst
            
            # format gates to be a list of alias lists
            # each of length group and appropriately named
            gates <- lapply(seq_len(length(alias)), function(y) {
              gates <- lapply(gates, function(x) {
                gts <- filters(list(x[[y]]))
              })
              names(gates) <- unique(pd$group_by)
              return(gates)
            })
            
            # Prepare gatingTemplate entries
            pop <- "+"
            
            # Prepare group_by
            if (all(is.character(group_by))) {
              group_by <- paste(group_by, collapse = ":")
            } else if (all(is.na(group_by))) {
              group_by <- NA
            }
            
            # Use add_pop to apply gates to GatingSet and construct gatingTemplate
            if (is.null(gatingTemplate)) {
              message("Writing gatingTemplate.csv to store gates.")
              
              # need to extract alias from gates list into new named list
              pops <- list()
              for (i in seq_len(length(alias))) {
                pops[[i]] <- suppressWarnings(add_pop(
                  gs = x,
                  alias = alias[i],
                  parent = parent,
                  pop = pop,
                  dims = paste(channels, collapse = ","),
                  gating_method = "gate_draw",
                  gating_args = list(gate = gates[[i]]),
                  groupBy = group_by,
                  collapseDataForGating = TRUE,
                  preprocessing_method = "pp_gate_draw"
                ))
              }
              pops <- do.call("rbind", pops)
              
              write.csv(pops, "gatingTemplate.csv", row.names = FALSE)
            } else if (.file_wd_check(gatingTemplate) == FALSE) {
              message(
                paste(
                  gatingTemplate, "not in this working directory. Writing",
                  gatingTemplate
                )
              )
              
              pops <- list()
              for (i in seq_len(length(alias))) {
                pops[[i]] <- suppressWarnings(add_pop(
                  gs = x,
                  alias = alias[i],
                  parent = parent,
                  pop = pop,
                  dims = paste(channels, collapse = ","),
                  gating_method = "gate_draw",
                  gating_args = list(gate = gates[[i]]),
                  groupBy = group_by,
                  collapseDataForGating = TRUE,
                  preprocessing_method = "pp_gate_draw"
                ))
              }
              pops <- do.call("rbind", pops)
              
              write.csv(pops, gatingTemplate, row.names = FALSE)
            } else if (.file_wd_check(gatingTemplate) == TRUE) {
              gt <- read.csv(gatingTemplate, header = TRUE)
              
              pops <- list()
              for (i in seq_len(length(alias))) {
                pops[[i]] <- suppressWarnings(add_pop(
                  gs = x,
                  alias = alias[i],
                  parent = parent,
                  pop = pop,
                  dims = paste(channels, collapse = ","),
                  gating_method = "gate_draw",
                  gating_args = list(gate = gates[[i]]),
                  groupBy = group_by,
                  collapseDataForGating = TRUE,
                  preprocessing_method = "pp_gate_draw"
                ))
              }
              pops <- do.call("rbind", pops)
              gt <- rbind(gt, pops)
              
              write.csv(gt, gatingTemplate, row.names = FALSE)
            }
            
            invisible(pops)
          }
)