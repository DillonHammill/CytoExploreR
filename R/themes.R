## OZETTE THEMES ---------------------------------------------------------------

#' Dark theme to match Ozette's platform
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @export
ozette_dark_theme <- function(point_shape = "hex",
                              point_col_scale = c(
                                "#131D2A",
                                "#1853C9",
                                "#EFF412",
                                "#FF8000",
                                "#FF1200"
                              ),
                              grid_line_alpha = 0.1,
                              border_fill = "#131B24",
                              border_line_col = "grey95",
                              axes_text_col = "grey70",
                              axes_label_text_col = "white",
                              axes_ticks_line_col = "grey95",
                              title_text_col = "white",
                              key_text_col = "grey70",
                              key_title_text_col = "white",
                              key_ticks_line_col = "grey95",
                              key_border_line_col = "grey95",
                              page_fill = "#0B1117",
                              contour_line_col = "white",
                              header_text_col = "white",
                              legend_text_col = "white",
                              gate_line_col = "white",
                              gate_line_width = 1,
                              gate_fill_alpha = 0.1,
                              gate_fill = "white",
                              label_fill_alpha = 0,
                              label_text_col = "white",
                              label_text_font = 1,
                              ...) {
  
  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )
  
}

#' Light theme to match Ozette's platform
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
ozette_light_theme <- function(point_shape = "hex",
                               point_col_scale = c(
                                 "#131D2A",
                                 "#1853C9",
                                 "#EFF412",
                                 "#FF8000",
                                 "#FF1200"
                               ),
                               axes_text_col = "grey40",
                               key_text_col = "grey40",
                               gate_line_col = "black",
                               gate_line_width = 1,
                               gate_fill_alpha = 0.1,
                               gate_fill = "grey40",
                               label_fill_alpha = 0,
                               label_text_col = "black",
                               label_text_font = 1,
                               ...) {
  
  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )
  
}

#' Combination of light and dark themes for Ozette's slide decks
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#' 
#' @export
ozette_theme <- function(point_shape = "hex",
                         point_col_scale = c(
                           "#131D2A",
                           "#1853C9",
                           "#EFF412",
                           "#FF8000",
                           "#FF1200"
                         ),
                         grid_line_col = "grey40",
                         grid_line_alpha = 0.1,
                         border_fill = "white",
                         border_line_col = "grey95",
                         axes_text_col = "grey70",
                         axes_label_text_col = "white",
                         axes_ticks_line_col = "grey95",
                         title_text_col = "white",
                         key_text_col = "grey70",
                         key_title_text_col = "white",
                         key_ticks_line_col = "grey95",
                         key_border_line_col = "grey95",
                         page_fill = "#0B1117",
                         contour_line_col = "black",
                         header_text_col = "white",
                         legend_text_col = "white",
                         gate_line_col = "black",
                         gate_line_width = 1,
                         gate_fill_alpha = 0.1,
                         gate_fill = "grey40",
                         label_fill_alpha = 0,
                         label_text_col = "black",
                         label_text_font = 1,
                         ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

## GENERIC THEMES --------------------------------------------------------------

#' Minimal clean theme with reduced visual clutter
#'
#' A light, minimal theme that removes grid lines and uses subtle colors to
#' keep the focus on the data. Ideal for clean, distraction-free plots.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_minimal()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_minimal <- function(point_shape = "hex",
                                    point_col_scale = c(
                                      "grey90",
                                      "#4292C6",
                                      "#08519C",
                                      "#08306B"
                                    ),
                                    grid = FALSE,
                                    border_fill = "white",
                                    border_line_col = "grey80",
                                    border_line_width = 0.5,
                                    axes_text_col = "grey30",
                                    axes_text_size = 1,
                                    axes_label_text_col = "grey20",
                                    axes_ticks_line_col = "grey60",
                                    title_text_col = "grey10",
                                    title_text_font = 2,
                                    key_text_col = "grey30",
                                    key_title_text_col = "grey20",
                                    key_border_line_col = "grey80",
                                    key_ticks_line_col = "grey60",
                                    page_fill = "white",
                                    contour_line_col = "grey40",
                                    header_text_col = "grey10",
                                    legend_text_col = "grey20",
                                    gate_line_col = "grey20",
                                    gate_line_width = 0.8,
                                    gate_fill_alpha = 0.05,
                                    gate_fill = "grey50",
                                    label_fill_alpha = 0,
                                    label_text_col = "grey20",
                                    label_text_font = 1,
                                    ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

#' High contrast dark theme for presentations
#'
#' A bold dark theme with vibrant colors designed for maximum visibility on
#' projectors and large screens. Uses bright accent colors against a deep black
#' background.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_presentation()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_presentation <- function(point_shape = "hex",
                                         point_col_scale = c(
                                           "#0D0D0D",
                                           "#2166AC",
                                           "#66BD63",
                                           "#FEE08B",
                                           "#F46D43",
                                           "#D73027"
                                         ),
                                         grid = FALSE,
                                         border_fill = "#1A1A2E",
                                         border_line_col = "#E0E0E0",
                                         border_line_width = 1,
                                         axes_text_col = "#E0E0E0",
                                         axes_text_size = 1.2,
                                         axes_label_text_col = "white",
                                         axes_label_text_size = 1.3,
                                         axes_label_text_font = 2,
                                         axes_ticks_line_col = "#E0E0E0",
                                         title_text_col = "white",
                                         title_text_size = 1.5,
                                         title_text_font = 2,
                                         key_text_col = "#E0E0E0",
                                         key_text_size = 1.1,
                                         key_title_text_col = "white",
                                         key_ticks_line_col = "#E0E0E0",
                                         key_border_line_col = "#E0E0E0",
                                         page_fill = "#0D0D0D",
                                         contour_line_col = "white",
                                         contour_line_width = 1,
                                         header_text_col = "white",
                                         header_text_size = 1.3,
                                         legend_text_col = "white",
                                         legend_text_size = 1.1,
                                         gate_line_col = "#00E5FF",
                                         gate_line_width = 1.5,
                                         gate_fill_alpha = 0.1,
                                         gate_fill = "#00E5FF",
                                         label_fill_alpha = 0,
                                         label_text_col = "#00E5FF",
                                         label_text_size = 1.2,
                                         label_text_font = 2,
                                         ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

#' Publication-ready grayscale theme
#'
#' A grayscale theme optimized for journal publications and print. Uses only
#' black, white, and grey tones to ensure plots look good in both color and
#' black-and-white print.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_publication()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_publication <- function(point_shape = "hex",
                                        point_col_scale = c(
                                          "grey95",
                                          "grey70",
                                          "grey40",
                                          "grey10"
                                        ),
                                        grid = FALSE,
                                        border_fill = "white",
                                        border_line_col = "black",
                                        border_line_width = 1,
                                        axes_text_col = "black",
                                        axes_text_size = 1,
                                        axes_text_font = 1,
                                        axes_label_text_col = "black",
                                        axes_label_text_size = 1.1,
                                        axes_label_text_font = 2,
                                        axes_ticks_line_col = "black",
                                        axes_ticks_line_width = 1,
                                        title_text_col = "black",
                                        title_text_size = 1.2,
                                        title_text_font = 2,
                                        key_text_col = "black",
                                        key_title_text_col = "black",
                                        key_ticks_line_col = "black",
                                        key_border_line_col = "black",
                                        page_fill = "white",
                                        contour_line_col = "grey30",
                                        contour_line_width = 0.8,
                                        header_text_col = "black",
                                        legend_text_col = "black",
                                        gate_line_col = "black",
                                        gate_line_width = 1,
                                        gate_line_type = 2,
                                        gate_fill_alpha = 0,
                                        gate_fill = "black",
                                        label_fill_alpha = 0,
                                        label_text_col = "black",
                                        label_text_font = 2,
                                        label_text_size = 1,
                                        ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

#' Colorblind-friendly theme
#'
#' A theme using a colorblind-safe palette based on the Viridis color scale.
#' Ensures plots are accessible to viewers with color vision deficiencies while
#' maintaining clear visual distinction between density levels.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_colorblind()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_colorblind <- function(point_shape = "hex",
                                       point_col_scale = c(
                                         "#440154",
                                         "#31688E",
                                         "#35B779",
                                         "#FDE725"
                                       ),
                                       grid_line_col = "grey85",
                                       grid_line_alpha = 0.5,
                                       border_fill = "white",
                                       border_line_col = "grey40",
                                       border_line_width = 0.8,
                                       axes_text_col = "grey20",
                                       axes_label_text_col = "black",
                                       axes_ticks_line_col = "grey40",
                                       title_text_col = "black",
                                       title_text_font = 2,
                                       key_text_col = "grey20",
                                       key_title_text_col = "black",
                                       key_ticks_line_col = "grey40",
                                       key_border_line_col = "grey40",
                                       page_fill = "white",
                                       contour_line_col = "grey30",
                                       header_text_col = "black",
                                       legend_text_col = "black",
                                       gate_line_col = "#D55E00",
                                       gate_line_width = 1.2,
                                       gate_fill_alpha = 0.08,
                                       gate_fill = "#D55E00",
                                       label_fill_alpha = 0,
                                       label_text_col = "#D55E00",
                                       label_text_font = 2,
                                       ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

#' Soft pastel theme
#'
#' A gentle, pastel-colored theme with a warm off-white background. Provides a
#' softer appearance suitable for reports and dashboards where a less intense
#' visual style is preferred.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_pastel()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_pastel <- function(point_shape = "hex",
                                   point_col_scale = c(
                                     "#F7F7F7",
                                     "#92C5DE",
                                     "#F4A582",
                                     "#CA0020"
                                   ),
                                   grid_line_col = "#E8E0D8",
                                   grid_line_alpha = 0.6,
                                   border_fill = "#FAF8F5",
                                   border_line_col = "#C8B8A8",
                                   border_line_width = 0.5,
                                   axes_text_col = "#6B5B4F",
                                   axes_label_text_col = "#4A3728",
                                   axes_ticks_line_col = "#C8B8A8",
                                   title_text_col = "#3A2718",
                                   title_text_font = 2,
                                   key_text_col = "#6B5B4F",
                                   key_title_text_col = "#4A3728",
                                   key_ticks_line_col = "#C8B8A8",
                                   key_border_line_col = "#C8B8A8",
                                   page_fill = "#F5F0EB",
                                   contour_line_col = "#8B7B6B",
                                   header_text_col = "#3A2718",
                                   legend_text_col = "#4A3728",
                                   gate_line_col = "#7570B3",
                                   gate_line_width = 1,
                                   gate_fill_alpha = 0.08,
                                   gate_fill = "#7570B3",
                                   label_fill_alpha = 0,
                                   label_text_col = "#7570B3",
                                   label_text_font = 1,
                                   ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}

#' Nature-inspired blue-green theme
#'
#' A theme using ocean-inspired blue-green tones with a cool color palette.
#' Provides a professional, calming appearance with good density visualization
#' through a blue-to-green-to-yellow color scale.
#'
#' @param ... additional arguments to customize the \code{cyto_plot()} theme,
#'   see \code{cyto_plot_theme_args()} for details.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' cyto_plot_theme_ocean()
#' cyto_plot(gs[[1]], parent = "root", channels = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
cyto_plot_theme_ocean <- function(point_shape = "hex",
                                  point_col_scale = c(
                                    "#F7FBFF",
                                    "#6BAED6",
                                    "#2171B5",
                                    "#08306B"
                                  ),
                                  grid_line_col = "#D4E6F1",
                                  grid_line_alpha = 0.4,
                                  border_fill = "#EBF5FB",
                                  border_line_col = "#85C1E9",
                                  border_line_width = 0.8,
                                  axes_text_col = "#2C3E50",
                                  axes_label_text_col = "#1B2631",
                                  axes_ticks_line_col = "#85C1E9",
                                  title_text_col = "#1B2631",
                                  title_text_font = 2,
                                  key_text_col = "#2C3E50",
                                  key_title_text_col = "#1B2631",
                                  key_ticks_line_col = "#85C1E9",
                                  key_border_line_col = "#85C1E9",
                                  page_fill = "#D6EAF8",
                                  contour_line_col = "#2C3E50",
                                  header_text_col = "#1B2631",
                                  legend_text_col = "#1B2631",
                                  gate_line_col = "#E74C3C",
                                  gate_line_width = 1.2,
                                  gate_fill_alpha = 0.08,
                                  gate_fill = "#E74C3C",
                                  label_fill_alpha = 0,
                                  label_text_col = "#E74C3C",
                                  label_text_font = 2,
                                  ...) {

  cyto_plot_theme_reset()
  args <- .args_list(...)
  cyto_func_call(
    "cyto_plot_theme",
    args
  )

}
