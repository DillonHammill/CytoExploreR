## CUSTOM CYTO_PLOT() THEMES ---------------------------------------------------

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
