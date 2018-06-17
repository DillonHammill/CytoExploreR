#' Graphics Theme Used for editSpillover.
#' 
#' @export
editSpillover_theme <- function(){
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, color = "black", face = "bold"))
}