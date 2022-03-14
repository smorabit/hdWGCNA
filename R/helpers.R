

#' scale01
#'
#' Function to scale a numeric vector between 0 and 1.
#'
#' @param x numeric vector
#' @keywords scRNA-seq
#' @export
#' @examples scale01
scale01 <- function(x){
  y <- min(x); z <- max(x);
  (x-y) / (z-y)
}



#' umap_theme
#'
#' ggplot theme to remove axes etc.
#'
#' @keywords scRNA-seq
#' @export
#' @examples umap_theme
umap_theme <- function(){
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
}



#' shuffle_points
#'
#' Function to shuffle the rows of a dataframe.
#'
#' @param df dataframe
#' @keywords scRNA-seq
#' @export
#' @examples shuffle_points
shuffle_points <- function(df){
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}
