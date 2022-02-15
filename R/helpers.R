# scale between 0 & 1
scale01 <- function(x){
  y <- min(x); z <- max(x);
  (x-y) / (z-y)
}

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
