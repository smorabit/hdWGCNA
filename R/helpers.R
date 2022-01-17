# scale between 0 & 1
scale01 <- function(x){
  y <- min(x); z <- max(x);
  (x-y) / (z-y)
}
