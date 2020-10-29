#' This gives descriptive stat
#' @param x numeric vector
#' @example
#' summaryConet(x)
#' @export
summaryConet <- function(x){
  data.frame(min = min(x),
             max = max(x),
             med = median(X))
}
