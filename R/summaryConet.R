#' This gives descriptive stat
#' @param x numeric vector
#' @example
#' summaryConet(x)
#' @export
summaryConet <- function(x){
  tibble(min = min(x),
             max = max(x),
             med = median(x))
}
