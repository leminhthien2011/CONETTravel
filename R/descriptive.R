#' This gives descriptive stat as tibble
#' @param x numeric vector
#' @import tibble
#' @importFrom stats median quantile sd
#' @examples
#' x = 1:9
#' summaryConet(x)
#' \dontrun{A tibble: 1 x 9
#' len  mean    sd   min   max   med firstpt thirdpt   IQR
#' <int> <dbl> <dbl> <int> <int> <int>   <dbl>   <dbl> <dbl>
#'  9     5  2.74     1     9     5       3       7     4
#'  }
#' summaryConet(x)
#' @export
summaryConet <- function(x){
  tibble(len = length(x),
         mean = mean(x),
         sd = sd(x),
         min = min(x),
         max = max(x),
         med = median(x),
         firstpt = quantile(x,.25),
         thirdpt = quantile(x, .75),
         IQR = thirdpt - firstpt)

}
