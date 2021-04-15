#' This function estimate number infected in each country for a given data set
#'  under our model assumption.
#' @param data A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased)
#' @param x initial corona
#' @export
infectfunc = function(data,x){
  U = t(t(rowSums(data)))
  diff1 = diff(U, differences = 1)
  diffa = diff1[-1]
  diffb = diff1[-length(diff1)]
  infect = rep(0,length(diff1))

  ###############
  infect[1] = x[2]

  for (i in 2:length(diff1)){
    i1 = i-1
    tmp = infect[i1]


    if (diffb[i1]>.4){
      infect[i] = tmp*diffa[i1]/diffb[i1]

    } else{
      infect[i] = diff1[i]


    }

  }

  return (infect)

}

