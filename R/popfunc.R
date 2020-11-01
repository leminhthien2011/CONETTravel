#' This gives estimation for the dynamic population
#' @param country list include total travel in and out each day total_in, total_out, duration nrep, population P
#' @export

popfunc = function(country){
  fluc1 = country$total_in - country$total_out
  pop = rep(0,country$nrep)
  pop[1] = country$P

  for(i in 2:country$nrep){
    x = pop[i-1]
    pop[i] = x+ fluc1[i] #update total population
  }
  return(pop)

}
