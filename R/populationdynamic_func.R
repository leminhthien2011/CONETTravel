#' This function gives the estimation for the dynamic population.
#' @param country list include total travel in: total_in, total travel out each day,
#' total_out, durationtravel, population P
#' @export

populationdynamicfunc = function(country){
  fluc1 = country$total_in - country$total_out
  pop = rep(0,country$durationtravel)
  pop[1] = country$P

  for(i in 2:country$durationtravel){
    x = pop[i-1]
    pop[i] = x+ fluc1[i] #update total population
  }
  return(pop)

}
