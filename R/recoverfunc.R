#' This gives the average realization of 1 given country for a given A,R,D sequence
#' @param country is a list include data of A,R,D, pop,and infect
#' @param theta parameter
#' @export
recoverfunc = function(country, theta){

  data = country$data
  pop = country$pop

  beta = theta[3]

  Ut = data[,1] + data[,2] + data[,3]
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########

  recovermat = matrix(0,nrow = nrow(data),ncol=6)

  infected = country$infect

  dRus = rep(0,nrow(data))

  for (i in 1:nrow(data)){

    dRus[i] = dUt[i]*beta

  }

  Ruseq1 = cumsum(dRus)

  Ruseq = c(0,Ruseq1[-length(Ruseq1)])


  Sseq = pop - Ut - Ruseq - infected

  recovermat = data.frame(Sseq, infected, data[,1], data[,2], data[,3], Ruseq)



  return(recovermat)

}
