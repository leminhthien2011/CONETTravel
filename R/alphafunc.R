#' This function gives a rough estimation of of the transmission rate alpha
#'  and used to improve the estimation of alpha.
#' @param data data of A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased)
#' @param country a list include betas estimation based on data, infect,
#' populationdynamic, travelin_compartments, travelout_compartments
#' @export
alphafunc = function(data,country){


  beta = median(country$betas)
  Ut = rowSums(data)
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


  Sseq = country$populationdynamic - Ut - Ruseq - infected


  SIseq = Sseq*infected

  #########Construct dSseq
  Sseq_p = rep(0, nrow(data))

  for (i in 1:nrow(data)){
    Sseq_p[i] = Sseq[i] - country$travelin_compartments[i,1] + country$travelout_compartments[i,1]
  }
  ######S+ and Spre sequence for alpha
  Splus_al = Sseq[-c(nrow(data)-1, nrow(data))]

  Spre_al = Sseq_p[-c(1,nrow(data))]

  dSseq_al = Splus_al - Spre_al

  SIseq_al = SIseq[-c(nrow(data)-1, nrow(data))]

  pop_al = country$populationdynamic[-c(nrow(data)-1, nrow(data))]


  #####Remove 0s in the denom
  index1 = which(SIseq_al<=.5)
  ####Remove negative difference
  index2 = which(dSseq_al<=0.5)
  ###########
  index = c(index1,index2)
  index = unique(index)
  index = sort(index)

  ######
  if (length(index!=0)){
    dSseq_al = dSseq_al[-index]
    SIseq_al = SIseq_al[-index]
    pop_al = pop_al[-index]
  }


  alphas_sim = dSseq_al/SIseq_al*pop_al

  ########This help to avoid empty array
  if (length(alphas_sim) == 0){alphas_sim = 10^8}
  ###########


  return(alphas_sim)

}


