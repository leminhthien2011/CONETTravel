#' This gives estimation of alpha
#' @param data data of A,R,D
#' @param inp a list include beta estimation based on data, infect, pop, travelin, travelout
#' @export
alphafunc = function(data,inp){

  nr = nrow(data)
  beta = median(inp$beta)
  Ut = rowSums(data)
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########
  recovermat = matrix(0,nrow = nr,ncol=6)

  infected = inp$infect

  dRus = rep(0,nr)

  for (i in 1:nr){

    dRus[i] = dUt[i]*beta

  }


  Ruseq1 = cumsum(dRus)

  Ruseq = c(0,Ruseq1[-length(Ruseq1)])


  Sseq = inp$pop - Ut - Ruseq - infected


  SIseq = Sseq*infected

  #########Construct dSseq
  Sseq_p = rep(0, nr)

  for (i in 1:nr){
    Sseq_p[i] = Sseq[i] - inp$travelin[i,1] + inp$travelout[i,1]
  }
  ######S+ and Spre sequence for alpha
  Splus_al = Sseq[-c(nr-1, nr)]

  Spre_al = Sseq_p[-c(1,nr)]

  dSseq_al = Splus_al - Spre_al

  SIseq_al = SIseq[-c(nr-1, nr)]

  pop_al = inp$pop[-c(nr-1, nr)]


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


