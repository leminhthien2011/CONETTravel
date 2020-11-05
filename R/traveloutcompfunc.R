#' This gives travel out compartment in 6 categories from 1 country to another
#' @param total total travelers out
#' @param country a list include data of A,R,D sequence, population pop, infect
#' @export

traveloutcompfunc = function(total,country)
{


  mydat = country$data
  pop = country$pop
  infect = country$infect


  #######

  drecover = diff(mydat[,2], differences = 1)

  active = mydat[,1][-length(mydat[,1])]

  index0 = which(active<=.1)

  if (length(index0)!=0){
    drecover = drecover[-index0]
    active = active[-index0]
  }

  betas = drecover/active


  Ut = rowSums(mydat)
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########
  dRus = rep(0,nrow(mydat))

  for (i in 1:(nrow(mydat)-1)){

    dRus[i] = infect[i]*median(betas)

  }


  Ruseq1 = cumsum(dRus)

  Ruseq = c(0,Ruseq1[-length(Ruseq1)])


  Sseq = pop - Ut - Ruseq - infect



  compartments = cbind(Sseq, infect, Ruseq)


  update = matrix(0,ncol=6,nrow=nrow(mydat))

  for (i in 2:nrow(mydat)){
    i1= i-1
    x = compartments[i1,]

    ####
    if(min(x)>0){
      update[i,] = c(total[i]*x[1]/(x[1]+x[2]), total[i]*x[2]/(x[1]+x[2]), 0, 0, 0, 0)
    } else{
      update[i,] = c(total[i], 0, 0, 0, 0, 0)
    }
  }
  return(update)
}

