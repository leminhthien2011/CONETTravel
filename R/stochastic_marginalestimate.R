#' This gives is an marginal estimate of one country under travel
#' @param theta parameter
#' @param country is a list include durationtravel : number of days, travel in compartments each day,
#' total_out each day, and x_ini: initial compartments of the country
#' @return  The realization of the country during the period take into account of travel
#' @importFrom stats rpois
#' @export

###
stochastic_marginalestimate =  function(theta,country){

  ##################Defining harzard functions
  #New infected rate, alphas,c(alpha0,alpha, beta, delta, eta, gamma)
  harzard1 = function(x,theta){
    h1 = (theta[1] + theta[2])*x[1]*x[2]/sum(x)
    names(h1)=c("hazard1")
    return(h1)
  }

  #New confirmed rate, gamma, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard2 = function(x,theta){
    h2 = theta[6]*x[2]
    names(h2)=c("hazard2")
    return(h2)
  }
  #New confirmed recover, beta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard3 = function(x,theta){
    h3 = theta[3]*x[3]
    names(h3)=c("hazard3")
    return(h3)
  }
  #New confirmed death, delta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard4 = function(x,theta){
    h4 = theta[4]*x[3]
    names(h4)=c("hazard4")
    return(h4)
  }
  #New unconfirmed recover, eta*beta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard5 = function(x,theta){
    h5 = theta[5]*theta[3]*x[2]
    names(h5)=c("hazard5")
    return(h5)
  }
  ##############
  status_matrix = matrix(0, nrow = country$durationtravel, ncol = 6)
  status_matrix[1,] = country$x_ini

  f_in = matrix(0, nrow = country$durationtravel, ncol = 6)

  f_out = matrix(0, nrow = country$durationtravel, ncol = 6)

  for (i in 2:country$durationtravel){
    x = status_matrix[(i-1),]
    ##Updated traveling flow
    #Number out from country 2
    out = country$total_out[i]
    if (x[1]+x[2] > 0){
      out_compartments = c(round(out*x[1]/(x[1]+x[2]),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
    } else{

      out_compartments = c(out, 0,0,0,0,0)
    }

    ##Generating Poisson values based on hazard functions
    y1 = rpois(1, harzard1(x,theta))
    #
    y2 =  rpois(1, harzard2(x,theta))
    #
    y3 = rpois(1, harzard3(x,theta))
    #
    y4 =  rpois(1, harzard4(x,theta))
    #
    y5 = rpois(1, harzard5(x,theta))


    ##Susceptible
    if(y1<= x[1]){
      x[1] = x[1] - y1} else{
        y1 = x[1]
        x[1] =0
      }



    #######Infect
    if(y1-y2-y5+x[2]>= 0){
      x[2] = x[2] + y1 - y2 - y5} else{
        y2 =  x[2] + y1 - y5
        x[2] =0
      }

    if(y2 <0){
      y2 = 0
      y5 = x[2] + y1
    }


    #######Active
    if(y2-y3-y4+x[3] >= 0){
      x[3] = x[3] + y2 - y3 - y4} else{
        y3 = y2-y4+x[3]
        x[3] =0
      }

    if(y3 <0){
      y3 = 0
      y4 = x[3] + y2
    }

    #Recover
    x[4] =x[4] + y3
    #Death
    x[5]= x[5] + y4
    #Recover Unconfirmed
    x[6]= x[6] + y5

    #update the flow traffic

    update = x  + country$travelin_compartments[i,] - out_compartments

    update[update<0.1]=0

    status_matrix[i,] = update


  }
  return(status_matrix)
}