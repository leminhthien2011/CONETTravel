#' This gives stochastic realization for 3 countries
#' @param ini initial conditions
#' @param thetas parameter
#' @param travel travel data
#' @param d duration travel
#' @importFrom stats rpois
#' @export
stofunc_travel3 =  function(ini, thetas , travel, d){

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




  status_matrix = matrix(0, nrow = d,ncol = length(ini))

  f_in = matrix(0, nrow = d,ncol = length(ini))
  f_in
  f_out = matrix(0, nrow = d, ncol = length(ini))
  status_matrix[1,] = ini



  for (i in 2:d){
    for (j in 1:3){

      c1 = (j-1)*6 + 1
      c2 = j*6
      x = status_matrix[(i-1),c1:c2]
      ##Updated traveling flow
      #Number out from country j
      out = travel[i,j]
      if (x[1]+x[2] > 0){
        outj = c(round(out*x[1]/(x[1]+x[2]),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
        #########
        f_out[i,c1:c2] = outj
      }
      #travel out from country i with sick and susceptible

      theta = thetas[c1:c2]
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



      status_matrix[i,c1:c2] = x

    }



    ######################No regulation ####
    P1 = sum(ini[1:6])
    P2 = sum(ini[7:12])
    P3 = sum(ini[13:18])


    f_in1 = f_out[i,7:12]*P1/(P1 + P3) + f_out[i,13:18]*P1/(P1+P2)
    f_in2 = f_out[i,1:6]*P2/(P2+P3) + f_out[i,13:18]*P2/(P1+P2)
    f_in3 = f_out[i,1:6]*P3/(P2+P3) + f_out[i,7:12]*P3/(P1+P3)


    f_in[i,] = c(round(f_in1,digits = 0), round(f_in2,digits=0), round(f_in3,digits=0))
    ########
    update = status_matrix[i,]  +  f_in[i,] - f_out[i,]
    update[update<0.1]=0
    status_matrix[i,] = update

  }
  return(status_matrix)
}



