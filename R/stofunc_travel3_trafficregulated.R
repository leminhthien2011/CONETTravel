#' This function gives stochastic realization for 3 countries for a given regulated strategy
#' @param combinetheta parameter
#' @param inp is a list include duration : durationtravel (number of days),
#' travelregulated: a list of matrix travel allowed from 1 country to another during the duration,
#' ini: initial compartments of countries)
#' @return  The average realization of 3 countries with travel data regulated
#' @importFrom stats rpois
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' library(matrixcalc)
#' P1 = 10^7
#' I1 = 250
#' A1 = 130
#' S1 = P1 - I1 - A1
#' x1 = c(S1,I1,A1,0,0,0) # State corresponding S,I,A,R,D,Ru country 1
#' P2 = 3*10^6
#' I2 = 20
#' A2 = 10
#' S2 = P2 - I2 - A2
#' x2 = c(S2,I2,A2,0,0,0) # State corresponding S,I,A,R,D,Ru country 1
#' P3 = 2*10^6
#' I3 = 15
#' A3 = 15
#' S3 = P3 - I3 - A3
#' x3 = c(S3,I3,A3,0,0,0) # State corresponding S,I,A,R,D,Ru country 1
#' initial_corona = c(x1,x2,x3) #initial condition of 3 countries
#' P = c(P1, P2, P3) #population 3 countries
#' k = 13
#' theta0 = as.numeric(thetas_3travel[[k]]) #choose a combo of theta
#' traveloutdat = travelout_3dat
#' ratein = 1 # policy that allows full rate of travel in
#' r12= rep(ratein,nrow(traveloutdat))
#' r13= rep(ratein,nrow(traveloutdat))
#' r21= rep(ratein,nrow(traveloutdat))
#' r23= rep(ratein,nrow(traveloutdat))
#' r31= rep(ratein,nrow(traveloutdat))
#' r32= rep(ratein,nrow(traveloutdat))
#' policy = list(rate12 = r12, rate13 = r13, rate21 = r21, rate23 = r23, rate31 = r31, rate32 = r32)
#' traveloutDivideRegulated =  totaltravelout_regulated(traveloutdat, policy, P) #total travelers regulated
#' inp = list(duration = nrow(traveloutdat), travelregulated = traveloutDivideRegulated, ini = initial_corona )
#' stofunc_travel3_trafficregulated(theta0,inp)}
#' @export

stofunc_travel3_trafficregulated =  function(combinetheta, inp){

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


  status_matrix = matrix(0, nrow = inp$duration,ncol = length(inp$ini))
  f_in = matrix(0, nrow =  inp$duration,ncol = length(inp$ini))
  f_out = matrix(0, nrow =  inp$duration, ncol = length(inp$ini))
  status_matrix[1,] = inp$ini



  for (i in 2: inp$duration){
    traveloutregulated = as.matrix(inp$travelregulated[[i]])
    totaltravelout = rowSums(traveloutregulated)

    for (j in 1:nrow(traveloutregulated)){

      c1 = (j-1)*6 + 1
      c2 = j*6
      x = status_matrix[(i-1),c1:c2]
      ##Updated traveling flow
      #Number out from country j
      out = totaltravelout[j]
      if (x[1]+x[2] > 0){
        outj = c(round(out*x[1]/(x[1]+x[2]),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
        #########
        f_out[i,c1:c2] = outj
      }else{

        f_out[i,c1:c2] = c(out, 0,0,0,0,0)
      }

      #travel out from country i with sick and susceptible

      theta = combinetheta[c1:c2]
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

    ##Construct matrix out of compartments for 3 countries

    f_outmat = matrix(0, nrow = nrow(traveloutregulated), ncol = length(inp$ini) )

    for (val in 1:nrow(traveloutregulated)){
      d1 = (val-1)*6 + 1
      d2 = val*6
      f_outtotal = f_out[i,][d1:d2]

      for (val1 in 1:nrow(traveloutregulated)){
        e1 = (val1 -1)*6 + 1
        e2 = val1*6
        if(sum(traveloutregulated[val,])>0){
          f_outmat[val,e1:e2] = f_outtotal*traveloutregulated[val,val1]/sum(traveloutregulated[val,])
        }else{
          f_outmat[val,e1:e2] = rep(0,6)
        }

        }
    }

    ##Construct matrix in of compartments for 3 countries
    for (val2 in 1:nrow(traveloutregulated)){
      f1 = (val2 -1)*6 + 1
      f2 = val2*6
      f_in[i, f1:f2] = colSums(f_outmat[,f1:f2])

    }

    f_in[i,] = round(f_in[i,], digits=0)
    ########
    update = status_matrix[i,]  +  f_in[i,] - f_out[i,]
    update[update<0.1]=0
    status_matrix[i,] = update

  }
  return(status_matrix)
}


