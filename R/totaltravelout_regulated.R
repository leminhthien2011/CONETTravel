#' This function gives a list of total travel data allowed from 1 country to another for a given travel policy
#' @param traveloutdat is the data of travel out from each country without pandemic
#' @param P vector population of countries in order as travel data
#' @param policy is a list of policy of rate 12 country 1 let country 2 in during the duration travel
#' rate 13 country 1 let country 3 in, rate 21 country 2 let country 1 in, rate 23 country 2 let country 3 in,
#' rate 31 country 3 let country 1 in, rate 32 country 3 let country 2 in
#' @return  A list of number travelers can go from 1 country to another each day during the duration
#' @examples
#' \dontrun{
#' P1 = 10^7
#' P2 = 3*10^6
#' P3 = 2*10^6
#' P = c(P1, P2, P3) #population of 3 countries
#' traveloutdat = travelout_3dat
#' ratein = 1 # policy that allows full rate of travel in
#' r12= rep(ratein,nrow(traveloutdat))
#' r13= rep(ratein,nrow(traveloutdat))
#' r21= rep(ratein,nrow(traveloutdat))
#' r23= rep(ratein,nrow(traveloutdat))
#' r31= rep(ratein,nrow(traveloutdat))
#' r32= rep(ratein,nrow(traveloutdat))
#' policy = list(rate12 = r12, rate13 = r13, rate21 = r21, rate23 = r23, rate31 = r31, rate32 = r32)
#' totaltravelout_regulated(traveloutdat, policy, P)
#' }


#' @import matrixcalc
#' @export


totaltravelout_regulated =  function(traveloutdat,policy, P){
##Construct a divide list travel out
numberCountries = ncol(traveloutdat)

traveloutDivide = list()

for (j in 1:nrow(traveloutdat)){
  tmp = traveloutdat[j,]
  ##CREATE MATRIX TRAVELOUT EACH ROW
  travelmat = matrix(0,nrow = numberCountries, ncol = numberCountries)
  for(i in 1:numberCountries){
    ones = rep(1,numberCountries)
    travelmat[i,] = kronecker(ones,as.numeric(tmp[i]) )
  }
  ##CREATE THE PROPORTION OF POPULATION MATRIX
  ones = t(t(rep(1,numberCountries)))
  popmat = kronecker(ones,t(P))
  ##
  denommat = matrix(0,nrow = numberCountries, ncol = numberCountries)
  for(i in 1:numberCountries){
    tmp = sum(P[-i])
    ones = rep(1,numberCountries)
    denommat[i,] = kronecker(ones,as.numeric(tmp) )
  }
  ##Take product for proportion
  promat = popmat*denommat^(-1)
  diag(promat) = 0
  # Take Hadamard product for actual number
  traveloutDivide[[j]] = hadamard.prod(promat, travelmat)
}


##Matrix list of traveling
durationtravel = nrow(traveloutdat)


c11 = rep(0, durationtravel)
c12 = policy$rate12 #country 1 let percent country 2 in
c13 = policy$rate13 ##country 1 let percent country 3 in

c21 = policy$rate21 #country 2 let percent country 1 in
c22 = rep(0, durationtravel)
c23 = policy$rate23 ##country 2 let percent country 3 in

c31 = policy$rate31 #country 3 let percent country 1 in
c32 = policy$rate32 ##country 3 let percent country 2 in
c33 = rep(0,durationtravel)

policymat = rbind(c11, c12, c13, c21, c22, c23, c31, c32, c33)

policymatlist = list()

for(i in 1:durationtravel){
  tmp =  matrix(0,nrow = numberCountries, ncol = numberCountries)
  tmp = matrix(policymat[,i], nrow = numberCountries, ncol = numberCountries)
  policymatlist[[i]] = tmp

}

############Matrix travel when percentage of travel in regulated
traveloutDivideRegulated = list()
for (i in 1:nrow(traveloutdat)){

  tmp1 = as.matrix(traveloutDivide[[i]])
  tmp2 = as.matrix(policymatlist[[i]])

  traveloutDivideRegulated[[i]] = round(hadamard.prod(tmp1, tmp2), digits=0)
}

return(traveloutDivideRegulated)
}
######
