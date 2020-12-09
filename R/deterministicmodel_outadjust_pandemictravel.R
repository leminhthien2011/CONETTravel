#' This function gives deterministic realization for n countries with a given regulated
#' strategy and quarantine duration required by the destination countries
#' for each travel out country depends on the situation of the travel out country
#' It also returns number new cases active confirmed each day during the pandemic if quarantine
#' exist and traveler status after done quarantine if any quarantine need.
#' @param thetamatrix is a matrix of parameters, parameters of each country is on 1 row
#' @param inp is a list include durationtravel : durationtravel (days),
#'  durationquarantine_adjustedin : number of days people travel in have to quarantine based on each country policy,
#' travelregulated: a list of travel allowed from 1 country to another during the duration,
#' initialmatrix is a matrix of initial compartments of countries, each country is on 1 row, and
#' quarantinerate is the rate people follow quarantine
#'
#' @return  The average realization of n countries with travel data regulated, also return
#' travelers status after quarantine, and active confirm update during the quarantine time

#' @examples
#' \dontrun{
#' library(CONETTravel)

#' thetagenerating = function(lowerbound, upperbound){
#'  tmp2 = 1 # need for kick off
#'  while(tmp2 >0){
#'   theta = c( alpha0 = 0,alpha = runif(1,0,1),beta = runif(1,0,.25), delta=runif(1,0,.25),
#'              eta=1, gamma=runif(1,0,1) )
#'   tmp1 = theta[2]/(theta[3] + theta[6])
#'   tmp2 = (tmp1 - lowerbound)*(tmp1 - upperbound)
#' }
#' return(theta)
#'}

##############################
#' initialmatrix_func =  function(numbercountries){
#'  initialmatrix = matrix(0, numbercountries, 6)
#'  for (country in 1:numbercountries){
#'    P = round(runif(1, 50000, 20000000000), digits=0)
#'    I = round(runif(1, 100000, 200000), digits=0)
#'    S = P - I
#'   initialmatrix[country,] = c(S, I, 0, 0,0,0)
#'  }

#' return(initialmatrix)
#' }

#' travelout_data = travelout_3dat

#' initial_corona =  initialmatrix_func(3)
#' P = rowSums(initial_corona)
#' theta0 = rbind(thetagenerating(1,6.5), thetagenerating(.5,6.5), thetagenerating(.5,6.5))
#' thetamatrix = theta0
#' ratein = 1 # policy that allows full rate of travel in
#' traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout_data, ratein, P)
#'inp = list(durationtravel = 84, travelregulated = traveloutDivideRegulated,
#'          initialmatrix = initial_corona, quarantinerate = 1, durationquarantine_adjustedout = c(0,7,14))

#'  deterministicmodel_outadjust_pandemictravel(theta0, inp)

#' }

#' @export




deterministicmodel_outadjust_pandemictravel = function (thetamatrix, inp)
{

  harzard1 = function(x, theta) {
    h1 = (theta[1] + theta[2]) * x[1] * x[2]/sum(x)
    names(h1) = c("hazard1")
    return(h1)
  }
  harzard2 = function(x, theta) {
    h2 = theta[6] * x[2]
    names(h2) = c("hazard2")
    return(h2)
  }
  harzard3 = function(x, theta) {
    h3 = theta[3] * x[3]
    names(h3) = c("hazard3")
    return(h3)
  }
  harzard4 = function(x, theta) {
    h4 = theta[4] * x[3]
    names(h4) = c("hazard4")
    return(h4)
  }
  harzard5 = function(x, theta) {
    h5 = theta[5] * theta[3] * x[2]
    names(h5) = c("hazard5")
    return(h5)
  }
  numbercountries = nrow(inp$initialmatrix)
  compartments = ncol(inp$initialmatrix)
  status_matrix = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                           compartments)
  f_in = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                  compartments)
  f_out = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                   compartments)
  totalduration = inp$durationtravel + max(inp$durationquarantine_adjustedout)
  f_out_donequarantine_list = list()
  f_out_activequarantine_addlist = list()

  for (outdone in 1:numbercountries) {
    f_out_donequarantine_list[[outdone]] = matrix(0, nrow = totalduration,
                                                  ncol = numbercountries * compartments)
    f_out_activequarantine_addlist[[outdone]] = matrix(0, nrow = totalduration,
                                                       ncol = numbercountries * compartments)


  }



  initial = rep(0, numbercountries * compartments)
  for (val3 in 1:numbercountries) {
    h1 = 1 + (val3 - 1) * 6
    h2 = 6 + (val3 - 1) * 6
    initial[h1:h2] = inp$initialmatrix[val3, ]
  }
  status_matrix[1, ] = initial


  for (i in 2:inp$durationtravel) {

    traveloutregulated = as.matrix(inp$travelregulated[[i]])
    totaltravelout = rowSums(traveloutregulated)

    for (j in 1:numbercountries) {
      c1 = (j - 1) * 6 + 1
      c2 = j * 6
      x = status_matrix[(i - 1), c1:c2]
      out = totaltravelout[j]
      if (x[1] + x[2] > 0) {
        outj = c(round(out * x[1]/(x[1] + x[2]), digits = 0),
                 round(out * x[2]/(x[1] + x[2]), digits = 0),
                 0, 0, 0, 0)
        f_out[i, c1:c2] = outj
      }
      else {
        f_out[i, c1:c2] = c(0, 0, 0, 0, 0, 0)
      }
      theta = thetamatrix[j, ]
      y1 = harzard1(x, theta)
      y2 = harzard2(x, theta)
      y3 = harzard3(x, theta)
      y4 = harzard4(x, theta)
      y5 = harzard5(x, theta)
      if (y1 <= x[1]) {
        x[1] = x[1] - y1
      }
      else {
        y1 = x[1]
        x[1] = 0
      }
      if (y1 - y2 - y5 + x[2] >= 0) {
        x[2] = x[2] + y1 - y2 - y5
      }
      else {
        y2 = x[2] + y1 - y5
        x[2] = 0
      }
      if (y2 < 0) {
        y2 = 0
        y5 = x[2] + y1
      }
      if (y2 - y3 - y4 + x[3] >= 0) {
        x[3] = x[3] + y2 - y3 - y4
      }
      else {
        y3 = y2 - y4 + x[3]
        x[3] = 0
      }
      if (y3 < 0) {
        y3 = 0
        y4 = x[3] + y2
      }
      x[4] = x[4] + y3
      x[5] = x[5] + y4
      x[6] = x[6] + y5
      status_matrix[i, c1:c2] = x
    } #end j loop


    f_outmat = matrix(0, nrow = numbercountries, ncol = numbercountries *
                        compartments)
    for (val in 1:numbercountries) {
      d1 = (val - 1) * 6 + 1
      d2 = val * 6
      f_outtotal = f_out[i, ][d1:d2]
      d3 = d1 + 1
      infect_outtotal = f_out[i, ][d3]
      probdistribute = rep(0, numbercountries)
      for (val6 in 1:numbercountries) {
        if (sum(traveloutregulated[val, ]) > 0) {
          probdistribute[val6] = traveloutregulated[val,
                                                    val6]/sum(traveloutregulated[val, ])
        }
        else {
          probdistribute[val6] = 0
        }
      }
      if (sum(traveloutregulated[val, ]) > 0) {
        infect_outdistribute = rmultinom(1, size = infect_outtotal,
                                         prob = probdistribute)
      }
      else {
        infect_outdistribute = rep(0, numbercountries)
      }
      for (val1 in 1:numbercountries) {
        e1 = (val1 - 1) * 6 + 1
        e2 = val1 * 6
        if (sum(traveloutregulated[val, ]) > 0) {
          tmp = round(f_outtotal * traveloutregulated[val,
                                                      val1]/sum(traveloutregulated[val, ]), digits = 0)
          suseptible = tmp[1] + tmp[2] - infect_outdistribute[val1]
          f_outmat[val, e1:e2] = tmp
        }
        else {
          f_outmat[val, e1:e2] = rep(0, 6)
        }
      }
    } # end val loop

    #####create matrix active out each time step
    f_out_activequarantine_list = list()

    for (outdone in 1:numbercountries) {

      f_out_activequarantine_list[[outdone]] = matrix(0, nrow = totalduration,
                                                      ncol = numbercountries * compartments)
    }


    for (countryout in 1:numbercountries) {


      durationquarantine_countryout = inp$durationquarantine_adjustedout[countryout]

      for (countryin in 1:numbercountries) {


        a1 = 1 + (countryin - 1) * 6
        a2 = 6 + (countryin - 1) * 6
        a3 = 3 + (countryin - 1) * 6
        quarantineinp = inp$quarantinerate * f_outmat[countryout,
                                                      a1:a2]
        inp1 = list(durationquarantine = durationquarantine_countryout,
                    ini = quarantineinp)
        theta1 = thetamatrix[countryin, ]


        i1 = i + durationquarantine_countryout

        tmp = deterministic_postquarantine_separate(theta1, inp1)

        if(durationquarantine_countryout > 0){
          f_out_donequarantine_list[[countryout]][i1, a1:a2] = tmp$donequarantine

          f_out_activequarantine_list[[countryout]][i:i1,a3] = tmp$activeconfirm_eachday[,3]
        } else{
          f_out_donequarantine_list[[countryout]][i, a1:a2] = tmp$donequarantine

          f_out_activequarantine_list[[countryout]][i,a3] = tmp$activeconfirm_eachday[3] # this reduced to a vector

        } #end loop that making sure a valid code when no quarantine required


      }
    } # end quarantine loop, S,I, Ru added; and active confirm each day




    ##########Accumulate ctive cases from travel
    for(country in 1:numbercountries){
      f_out_activequarantine_addlist[[country]] = f_out_activequarantine_addlist[[country]] + f_out_activequarantine_list[[country]]
    }


    ############ Done quarantine return for each country
    f_in_donequarantine = matrix(0, nrow = totalduration,
                                 ncol = numbercountries * compartments) # avoid overlap adding, set to 0

    for (countryout1 in 1:numbercountries) {
      for (countryin1 in 1:numbercountries) {
        aa1 = 1 + (countryin1 - 1) * 6
        aa2 = 6 + (countryin1 - 1) * 6
        f_in_donequarantine[, aa1:aa2] = f_in_donequarantine[,
                                                             aa1:aa2] + f_out_donequarantine_list[[countryout1]][,
                                                                                                                 aa1:aa2]
      }
    }


    #########Active confirmed import updated for each country
    f_in_activeupdated = matrix(0, nrow = totalduration,
                                ncol = numbercountries * compartments) # avoid overlap adding, set to 0

    for (countryout1 in 1:numbercountries) {
      for (countryin1 in 1:numbercountries) {
        aa1 = 1 + (countryin1 - 1) * 6
        aa2 = 6 + (countryin1 - 1) * 6
        f_in_activeupdated[, aa1:aa2] = f_in_activeupdated[,
                                                           aa1:aa2] + f_out_activequarantine_addlist[[countryout1]][,
                                                                                                                    aa1:aa2]
      }
    }




    ##############

    for (val2 in 1:numbercountries) {
      f1 = (val2 - 1) * 6 + 1
      f2 = val2 * 6
      f_in[i, f1:f2] = colSums(f_outmat[, f1:f2])
    }

    update = status_matrix[i, ] + f_in_donequarantine[i,
    ] + (1 - inp$quarantinerate) * f_in[i, ] - f_out[i,
    ] + f_in_activeupdated[i,]



    update[update < 0.5] = 0
    status_matrix[i, ] = update
  }


  return(list(model_output = round(status_matrix, digits = 0), activeconfirmtravel = f_in_activeupdated[1:inp$durationtravel,],
              travelarrival_postquarantine = f_in_donequarantine[1:inp$durationtravel,] ))
}
