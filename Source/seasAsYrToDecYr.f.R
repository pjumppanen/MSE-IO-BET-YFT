#function to convert a MSE year to a decimal year, eg. for labels 2013.125
mseYrToDecYr.f <- function(mseYr=1, nyearsMSE = 43, lastDecYear=2014){

  decYr <- mseYr + lastDecYear - nyearsMSE

  return(decYr)
}

mseYrSeasToDecYrSeas.f <- function(mseYrQtr=1, decYrQtr=1972.125){

  decYrQtr <- decYrQtr + (mseYrQtr -1) *0.25 

  return(decYrQtr)
}



#function to convert a Stock Synthesis season as year index to a decimal year, eg. for labels 2013.125
seasAsYrToDecYr.f <- function(seasAsYr=101, endSeasAsYr=272, numSeas=4, endYr=2014, endSeas=4){

  decYr <- (endYr + endSeas/numSeas - 1/(2*numSeas)) - (endSeasAsYr-seasAsYr)/numSeas

  return(decYr)
}

#function to convert a Stock Synthesis season as year index to a year and integer season, eg. 2013.125 -> 2013,1 
seasAsYrToYrSeas.f <- function(seasAsYr=101, endSeasAsYr=272, numSeas=4, endYr=2014, endSeas=4){

  decYr  <- (endYr + endSeas/numSeas - 1/(2*numSeas)) - (endSeasAsYr-seasAsYr)/numSeas
  calendarYear   <- floor(decYr)
  season <- floor((decYr - calendarYear)*numSeas)+1

  return(c(calendarYear,season))
}

#function to convert a Stock Synthesis season as year index to a tcMSE year and integer season, eg. 1972.125 -> 1 1
#assumes that nobody is going to define model years misaligned with calendar years  
seasAsYrToMSEYrSeas.f <- function(seasAsYr=101, endSeasAsYr=272, numSeas=4, endYr=2014, endSeas=4, firstSeasAsYr=101, firstSeas=1){

  decYr          <- (endYr + endSeas/numSeas - 1/(2*numSeas)) - (endSeasAsYr-seasAsYr)/numSeas
  calendarYear   <- floor(decYr)
  season         <- floor((decYr - calendarYear)*numSeas)+1
  MSEYear        <- floor((seasAsYr - firstSeasAsYr)/numSeas+1) 

  return(c(MSEYear,season))
}

#function to convert a Stock Synthesis season as year index to a tcMSE year and integer season, eg. 1972.125 -> 1 1
#assumes that nobody is going to define model years misaligned with calendar years  
seasAsYrToMSEYr.f <- function(seasAsYr=101, endSeasAsYr=272, numSeas=4, endYr=2014, endSeas=4, firstSeasAsYr=101, firstSeas=1){

  decYr          <- (endYr + endSeas/numSeas - 1/(2*numSeas)) - (endSeasAsYr-seasAsYr)/numSeas
  calendarYear   <- floor(decYr)
  season         <- floor((decYr - calendarYear)*numSeas)+1
  MSEYear        <- floor((seasAsYr - firstSeasAsYr)/numSeas+1) 

  return(c(MSEYear))
}
