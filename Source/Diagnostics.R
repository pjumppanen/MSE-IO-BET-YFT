# ======================================================================================================================
# ==== MSE Diagnostics =================================================================================================
# ======================================================================================================================

getperf <- function(object)
{
  MSEobj    <- object
  nsim      <- MSEobj@nsim
  proyears  <- MSEobj@proyears
  nMPs      <- MSEobj@nMPs
  ntargpop  <- length(MSEobj@targpop)

  prjd  <- (MSEobj@nyears + 1):(MSEobj@nyears + MSEobj@proyears)
  CIind <- prjd - MSEobj@nyears - 1

  GK    <- round(apply(MSEobj@F_FMSY[,,prjd] < 1 & MSEobj@B_BMSY[,,prjd] > 1, 1:2, sum) / (MSEobj@proyears), 3)
  Y     <- round(apply(MSEobj@CM[,,MSEobj@targpop,prjd], 1:2, mean), 0)
  Y5    <- round(apply(array(MSEobj@CM[,,MSEobj@targpop,prjd], c(nMPs,nsim,ntargpop,proyears)) * array(rep(0.95 ^ CIind, each=nMPs * nsim * ntargpop), c(nMPs,nsim,ntargpop,proyears)), 1:2, mean), 0)
  Y10   <- round(apply(array(MSEobj@CM[,,MSEobj@targpop,prjd], c(nMPs,nsim,ntargpop,proyears)) * array(rep(0.9 ^ CIind, each=nMPs * nsim * ntargpop), c(nMPs,nsim,ntargpop,proyears)), 1:2, mean), 0)
  Y_5   <- round(apply(array(MSEobj@CM[,,MSEobj@targpop,prjd], c(nMPs,nsim,ntargpop,proyears)) * array(rep(1.05 ^ CIind, each=nMPs * nsim * ntargpop), c(nMPs,nsim,ntargpop,proyears)), 1:2, mean), 0)
  ind   <- prjd[1:(length(prjd)-1)]
  ind2  <- prjd[2:length(prjd)]
  AAVY  <- apply(((MSEobj@CM[,,MSEobj@targpop,ind2] - MSEobj@CM[,,MSEobj@targpop,ind]) ^ 2) ^ 0.5, 1:2, mean)

  list(GK,Y,AAVY,Y5,Y10,Y_5)
}


createMsevizPerformanceData <- function(object, MPsSub=NA, nsimSub=NA)
{
  perfAddRows <- function(df, yr, MPs, MPnames, data, indicator, name)
  {
    for (mp in MPs)
    {
      nsims <- dim(data)[2]
      C1 <- rep(indicator, times=nsims)
      C2 <- rep(yr, times=nsims)
      C3 <- data[mp,]
      C4 <- 1:nsims
      C5 <- rep(name, times=nsims)
      C6 <- rep(MPnames[mp], times=nsims)
      df <- rbind(data.frame(df), data.frame(indicator=C1,year=C2,data=C3, iter=C4, name=C5, mp=C6))
    }

    return (df)
  }

  result.list <- list()

  MSEobj    <- object
  npop      <- MSEobj@npop
  proyears  <- MSEobj@proyears
  allyears  <- MSEobj@proyears + MSEobj@nyears
  MPs       <- if (any(is.na(MPsSub))) 1:MSEobj@nMPs else MPsSub
  Sims      <- if (any(is.na(nsimSub))) 1:MSEobj@nsim else nsimSub
  nsim      <- length(Sims)
  nMPs      <- length(MPs)
  ntargpop  <- length(MSEobj@targpop)

  firstMPy  <- MSEobj@nyears + MSEobj@firstMPYr - MSEobj@lastCalendarYr

  # projection period to report on (starting from the first MPyear and ignoring the bridging years)
  projPeriodList <- c(1,3,5,10,20,allyears-firstMPy,1001,1002) #1001-2 are YFT tuning years (temporary 'cause tunin spects numeric)

  for (pp in projPeriodList)
  {
    df <- NULL

    if(pp<1000)
    {
      ppnum       <- as.numeric(pp)
      projPeriod  <-(firstMPy):(firstMPy + ppnum - 1)
    }

    if(pp>1000)
    {
      if (pp == 1001) projPeriod  <- MSEobj@nyears - MSEobj@lastCalendarYr + 2019:2039
      if (pp == 1002) projPeriod  <- MSEobj@nyears - MSEobj@lastCalendarYr + 2024

      ppnum <- length(projPeriod)
    }

    projPeriodm1 <- projPeriod - 1
    MPnames      <- MSEobj@MPs

    #Performance statistics dim(nMPs, nsim)
    # F1 Pr(SB>0.2SB0)
    PrSBgtp2SB0   <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod]>0.2, MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, PrSBgtp2SB0, "F1", "Pr(SB>0.2SB0)")

    # F2 Pr(SB>SBlim) where SBlim = 0.4SSBMSY
    PrSBgtSBlim   <- round(apply(as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod]>MSEobj@SBlim, MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, PrSBgtSBlim, "F2", "Pr(SB>SBlim)")

    # S1 mean(SB/SB_0)
    SBoSB0        <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, SBoSB0, "S1", "mean(SB/SB_0)")

    # S2 min(SB/SB0)
    minSBoSB0     <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), min), digits=3)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, minSBoSB0, "S2", "min(SB/SB_0)")

    # S3 mean(SB/SB_MSY)
    SBoSBMSY      <- round(apply(as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, SBoSBMSY, "S3", "mean(SB/SB_MSY)")

    # S4 mean(F/F_target), in this case...Ftarget = FMSY
    FoFtarg       <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, FoFtarg, "S4", "mean(F/F_target)")

    # S5 mean(F/F_MSY)
    FoFMSY        <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, FoFMSY, "S5", "mean(F/F_MSY)")

    # S6 Pr(Green)
    PrGreen       <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod] < 1 & as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod] > 1, 1:2, sum) / (ppnum), 3)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, PrGreen, "S6", "Pr(Green)")

    # S7 Pr(Red)
    PrRed         <-round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod] > 1 & as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod] < 1, 1:2, sum) / (ppnum), 3)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, PrRed, "S7", "Pr(Red)")

    # S8 Pr(SB>SB_MSY)
    SBgtSBMSY     <- round(apply(as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod] > 1.0, MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, SBgtSBMSY, "S8", "Pr(SB>SB_MSY)")

    # T1 mean(C(t)/C(t-1))
    CtonCtm1      <- apply((as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod] / as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriodm1]) , 1:2, mean)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, CtonCtm1, "T1", "mean(C(t)/C(t-1))")

    # T2 var(C)
    varC          <-round(apply(as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod], MARGIN=c(1:2), var), 2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, varC, "T2", "var(C)")

    # T3 var(F)
#    varF          <-round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod] * FMSY1, MARGIN=c(1:2), var), 2)
#    df            <- perfAddRows(df, firstMPy, MPs, MPnames, varF, "T3", "var(F)")

    # T4 Pr(C<0.1MSY)
    PrCltp1MSY    <- round(apply(as.karray(MSEobj@C_MSY)[keep(MPs),keep(Sims),,projPeriod]<0.1, MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, PrCltp1MSY, "T4", "Pr(C<0.1MSY)")

    # Y1 mean(C)
    C             <-round(apply(as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod], MARGIN=c(1:2), mean), 0)/1000
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, C, "Y1", "mean(C)")

    # Y3 mean(C/MSY)
    CoMSY         <- round(apply(as.karray(MSEobj@C_MSY)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), mean), digits=2)
    df            <- perfAddRows(df, firstMPy, MPs, MPnames, CoMSY, "Y3", "mean(C/MSY)")

    result.list[[paste(pp)]] <- as.data.table(df)
  } #pp

  return(result.list)
}


createMsevizOmRunsData <- function(mseObj, plotByRF=T, param="", MPsSub=NA, nsimSub=NA)
{
  runsAddRows <- function(df, MP_SY, MPnames, YrLabels, data, name)
  {
    C1 <- MPnames[MP_SY[,1]]
    C2 <- YrLabels[MP_SY[,3]]
    C3 <- data[MP_SY]
    C4 <- MP_SY[,2]
    C5 <- rep(name, times=length(MP_SY[,1]))
    df <- rbind(data.frame(df), data.frame(mp=C1, year=C2, data=C3, iter=C4, qname=C5))

    return (df)
  }

  omsAddRows <- function(df, SY, YrLabels, data, name)
  {
    C1 <- YrLabels[SY[,2]]
    C2 <- data[SY]
    C3 <- SY[,1]
    C4 <- rep(name, times=length(SY[,1]))
    df <- rbind(data.frame(df), data.frame(year=C1, data=C2, iter=C3, qname=C4))

    return (df)
  }

  suffix    <- ""
  proyears  <- mseObj@proyears
  allyears  <- mseObj@proyears + mseObj@nyears
  MPs       <- if (any(is.na(MPsSub))) 1:mseObj@nMPs else MPsSub
  nSims     <- if (any(is.na(nsimSub))) 1:mseObj@nsim else nsimSub

  if (nchar(param) > 0)
  {
    if (!is.null(mseObj@sp_names[[param]]))
    {
      nSims <- which(mseObj@sp_idx == mseObj@sp_names[[param]])
    }
    else if (!is.null(mseObj@h_names[[param]]))
    {
      nSims <- which(mseObj@h_idx == mseObj@h_names[[param]])
    }
    else if (!is.null(mseObj@M_names[[param]]))
    {
      nSims <- which(mseObj@M_idx == mseObj@M_names[[param]])
    }
    else if (!is.null(mseObj@t_names[[param]]))
    {
      nSims <- which(mseObj@t_idx == mseObj@t_names[[param]])
    }
    else if (!is.null(mseObj@q_names[[param]]))
    {
      nSims <- which(mseObj@q_idx == mseObj@q_names[[param]])
    }

    suffix <- "(" %&% param %&% ")"
  }

  oms_df  <- NULL
  runs_df <- NULL
  SY      <- as.matrix(expand.grid(nsims=nSims, nyrs=(1:mseObj@nyears)))
  MP_SY   <- as.matrix(expand.grid(mps=MPs, nsims=nSims, nyrs=(mseObj@nyears:allyears)))
  MPnames <- mseObj@MPs

  # Aggregate CPUE series skips final year
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@IobsArchive, "CPUE(aggregate)")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@IobsArchive[1,,], "CPUE(aggregate)")

  # Recruitment
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@Rec, "Recruitment")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@Rec[1,,], "Recruitment")

  #B/B0
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, apply(mseObj@B_B0, MARGIN=c(1,2,4), FUN=sum), "B/B0")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, apply(mseObj@B_B0[1,,,], MARGIN=c(1,3), FUN=sum), "B/B0")

  #B/BMSY
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@B_BMSY, "B/BMSY")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@B_BMSY[1,,], "B/BMSY")

  #SSB/SSB0
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, apply(mseObj@SSB_SSB0, MARGIN=c(1,2,4), FUN=sum), "SSB/SSB0")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, apply(mseObj@SSB_SSB0[1,,,], MARGIN=c(1,3), FUN=sum), "SSB/SSB0")

  #SSB/SSBMSY
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@SSB_SSBMSY, "SSB/SSBMSY")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@SSB_SSBMSY[1,,], "SSB/SSBMSY")

  #F/FMSY
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@F_FMSY, "F/FMSY")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@F_FMSY[1,,], "F/FMSY")

  #Catch
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, apply(mseObj@CM, MARGIN=c(1,2,4), FUN=sum), "C")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, apply(mseObj@CM[1,,,], MARGIN=c(1,3), FUN=sum), "C")

  # Rec by Qtr
  QtrYrBases <- (0:(allyears-1)) * 4

  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@RecYrQtr[,,QtrYrBases + 1], "Recruitment Q1")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@RecYrQtr[1,,QtrYrBases + 1], "Recruitment Q1")
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@RecYrQtr[,,QtrYrBases + 2], "Recruitment Q2")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@RecYrQtr[1,,QtrYrBases + 2], "Recruitment Q2")
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@RecYrQtr[,,QtrYrBases + 3], "Recruitment Q3")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@RecYrQtr[1,,QtrYrBases + 3], "Recruitment Q3")
  runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@RecYrQtr[,,QtrYrBases + 4], "Recruitment Q4")
  oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@RecYrQtr[1,,QtrYrBases + 4], "Recruitment Q4")

  if (plotByRF)
  {
    #Catch by fishery
    if (mseObj@nfleets > 1)
    {
      for (fi in 1:mseObj@nfleets)
      {
        runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@CMbyF[,,mseObj@targpop,,fi], paste("C by Fleet",fi))
        oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@CMbyF[1,,mseObj@targpop,,fi], paste("C by Fleet",fi))
      }
    }

    #CPUE by region
    if (mseObj@nareas > 1)
    {
      for (ri in 1:mseObj@nareas)
      {
        runs_df <- runsAddRows(runs_df, MP_SY, MPnames, mseObj@yrLabels, mseObj@IobsRArchive[,,,ri], paste("CPU by Area",ri))
        oms_df  <- omsAddRows(oms_df, SY, mseObj@yrLabels, mseObj@IobsRArchive[1,,,ri], paste("CPU by Area",ri))
      }
    }
  }

  return (list(oms=as.data.table(oms_df), runs=as.data.table(runs_df)))
}


#DK modified summary for MWG statistics
tableMSE.f <- function(object, percentiles=c(0.1,0.25,0.5,0.75,0.9), MPsSub=NA, nsimSub=NA)
{
  stats     <- NULL
  MSEobj    <- object
  npop      <- MSEobj@npop
  proyears  <- MSEobj@proyears
  allyears  <- MSEobj@proyears + MSEobj@nyears
  MPs       <- if (any(is.na(MPsSub))) 1:MSEobj@nMPs else MPsSub
  Sims      <- if (any(is.na(nsimSub))) 1:MSEobj@nsim else nsimSub
  nsim      <- length(Sims)
  nMPs      <- length(MPs)
  ntargpop  <- length(MSEobj@targpop)

  firstMPy  <- MSEobj@nyears + MSEobj@firstMPYr - MSEobj@lastCalendarYr

  # projection period to report on (starting from the first MPyear and ignoring the bridging years)
#  projPeriodList <- c(1,3,5,10,20,allyears-firstMPy)
#  for (pp in projPeriodList)
#  {
#    projPeriod   <-(firstMPy):(firstMPy + pp -1)
#    projPeriodm1 <-(firstMPy-1):(firstMPy + pp -2)

  projPeriodList <- c(1,3,5,10,20,allyears-firstMPy,1001,1002) #1001-2 are YFT tuning years (temporary 'cause tunin spects numeric)

  for (pp in projPeriodList)
  {
    if(pp<1000){
      ppnum <- as.numeric(pp)
      projPeriod   <-(firstMPy):(firstMPy + ppnum -1)
    }
    if(pp>1000){
      if (pp == 1001) projPeriod  <- MSEobj@nyears - MSEobj@lastCalendarYr + 2019:2039
      if (pp == 1002) projPeriod  <- MSEobj@nyears - MSEobj@lastCalendarYr + 2024
      ppnum <- length(projPeriod)
    }

    projPeriodm1 <- projPeriod-1

    #Performance statistics dim(nMPs, nsim) #numbers from final report
    #1
    SBoSB0        <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), mean), digits=2)

    #2 - min SB relative to SB0
    minSBoSB0     <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), min), digits=3)

    #3
    SBoSBMSY      <- round(apply(as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)

    #4
    FoFMSY        <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)

    #5 = #4 in this case...Ftarget = FMSY
    FoFtarg       <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)

    #6 - Probability in green Kobe quadrant
    GK            <- round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod] < 1 & as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod] > 1, 1:2, sum) / (ppnum), 3)

    #7 - Probability in red Kobe quadrant
    RK            <-round(apply(as.karray(MSEobj@F_FMSY)[keep(MPs),keep(Sims),projPeriod] > 1 & as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod] < 1, 1:2, sum) / (ppnum), 3)

    # 8 Pr(SB > 0.2SB0)
    PrSBgt0.2SB0  <- round(apply(as.karray(MSEobj@SSB_SSB0)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod]>0.2, MARGIN=c(1:2), mean), digits=2)

    # 9 Pr(SB > SBlim) where SBlim = 0.4SSBMSY
    PrSBgtSBlim   <- round(apply(as.karray(MSEobj@SSB_SSBMSY)[keep(MPs),keep(Sims),projPeriod]>MSEobj@SBlim, MARGIN=c(1:2), mean), digits=2)

    #10 mean Catch
    Y             <-round(apply(as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod], MARGIN=c(1:2), mean), 0)/1000

    #11 mean Catch by fishery/area
    #see summaryByAF() for region and gear breakdown

    #12
    YoMSY         <- round(apply(as.karray(MSEobj@C_MSY)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), mean), digits=2)

    #13 mean catch rates relative to catch rates over last data year; see summaryByAF for region and gear breakdown
    relCPUE       <- round(apply(as.karray(MSEobj@IobsArchive)[keep(MPs),keep(Sims),projPeriod],   MARGIN=c(1:2), mean)/
                           (as.karray(MSEobj@IobsArchive)[keep(MPs),keep(Sims),MSEobj@nyears]), digits=2)

    #14 mean absolute proportional change in catch
    #Average proportional change in Yield between consecutve years: abs(1-Ct/C(t-1))
    APCY          <- apply(abs(1-(as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod] / as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriodm1])) , 1:2, mean)

    #15 Var(catch) ...changed to CV%
    YcvPct        <-round(apply(as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod], MARGIN=c(1:2), sd), 2)/Y

    #16 Pr(fishery collapse) = Pr(years) C<0.1MSY
    PrYlt0.1MSY   <- round(apply(as.karray(MSEobj@C_MSY)[keep(MPs),keep(Sims),,projPeriod]<0.1, MARGIN=c(1:2), mean), digits=2)

    # others not listed in MWP list
    BoBMSY        <- round(apply(as.karray(MSEobj@B_BMSY)[keep(MPs),keep(Sims),projPeriod], MARGIN=c(1:2), mean), digits=2)
    BoB0          <- round(apply(as.karray(MSEobj@B_B0)[keep(MPs),keep(Sims),,projPeriod], MARGIN=c(1:2), mean), digits=2)

    #Average Absolute value change in Yield
    AAVY          <- apply(((as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriod] - as.karray(MSEobj@CM)[keep(MPs),keep(Sims),MSEobj@targpop,projPeriodm1]) ^ 2) ^ 0.5, 1:2, mean)

    temp <- data.frame(cbind(apply(SBoSB0, MARGIN=1,mean),      t(apply(SBoSB0, MARGIN=1,quantile, probs=percentiles)),
                             apply(minSBoSB0, MARGIN=1,mean),   t(apply(minSBoSB0, MARGIN=1,quantile, probs=percentiles)),
                             apply(SBoSBMSY, MARGIN=1,mean),    t(apply(SBoSBMSY, MARGIN=1,quantile, probs=percentiles)),
                             apply(FoFMSY, MARGIN=1,mean),      t(apply(FoFMSY, MARGIN=1,quantile, probs=percentiles)),
                             apply(FoFtarg, MARGIN=1,mean),     t(apply(FoFtarg, MARGIN=1,quantile, probs=percentiles)),
                             apply(GK,MARGIN=1,mean),           t(apply(GK,MARGIN=1,quantile, probs=percentiles)),
                             apply(RK,MARGIN=1,mean),           t(apply(RK,MARGIN=1,quantile, probs=percentiles)),
                             apply(PrSBgt0.2SB0,MARGIN=1,mean), t(apply(PrSBgt0.2SB0,MARGIN=1,quantile, probs=percentiles)),
                             apply(PrSBgtSBlim,MARGIN=1,mean),  t(apply(PrSBgtSBlim,MARGIN=1,quantile, probs=percentiles)),
                             apply(Y, MARGIN=1,mean),           t(apply(Y, MARGIN=1,quantile, probs=percentiles)),
                             apply(relCPUE, MARGIN=1,mean),     t(apply(relCPUE, MARGIN=1,quantile, probs=percentiles)),
                             apply(YoMSY, MARGIN=1,mean),       t(apply(YoMSY, MARGIN=1,quantile, probs=percentiles)),
                             apply(APCY,MARGIN=1,mean),         t(apply(APCY,MARGIN=1,quantile, probs=percentiles)),
                             apply(YcvPct,MARGIN=1,mean),       t(apply(YcvPct,MARGIN=1,quantile, probs=percentiles,na.rm=T)),
                             apply(PrYlt0.1MSY,MARGIN=1,mean),  t(apply(PrYlt0.1MSY,MARGIN=1,quantile, probs=percentiles))))

    colnames <- c("SBoSB0" %&% c("mean",percentiles),
                  "minSBoSB0" %&% c("mean",percentiles),
                  "SBoSBMSY" %&% c("mean",percentiles),
                  "FoFMSY" %&% c("mean",percentiles),
                  "FoFtarg" %&% c("mean",percentiles),
                  "GK" %&% c("mean",percentiles),
                  "RK" %&% c("mean",percentiles),
                  "PrSBgt0.2SB0" %&% c("mean",percentiles),
                  "PrSBgtSBlim" %&% c("mean",percentiles),
                  # Y by area/fishery
                  "Y" %&% c("mean",percentiles),
                  "relCPUE" %&% c("mean",percentiles),
                  # CPUE by area/fishery
                  "YoMSY" %&% c("mean",percentiles),
                  "APCY" %&% c("mean",percentiles),
                  "YcvPct" %&% c("mean",percentiles),
                  "PrYlt0.1MSY" %&% c("mean",percentiles))

    if (any(is.na(nsimSub)))
    {
      sub_names <- (if (length(MSEobj@t_names) == 0) c("h_", "M_", "q_") else c("h_", "M_", "t_", "q_"))

      for (sub_name in sub_names)
      {
        sub_factor      <- factor(as.array(slot(MSEobj, paste(sub_name, "idx", sep=""))))
        factors         <- list(sub_factor)
        dim_mean        <- c(length(levels(sub_factor)), 1)
        dim_percentiles <- c(length(percentiles), length(levels(sub_factor)))

        ExtraResults <- NULL

        for (MP in MPs)
        {
          resultsByFactor <- data.frame(cbind(array(tapply(SBoSB0[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(SBoSB0[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(minSBoSB0[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(minSBoSB0[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(SBoSBMSY[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(SBoSBMSY[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(FoFMSY[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(FoFMSY[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(FoFtarg[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(FoFtarg[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(GK[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(GK[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(RK[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(RK[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(PrSBgt0.2SB0[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(PrSBgt0.2SB0[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(PrSBgtSBlim[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(PrSBgtSBlim[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(Y[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(Y[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(relCPUE[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(relCPUE[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(YoMSY[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(YoMSY[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(APCY[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(APCY[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles)),
                                              array(tapply(YcvPct[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(YcvPct[MP,], factors, quantile, probs=percentiles, names=FALSE, na.rm=T)), dim=dim_percentiles)),
                                              array(tapply(PrYlt0.1MSY[MP,], factors, mean), dim=dim_mean), t(array(unlist(tapply(PrYlt0.1MSY[MP,], factors, quantile, probs=percentiles, names=FALSE)), dim=dim_percentiles))))

          # results are organised one row per factor. Need to unwrap multiple rows into one.
          unwrapped <- NULL
          first     <- TRUE

          for (row_name in levels(sub_factor))
          {
            if (first)
            {
              unwrapped <- cbind(resultsByFactor[row_name,])
              first     <- FALSE
            }
            else
            {
              unwrapped <- cbind(unwrapped, resultsByFactor[row_name,])
            }
          }

          rownames(unwrapped) <- paste(MSEobj@MPs[MP],"y",pp,sep="")
          colnames(unwrapped) <- "x" %&% seq(1, dim_mean[1] * 15 * (1 + dim_percentiles[1]))

          ExtraResults <- rbind(data.frame(ExtraResults), unwrapped)
        }

        for (stat_name  in names(slot(MSEobj, paste(sub_name, "names", sep=""))))
        {
          colnames <- c(colnames,
                        "SBoSB0" %&% stat_name %&% c("mean",percentiles),
                        "minSBoSB0" %&% stat_name %&% c("mean",percentiles),
                        "SBoSBMSY" %&% stat_name %&% c("mean",percentiles),
                        "FoFMSY" %&% stat_name %&% c("mean",percentiles),
                        "FoFtarg" %&% stat_name %&% c("mean",percentiles),
                        "GK" %&% stat_name %&% c("mean",percentiles),
                        "RK" %&% stat_name %&% c("mean",percentiles),
                        "PrSBgt0.2SB0" %&% stat_name %&% c("mean",percentiles),
                        "PrSBgtSBlim" %&% stat_name %&% c("mean",percentiles),
                        # Y by area/fishery
                        "Y" %&% stat_name %&% c("mean",percentiles),
                        "relCPUE" %&% stat_name %&% c("mean",percentiles),
                        # CPUE by area/fishery
                        "YoMSY" %&% stat_name %&% c("mean",percentiles),
                        "APCY" %&% stat_name %&% c("mean",percentiles),
                        "YcvPct" %&% stat_name %&% c("mean",percentiles),
                        "PrYlt0.1MSY" %&% stat_name %&% c("mean",percentiles)
                        )
          }

        temp <- cbind(temp, ExtraResults)
      }
    }

    rownames(temp) <- paste(MSEobj@MPs[MPs],"y",pp,sep="")
    colnames(temp) <- colnames

    stats <- rbind(data.frame(stats), temp)
  } #pp

  colnames(stats) <- colnames

  return(stats)
}




# key summary stats by region and fishery
tableMSEbyRF.f <- function(object,percentiles=c(0.1,0.25,0.5,0.75,0.9))
{
  stats <- NULL
  MSEobj    <- object
  nsim      <- MSEobj@nsim
  npop      <- MSEobj@npop
  proyears  <- MSEobj@proyears
  allyears  <- MSEobj@proyears + MSEobj@nyears
  nMPs      <- MSEobj@nMPs
  ntargpop  <- length(MSEobj@targpop)
  firstMPy  <- MSEobj@nyears + MSEobj@firstMPYr - MSEobj@lastCalendarYr


  # projection period to report on (starting from the first MPyear and ignoring the bridging years)
  projPeriodList <- c(1,5,10,20,allyears-firstMPy) #WP Methods defined
  # projPeriodList <- c(1,3,5,10,20,allyears-firstMPy)

  YbyFStats <- NULL
  IbyRStats <- NULL

  for (pp in projPeriodList)
  {
    projPeriod   <-(firstMPy):(firstMPy + pp -1)
    projPeriodm1 <-(firstMPy-1):(firstMPy + pp -2)

    #11 mean Catch by fishery
    YbyF <-round(apply(as.karray(MSEobj@CMbyF)[,,keep(MSEobj@targpop),keep(projPeriod),], MARGIN=c(1:2,5), mean), 0)
    tmpF2 <- NULL

    for (fi in 1:MSEobj@nfleets)
    {
      tmpF1 <- cbind(apply(as.karray(YbyF)[,,fi], MARGIN=c(1),mean),       t(apply(as.karray(YbyF)[,,fi], MARGIN=c(1),quantile, probs=percentiles)))
      colnames(tmpF1) <- c("Yf" %&% fi %&% c("mean",percentiles))
      tmpF2 <- cbind(tmpF2,tmpF1)
    }

    rownames(tmpF2) <- paste(MSEobj@MPs,"y",pp,sep="")
    YbyFStats <- rbind(YbyFStats,tmpF2)

    #13 mean catch rates relative to catch rates over last data year by region and fishery
    IbyR <-round(apply(as.karray(MSEobj@IobsRArchive)[,,keep(projPeriod),], MARGIN=c(1:2,4), mean), 3)/
           round(apply(as.karray(MSEobj@IobsRArchive)[,,keep(MSEobj@nyears),], MARGIN=c(1:2,4), mean), 3)

    tmpI2 <- NULL

    for (ri in 1:MSEobj@nareas)
    {  # should replace nares with ncpue
      tmpI1 <- cbind(apply(as.karray(IbyR)[,,ri], MARGIN=c(1),mean),       t(apply(as.karray(IbyR)[,,ri], MARGIN=c(1),quantile, probs=percentiles)))
      colnames(tmpI1) <- c("CPUEr" %&% ri %&% c("mean",percentiles))
      tmpI2 <- cbind(tmpI2,tmpI1)
    }

    rownames(tmpI2) <- paste(MSEobj@MPs,"y",pp,sep="")

    IbyRStats <- rbind(IbyRStats,tmpI2)
  } #pp

    par(mfrow=c(2,2))

  for (MPi in 1:MSEobj@nMPs)
  {
      plot(MSEobj@IobsRArchive[MPi,1,,ri],type='l')

    for (ri in 1:MSEobj@nareas)
    {
        lines(MSEobj@IobsRArchive[MPi,1,,ri],type='l',col=ri)
      }
    }

  return(list(YbyFStats, IbyRStats))
}


#DK modified with karray for nMP=1
setMethod("summary",
  signature(object = "MSE"),
  function(object,bysim=F)
{
  MSEobj<-object
  nsim<-MSEobj@nsim
  proyears<-MSEobj@proyears
  allyears<-MSEobj@proyears+MSEobj@nyears
  nMPs<-MSEobj@nMPs
  ntargpop<-length(MSEobj@targpop)

  # projection period
    projPeriod <-(MSEobj@nyears+1):(MSEobj@nyears+MSEobj@proyears)
    CIind<-projPeriod-MSEobj@nyears-1

     #Probability in green Kobe
    GK  <-round(apply(as.karray(MSEobj@F_FMSY)[,,projPeriod]<1 & as.karray(MSEobj@B_BMSY)[,,projPeriod]>1,1:2,sum)/(MSEobj@proyears),3)
    GK10<-round(apply(as.karray(MSEobj@F_FMSY)[,,projPeriod]<1 & as.karray(MSEobj@B_BMSY)[,,projPeriod]>1,1:2,sum)/(MSEobj@proyears),3)

#Y<-round(apply(MSEobj@C[,,MSEobj@targpop,projPeriod],1:2,mean),0)
    Y<-round(apply(as.karray(MSEobj@CM)[,,MSEobj@targpop,projPeriod],1:2,mean),0)
#Dend <- MSEobj@B_BMSY[,,allyears,drop=F]
    Dend <- as.karray(MSEobj@B_BMSY)[,,allyears]

    Y5<-round(apply(array(MSEobj@CM[,,MSEobj@targpop,projPeriod],c(nMPs,nsim,ntargpop,proyears))*array(rep(0.95^CIind,each=nMPs*nsim*ntargpop),c(nMPs,nsim,ntargpop,proyears)),1:2,mean),0)
    Y10<-round(apply(array(MSEobj@CM[,,MSEobj@targpop,projPeriod],c(nMPs,nsim,ntargpop,proyears))*array(rep(0.9^CIind,each=nMPs*nsim*ntargpop),c(nMPs,nsim,ntargpop,proyears)),1:2,mean),0)
    Y_5<-round(apply(array(MSEobj@CM[,,MSEobj@targpop,projPeriod],c(nMPs,nsim,ntargpop,proyears))*array(rep(1.05^CIind,each=nMPs*nsim*ntargpop),c(nMPs,nsim,ntargpop,proyears)),1:2,mean),0)
    ind<-projPeriod[1:(length(projPeriod)-1)]
    ind2<-projPeriod[2:length(projPeriod)]
#AAVY<-apply(((MSEobj@C[,,MSEobj@targpop,ind2]-MSEobj@C[,,MSEobj@targpop,ind])^2)^0.5,1:2,mean)
    AAVY<-apply(((as.karray(MSEobj@CM)[,,MSEobj@targpop,ind2]-as.karray(MSEobj@C)[,,MSEobj@targpop,ind])^2)^0.5,1:2,mean)

  if (bysim)
  {
    list(GK,Dend,AAVY,Y,Y5,Y10,Y_5)
  }
  else
  {
    temp<-data.frame(cbind(apply(GK,1,mean),
                              apply(Dend,1,mean),
                              apply(AAVY,1,mean),
                              apply(Y,1,mean),
                              apply(Y5,1,mean),
                              apply(Y10,1,mean),
                              apply(Y_5,1,mean)),
                              row.names=MSEobj@MPs)

    names(temp)=c("Green Kobe","Final depletion","AAV Yield","Yield","Yield 5% DR", "Yield 10% DR", "Yield -5% DR")
    temp
  }
})


sumplot<-function(dat,field,adjv=c(1,1,1),pm=2,UB=10)
{
  perc=c(0.02,0.02,0.02)

  perc[3]<-perc[3]*pm
  col<-c("black","red","green","blue","orange","grey","purple","pink","brown")
  coln<-match(field,names(dat))
  levs<-unique(dat[,coln])

  mnam<-c("Yield (% Disc. Rate)","Prob. Green Kobe","Av Ann. Var. Yield")
  mind<-c(13,11,10)

  for (met in 1:length(mnam))
  {
    ymax<--1000
    xlim<-c(10000,-10000)

    for (i in 1:length(levs))
    {
      tdat<-dat[dat[,coln]==levs[i],mind[met]]
      tdat<-tdat[tdat<UB&tdat>-0.001]
      dd<-density(tdat,adj=adjv[met],from=0)
      xt<-quantile(tdat,c(perc[met]/2,1-perc[met]))
      xlim[1]<-min(xlim[1],xt[1])
      xlim[2]<-max(xlim[2],xt[2])
      ymax<-max(ymax, max(dd$y))
    }

    for (i in 1:length(levs))
    {
      tdat<-as.numeric(dat[dat[,coln]==levs[i],mind[met]])
      tdat<-tdat[tdat<UB&tdat>-0.001]

      if(i==1)plot(density(tdat,adj=adjv[met],from=0),ylim=c(0,ymax),xlim=xlim,col=col[1],type='l',main=mnam[met])

      if(i>1)lines(density(tdat,adj=adjv[met],from=0),col=col[i])
    }
  }

  legend('topright',legend=levs,text.col=col[1:length(levs)],bty='n')
}


sumplot2<-function(dat,fieldv,adjv=c(1,1,1),pm=2,UB=10,refMP="UMSY_PI")
{
  dat<-dat[dat$MP!=refMP,]
  perc=c(0.02,0.02,0.02)
  perc[3]<-perc[3]*pm
  col<-c("black","red","green","blue","orange","grey","purple","pink","brown")

  for (ff in 1:length(fieldv))
  {
  field<-fieldv[ff]
  coln<-match(field,names(dat))
  levs<-unique(dat[,coln])

  mnam<-c("Yield (5% Disc. Rate)","Av Ann. Var. Yield","Prob. Green Kobe")
  mind<-c(13,10,11)

    for (met in 1:length(mnam))
    {
    ymax<--1000
    xlim<-c(10000,-10000)

      for (i in 1:length(levs))
      {
      tdat<-dat[dat[,coln]==levs[i],mind[met]]
      tdat<-tdat[tdat<UB&tdat>-0.001]
      dd<-density(tdat,adj=adjv[met],from=0)
      xt<-quantile(tdat,c(perc[met]/2,1-perc[met]))
      xlim[1]<-min(xlim[1],xt[1])
      xlim[2]<-max(xlim[2],xt[2])
      ymax<-max(ymax, max(dd$y))
    }

      for (i in 1:length(levs))
      {
      tdat<-as.numeric(dat[dat[,coln]==levs[i],mind[met]])
      tdat<-tdat[tdat<UB&tdat>-0.001]

      if(i==1)plot(density(tdat,adj=adjv[met],from=0),ylim=c(0,ymax),xlim=xlim,xlab="",ylab="",col=col[1],type='l',main="")

      if(i>1)lines(density(tdat,adj=adjv[met],from=0),col=col[i])
    }

    if(ff==1)mtext(mnam[met],side=3,line=0.3)
  }

  legend('topright',legend=levs,text.col=col[1:length(levs)],bty='n')
  }

  mtext("Relative frequency",side=2,line=0.5,outer=T)
}


Tplot<-function(dat,field,refMP="UMSY_PI",legy=F,legpos='top')
{
  dat<-dat[dat$MP!=refMP,]
  col<-c("black","orange","blue","green")
  mnam<-c("Yield (no Disc. rate)","Yield (5% Disc. rate)",
          "Yield (10% Disc. rate)","Prob. Green Kobe",
          "Av Ann. Var. Yield")

  MPs<-unique(dat$MP)
  nMP<-length(MPs)

  coln<-match(field,names(dat))
  levs<-unique(dat[,coln])
  nlev<-length(levs)

  Y<-aggregate(dat$Y,by=list(dat$MP,dat[,coln]),FUN="mean")$x
  Y5<-aggregate(dat$Y5,by=list(dat$MP,dat[,coln]),FUN="mean")$x
  Y10<-aggregate(dat$Y10,by=list(dat$MP,dat[,coln]),FUN="mean")$x
  PGK<-aggregate(dat$PGK,by=list(dat$MP,dat[,coln]),FUN="mean")$x
  AAV<-aggregate(dat$AAV,by=list(dat$MP,dat[,coln]),FUN="mean")$x

  Yinc<-(max(Y)-min(Y))/25
  Y5inc<-(max(Y5)-min(Y5))/25
  Y10inc<-(max(Y10)-min(Y10))/25
  PGKinc<-(max(PGK)-min(PGK))/25
  AAVinc<-(max(AAV)-min(AAV))/25

  cols<-rep(col,each=nMP)

  plot(range(PGK)+c(-PGKinc,PGKinc),range(Y5)+c(-Y5inc,Y5inc),col='white',xlab="",ylab="")
  addgg(PGK,Y5)
  textplot(PGK,Y5,rep(MPs,nlev),col=cols,new=F)

  plot(range(AAV)+c(-AAVinc,AAVinc),range(Y5)+c(-Y5inc,Y5inc),col='white',xlab="",ylab="")
  addgg(AAV,Y5)
  textplot(AAV,Y5,rep(MPs,nlev),col=cols,new=F)

  legend(legpos,legend=levs,text.col=col[1:length(levs)],bg=makeTrans('white',97),box.col='white')
}


Tplot2<-function(dat,fieldv,legpos='top')
{
  for(ll in 1:length(fieldv))Tplot(dat,fieldv[ll],legpos)

  mtext("Yield relative to MSY (5% Disc. rate)",side=2,line=0.5,outer=T)
  mtext(c("Prob. green Kobe","AAVY"),side=1,at=c(0.25,0.75),line=0.8,outer=T)
}


addgg<-function(x,y,pcol='azure2')
{
  resx<-(max(x)-min(x))/10
  resy<-(max(y)-min(y))/10
  xlim<-c(min(x)-(2*resx),max(x)+(2*resx))
  ylim<-c(min(y)-(2*resy),max(y)+(2*resy))

  divx<-pretty(seq(xlim[1],xlim[2],length.out=20))
  divy<-pretty(seq(ylim[1],ylim[2],length.out=20))

  polygon(c(xlim,xlim[2:1]),rep(ylim,each=2),col=pcol)
  abline(v=divx,col='white')
  abline(h=divy,col='white')
}


makeTrans<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# Plot performance summary of the mse object
#setMethod("summary",
#  signature(object = "MSE"),
stats<-function(object)
{
  nsim<-object@nsim
  nyears<-object@nyears
  proyears<-object@proyears
  nMPs<-object@nMPs
  targpop<-object@targpop

  CM<-apply(array(object@CM[,,targpop,(nyears+1):(nyears+proyears)],
                 c(nMPs,nsim,length(targpop),proyears)),c(1,2,4),sum)

  AAVY<-array(NA,c(nMPs,nsim,proyears-1))
  ind1<-as.matrix(expand.grid(1:nMPs,1:nsim,1:(proyears-1)))
  ind2<-as.matrix(expand.grid(1:nMPs,1:nsim,2:proyears))
  AAVY[ind1]<-((CM[ind1]-C[ind2])^2)^0.5
  AAVY<-apply(AAVY,1:2,mean)

  Y<-apply(CM[,,(proyears-4):proyears],1:2,mean)

  F_FMSY<-object@F_FMSY[,,(nyears+1):(nyears+proyears)]
  B_BMSY<-object@B_BMSY[,,(nyears+1):(nyears+proyears)]
  Pgreen<-apply(array(as.integer(F_FMSY<1&B_BMSY>1),dim(B_BMSY)),1:2,mean)

  list("Y"=Y,"AAVY"=AAVY,"Pgreen"=Pgreen,"Dep"=B_BMSY[,,proyears],"CM"=CM,"F_FMSY"=F_FMSY,"B_BMSY"=B_BMSY)
}#)


anim8mov<-function(.Object)
{
  pref<-.Object@npop:1
  mref<-c(2,3,4,1)
  Pdist<-array(NA,c(.Object@npop,2,.Object@nsubyears,.Object@nareas))

  temp<-.Object@mov[1,,1,,,]
  Mtemp<-.Object@Mmov[1,,20,,,]
  MIdist<-Idist<-.Object@excl/apply(.Object@excl,1,sum)

  for (i in 1:100)
  {
    for(mm in 1:4)
    {
    Idist<-domov3(Idist,temp[,mm,,])
    MIdist<-domov3(MIdist,Mtemp[,mm,,])

    if(i==100)Pdist[,1,mref[mm],]<-Idist
    if(i==100)Pdist[,2,mref[mm],]<-MIdist
    }
  }

  fac<-40
  nint<-(4-1)*(fac+1)+1
  Pint<-array(NA,c(.Object@npop,2,nint,.Object@nareas))

  for(pp in 1:.Object@npop)for(mat in 1:2)Pint[pp,mat,,]<-tinter(Pdist[pp,mat,,],fac)

  for (i in 1:nint)
  {
    par(mfrow=c(2,.Object@npop),mai=c(0.1,0.1,0.01,0.01),omi=c(0.01,0.4,0.6,0.05))

    for (mat in 1:2)
    {
      for (pp in 1:.Object@npop)
      {
        sdensplot(Pint[pref[pp],mat,i,],.Object@Area_defs)
        if(pp==.Object@npop&mat==1)mtext(paste("DOY =" ,ceiling(i/nint*365)),3,line=1,adj=0.95)
        if(pp==1&mat==1)mtext("Juvenile",2,line=1)
        if(pp==1&mat==2)mtext("Mature",2,line=0.3)
        if(mat==1)mtext(paste("Population",pp),3,line=1)
        ani.pause()
      }
    }
  }
}


domov3<-function(Ntemp,movtemp)
{ # P R  x  P R R
  #Ntemp<-Idist
  #movtemp<-.Object@mov[,,,1,1,,]
  nareas<-dim(movtemp)[3]
  apply(array(Ntemp, c(dim(Ntemp),nareas))*movtemp,c(1,3),sum)
}


tinter<-function(grid,fac=2)
{
  nr<-nrow(grid)
  nr2<-dim(grid)[1]*(fac+2)-fac
  nc<-ncol(grid)
  multi1<-c(rep(seq(1,0,length.out=fac+2)[1:(fac+1)],nr-1),0)
  multi2<-1-multi1
  ind<-rep(1:nr,each=fac+1)
  ind1<-rep(c(ind[1:(length(ind)-fac-1)],nr-1),nc)
  ind2<-rep(c(ind[(2+fac):length(ind)],nr),nc)
  indy<-rep(1:nc,each=length(ind1)/nc)

  array(grid[cbind(ind1,indy)]*multi1+
          grid[cbind(ind2,indy)]*multi2,c(length(ind)-fac,nc))
}
