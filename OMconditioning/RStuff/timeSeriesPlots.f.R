# Does some aggregate weighted plots
# Kobe plots (and matrix) removed, because they are broken in the demo context


timeSeriesPlots.f <- function (mList = gridY1List, doProj = F, doLegend = F, plotRefCase=F, opt = c("h55",
    "h75", "h95", "Mest", "L69", "U1", "r0", "rqs", "t1", "t2",
    "t3", "od02", "od20", "od70", "ss", "sa", "rttp", "rtss",
    "CL2"), optWt = c(0.33, 0.33, 0.33, 1, 1, 1, 0.5, 0.5, 0.33,
    0.33, 0.33, 0.33, 0.33, 0.33, 0.5, 0.5, 0.5, 0.5, 1), SSYr=0.25, yrSeasY=272,
    CIupper=0.05,CIlower=0.95, CI,
    seasAsYr=101, endSeasAsYr=272, numSeas=4, endYr=2014, endSeas=4, xlims=c(1950,2014))
{
    factorList <- array(NA, dim = c(length(mList), length(unlist(strsplit(mList[1],
        split = "_")))))
    modWtList <- c(1:length(mList)) * 0
    for (mi in 1:length(mList)) {
        factors <- strsplit(mList[mi], split = "_")
        factorList[mi, ] <- unlist(factors)
        modWt <- 1
        for (o in 1:length(opt)) {
            if (opt[o] %in% factorList[mi, ])
                modWt <- modWt * optWt[o]
        }
        modWtList[mi] <- modWt
    }
    print(cbind(unlist(mList), modWtList))
    m <- get(mList[1])

#SSBY <- m$derived_quants[m$derived_quants$LABEL == 'SPB_272',]$Value

    SPBSeries <- m$derived_quants[substring(m$derived_quants$LABEL, 1, 3) == "SPB", ]
    SPBSeries <- SPBSeries[3:nrow(SPBSeries), ]
    numTime <- nrow(SPBSeries)
# translate seasons as years back to readable time
# m[[7]] extraction probably not intended for multiple areas...
#year = 1947 + SSYr/4 - 0.125 (maybe off by one?)
    firstQ <- m[[7]]$Seas == 1
    TBYear <- seasAsYrToDecYr.f(m[[7]]$Yr[firstQ],endSeasAsYr=endSeasAsYr, numSeas=numSeas, endYr=endYr, endSeas=endSeas)
    TBSeries <- m[[7]]$Bio_all[firstQ]
    #recYrQtr <- m[[7]]$Yr + (m[[7]]$Seas - 1)/length(unique(get(mList[1])[[7]]$Seas))
    #recYrQtr <- seasAsYrToDecYr.f(m[[7]]$Yr)
    #recYrQtr <- seasAsYrToDecYr.f(seasAsYr=m[[7]]$Yr, endSeasAsYr=384, numSeas=4, endYr=2022, endSeas=4)
    recYrQtr <- seasAsYrToDecYr.f(seasAsYr=m[[7]]$Yr, endSeasAsYr=endSeasAsYr, numSeas=numSeas, endYr=endYr, endSeas=endSeas)
browser()
    recSeries <- m[[7]]$Recruit_0
    plot(TBYear, TBSeries)
    plot(recYrQtr, recSeries)
    SPBMat <- array(NA, dim = c(length(mList), numTime))
    SPBMatSeries <- array(NA, dim = c(numTime, 6))
    TBMat <- array(NA, dim = c(length(mList), length(TBYear)))
    TBMatSeries <- array(NA, dim = c(length(TBYear), 6))
    recMat <- array(NA, dim = c(length(mList), length(recYrQtr)))
    recMatSeries <- array(NA, dim = c(length(recYrQtr), 6))
    SPBoSPBMSYSeries <- SPBSeries
    SPBoSPBMSYMat <- array(NA, dim = c(length(mList), numTime))
    SPBoSPBMSYMatSeries <- array(NA, dim = c(numTime, 6))
    SPBoSPB0Series <- SPBSeries
    SPBoSPB0Mat <- array(NA, dim = c(length(mList), numTime))
    SPBoSPB0MatSeries <- array(NA, dim = c(numTime, 6))
    FMat <- SPBMat
    FMatSeries <- SPBMatSeries
    MSYMat <- array(NA, dim = c(length(mList)))
    SPBY <- FY <- NULL
    for (mi in 1:length(mList)) {
        print(mList[mi])
        m <- get(mList[mi])
        SPBSeries <- m$derived_quants[substring(m$derived_quants$LABEL, 1, 3) == "SPB", ]

        SPBSeries <- SPBSeries[3:nrow(SPBSeries), ]
        yearLab <- seasAsYrToDecYr.f(as.numeric(substring(SPBSeries$LABEL, 5)),endSeasAsYr=endSeasAsYr, numSeas=numSeas, endYr=endYr, endSeas=endSeas)
        TBSeries <- m[[7]]$Bio_all[firstQ]
        recSeries <- m[[7]]$Recruit_0
        SPBoSPBMSYSeries[, 2] <- SPBSeries[, 2]/as.numeric(subset(m$derived_quants,
            m$derived_quants$LABEL == "SSB_MSY")$Value)
        SPBoSPB0Series[, 2] <- SPBSeries[, 2]/as.numeric(subset(m$derived_quants,
            m$derived_quants$LABEL == "SPB_Virgin")$Value)
        FSeries <- m$derived_quants[substring(m$derived_quants$LABEL, 1, 2) == "F_", ]
        # note F_std is already relative to FMSY
        #FSeries[, 2] <- FSeries[, 2]/as.numeric(subset(m$derived_quants, m$derived_quants$LABEL == "Fstd_MSY")$Value)
        TBMat[mi, ] <- TBSeries
        recMat[mi, ] <- recSeries
        SPBMat[mi, ] <- SPBSeries[, 2]
        SPBoSPBMSYMat[mi, ] <- SPBoSPBMSYSeries[, 2]
        SPBoSPB0Mat[mi, ] <- SPBoSPB0Series[, 2]
        FMat[mi, ] <- FSeries[, 2]
        #MSY multipled by 4 quarters
        MSYMat[mi] <- (m$derived_quants[m$derived_quants$LABEL ==
            "TotYield_MSY", ]$Value)/SSYr
        SPBY <- c(SPBY, as.numeric(m$derived_quants[substring(m$derived_quants$LABEL,
            1, 8) == "SPB_" %&% yrSeasY, 2]))
        #F averaged over 4 quarters
        if(SSYr==0.25){
          FY <- c(FY, mean(as.numeric(m$derived_quants[substring(m$derived_quants$LABEL,1, 8) %in% c("F_269","F_270","F_271","F_272"), 2])))
        }
        if(SSYr!=0.25){
          print("not set up for this time structure")
          break
        }



    }
    par(mfrow = c(1, 1))
    hist(SPBY)
    hist(log(SPBY))
    hist(FY)
    hist(log(FY))
    print(sd(SPBY))
    print(sd(log(SPBY)))
    print(sd(SPBY))
    print(sd(log(FY)))
    for (t in 1:numTime) {
        SPBoSPBMSYMatSeries[t, 1] <- weighted.mean(SPBoSPBMSYMat[,
            t], w = modWtList)
        SPBoSPBMSYMatSeries[t, 2] <- min(SPBoSPBMSYMat[, t])
        SPBoSPBMSYMatSeries[t, 3] <- max(SPBoSPBMSYMat[, t])
        tmpSort <- sort(SPBoSPBMSYMat[, t])
        modWtSort <- modWtList[order(SPBoSPBMSYMat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        SPBoSPBMSYMatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        SPBoSPBMSYMatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        SPBoSPBMSYMatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
        SPBoSPB0MatSeries[t, 1] <- weighted.mean(SPBoSPB0Mat[,
            t], w = modWtList)
        SPBoSPB0MatSeries[t, 2] <- min(SPBoSPB0Mat[, t])
        SPBoSPB0MatSeries[t, 3] <- max(SPBoSPB0Mat[, t])
        tmpSort <- sort(SPBoSPB0Mat[, t])
        modWtSort <- modWtList[order(SPBoSPB0Mat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        SPBoSPB0MatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        SPBoSPB0MatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        SPBoSPB0MatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
        SPBMatSeries[t, 1] <- weighted.mean(SPBMat[, t], w = modWtList)
        SPBMatSeries[t, 2] <- min(SPBMat[, t])
        SPBMatSeries[t, 3] <- max(SPBMat[, t])
        tmpSort <- sort(SPBMat[, t])
        modWtSort <- modWtList[order(SPBMat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        SPBMatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        SPBMatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        SPBMatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
        FMatSeries[t, 1] <- weighted.mean(FMat[, t], w = modWtList)
        FMatSeries[t, 2] <- min(FMat[, t])
        FMatSeries[t, 3] <- max(FMat[, t])
        tmpSort <- sort(FMat[, t])
        modWtSort <- modWtList[order(FMat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        FMatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        FMatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        FMatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
        if (yearLab[t] == yrSeasY) {
            print("BY/BMSY, BY/B0, and FY / FMSY  mean, range, 5%, 50%, 95%")
            print(SPBoSPBMSYMatSeries[t, ])
            print(SPBoSPB0MatSeries[t, ])
            print(FMatSeries[t, ])
        }
    }
    print("SPB/SPBMSY weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(yearLab, SPBoSPBMSYMatSeries))
    print("SPB/SPB0 weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(yearLab, SPBoSPB0MatSeries))
    print("SP Biomass weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(yearLab, SPBMatSeries))
    print("F weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(yearLab, FMatSeries))
    for (t in 1:length(TBYear)) {
        TBMatSeries[t, 1] <- weighted.mean(TBMat[, t], w = modWtList)
        TBMatSeries[t, 2] <- min(TBMat[, t])
        TBMatSeries[t, 3] <- max(TBMat[, t])
        tmpSort <- sort(TBMat[, t])
        modWtSort <- modWtList[order(TBMat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        TBMatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        TBMatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        TBMatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
        if (TBYear[t] == yrSeasY) {
        }
    }
    print("Total Biomass weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(TBYear, TBMatSeries))
    for (t in 1:length(recYrQtr)) {
        recMatSeries[t, 1] <- weighted.mean(recMat[, t], w = modWtList)
        recMatSeries[t, 2] <- min(recMat[, t])
        recMatSeries[t, 3] <- max(recMat[, t])
        tmpSort <- sort(recMat[, t])
        modWtSort <- modWtList[order(recMat[, t])]
        tmpcum <- cumsum(modWtSort)/sum(modWtSort)
        recMatSeries[t, 4] <- tmpSort[tmpcum > CIlower][1]
        recMatSeries[t, 5] <- tmpSort[tmpcum > 0.5][1]
        recMatSeries[t, 6] <- tmpSort[tmpcum > CIupper][1]
    }
    print("Rec weighted.mean, min, max, 5%, 50%, 95%")
    print(cbind(recYrQtr, recMatSeries))
    print("MSY")
    print(c(weighted.mean(MSYMat, w = modWtList), range(MSYMat)))
    MSYSort <- sort(MSYMat)
    modWtSort <- modWtList[order(MSYMat)]
    MSYcum <- cumsum(modWtSort)/sum(modWtSort)
    print(c(MSYSort[MSYcum > 0.05][1], MSYSort[MSYcum > 0.5][1],
        MSYSort[MSYcum > 0.95][1]))

    if(plotRefCase){
      #reference series for plotting
      refSPBSeries        <- ref$derived_quants[substring(ref$derived_quants$LABEL, 1, 3) == "SPB", ]
      refSPBSeries        <- refSPBSeries[3:nrow(refSPBSeries), ]
      refSPBoSPBMSYSeries <- refSPBSeries[, 2]/as.numeric(subset(ref$derived_quants, ref$derived_quants$LABEL == "SSB_MSY")$Value)
      refSPBoSPB0Series   <- refSPBSeries[, 2]/as.numeric(subset(ref$derived_quants, ref$derived_quants$LABEL == "SPB_Virgin")$Value)
      refFSeries          <- ref$derived_quants[substring(ref$derived_quants$LABEL, 1, 2) == "F_", ]
    }

    plot(yearLab, SPBoSPBMSYMatSeries[, 5], ylim = c(0, 7), main = "SSB/SSBMSY", xlim=xlims,
        ylab = "SSB/SSBMSY", xlab = "", type = "l", col = 1,
        lwd = 2)
    lines(yearLab, SPBoSPBMSYMatSeries[, 4], col = 1, lty = 1,
        lwd = 1)
    lines(yearLab, SPBoSPBMSYMatSeries[, 6], col = 1, lty = 1, lwd = 1)
    lines(yearLab, SPBoSPBMSYMatSeries[, 2], col = 1, lty = 2, lwd = 1)
    lines(yearLab, SPBoSPBMSYMatSeries[, 3], col = 1, lty = 2, lwd = 1)
    abline(h=1, col = 2, lty = 1)
    if(plotRefCase)lines(yearLab, refSPBoSPBMSYSeries,      col = 3, lty = 1, lwd = 2)


    plot(yearLab, SPBoSPB0MatSeries[, 5], ylim = c(0, 2), main = "SSB/SSB0",  xlim=xlims,
        ylab = "SSB/SSB0", xlab = "", type = "l", col = 1, lwd = 2)
    lines(yearLab, SPBoSPB0MatSeries[, 4], col = 1, lty = 1,
        lwd = 1)
    lines(yearLab, SPBoSPB0MatSeries[, 6], col = 1, lty = 1, lwd = 1)
    lines(yearLab, SPBoSPB0MatSeries[, 2], col = 1, lty = 2, lwd = 1)
    lines(yearLab, SPBoSPB0MatSeries[, 3], col = 1, lty = 2, lwd = 1)
    abline(h=1, col = 2, lty = 1)
    abline(h=0.2, col = 2, lty = 2)
    abline(h=0.4, col = 2, lty = 2)
    if(plotRefCase)lines(yearLab, refSPBoSPB0Series,      col = 3, lty = 1, lwd = 2)


    plot(yearLab, FMatSeries[, 5], ylim = c(0, 3), main = "F/FMSY", xlim=xlims,
        ylab = "F/FMSY", xlab = "", type = "l", col = 1, lwd = 2)
    lines(yearLab, FMatSeries[, 4], col = 1, lty = 1, lwd = 1)
    lines(yearLab, FMatSeries[, 6], col = 1, lty = 1, lwd = 1)
    lines(yearLab, FMatSeries[, 2], col = 1, lty = 2, lwd = 1)
    lines(yearLab, FMatSeries[, 3], col = 1, lty = 2, lwd = 1)
    abline(h=1, col = 2, lty = 1)
    if(plotRefCase)lines(yearLab, refFSeries[,2],      col = 3, lty = 1, lwd = 2)

return(-999)
}

