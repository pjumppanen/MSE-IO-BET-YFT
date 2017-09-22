# plot.f.R - plotting routine for IOTC bigeye and yellowfin MSE results
# the essence of the stockPlot() function is derived from:

# ggplotFL/R/plot.R
# Copyright 2003-2014 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC, Laurie Kell, ICCAT

# Dale Kolody modified to conform with IOTC WP Methods graphics specifications, including:
# -identification and plotting of worms corresponding to 25th, 50th and 75th %ile of terminal SB/SBMSY
# - fishery and area-specific plots
# - multi-MP trade-off plots
# - udon-soba-style boxplots

library(ggplot2)
library(reshape2)


#Confidence interval plots (double whisker aka udon-soba plots).
# median, thick confidence interval whiskers for the 25-75th percentiles, and thin whiskers indicating 10-90th percentiles
plotUS.f <- function(mseObj, MPs=c("CE1.0n","CC400n","ITarg2.50","PT4010.1.0"), plotStats="def", param="")
{
  suffix <- if (nchar(param) > 0) "(" %&% param %&% ")" else ""

  fullPIList <- c('SBoSB0',
                  'minSBoSB0',
                  'SBoSBMSY',
                  'FoFMSY',
                  'FoFtarg',
                  'GK',
                  'RK',
                  'PrSBgt0.2SB0',
                  'PrSBgtSBlim',
                  'Y',
                  'relCPUE',
                  'YoMSY',
                  'APCY',
                  'YcvPct',
                  'PrYlt0.1MSY')

  fullPILabels <- c('SB/SB0',
                    'min(SB/SB0)',
                    'SB/SB(MSY)',
                    'F/FMSY',
                    'F/F(target)',
                    'Prob. in Green Kobe',
                    'Prob. in Red Kobe',
                    'Pr( SB > 0.2 SB0)',
                    'Pr( SB > SB(limit))',
                    'Catch (1000t)',
                    'CPUE/CPUE(YSA)',
                    'Catch / MSY',
                    'Average Percent Change in Catch',
                    'Catch CV (%)',
                    'Pr(Catch < 0.1 MSY)')

  if (plotStats == 'def')
  {
    PIList <- fullPIList %&% param
    PILabels <- fullPILabels %&% suffix
  }
  else
  {
    PIList <- plotStats %&% param
    PILabels <- fullPILabels[fullPIList %in% plotStats] %&% suffix
  }

  ms <- tableMSE.f(mseObj)

  t <- c(1,5,10,20) #summary periods (years)

  par(mfrow=c(4,3))
  par(mar=c(6, 3, 2.5, 0) + 0.1)

  for (ti in t)
  {
    for (paneli in 1:length(PIList))
    {
      PI50 <- PI10 <- PI25 <- PI75 <- PI90 <- NULL

      for (mpi in 1:length(MPs))
      {
        PI50 <- c(PI50,(ms[MPs[mpi] %&% "y" %&% ti, PIList[paneli] %&% "0.5"]))
        PI10 <- c(PI10,(ms[MPs[mpi] %&% "y" %&% ti, PIList[paneli] %&% "0.1"]))
        PI25 <- c(PI25,(ms[MPs[mpi] %&% "y" %&% ti, PIList[paneli] %&% "0.25"]))
        PI75 <- c(PI75,(ms[MPs[mpi] %&% "y" %&% ti, PIList[paneli] %&% "0.75"]))
        PI90 <- c(PI90,(ms[MPs[mpi] %&% "y" %&% ti, PIList[paneli] %&% "0.9"]))
      }

      if (sum(is.na(c(PI50, PI10, PI25, PI75,PI90))) == 0)
      {
        plot(PI50, xlab='', ylab='',pch=3, xaxt='n',col=1:length(PI50),cex=2, cex.main=0.8, ylim=c(0,max(PI90)),
             main=PILabels[paneli] %&% " over " %&% ti %&% " y \n" %&% 'OM ' %&% mseObj@Label)
        axis(side=1, labels=MPs, at=c(1:length(MPs)),las=2)

        for (mpi in 1:length(MPs))
        {
          lines(c(mpi,mpi), c(PI25[mpi],PI75[mpi]), lwd=3,col=mpi)
          lines(c(mpi,mpi), c(PI10[mpi],PI90[mpi]), lwd=1,col=mpi)
        }
      }
    }
  }
}


# WPMethods bivariate trade-off plots
#Four core trade-off plots, computed for each of 10 and 20 years of projection (i.e. year 1 = first TAC/TAE implementation)
# 1. SB/SBMSY (or SB/SB0 for skipjack) vs. Yield
# 2. Pr(Green Kobe) vs. Yield
# 3. Pr(SB > BLim) vs. Yield
# 4. mean(1 – Cy/Cy-1) vs. Yield
plotTO.f <- function(mseObj=MSEList, MPs=MPList, ylims='def', xlims='def', param='')
{
  numOMs <- length(mseObj)
  ms <- tableMSE.f(mseObj[[1]])

  suffix <- if (nchar(param) > 0) "(" %&% param %&% ")" else ""
  y <- c("SBoSBMSY","GK","PrSBgtSBlim","APCY") %&% param
  ylabs  <- c("SB / SB(MSY)","Probability","Probability","mean(|1-C(y)/C(y-1|)")
  titles <- c("Spawning Biomass","Prob. in Green Kobe Quad.","Prob. SB > SB(limit)","Mean Change in Yield ratio") %&% suffix

  if (ylims=='def') ylims <- c(4,1,1,10)

  x <- rep("Y" %&% param, length(y))

  sel <- "Y" %&% param %&% "0.9"

  if (xlims=='def') xlims <- rep(max(ms[[sel]]),length(x))

  t <- c(10,20) #summary periods (years)

  for(ti in t)
  {
#    dev.new(width=6, height=4)
    par(mfrow=c(1,1),mar=c(4, 3, 2, 0) + 0.1, mgp=c(2,1,0))
#    par(mfrow=c(2,2),mar=c(4, 3, 2, 0) + 0.1, mgp=c(2,1,0))

    for(paneli in 1: length(y))
    {
      if (numOMs==1)
      {
        mainTitle <- titles[paneli] %&% " over " %&% ti %&% " y" %&% suffix %&% " \n" %&% 'OM ' %&% mseObj[[1]]@Label
      }
      else
      {
        mainTitle <- titles[paneli] %&% " over " %&% ti %&% " y" %&% suffix
      }

      plot(ms[MPs[1] %&% "y" %&% ti, x[paneli] %&% "0.5"], ms[MPs[1] %&% "y" %&% ti,y[paneli] %&% "0.5"],type="n",
           xlim=c(0,xlims[paneli]), ylim=c(0,ylims[paneli]),
           ylab=ylabs[paneli], xlab='Yield (1000t)',cex.main=0.8,
           main=mainTitle)
           #main=titles[paneli] %&% " over " %&% ti %&% " y \n" %&% 'OM ' %&% mseObj@Name)

      OMNames <- NULL

      for (msei in 1:length(mseObj))
      {
        ms <- tableMSE.f(mseObj[[msei]])
        OMNames <- c(OMNames, mseObj[[msei]]@Label)

        for (mpi in 1:length(MPs))
        {
          #print(c(ti,paneli,mpi))
          #print(c(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"], ms[MPs[mpi] %&% "y" %&% ti, y[paneli]] ))
          #print(c(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.1"], ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.9"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.5"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.5"]))
          #print(c(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"], ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.1"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.9"]))

          if (numOMs == 1)
          {
            colour <- mpi
          }
          else
          {
            colour <- msei
          }

          if (numOMs > 1)
          {
            lineWid   <- numOMs+1 - msei
            pointSize <- 1+(numOMs+1 - msei)/5
          }
          else
          {
            lineWid   <- 1
            pointSize <- 1
          }

          points(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"], ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.5"], pch=mpi,col=colour)
          lines(c(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.1"],ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.9"]), c(ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.5"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.5"]),col=colour, lwd=lineWid, cex=pointSize)
          lines(c(ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"],ms[MPs[mpi] %&% "y" %&% ti, x[paneli] %&% "0.5"]), c(ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.1"],ms[MPs[mpi] %&% "y" %&% ti,y[paneli] %&% "0.9"]),col=colour, lwd=lineWid, cex=pointSize)
        }
      }
    }

    dev.new(width=6, height=1.75)        # set the graphics window size (inches)
    par(mfrow=c(1,1), mai=c(0,0,0,0))
    plot(1,1, type="n", axes=F, xlab="", ylab="",xlim=c(0,1),ylim=c(0,1))

    if (numOMs == 1)
    {
      markers <- colour <- 1:length(MPs)
      legLabs <- MPs
    }
    else
    {
      colour  <- rep(1:numOMs, each=length(MPs))
      markers <- rep(1:length(MPs), times=numOMs)
      legLabs <- rep(OMNames, each=length(MPs)) %&% " " %&% rep(MPs, times=length(numOMs))
    }

    legend(x=0,y=1, legend=rep(legLabs), col=colour,pch=markers, ncol=3,cex=0.8, bty='n')
    par(mar=c(5, 4, 4, 2) + 0.1)
  }
}


# wrapper to provide labels etc.
plotTS.f <- function (mseObj = tmse, doWorms=T, wormProbs=c(0.25, 0.50, 0.75), plotByRF=T, mwgPlots=F, param="")
{
  suffix <- ""
  nSims  <- 1:mseObj@nsim

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

  par(mfrow=c(5,1))

  for (iMP in 1:mseObj@nMPs)
  {
    projRefLine <- c(mseObj@lastCalendarYr, mseObj@firstMPYr) # + 1

    #SSB/SSBMSY (includes finding reference worms)
    rawPlotDat <- mseObj@SSB_SSBMSY[iMP,nSims,] #nMPs,nsim,years
    plotDat <- as.data.frame(array(NA, dim=c(prod(dim(rawPlotDat)),3)))
    names(plotDat) <- c("year","iter","data")
    plotDat$year <- rep(mseObj@yrLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
    plotDat$iter <- rep(1:dim(rawPlotDat)[1], times = dim(rawPlotDat)[2]) #nsim nyears
    plotDat$data <- as.vector(rawPlotDat)

    #find sim numbers based on terminal SSB/SSBMSY percntiles to plot as worms
    if (doWorms)
    {
      xlast <- plotDat[plotDat$year==max(plotDat$year),]
      xlast$data <- rank(xlast$data)/length(xlast$data)

      wormSims <- NULL

      for (iprobs in wormProbs)
      {
        wormSims <- c(wormSims, which(abs(xlast$data-iprobs)==min(abs(xlast$data-iprobs)))[1])
      }
    }

    #SSB/SSBMSY
    p2b <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&% ": SSB/SSBMSY" %&% suffix, xlab="", ylab="SSB/SSBMSY",xref=projRefLine,yreflim=mseObj@SBlim,yreftarget=1., doWorms=doWorms, wormSims=wormSims)
    plot(p2b)

    #aggregate CPUE series skips final year
    rawPlotDat <- mseObj@IobsArchive[iMP,nSims,] #nMPs,nsim,years
    plotDat <- as.data.frame(array(NA, dim=c(prod(dim(rawPlotDat)),3)))
    names(plotDat) <- c("year","iter","data")
    #plotDat$year <- rep(1:dim(rawPlotDat)[2], each  = dim(rawPlotDat)[1]) #nyears, nsim
    #plotDat$year <- mseYrToDecYr.f(plotDat$year)
    plotDat$year <- rep(mseObj@yrLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
    plotDat$iter <- rep(1:dim(rawPlotDat)[1], times = dim(rawPlotDat)[2]) #nsim nyears
    plotDat$data <- as.vector(rawPlotDat)
    p6 <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&%": CPUE (aggregate)" %&% suffix, xlab="Year", ylab="CPUE",xref=projRefLine, yref=0, doWorms=doWorms, wormSims=wormSims)
    #plot(p6)

    rawPlotDat <- mseObj@Rec[iMP,nSims,] #nMPs,nsim,years
    plotDat <- as.data.frame(array(NA, dim=c(prod(dim(rawPlotDat)),3)))
    names(plotDat) <- c("year","iter","data")
#    plotDat$year <- rep(1:dim(rawPlotDat)[2], each  = dim(rawPlotDat)[1]) #nyears, nsim
#    plotDat$year <- mseYrToDecYr.f(plotDat$year)
    plotDat$year <- rep(mseObj@yrLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
    plotDat$iter <- rep(1:dim(rawPlotDat)[1], times = dim(rawPlotDat)[2]) #nsim nyears

    #Rec
    plotDat$data <- as.vector(rawPlotDat)/1E+3
    p1 <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&% ": Recruitment" %&% suffix, xlab="", ylab="Millions",xref=projRefLine, doWorms=doWorms, wormSims=wormSims)
    #plot(p1)

    #B/B0
    rawPlotDat <- mseObj@B_B0[iMP,nSims,1,] #nMPs,nsim,npop,years
    plotDat$data <- as.vector(rawPlotDat)
    p2 <- plotStock.f(plotDat, main="MP " %&%mseObj@MPs[iMP]  %&%": B/B0" %&% suffix, xlab="", ylab="B/B0",xref=projRefLine, ,yreflim=0.2, doWorms=doWorms, wormSims=wormSims)

    #B/BMSY
    rawPlotDat <- mseObj@B_BMSY[iMP,nSims,] #nMPs,nsim,years
    plotDat$data <- as.vector(rawPlotDat)
    p2.2 <- plotStock.f(plotDat, main="MP " %&%mseObj@MPs[iMP]  %&%": B/BMSY" %&% suffix, xlab="", ylab="B/BMSY",xref=projRefLine, doWorms=doWorms, wormSims=wormSims)

    #SSB/SSB0
    rawPlotDat <- mseObj@SSB_SSB0[iMP,nSims,1,] #nMPs,nsim,npop,years
    plotDat$data <- as.vector(rawPlotDat)
    p2.1 <- plotStock.f(plotDat, main="MP " %&%mseObj@MPs[iMP]  %&%": SSB/SSB0" %&% suffix, xlab="", ylab="SSB/SSB0",xref=projRefLine, ,yreflim=0.2, doWorms=doWorms, wormSims=wormSims)

    #SSB/SSBMSY
    #rawPlotDat <- mseObj@SSB_SSBMSY[iMP,nSims,] #nMPs,nsim,years
    #plotDat$data <- as.vector(rawPlotDat)
    #p2b <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&% ": SSB/SSBMSY" %&% suffix, xlab="", ylab="SSB/SSBMSY",xref=projRefLine,yreflim=0.4,yreftarget=1., doWorms=doWorms, wormSims=wormSims)
    #plot(p2b)

    #F/FMSY
    rawPlotDat <- mseObj@F_FMSY[iMP,nSims,] #nMPs,nsim,years
    rawPlotDat[rawPlotDat>5] <- 5
    plotDat$data <- as.vector(rawPlotDat)
    p3 <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&% ": F/FMSY" %&% suffix, xlab="", ylab="F/FMSY",xref=projRefLine,yreflim=mseObj@Flim,yreftarget=1., doWorms=doWorms, wormSims=wormSims)
    #plot(p3)

    #Catch
    rawPlotDat <- mseObj@CM[iMP,nSims,1,] #nMPs,nsim,years
    plotDat$data <- as.vector(rawPlotDat)/1000
    p4 <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&%": Catch" %&% suffix, xlab="Year", ylab="1000 t",xref=projRefLine, doWorms=doWorms, wormSims=wormSims)
    #plot(p4)

    #multiplot(p1, p2, p3, p4, cols=1)

    # Rec by Qtr
    qPlotDat <- mseObj@RecYrQtr[iMP,nSims,] #nMPs,nsim,years
    plotDat <- as.data.frame(array(NA, dim=c(prod(dim(qPlotDat)),3)))
    names(plotDat) <- c("year","iter","data")
    #plotDat$year <- rep(1:dim(qPlotDat)[2], each  = dim(qPlotDat)[1])
    #plotDat$year <- mseYrSeasToDecYrSeas.f(plotDat$year)
    plotDat$year <- rep(mseObj@yrSeasLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
    plotDat$iter <- rep(1:dim(qPlotDat)[1], times = dim(qPlotDat)[2]) #nsim nyears*nsubyears
    plotDat$data <- as.vector(qPlotDat)/1E+3
    p5 <- plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&% ": Recruitment" %&% suffix, xlab="", ylab="Millions",xref=projRefLine, doWorms=doWorms, wormSims=wormSims)
    #plot(p5)

    if (mwgPlots)
    {
      multiplot(list(p2b, p6, p3, p4),cols=2)
    }
    else
    {
      multiplot(list(p2, p2b, p6, p5, p3, p4),cols=2)
    }

    if (plotByRF)
    {
      #Catch by fishery
      if (mseObj@nfleets > 1)
      {
        plotsPerPage <- 12
        pCFList <- list(NULL)
        paneli <- 1

        for (fi in 1:mseObj@nfleets)
        {
          rawPlotDat <- mseObj@CMbyF[iMP,nSims,mseObj@targpop,,fi] #nMPs,nsim,years
          plotDat$data <- as.vector(rawPlotDat)
          names(plotDat) <- c("year","iter","data")
          plotDat$year <- rep(mseObj@yrLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
          plotDat$iter <- rep(1:dim(rawPlotDat)[1], times = dim(rawPlotDat)[2]) #nsim nyears
          plotDat$data <- as.vector(rawPlotDat)
          pName <- "pCF" %&% fi
          assign(pName, plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&%": Catch " %&% "F" %&% fi %&% suffix, xlab="Year", ylab="Catch",xref=projRefLine, yref=0, doWorms=doWorms, wormSims=wormSims))
          #plot(get(pName))
          pCFList[[paneli]] <- get(pName)

          if (paneli%%plotsPerPage == 0 | fi==mseObj@nfleets)
          {
            multiplot(pCFList,cols=3)
            paneli <- 0
            pCFList <- list(NULL)
          }

          paneli <- paneli+1
        }
      }

      #CPUE by region
      if (mseObj@nareas > 1)
      {
        pIRList <- list()

        for (ri in 1:mseObj@nareas)
        {
          rawPlotDat <- mseObj@IobsRArchive[iMP,nSims,,ri] #nMPs,nsim,years
          plotDat$data <- as.vector(rawPlotDat)
          names(plotDat) <- c("year","iter","data")
          plotDat$year <- rep(mseObj@yrLabels, each  = dim(rawPlotDat)[1]) #nyears, nsim
          plotDat$iter <- rep(1:dim(rawPlotDat)[1], times = dim(rawPlotDat)[2]) #nsim nyears
          plotDat$data <- as.vector(rawPlotDat)
          pName <- "pIR" %&% ri
          assign(pName, plotStock.f(plotDat, main="MP " %&% mseObj@MPs[iMP]  %&%": CPUE " %&% "R" %&% ri %&% suffix, xlab="Year", ylab="CPUE",xref=projRefLine, yref=0, doWorms=doWorms, wormSims=wormSims))
          #plot(get(pName))
          pIRList[[ri]] <- get(pName)
        }

        multiplot(pIRList,cols=1)
        #multiplot(unlist(pIRList),cols=1)
      }
    }
  }
}


# data.frame with
# year, iter, data
plotStock.f <- function(x, main="", xlab="", ylab="", na.rm=TRUE, xref=2014, yref=0, yreflim=0, yreftarget=0, doWorms=F,
    probs=c(0.10, 0.25, 0.50, 0.75, 0.90), wormSims=NA, type=7)
{
    # check probs length is odd
    if((length(probs) %% 2 ==0)) stop("quantile probs can only be a vector of odd length")
    quans <- paste0(probs * 100, "%")
    mid <- ceiling(length(quans)/2)
    mquan <- quans[mid]

    # compute quantiles
    df <- as.data.frame(do.call(rbind, as.list(by(x$data, as.factor(x$year), quantile, probs=probs, type=7,na.rm=T))))
    df$year <- as.numeric(rownames(df))


    # basic plot data vs. year
    p <- ggplot(data=df, aes_q(x=quote(year), y=as.name(mquan))) +
      # line + xlab + ylab +
      geom_line(colour="black") + xlab(xlab) + ylab(ylab) +
      # limits to include 0 +
      expand_limits(y=0) +
      theme(plot.title = element_text(size=10))+
      theme(axis.title=element_text(size=10)) +
       ggtitle(main) +
      # no legend +
      theme(legend.title = element_blank())

    #reference lines
    p <- p + geom_hline(yintercept=0) +
             geom_hline(yintercept=yreftarget, colour="green") +
             geom_hline(yintercept=yreflim, colour="red") +
             geom_hline(yintercept=yref) +
             geom_vline(xintercept=xref)

    #worms
  if (doWorms)
  {
      #duplicate percentiles due to crashing ignored
      worms <- subset(x[,], x$iter %in% wormSims)
      worms <- reshape(worms,timevar="iter", direction='wide',idvar='year')

     for (iprobs in 1:(ncol(worms)-1))
     {
         p <- p +  geom_line(aes_q(x=worms$year), y = worms[,iprobs+1],
            colour=iprobs+3, alpha = 0.5, linetype=1)
       }
    }

    #p <- p +
  #    # extreme probs as dotted line
  #    geom_line(aes_q(x=quote(year), y = as.name(quans[1])),
  #      colour="red", alpha = .50, linetype=3) +
  #    geom_line(aes_q(x=quote(year), y = as.name(quans[length(quans)])),
  #      colour="red", alpha = .50, linetype=3)

      # all others as ribbons of changing alpha
  if (length(quans) > 3)
  {
        #ids <- seq(2, mid-1)
        ids <- seq(1, mid-1)
        for(i in ids)
          p <- p + geom_ribbon(aes_q(x=quote(year),
            ymin = as.name(quans[i]),
            ymax = as.name(quans[length(quans)-i+1])),
            #fill="red", alpha = probs[i])
            fill=rgb(0.01,0.01,0.01), alpha = probs[i])
      }

#    plot(p)
    return(p)
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL)
{
  library(grid)
  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)

  plots <- c(..., plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout))
  {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1)
  {
    print(plots[[1]])
  }
  else
  {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots)
    {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
