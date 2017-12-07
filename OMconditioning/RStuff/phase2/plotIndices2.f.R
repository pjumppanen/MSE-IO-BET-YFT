plotIndices2.f <- function (modList = gridZList, MSYyLim=c(0,300), mFile = F, mfrowLayout = c(3,2), SPB_Yr="SPB_276",
   inputRefLines=T, #override list below which are extracted from reports etc
   refSSB0=c(3421,2492,4628),refSSBY=c(786,457,1009), refSSBMSY=c(1200,751,1483), refSSBYoSSBMSY=c(0.66,0.48,0.72), refSSBYoSSB0=c(0.23,0.17,0.32),refMSY=c(423))
{
    if(inputRefLines){
      refMSY           <- ref$derived_quants[ref$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
      refSSB0          <- ref$derived_quants[ref$derived_quants$LABEL == 'SPB_Virgin',]$Value/1000
      refSSBY          <- ref$derived_quants[ref$derived_quants$LABEL == SPB_Yr,]$Value/1000
      refSSBMSY        <- ref$derived_quants[ref$derived_quants$LABEL == 'SSB_MSY',]$Value/1000
      refSSBYoSSBMSY   <- refSSBY/refSSBMSY
      refSSBYoSSB0     <- refSSBY/refSSB0
    }

    par(mar = c(4, 4, 3, 1) + 0.1)
    factorList <- array(NA, dim = c(length(modList), length(unlist(strsplit(modList[1], 
        split = "_")))))
    for (mi in 1:length(modList)) {
        factors <- strsplit(modList[mi], split = "_")
        #one-off kludge:
        #if(length(unlist(factors))==8) factors[[9]] <- "ess5"

        factorList[mi, ] <- unlist(factors)
    }
    print(factorList)
    numMod <- length(modList)
    numFleets <- length(get(modList[1])$Length_comp_Eff_N_tuning_check$mean_effN)
    numCPUE <- length(levels(as.factor(get(modList[1])$cpue$FleetName)))
    pBoundall <- NULL
    MSYall <- NULL
    SSBYall <- NULL
    BMinall <- NULL
    BMinoSSBMSYall <- NULL
    SSBYoSSBMSYall <- NULL
    SSBMSYall <- NULL
    SSBYoSSB0all <- NULL
    SSBMSYoSSB0all <- NULL
    FYoFMSYall <- NULL
    CYoMSYall <- NULL
    essall <- array(NA, dim = c(numMod, numFleets))
    cpueRMSEall <- array(NA, dim = c(numMod, numCPUE))
    gradall <- NULL
    boundCheckall <- NULL
    LLHall <- NULL
    tagFit <- NULL
    recFit <- NULL     #rec LLH from SS3, including whatever treatment of constants and bias corrections
    recFit2 <- NULL    #
    rTrendall <- NULL
    for (i in 1:length(modList)) {
        print(modList[i])
        if (modList[i] != "nullMod") {
            m <- get(modList[i])
            #MSY <- max(m$equil_yield$Catch)

            MSY <- m$derived_quants[m$derived_quants$LABEL == 'TotYield_MSY',]$Value
            MSYall <- c(MSYall, MSY)
            #BMin <- min(m$equil_yield$Depletion)
            #BMinall <- c(BMinall, BMin)
            #SSBY <- as.numeric(m$derived_quants[substring(m$derived_quants$LABEL, 1, 8) == SPB_Yr, ]$Value)
            SSBY <- m$derived_quants[m$derived_quants$LABEL == SPB_Yr,]$Value
            SSBYall <- c(SSBYall,SSBY)
            #SSB0 <- as.numeric(subset(m$derived_quants, m$derived_quants$LABEL == "SPB_Virgin")$Value)
            SSB0 <- m$derived_quants[m$derived_quants$LABEL == 'SPB_Virgin',]$Value

            #SSBMSY <- m$equil_yield$Depletion[m$equil_yield$Catch == MSY][1] * SSB0
            SSBMSY <- m$derived_quants[m$derived_quants$LABEL == 'SSB_MSY',]$Value
            SSBMSYall <- c(SSBMSYall,SSBMSY)
          
            #BMinoSSBMSY <- BMin/m$equil_yield$Depletion[m$equil_yield$Catch == MSY][1]
            #BMinoSSBMSYall <- c(BMinoSSBMSYall, BMinoSSBMSY)
            SSBYoSSBMSYall <- c(SSBYoSSBMSYall, SSBY/SSBMSY)
            SSBMSYoSSB0all <- c(SSBMSYoSSB0all, SSBMSY/SSB0)
            SSBYoSSB0all <- c(SSBYoSSB0all, SSBY/SSB0)
            #CY = sum(m$catchc1.0[m$timeseries$Yr == 2009])
            #CYoMSYall <- c(CYoMSYall, CY/MSY)
            #print(c("BMin", BMin))
            print(c("MSY", MSY))
            print(c("SSBY", SSBY))
            print(c("SSB0", SSB0))
            #print(c("dep", m$equil_yield$Depletion[m$equil_yield$Catch ==  MSY]))
            print(c("SSBMSY", SSBMSY))
            #print(c("CY", CY))
            #pBound <- nrow(subset(m$parameters, (m$parameters[9] ==     "HI" | m$parameters[9] == "LO"))[1:9])
            #pBoundall <- c(pBound, pBoundall)
            ess <- m$Length_comp_Eff_N_tuning_check$mean_effN
            names(ess) <- m$Length_comp_Eff_N_tuning_check$Fleet
            colnames(essall) <- m$Length_comp_Eff_N_tuning_check$Fleet
            essall[i, ] <- as.numeric(ess)
            tmp <- aggregate(m$cpue$Dev^2, by = list(fleetName = m$cpue$FleetName), 
                FUN = mean)
            cpueRMSE <- sqrt(tmp[, 2])
            names(cpueRMSE) <- tmp[, 1]
            #adapt this and include autocorrelation
            #if ("SURVEY1" %in% names(cpueRMSE)) {
            #    tmp <- m$cpue
            #    tmp <- subset(tmp, tmp$FleetName %in% c("SURVEY1", 
            #      "SURVEY2"))
            #    tmp <- sqrt(mean(tmp$Dev^2))
            #    cpueRMSE[names(cpueRMSE) == "SURVEY1"] <- tmp
            #}
            #cpueRMSE <- cpueRMSE[names(cpueRMSE) != "UJ91p_SW"]
            if ("JS" %in% factorList[i, ]) {
            }
            cpueRMSEall[i, ] <- as.numeric(cpueRMSE)
            gradall <- c(gradall, m$maximum_gradient_component)
            boundCheckall <- c(boundCheckall, sum(m$ref[[131]]$Status=="HI" | ref[[131]]$Status=="LO"))

            LLHall <- c(LLHall, m$likelihoods_used$values[1])

            tmp <- m$likelihoods_raw_by_fleet[m$likelihoods_raw_by_fleet[,1]=="Tag_comp_Like" | m$likelihoods_raw_by_fleet[,1]=="Tag_negbin_Like",]
            tagFit <- c(tagFit, sum(as.numeric(tmp[1,3:length(tmp)])) + sum(as.numeric(tmp[2,3:length(tmp)])))
            recFit <- c(recFit, m$likelihoods_used['Recruitment','values'])

            rdev <- m$recruit$dev[m$recruit$era == "Main"]
            nrdev <- length(rdev)
            rSlope <- coefficients(lm(rdev ~ c(1:nrdev)))[2]
            #rTrend <- length(m$recruit$dev[m$recruit$era == "Main"]) * 
            #    rSlope * 100

            sigR <- ref$sigma_R_in
            #recFit2 <- c(recFit2, sum(log(sqrt(2*pi*sigR^2)) + ((rdev)^2 - 0.5*sigR^2)/(2*sigR^2)) ) #exactly the same as SS output except for a constant
            recFit2 <- c(recFit2,sqrt(sum(rdev^2)/nrdev))

            rTrendall <- c(rTrendall, rSlope)
        }
        else {
            MSYall <- c(MSYall, NA)
            SSBMSYall <- c(SSBMSYall, NA)
            SSBYoSSBMSYall <- c(SSBYoSSBMSYall, NA)
            SSBYoSSB0all <- c(SSBYoSSB0all, NA)
            SSBMSYoSSB0all <- c(SSBMSYoSSB0all, NA)
            FYoFMSYall <- c(FYoFMSYall, NA)
            CYoMSYall <- c(CYoMSYall, NA)
            SSBYall <- c(SSBYall, NA)
            essall[i, ] <- NA
            cpueRMSEall[i, ] <- NA
            gradall <- c(gradall, NA)
            boundCheckall <- c(boundCheckall, NA)
            LLHall <- c(LLHall, NA)
            tagFit <- c(tagFit, NA)
            recFit <- c(recFit, NA)
            recFit2 <- c(recFit2, NA)
            rTrendall <- c(rTrendall, NA)
        }
    }   # modList i

    print(LLHall)
    LLHall <- LLHall - min(LLHall, na.rm = T)
    fList <- NULL
    for (i in 1:ncol(factorList)) {
        if (length(unique(factorList[, i])) > 1) 
            fList <- cbind(fList, factorList[, i])
    }
    print(LLHall)
    print(dim(fList))
    par(mfrow = mfrowLayout)
    boxplot(LLHall, ylab = "Likelihood units")
    if(!is.null(fList)) for (i in 1:ncol(fList)) {
        boxplot(LLHall ~ fList[, i], ylab = "Likelihood units")
    }
    par(mfrow = c(1, 1))
    boxplot(LLHall ~ fList[, 1]:fList[, 2], ylab = "Likelihood units")
    boxplot(LLHall ~ fList[, 1]:fList[, 3], ylab = "Likelihood units")

    #par(mfrow = c(1,1))
    #boxplot(tagFit ~ fList[,], ylab = "Tag likelihood (no lambda)", col=8)
    #first <- length(unique(fList[,1]))+1
    #for (i in 2:ncol(fList)) {
    #  last <- first + length(unique(fList[,i]))-1
    #  boxplot(tagFit ~ fList[,i], col=i, add=T, at=first:last)
    #  first <- last+1
    #}



    par(mfrow = c(1,1))
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    boxplot(cbind(tagFit,spacer), ylab = "Tag likelihood (no lambda)", names=c("all",rep(" ",length(unique(as.vector(fList))))))

    first <- 2
    colList=c(8,2:100)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(tagFit ~ fList[,i], col=colList[i], add=T, at=first:last)
      first <- last+1
    }

    #tagLLHs for full mixing, t10 only - one off code not robust
    tagPartition <- F
    if(tagPartition){

      # partition by x3 and x8
      par(mfrow = c(1,1))
      dropCol <- which.max(colSums(fList[] == "x3"))
      tagFit2 <- tagFit[fList[,dropCol]=="x3"]
      fList2  <- fList[fList[,dropCol]=="x3",]
      fList2  <- fList2[,fList2[1,] != "x3"]
      boxplot(tagFit2 ~ fList2[,], ylab = "Tag likelihood (no lambda)", main="x3 - 3 qtr tag mixing", names=rep(" ",length(unique(as.vector(fList2)))))
      first <- 1
      colList=c(8,2,3,4,5,6,7,8,9)
      for (i in 1:ncol(fList2)) {
        last <- first + length(unique(fList2[,i]))-1
        boxplot(tagFit2 ~ fList2[,i], col=colList[i], add=T, at=first:last)
        first <- last+1
      }
      dropCol <- which.max(colSums(fList[] == "x8"))
      tagFit2 <- tagFit[fList[,dropCol]=="x8"]
      fList2  <- fList[fList[,dropCol]=="x8",]
      fList2  <- fList2[,fList2[1,] != "x8"]
      boxplot(tagFit2 ~ fList2[,], ylab = "Tag likelihood (no lambda)", main="x8 - 8 qtr tag mixing", names=rep(" ",length(unique(as.vector(fList2)))))
      first <- 1
      for (i in 1:ncol(fList2)) {
        last <- first + length(unique(fList2[,i]))-1
        boxplot(tagFit2 ~ fList2[,i], col=colList[i], add=T, at=first:last)
        first <- last+1
      }

      par(mfrow = c(1,1))
      dropCol <- which.max(colSums(fList[] == c("x3")))
      tagFit2 <- tagFit[fList[,dropCol]=="x3"]
      fList2  <- fList[fList[,dropCol]=="x3",]
      fList2  <- fList2[,fList2[1,] != "x3"]

      dropCol <- which.max(colSums(fList2[] == c("t10")))
      tagFit3 <- tagFit2[fList2[,dropCol]=="t10"]
      fList3  <- fList2[fList2[,dropCol]=="t10",]
      fList3  <- fList3[,fList3[1,] != "t10"]


      boxplot(tagFit3 ~ fList3[,], ylab = "Tag likelihood (no lambda)", main="(x3, t10) 3 qtr tag mixing and lambda=1.0", names=rep(" ",length(unique(as.vector(fList3)))))
      first <- 1
      colList=c(8,2,3,4,5,6,7,8,9)
      for (i in 1:ncol(fList3)) {
        last <- first + length(unique(fList3[,i]))-1
        boxplot(tagFit3 ~ fList3[,i], col=colList[i], add=T, at=first:last)
        first <- last+1
      }
    } #tagPartiion










    par(mfrow = mfrowLayout)
    boxplot(gradall, ylab = "max. gradient")
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(gradall ~ fList[, i], ylab = "max. gradient xxx")
    }
    par(mfrow = mfrowLayout)
    boxplot(log10(gradall), ylab = "max. gradient")
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(log10(gradall) ~ fList[, i], ylab = "log10(max Grad)")
    }

    #par(mfrow = c(1,1))
    #boxplot(gradall ~ fList[,], ylab = "log10(max. gradient)", col=8)
    #3first <- length(unique(fList[,1]))+1
    #for (i in 2:ncol(fList)) {
    #  last <- first + length(unique(fList[,i]))-1
    #  boxplot(log10(gradall) ~ fList[,i], col=i, add=T, at=first:last)
    #      first <- last+1
    #}

    par(mfrow = c(1,1))
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    boxplot(cbind(log10(gradall),spacer), ylab = "log10(max. gradient)", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(log10(gradall) ~ fList[,], ylab = "log10(max. gradient)", names=rep(" ",length(unique(as.vector(fList)))))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(log10(gradall) ~ fList[,i], col=colList[i], add=T, at=first:last, outcol=colList[i])
      first <- last+1
    }

    #check status on bounds - relies on r4ss definition of close to bounds Status = HI or LO(soft limit presumably)
    par(mfrow = c(1,1))
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    boxplot(cbind(boundCheckall,spacer), ylab = "Par bounds issues", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(boundCheckall ~ fList[,], ylab = "Par bounds issues", names=rep(" ",length(unique(as.vector(fList)))))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(boundCheckall ~ fList[,i], col=colList[i], add=T, at=first:last, outcol=colList[i])
      first <- last+1
    }



    #par(mfrow = mfrowLayout)
    #boxplot(pBoundall, ylab = "N par near bounds")
    #for (i in 1:ncol(fList)) {
    #    boxplot(pBoundall ~ fList[, i], ylab = "N par near bounds")
    #}
    #par(mfrow = mfrowLayout)
    #boxplot(BMinall, ylab = "Minimum depletion at high F (B/B0)", 
    #    ylim = c(0, 0.25))
    #for (i in 1:ncol(fList)) {
    #    boxplot(BMinall ~ fList[, i], ylab = "B/B0 (at high F)", 
    #        ylim = c(0, 0.25))
    #}
    #par(mfrow = mfrowLayout)
    #boxplot(BMinoSSBMSYall, ylab = "BMin/BMSY)", ylim = c(0, 
    #    1))
    #for (i in 1:ncol(fList)) {
    #    boxplot((BMinoSSBMSYall) ~ fList[, i], ylab = "BMin/BMSY ", 
    #        ylim = c(0, 1))
    #}
    for (f in 1:numFleets) {
        par(mfrow = mfrowLayout)
        boxplot(essall[, f], ylab = "mean ESS", main = "ESS " %&%
            colnames(essall)[f])
        if(!is.null(fList))for (i in 1:ncol(fList)) {
            boxplot(essall[, f] ~ fList[, i], ylab = "mean ESS",
                main = "ESS " %&% colnames(essall)[f])
        }
    }
    par(mfrow = c(1,1))
    boxplot(essall[,],ylim=c(0,150),xlab="Fishery", ylab="ESS")
    abline(h=5)
    boxplot(essall[,],ylim=c(0,600),xlab="Fishery", ylab="ESS")
    abline(h=5)
    boxplot(log(essall[,]),ylim=c(0,7),xlab="Fishery", ylab="log(ESS)")
    abline(h=log(5))

    par(mfrow = mfrowLayout)
    boxplot(4*MSYall/1000, ylab = "MSY (1000 t) (X4)", ylim=MSYyLim)
    abline(h=refMSY,col=3)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(4*MSYall/1000 ~ fList[, i], ylab = "MSY (1000 t) (X4)", ylim=MSYyLim)
        abline(h=refMSY,col=3)
    }

    #par(mfrow = c(1,1))
    #boxplot(4*MSYall/1000 ~ fList[,], ylab = "MSY", col=8, ylim=c(0,1000))
    #first <- length(unique(fList[,1]))+1
    #for (i in 2:ncol(fList)) {
    #  last <- first + length(unique(fList[,i]))-1
    #  boxplot(4*MSYall/1000 ~ fList[,i], col=i, add=T, at=first:last)
    #  first <- last+1
    #}
    #abline(h=refMSY,col=1)

    par(mfrow = c(1,1))
    hist(4*MSYall/1000, breaks=20, main = "MSY (1000 t)",xlab="")
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    plotDat <- cbind(4*MSYall/1000,spacer)
    boxplot(plotDat, ylab = "MSY (1000 t)", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(4*MSYall/1000 ~ fList[,], ylab = "MSY (1000 t)", names=rep(" ",length(unique(as.vector(fList)))), ylim=c(0,1000))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(4*MSYall/1000 ~ fList[,i], col=colList[i], add=T, at=first:last)
      first <- last+1
    }
    abline(h=refMSY,col=1)


    par(mfrow = mfrowLayout)
    boxplot(SSBYall/1000, ylab = "SSBY (1000 t)", ylim = c(0, 2500))
    abline(h=refSSBY,col=3)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(SSBYall/1000 ~ fList[, i], ylab = "SSBY (1000 t)", ylim = c(0, 2500))
      abline(h=refSSBY,col=3)
    }

    par(mfrow = mfrowLayout)
    boxplot(SSBMSYall/1000, ylab = "SSBMSY (1000 t)", ylim = c(0, 3000))
    abline(h=refSSBMSY,col=3)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(SSBMSYall/1000 ~ fList[, i], ylab = "SSBMSY (1000 t)", ylim = c(0, 3000))
        abline(h=refSSBMSY,col=3)
    }
    
    par(mfrow = mfrowLayout)
    boxplot(SSBYoSSBMSYall, ylab = "SSBY / SSBMSY", ylim=c(0,3.5))
    abline(h=refSSBYoSSBMSY,col=3)
    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(SSBYoSSBMSYall ~ fList[, i], ylab = "SSBY / SSBMSY", ylim=c(0,3.5))
        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    abline(h=refSSBYoSSBMSY,col=3)
    }


    #par(mfrow = c(1,1))
    #boxplot(SSBYoSSBMSYall ~ fList[,], ylab = "SSB(2014)/SSBMSY", col=8, ylim=c(0,1.5))
    #first <- length(unique(fList[,1]))+1
    #for (i in 2:ncol(fList)) {
    #  last <- first + length(unique(fList[,i]))-1
    #  boxplot(SSBYoSSBMSYall ~ fList[,i], col=i, add=T, at=first:last)
    #  first <- last+1
    #}
    #abline(h=refSSBYoSSBMSY,col=1)

    par(mfrow = c(1,1))
    hist(SSBYoSSBMSYall, breaks=20, main = "SSB(Current) / SSBMSY",xlab="")
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    plotDat <- cbind(SSBYoSSBMSYall,spacer)
    boxplot(plotDat, ylab = "SSB(2015)/SSBMSY", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(SSBYoSSBMSYall ~ fList[,], ylab = "SSB(2014)/SSBMSY", names=rep(" ",length(unique(as.vector(fList)))), ylim=c(0,1.5))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(SSBYoSSBMSYall ~ fList[,i], col=colList[i], add=T, at=first:last)
      first <- last+1
    }
    abline(h=refSSBYoSSBMSY,col=1)




    par(mfrow = mfrowLayout)
    boxplot(SSBMSYoSSB0all, ylab = "SSBMSY / SSB0", ylim = c(0,
        1))
    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(SSBMSYoSSB0all ~ fList[, i], ylab = "SSBMSY / SSB0", 
            ylim = c(0, 1))
        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    }

    par(mfrow = mfrowLayout)
    boxplot(SSBYoSSB0all, ylab = "SSBY / SSB0", ylim = c(0,
        1))
    abline(h=refSSBYoSSB0,col=3)
    lines(c(0, 10000), c(0.2, 0.2), lty = 2, col = 2)
    lines(c(0, 10000), c(0.4, 0.4), lty = 2, col = 2)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(SSBYoSSB0all ~ fList[, i], ylab = "SSBY / SSB0", 
            ylim = c(0, 1))
        lines(c(0, 10000), c(0.2, 0.2), lty = 2, col = 2)
        lines(c(0, 10000), c(0.4, 0.4), lty = 2, col = 2)
        abline(h=refSSBYoSSB0,col=3)
    }

    #par(mfrow = c(1,1))
    #boxplot(SSBYoSSB0all ~ fList[,], ylab = "SSB(2014)/SSB0", col=8, ylim=c(0,1))
    #first <- length(unique(fList[,1]))+1
    #for (i in 2:ncol(fList)) {
    #  last <- first + length(unique(fList[,i]))-1
    #  boxplot(SSBYoSSB0all ~ fList[,i], col=i, add=T, at=first:last)
    #  first <- last+1
    #}
    #abline(h=refSSBYoSSB0,col=1)


    par(mfrow = c(1,1))
    hist(SSBYoSSB0all, breaks=20, main = "SSB(current) / SSB0",xlab="")
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    plotDat <- cbind(SSBYoSSB0all,spacer)
    boxplot(plotDat, ylab = "SSB(current)/SSB0", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(SSBYoSSB0all ~ fList[,], ylab = "SSB(current)/SSB0", names=rep(" ",length(unique(as.vector(fList)))), ylim=c(0,1))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(SSBYoSSB0all ~ fList[,i], col=colList[i], add=T, at=first:last)
      first <- last+1
    }
    abline(h=refSSBYoSSB0,col=1)



    #par(mfrow = mfrowLayout)
    #boxplot(CYoMSYall, ylab = "C(2009) / MSY")
    #lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    #for (i in 1:ncol(fList)) {
    #    boxplot(CYoMSYall ~ fList[, i], ylab = "C(2009) / MSY")
    #    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    #}
    par(mfrow = mfrowLayout)
    boxplot(cpueRMSEall, ylab = "CPUE RMSE All", ylim = c(0, 0.8))
    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(cpueRMSEall ~ fList[, i], ylab = "CPUE RMSE All",
            ylim = c(0, 0.8))
        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    }

    par(mfrow = c(2,1))
    for(r in 1:4){  #region
      boxplot(cpueRMSEall[,r] ~ fList[,], ylab = "CPUE RMSE", main="Region " %&% r, col=8)
      first <- length(unique(fList[,1]))+1
      for (i in 2:ncol(fList)) {
        last <- first + length(unique(fList[,i]))-1
        boxplot(cpueRMSEall[,r] ~ fList[,i], col=i, add=T, at=first:last)
        first <- last+1
      }
    }

    par(mfrow = c(1,1))
    spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
    plotDat <- cbind(rowMeans(cpueRMSEall),spacer)
    boxplot(plotDat, ylab = "CPUE RMSE", main="mean(all regions) ", names=c("all",rep(" ",length(unique(as.vector(fList))))))
    first <- 2
    #boxplot(cpueRMSEall[,r] ~ fList[,], ylab = "CPUE RMSE", main="Region " %&% r, names=rep(" ",length(unique(as.vector(fList)))))
    #first <- 1
    colList=c(8,2,3,4,5,6,7,8,9)
    for (i in 1:ncol(fList)) {
      last <- first + length(unique(fList[,i]))-1
      boxplot(rowMeans(cpueRMSEall) ~ fList[,i], col=colList[i], add=T, at=first:last)
      text((first:last),1.01*max(rowMeans(cpueRMSEall)), round(table(fList[,i])/sum(table(fList[,i])),2))
      first <- last+1
    }

    par(mfrow = c(1,1))
    for(r in 1:4){  #region
      spacer <- array(NA, dim=c(dim(fList)[1], length(unique(as.vector(fList)))))
      plotDat <- cbind(cpueRMSEall[,r],spacer)
      boxplot(plotDat, ylab = "CPUE RMSE", main="Region " %&% r, names=c("all",rep(" ",length(unique(as.vector(fList))))))
      first <- 2
      #boxplot(cpueRMSEall[,r] ~ fList[,], ylab = "CPUE RMSE", main="Region " %&% r, names=rep(" ",length(unique(as.vector(fList)))))
      #first <- 1
      colList=c(8,2,3,4,5,6,7,8,9)
      for (i in 1:ncol(fList)) {
        last <- first + length(unique(fList[,i]))-1
        boxplot(cpueRMSEall[,r] ~ fList[,i], col=colList[i], add=T, at=first:last)
        first <- last+1
      }
    }


    #par(mfrow = mfrowLayout)
    #for (c in 1:numCPUE) {
    #    boxplot(cpueRMSEall[, c], ylab = "CPUE RMSE", ylim = c(0, 
    #        0.8))
    #    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    #    for (i in 1:ncol(fList)) {
    #        boxplot(cpueRMSEall[, c] ~ fList[, i], ylab = "CPUE RMSE", 
    #            ylim = c(0, 0.8), main = "CPUE RMSE " %&% names(cpueRMSE)[c])
    #        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    #    }
    #}
    #par(mfrow = c(4, 2))
    #for (f in 1:numFleets) {
    #    boxplot(essall[, f] ~ fList[, 5], ylab = "mean ESS", 
    #        main = "ESS " %&% colnames(essall)[f])
    #}
    #par(mfrow = c(4, 2)) #what's this?
    #for (f in 1:numFleets) {
    #    boxplot(essall[, f] ~ fList[, 2], ylab = "mean ESS", 
    #        main = "ESS " %&% colnames(essall)[f])
    #}
    
    par(mfrow = mfrowLayout)
    boxplot(rTrendall, ylab = " Slope in rec devs (per qtr)")
    if(!is.null(fList))for (i in 1:ncol(fList)) {
        boxplot(rTrendall ~ fList[, i], ylab = " Slope in rec devs (per qtr)")
        lines(c(0, 10000), c(0, 0), lty = 2, col = 2)
    }
    par(mar = c(5, 4, 4, 2) + 0.1)


    CPUE.fit  <- apply(cpueRMSEall[,], FUN=mean, MARGIN=1) # mean CPUE fit
    CL.Fit    <- apply(essall[,], FUN=mean, MARGIN=1)     # mean CL ESS
    tag.Fit   <- tagFit
    rec.Fit   <- recFit
    rec.RMSE  <- recFit2
    #tag.Fit3  <- tagFit[fList[,5]=="x3"]
    #tag.Fit8  <- tagFit[fList[,5]=="x8"]

    max.Grad  <- log10(gradall)
    rec.Trend <- rTrendall
    MSY       <- 4.*MSYall/1000.
    B_B.MSY   <- SSBYoSSBMSYall
    B_B0      <- SSBYoSSB0all
    F_F.MSY   <- FYoFMSYall
    corrDat<- cbind(jitter(CPUE.fit, amount=diff(range(CPUE.fit))/50),
                    jitter(CL.Fit, amount=diff(range(CL.Fit))/50),
                    jitter(tag.Fit, amount=diff(range(tag.Fit))/50),
                    jitter(rec.Trend, amount=diff(range(rec.Trend))/50),
                    jitter(rec.Fit, amount=diff(range(rec.Fit))/50),
                    jitter(MSY, amount=diff(range(MSY)/50)),
                    jitter(B_B.MSY, amount=diff(range(B_B.MSY))/50))
    colnames(corrDat) <- c('CPUE.fit', 'CL.Fit', 'tag.Fit', 'rec.Trend', 'rec.Fit', 'MSY', 'B_B.MSY')
    rownames(corrDat) <- NULL

    for (i in 1:ncol(fList)){
      #dev.new() - need to do something to graphics pars here of legend layout screws up
      reset_graphics_par()

      f <- fList[,i]
      ff <- as.factor(f)
      chart.Correlation(corrDat,bg=seq(1:length(levels(ff)))[ff],pch=21)
      par(xpd=T)
      legend( .85,.95,sort(as.vector(unique(ff))), fill=seq(1:length(levels(ff))),bg='white')

    }

    corrDat<- cbind(rec.Trend, rec.Fit, tag.Fit, MSY)
    colnames(corrDat) <- c('rec.Trend', 'rec.Fit', 'tag.Fit','MSY')
    rownames(corrDat) <- NULL
    for (i in 1:ncol(fList)){
      #dev.new() - need to do something to graphics pars here of legend layout screws up
      reset_graphics_par()

      f <- fList[,i]
      ff <- as.factor(f)
      chart.Correlation(corrDat,bg=seq(1:length(levels(ff)))[ff],pch=21)
      par(xpd=T)
      legend( .85,.95,sort(as.vector(unique(ff))), fill=seq(1:length(levels(ff))),bg='white')

    }
      reset_graphics_par()


#    reset_graphics_par()
    #look at ESS by fishery vs MSY
    #pairs(cbind(MSY, essall[,1:5]))
    #pairs(cbind(MSY, essall[,6:10]))
    #pairs(cbind(MSY, essall[,11:15]))
    #pairs(cbind(MSY, essall[,16:20]))
    #pairs(cbind(MSY, essall[,21:23]))


    #look at CPUE by area vs MSY
    #corrDat<- cbind(MSY, cpueRMSEall)
    #pairs(corrDat)

    plot(MSY, rec.Trend, ylab='slope (dev/qtr)', main="Recruitment Deviation Trend" )

    plot(x=jitter(MSY, amount=mean(MSY)/100),y=jitter(rec.Trend, amount=max(rec.Trend)/25), ylab='slope (dev/qtr)', main="Recruitment Deviation Trend", cex=0.1)
    par(mfrow=c(2,2))

    par(mfrow=c(2,2))
    for(r in 1:4){
      plot(MSY, cpueRMSEall[,r], ylab='RMSE', main="CPUE Region " %&% r)
    }

    par(mfrow=c(3,3))
    for(f in 1:numFleets){
      plot(MSY, essall[,f], ylab='post-fit ESS', main="Fishery " %&% f)
    }

    BY   <- SSBYall
    BMSY <- SSBMSYall
    return(cbind(CPUE.fit, CL.Fit, tag.Fit, max.Grad, rec.Trend, rec.Fit, rec.RMSE, MSY, B_B.MSY, B_B0, BY, BMSY))
}

# this fixes the problem arising from chart.Correlation
reset_graphics_par <- function(){
    op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                         ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                         cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                         col.axis = "black", col.lab = "black", col.main = "black",
                         col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                         fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                         ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                         font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                         lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                         lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                           4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                          1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                         oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                           0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                         #xpd = F,
                         plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                                 0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                         tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                             1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                         yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                  "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                  "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                  "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                  "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                  "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                  "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                  "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                  "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                  "yaxt", "ylbias"))
    par(op)
}


