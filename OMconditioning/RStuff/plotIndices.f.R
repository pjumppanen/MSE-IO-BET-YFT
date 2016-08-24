plotIndices.f <- function (modList = gridZList, MSYyLim=c(0,300), mFile = F, mfrowLayout = c(3,2), SPB_Yr="SPB_272",
  refSSB0=c(3421,2492,4628),refSSBY=c(786,457,1009), refSSBMSY=c(1200,751,1483), refSSBYoSSBMSY=c(0.66,0.48,0.72), refSSBYoSSB0=c(0.23,0.17,0.32),refMSY=c(402,309,530)) 
{
    par(mar = c(4, 4, 3, 1) + 0.1)
    factorList <- array(NA, dim = c(length(modList), length(unlist(strsplit(modList[1], 
        split = "_")))))
    for (mi in 1:length(modList)) {
        factors <- strsplit(modList[mi], split = "_")
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
    LLHall <- NULL
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
            LLHall <- c(LLHall, m$likelihoods_used$values[1])
            rSlope <- coefficients(lm(exp(m$recruit$dev[m$recruit$era == 
                "Main"]) ~ c(1:length(m$recruit$dev[m$recruit$era == 
                "Main"]))))[2]
            #rTrend <- length(m$recruit$dev[m$recruit$era == "Main"]) * 
            #    rSlope * 100
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
            LLHall <- c(LLHall, NA)
            rTrendall <- c(rTrendall, NA)
        }
    }

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
    for (i in 1:ncol(fList)) {
        boxplot(LLHall ~ fList[, i], ylab = "Likelihood units")
    }
    par(mfrow = c(1, 1))
    boxplot(LLHall, ylab = "Likelihood units")
    for (i in 1:ncol(fList)) {
        boxplot(LLHall ~ fList[, 1]:fList[, 2], ylab = "Likelihood units")
    }
    par(mfrow = mfrowLayout)
    boxplot(gradall, ylab = "max. gradient")
    for (i in 1:ncol(fList)) {
        boxplot(gradall ~ fList[, i], ylab = "max. gradient xxx")
    }
    par(mfrow = mfrowLayout)
    boxplot(log10(gradall), ylab = "max. gradient")
    for (i in 1:ncol(fList)) {
        boxplot(log10(gradall) ~ fList[, i], ylab = "log10(max Grad)")
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
        for (i in 1:ncol(fList)) {
            boxplot(essall[, f] ~ fList[, i], ylab = "mean ESS", 
                main = "ESS " %&% colnames(essall)[f])
        }
    }

    par(mfrow = mfrowLayout)
    boxplot(4*MSYall/1000, ylab = "MSY (1000 t) (X4)", ylim=MSYyLim)
    abline(h=refMSY,col=3)
    for (i in 1:ncol(fList)) {
        boxplot(4*MSYall/1000 ~ fList[, i], ylab = "MSY (1000 t) (X4)", ylim=MSYyLim)
        abline(h=refMSY,col=3)
    }
    
    par(mfrow = mfrowLayout)
    boxplot(SSBYall/1000, ylab = "SSBY (1000 t)", ylim = c(0, 2500))
    abline(h=refSSBY,col=3)
    for (i in 1:ncol(fList)) {
        boxplot(SSBYall/1000 ~ fList[, i], ylab = "SSBY (1000 t)", ylim = c(0, 2500))
      abline(h=refSSBY,col=3)
    }

    par(mfrow = mfrowLayout)
    boxplot(SSBMSYall/1000, ylab = "SSBMSY (1000 t)", ylim = c(0, 3000))
    abline(h=refSSBMSY,col=3)
    for (i in 1:ncol(fList)) {
        boxplot(SSBMSYall/1000 ~ fList[, i], ylab = "SSBMSY (1000 t)", ylim = c(0, 3000))
        abline(h=refSSBMSY,col=3)
    }
    
    par(mfrow = mfrowLayout)
    boxplot(SSBYoSSBMSYall, ylab = "SSBY / SSBMSY", ylim=c(0,3.5))
    abline(h=refSSBYoSSBMSY,col=3)
    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    for (i in 1:ncol(fList)) {
        boxplot(SSBYoSSBMSYall ~ fList[, i], ylab = "SSBY / SSBMSY", ylim=c(0,3.5))
        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    abline(h=refSSBYoSSBMSY,col=3)
    }


    par(mfrow = mfrowLayout)
    boxplot(SSBMSYoSSB0all, ylab = "SSBMSY / SSB0", ylim = c(0, 
        1))
    lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
    for (i in 1:ncol(fList)) {
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
    for (i in 1:ncol(fList)) {
        boxplot(SSBYoSSB0all ~ fList[, i], ylab = "SSBY / SSB0", 
            ylim = c(0, 1))
        lines(c(0, 10000), c(0.2, 0.2), lty = 2, col = 2)
        lines(c(0, 10000), c(0.4, 0.4), lty = 2, col = 2)
        abline(h=refSSBYoSSB0,col=3)
    }
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
    for (i in 1:ncol(fList)) {
        boxplot(cpueRMSEall ~ fList[, i], ylab = "CPUE RMSE All",
            ylim = c(0, 0.8))
        lines(c(0, 10000), c(1, 1), lty = 2, col = 2)
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
    for (i in 1:ncol(fList)) {
        boxplot(rTrendall ~ fList[, i], ylab = " Slope in rec devs (per qtr)")
        lines(c(0, 10000), c(0, 0), lty = 2, col = 2)
    }
    par(mar = c(5, 4, 4, 2) + 0.1)
    return(1)
}