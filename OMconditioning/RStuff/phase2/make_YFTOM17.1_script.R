#=============================================================================================================================
# Basic R script for making and inspecting the conditioned YFT Operating Model OMyft17.1 and associated robustness scenarios
#=============================================================================================================================
#online
rootDir <- "H:\\C-offline\\MSE-IO-BET-YFT\\"  #modify for local path
#offline
#rootDir <- "C:\\MSE-IO-BET-YFT\\"  #modify for local path

library(TinnRcom)                            #if using TinnR IDE
library (r4ss)                                #R package supporting Stock Synthesis - if something's broken it might be becuase this changed
library(PerformanceAnalytics)                 # for chart.Correlation
source(paste(rootDir,"OMconditioning\\RStuff\\phase2\\pasteOperator.R",sep=""))

#source(paste(rootDir,"Source\\pasteOperator.R",sep=""))
source(rootDir %&% "OMconditioning\\RStuff\\seasAsYrToDecYr.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.1.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.1tagWt.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.1selTrend.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.1I10LLdev.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\importGrid.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\plotIndices.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\plotIndices2.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\timeSeriesPlots.f.R")


#extract some reference case assessment values
#For comparison: Langley 2016 YFT assessment ref case (provided by Dan Fu with nod from Adam)
ref2017 <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\AssessmentFiles2017\\refYFT2016", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=T, forecast=F)
ref <- ref2017
refMSY   <- 422
refMSYcv <- 0.035  #STD TotYield_MSY 3.67914e+03 (times 4 = 14716.56 )... CV = 3.5% ?!
refB_B.MSY   <- 0.89 # from Langley report
refB_B0 <- 0.29
ref2017SSBCurrent    <- ref2017[['derived_quants']][ref2017[['derived_quants']]$"LABEL"=="SPB_276","Value"]  #
ref2017SSBCurSD      <- ref2017[['derived_quants']][ref2017[['derived_quants']]$"LABEL"=="SPB_276","StdDev"]  #
ref2017SSBMSY        <- ref2017[['derived_quants']][ref2017[['derived_quants']]$"LABEL"=="SSB_MSY","Value"]  #
ref2017SSBMSYSD      <- ref2017[['derived_quants']][ref2017[['derived_quants']]$"LABEL"=="SSB_MSY","StdDev"]  #
# SPB_MSY is not included in covar
ref2017SSBCurMSYcorr <- ref2017[["CoVar"]][ref2017[['CoVar']]$"label.i"=="SPB_276" & ref2017[['CoVar']]$"label.j"=="SPB_MSY","corr"]
# Use correlation from MSY filtered grid (SPB_276, SPB_MSY) = 0.59
# then using cavar calcns (at bottom) of BY/BMSY
refBYoBMSYcv <- 0.052  # assuming 0 correlation yields cv = 0.075

SS_plots(ref2017, uncertainty=F)
ref$derived_quants[ref$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
max(ref$equil_yield$Catch)




# calculate the lowest 3 year (12qtr) rec period as a recruitment shock robustness scenario
rd <- ref2017[[35]]$dev
rd <- rd[!is.na(rd)]
cx <- cumsum(rd)
n <- 12
rsum <- (cx[(n+1):length(rd)] - cx[1:(length(rd) - n)]) / n
min(rsum) # = -0.5497042


# reference case Y17.1
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1'
gridY17.1List <- makeGridY17.1.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1List)
gridY17.1List <- makeGridY17.1.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file

# exit this R script and run the SS batch file(s)

# import summary results & diagnostics from the gridY17.1List of models
gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1//'
importGrid.f(gridList=gridY17.1List, gridDir=gridDir, covar=F)

#use boxplots to plot some key model results and diagnostics, marginalized by grid Options for the whole ensemble
plotIndices2.f(modList=gridY17.1List, mfrowLayout = c(4,2), MSYyLim=c(0,1000))

timeSeriesPlots.f(mList = gridY17.1List, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h70", "h80", "h90",
             "M10", "M08","M06",
             "t10","t01","t00",
             "q0","q1",
             "x3","x8",
             "iC","iH",
             "SS"),
             CIupper=0.25,CIlower=0.75,
     optWt = rep(1,100))





#Various non-grid one-off explorations below here...

#check iterative reweighting
tmp <- refESSreWeight <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refESSreWeight", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tmp <- refESSreWeightt00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refESSreWeightt00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tmp <- reft00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\R4MvEst_h80_M10_t00_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000


SS_plots(refESSreWeight, uncertainty=F)
SS_plots(refESSreWeightt00, uncertainty=F)

testGridList <- c("R4MvEst_h80_M10_t00_q0_x3_iH_SS","R4MvEst_h80_M10_t10_q0_x3_iH_SS","R4MvEst_h80_M10_t00_q0_x3_iH_SSir","R4MvEst_h80_M10_t10_q0_x3_iH_SSir")
gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1//'
importGrid.f(gridList=testGridList, gridDir=gridDir, covar=F)

#fudge around min 3 hard-coded dimensons
R4MvEst1_h80_M10_t00_q0_x3_iH_SS <- R4MvEst_h80_M10_t00_q0_x3_iH_SS
testGridList <- c("R4MvEst1_h80_M10_t00_q0_x3_iH_SS","R4MvEst_h80_M10_t10_q0_x3_iH_SS","R4MvEst_h80_M10_t00_q0_x3_iH_SSir","R4MvEst_h80_M10_t10_q0_x3_iH_SSir")
plotIndices2.f(modList=testGridList, mfrowLayout = c(4,2), MSYyLim=c(0,1000))

#test
#importGrid.f(gridList=gridY17.1List[1:5], gridDir=gridDir, covar=F)
#plotIndices2.f(modList=gridY17.1List[1:5], mfrowLayout = c(4,2), MSYyLim=c(0,1000))

#check PSLS and FS time series blocks - 3 blocks, 21 fisheries - equivalent (identical?) to 25 fishery set up (x8, h9 for some reason...)
tmp <- refPSselTimeSeries1 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries1", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries1, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000


# as refPSselTimeSeries1 except more (6?) TS blocks
tmp <- refPSselTimeSeries2 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries2", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries2, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000

# as refPSselTimeSeries2 except tag weight t00
tmp <- refPSselTimeSeries2t00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries2t00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries2t00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1058 K

# as refPSselTimeSeries2 except tag weight t00 and ESS reWeighting
tmp <- refPSselTimeSeries2t00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries2t00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries2t00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1066 K

# as refPSselTimeSeries1 except PS sel devs
tmp <- refPSselTimeSeries3devs <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries3devs", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries3devs, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 578 K


# as refPSselTimeSeries1 except PS sel devs + t00
tmp <- refPSselTimeSeries3devst00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries3devst00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries3devst00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 981K


# as refPSselTimeSeries1 except PS sel devs + t00 + ESS reweighting
tmp <- refPSselTimeSeries3devst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries3devst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries3devst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 778 K

# as refPSselTimeSeries1 except PS sel devs more relaxed
tmp <- refPSselTimeSeries4devs <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries4devs", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries4devs, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 622 K

# as refPSselTimeSeries1 except PS sel devs more relaxed + t00
tmp <- refPSselTimeSeries4devst00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries4devst00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries4devst00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 922K

# as refPSselTimeSeries1 except PS sel devs more relaxed + t00 + ESS reweighting
tmp <- refPSselTimeSeries4devst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries4devst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries4devst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 691 K

# as refPSselTimeSeries1 except PS sel devs more relaxed + t00 + ESS reweighting PS ESS=5 for stationary sel years
tmp <- refPSselTimeSeries5devst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries5devst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries5devst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 912 K

# as refPSselTimeSeries1 except PS sel devs more relaxed + PS sel split between R1 and R2 + t00 + ESS reweighting PS ESS=5 for stationary sel years
tmp <- refPSselTimeSeries6devst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries6devst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries6devst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1001 K



# as refPSselTimeSeries4devst00I1  + sigmaCPUE = 0.1
tmp <- refPSselTimeSeries4devst00I1 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries4devst00I1", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries4devst00I1, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 554 K

# as refPSselTimeSeries6devst00rW except sel devs applied to most fleets for last ~15 years
tmp <- refselTimeSeries7Mostdevst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refselTimeSeries7Mostdevst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refselTimeSeries7Mostdevst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 425 K

# as refPSselTimeSeries6devst00rW except sel devs applied to miscellaneous fleets for last ~15 years
tmp <- refselTimeSeries7Miscdevst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refselTimeSeries7Miscdevst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refselTimeSeries7Miscdevst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY =  985 K

# as refPSselTimeSeries6devst00rW except sel devs applied to all LL fleets for last ~15 years
tmp <- refPSselTimeSeries7LLdevst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries7LLdevst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries7LLdevst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 419 K


# as refPSselTimeSeries6devst00 except sel devs applied to all LL fleets for last ~15 years + no reWeight
tmp <- refPSselTimeSeries7LLdevst00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refPSselTimeSeries7LLdevst00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refPSselTimeSeries7LLdevst00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 784 K

# as refPSselTimeSeries7LLdevst00rW except sel devs applied to LL fleets for last ~30 years
tmp <- refselTimeSeries8LLdevst00rW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refselTimeSeries8LLdevst00rW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refselTimeSeries8LLdevst00rW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 785 K

# as refselTimeSeries8LLdevst00rW except no CL reWeight
tmp <- refselTimeSeries8LLdevst00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refselTimeSeries8LLdevst00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refselTimeSeries8LLdevst00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1082 K

# as refPSselTimeSeries8LLdevst00 except LL sel split by region
tmp <- refselTimeSeries9LLdevst00 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refselTimeSeries9LLdevst00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refselTimeSeries9LLdevst00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 873 K




# for CAPAM paper

#create plots for OM-ref specifications yielding a range of MSY estimates
# refAnalogue <- "R4MvEst_h80_M10_t10_q0_x3_iH_SS"
tmp <- refAnalogue <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\refAnalogue", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnalogue, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 460 K

tmp <- refAnaloguehEst <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\refAnaloguehEst", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguehEst, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 481 K

tmp <- refAnaloguet00 <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\refAnaloguet00", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1275 K

tmp <- refAnaloguet00hEst <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\refAnaloguet00hEst", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00hEst, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 99 K

# as refAnalogue + sigmaCPUE = 0.1
tmp <- refAnalogueI1 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnalogueI1", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnalogueI1, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 369 K

# as refAnalogue + t00 + sigmaCPUE = 0.1
tmp <- refAnaloguet00I1 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00I1", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00I1, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 573 K

# as refAnalogue + t00 + 15y of LL devs (5 shared fleets only)
tmp <- refAnaloguet00LLdev <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LLdev", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLdev, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 477 K

# as refAnalogue + t00 + 3y of LL devs around 2004-6 strange productivity period (LL sel shared)
tmp <- refAnaloguet00LLdevIOI <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LLdevIOI", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLdevIOI, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1146 K

# as refAnalogue + t00 + 5y of LL devs around 2005-9 tag recovery period (LL sel shared)
tmp <- refAnaloguet00LLdevTP <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LLdevTP", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLdevTP, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1121 K

# as refAnalogue + t00 + LL sel independent by area
tmp <- refAnaloguet00LL5 <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LL5", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LL5, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 1277 K

# as refAnalogue + t00 + LL cubic spline
tmp <- refAnaloguet00LLspl <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LLspl", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLspl, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 2287 K

# as refAnalogue + t00 + LL cubic spline +CL reWeight
tmp <- refAnaloguet00LLsplrW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet00LLsplrW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLsplrW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 7711 K

# as refAnalogue (I think) + t00 + 15y LL devs + CL reWeight (from grid gridY17.1LLdevRW)
tmp <- refAnaloguet00LLdrW <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1LLdevRW\\R4MvEst_h80_M10_t0001_q0_x3_iH_LLd_CLRW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LLdrW, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 436 K

# pessimistic run from grid gridY17.1LLdevRW
tmp <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1LLdevRW\\R4MvEst_h80_M06_t0001_q0_x3_iH_LLd_CLRW", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(tmp, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 207 K




# as refAnalogue + t00 + 16 qtr LL sel blocks (5 shared fleets only)
tmp <- refAnaloguet00LL4yBlk <- SS_output(dir = "H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridY17.1\\refAnaloguet004yBlk", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnaloguet00LL4yBlk, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
#MSY = 2070 K







maxMSYmod <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\R4MvEst_h90_M10_t00_q0_x8_iC_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(maxMSYmod, uncertainty=F)
maxMSYmod$derived_quants[maxMSYmod$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000


MSY813mod <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\R4MvEst_h80_M08_t01_q0_x8_iC_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(MSY813mod, uncertainty=F)



# hiMSY - estimate steepness
maxMSYmodhEst  ... fit but maybe not used...




#make some tag fit plots
tt00   <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t00_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tt0001 <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t0001_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tt001  <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t001_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tt01   <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t01_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tt10   <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t10_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
tt15   <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1tagWt\\R4MvEst_h80_M10_t15_q0_x3_iH_SS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(tt00, uncertainty=F)
SS_plots(tt0001, uncertainty=F)
SS_plots(tt001, uncertainty=F)
SS_plots(tt01, uncertainty=F)
SS_plots(tt10, uncertainty=F)
SS_plots(tt15, uncertainty=F)



sum(as.numeric(tt00$likelihoods_raw_by_fleet[10,3:133]))
sum(as.numeric(tt0001$likelihoods_raw_by_fleet[10,3:133]))
sum(as.numeric(tt001$likelihoods_raw_by_fleet[10,3:133]))
sum(as.numeric(tt01$likelihoods_raw_by_fleet[10,3:133]))
sum(as.numeric(tt10$likelihoods_raw_by_fleet[10,3:133]))





#compare the reference case OM with the closest grid analogue  "R4MvEst_h80_M10_t10_q0_x3_iH"
#i.e. same assumptions except for environemental migration
refAnalogue <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\refAnalogue", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(refAnalogue, uncertainty=F)

refAnalogueList  <- "R4MvEst_h80_M10_t10_q0_x3_iH"
refAnalogue <- get(gridY17.1List[gridY17.1List==refAnalogueList])

SS_plots(refAnalogue, uncertainty=F)


plotIndices.f(modList=refAnalogueList, mfrowLayout = c(4,2), MSYyLim=c(0,700))

timeSeriesPlots.f(mList = refAnalogueList, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h70", "h80", "h90",
             "M10", "M08","M06",
             "t10","t01","t00",
             "q0","q1",
             "x3","x8",
             "iC","iH"),
     optWt = rep(1,100))



# confirm non-stationary LL selectivity implementation
# blocks look right in report.SSO (except 110 and 277 are defined for some reasons); SS_plots can't handle results
NSLLTest <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\NSLLTest", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(NSLLTest, uncertainty=F)  #time-varying plots fail

NSLL <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1NS\\R4MvEst_h80_M10_t10_q0_x3_iH_NS", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(NSLL, uncertainty=F)  #time-varying plots fail




selTrendTest <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\selTrendTest", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(selTrendTest, uncertainty=F)

seleLog <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\selSpline", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(seleLog, uncertainty=F)
tmp <- seleLog$parameters[ (seleLog$parameters$Value !="NA") & !is.na(seleLog$parameters$Value),]
tmp <- tmp[(tmp$Value >= tmp$Max) | (tmp$Value <= tmp$Min),c(2,3,6,7)]

selLeneLog <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\selLeneLogistic", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(selLeneLog, uncertainty=F)  #
#not sure which bound errors noted in warnings...
tmp <- selLeneLog$parameters[ (selLeneLog$parameters$Value !="NA") & !is.na(selLeneLog$parameters$Value),]
#tmp <- tmp[tmp$parameters$Value >= tmp$parameters$Max | tmp$parameters$Value <= tmp$parameters$Min,1:3]
tmp <- tmp[(tmp$Value >= tmp$Max) | (tmp$Value <= tmp$Min),c(2,3,6,7)]


selLenLog <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\selLenLogistic", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(selLenLog, uncertainty=F)  #



mvHi <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\mvHi", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(mvHi, uncertainty=F)

mvLo <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\mvLo", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(mvLo, uncertainty=F)

mvHiTrop <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1\\mvHiTrop", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(mvHiTrop, uncertainty=F)


plot(ref2017$likelihoods_used[,1],ref2017$likelihoods_used[,1],col=1, type='l', xlim=c(-200,20000), ylim=c(-200,20000))
points(ref2017$likelihoods_used[,1],refAnalogue$likelihoods_used[,1],col=2)
points(ref2017$likelihoods_used[,1],NSLL$likelihoods_used[,1],col=3,pch=3)
points(ref2017$likelihoods_used[,1],seleLog$likelihoods_used[,1],col=4)
points(ref2017$likelihoods_used[,1],selLeneLog$likelihoods_used[,1],col=5)
points(ref2017$likelihoods_used[,1],refAnalogue$likelihoods_used[,1],col=6,pch=15)
points(ref2017$likelihoods_used[,1],mvHi$likelihoods_used[,1],col=7,pch=4)
points(ref2017$likelihoods_used[,1],mvHiTrop$likelihoods_used[,1],col=8,pch=4)






tmp <- NSLLTest$ageselex[NSLLTest$ageselex$fleet==3 & NSLLTest$ageselex$factor=="Asel",]  #3,7,10,11,13
plot(0:27,as.numeric(tmp[tmp$year==13 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), type='l', col=1,lwd=6)
lines(0:27,as.numeric(tmp[tmp$year==77 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), col=2,lwd=5)
lines(0:27,as.numeric(tmp[tmp$year==117 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), col=3,lwd=4)
lines(0:27,as.numeric(tmp[tmp$year==157 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), col=4,lwd=3)
lines(0:27,as.numeric(tmp[tmp$year==197 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), col=5,lwd=2)
lines(0:27,as.numeric(tmp[tmp$year==237 & tmp$factor=="Asel" & !is.na(tmp[,1]),9:36]), col=6,lwd=1)

for(i in 13:117){
 print(c(i, as.numeric(tmp[tmp$year==i & tmp$factor=="Asel" & !is.na(tmp[,1]),19:22])))
}



##################################################################################
# robustness case Y17.1tagWt - includes upweighted tags
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1tagWt'
gridY17.1tagWtList <- makeGridY17.1tagWt.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1tagWtList)
gridY17.1tagWtList <- makeGridY17.1tagWt.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1tagWt//'
importGrid.f(gridList=gridY17.1tagWtList, gridDir=gridDir, covar=F)

ref <- ref2017
plotIndices2.f(modList=gridY17.1tagWtList, mfrowLayout = c(2,2), MSYyLim=c(0,700),inputRefLines=T)

timeSeriesPlots.f(mList = gridY17.1tagWtList, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h80",
             "M10", "M08","M06",
             "t0001","t001","t01","t10","t15",
             "q0",
             "x3","x8",
             "iH",
             "SS"),
     optWt = rep(1,100))




##################################################################################
# robustness case Y17.1selTrend - includes two forms of time-varying LL selectivity
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1selTrend'
gridY17.1selTrendList <- makeGridY17.1selTrend.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1selTrendList)
gridY17.1selTrendList <- makeGridY17.1selTrend.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1selTrend//'
importGrid.f(gridList=gridY17.1selTrendList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1selTrendList, mfrowLayout = c(4,2), MSYyLim=c(0,700))

timeSeriesPlots.f(mList = gridY17.1selTrendList, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h80",
             "M10", "M08","M06",
             "t01","t10",
             "q0",
             "x3","x8",
             "iH",
             "SS","NS","ST"),
     optWt = rep(1,100))




##################################################################################
# robustness case Y17.1selTrend - includes two forms of time-varying LL selectivity
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1I10LLdev'
gridY17.1I10LLdevList <- makeGridY17.1I10LLdev.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1I10LLdevList)
gridY17.1I10LLdevList <- makeGridY17.1I10LLdev.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1I10LLdev//'
importGrid.f(gridList=gridY17.1I10LLdevList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1I10LLdevList, mfrowLayout = c(4,2), MSYyLim=c(0,700))

timeSeriesPlots.f(mList = gridY17.1I10LLdevList, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h80",
             "M10", "M08","M06",
             "t01","t10",
             "q0",
             "x3","x8",
             "iH","i10H",
             "SS","NSd"),
     optWt = rep(1,100))


##################################################################################
# robustness case with  LL selectivity devs and iterative reweighting of size comp dev CV=0.3
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1LLdevRW'
gridY17.1LLdevRWList <- makeGridY17.1LLdevRW.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1LLdevRWList)
gridY17.1LLdevRWList <- makeGridY17.1LLdevRW.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1LLdevRW//'
importGrid.f(gridList=gridY17.1LLdevRWList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1LLdevRWList, mfrowLayout = c(4,2), MSYyLim=c(0,700))

timeSeriesPlots.f(mList = gridY17.1LLdevRWList, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h80",
             "M10", "M08","M06",
             "t01","t10",
             "q0",
             "x3","x8",
             "iH","i10H",
             "SS","NSd"),
     optWt = rep(1,100))

   # %&% gridY17.1LLdevRWList[1]
tmp <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1LLdevRW\\R4MvEst_h80_M06_t0001_q0_x3_iH_LLd_ess2" , repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(tmp, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
max(tmp$equil_yield$Catch)

tmp <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\gridY17.1LLdevRW\\R4MvEst_h80_M06_t0001_q0_x3_iH_LLd_CLRW" , repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F, forecast=F)
SS_plots(tmp, uncertainty=F)
tmp$derived_quants[tmp$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
max(tmp$equil_yield$Catch)

##################################################################################
# robustness case with  LL selectivity devs and iterative reweighting of size comp dev CV=0.1

path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1LLdevRW'
gridY17.1LLdevRWList <- makeGridY17.1LLdevRW.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1LLdevRWList)
gridY17.1LLdevRWList <- makeGridY17.1LLdevRW.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

#get rid of ESS2
gridY17.1LLdevRWList <- makeGridY17.1LLdevRW.f(path=path, makeGrid=F,
  sp.val  = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  h.val   = c( "h80"),  #SR steepness
  M.val   = c( "M10", "M08","M06"), #mort
  t.val   = c("t10","t0001"),   #tag weight
  q.val   = c("q0"),           #CPUE q trend % per y
  mix.val = c("x3","x8"),            #tag mixing period (qtrs)
  ind.val = c("iH"),            #CPUE standardization method: cluster vs HBF
  sel.val = c("LLd"),           #CPUE standardization method: cluster vs HBF
  ess.val = c("CLRW"))

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1LLdevRW//'
importGrid.f(gridList=gridY17.1LLdevRWList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1LLdevRWList, mfrowLayout = c(4,2), MSYyLim=c(0,700))


##################################################################################
# robustness case with mix of options already shown to reduce MSY in absence of tags
# LL selectivity devs, iterative reweighting of size comp dev, and CPUE CV=0.1

path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1MSY'
gridY17.1MSYList <- makeGridY17.1MSY.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1MSYList)
gridY17.1MSYList <- makeGridY17.1MSY.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

#subset
gridY17.1MSYList <- makeGridY17.1MSY.f(path=path, makeGrid=F,
  sp.val  = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  h.val   = c( "h80"),  #SR steepness
  M.val   = c( "M10", "M08","M06"), #mort
  t.val   = c("t0001"),   #tag weight
  q.val   = c("q0"),           #CPUE q trend % per y
  mix.val = c("x3"),            #tag mixing period (qtrs)
  ind.val = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  sel.val = c("SS","LLd"),           #CPUE standardization method: cluster vs HBF
  ess.val = c("CLRW","ess5"))          #size comp weighting

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1MSY//'
importGrid.f(gridList=gridY17.1MSYList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1MSYList, mfrowLayout = c(4,2), MSYyLim=c(0,700))



##################################################################################
# robustness case to test implication of sigmaR

path=rootDir %&% 'OMconditioning\\YFT\\gridY17.1sigmaR'
gridY17.1sigmaRList <- makeGridY17.1sigmaR.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.1sigmaRList)
gridY17.1sigmaRList <- makeGridY17.1sigmaR.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

#subset
gridY17.1sigmaRList <- makeGridY17.1sigmaR.f(path=path, makeGrid=F,
  sp.val  = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  h.val   = c( "h80"),  #SR steepness
  M.val   = c( "M10", "M08","M06"), #mort
  t.val   = c("t0001"),   #tag weight
  q.val   = c("q0"),           #CPUE q trend % per y
  mix.val = c("x3"),            #tag mixing period (qtrs)
  ind.val = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  sel.val = c("SS","LLd"),           #CPUE standardization method: cluster vs HBF
  ess.val = c("CLRW","ess5"))          #size comp weighting

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.1sigmaR//'
importGrid.f(gridList=gridY17.1sigmaRList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.1sigmaRList, mfrowLayout = c(4,2), MSYyLim=c(0,700))


##################################################################################
# gridY17.2noTag - suite of models without tags for OM-ref2
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.2noTag'
gridY17.2noTagList <- makeGridY17.2noTag.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.2noTagList)
gridY17.2noTagList <- makeGridY17.2noTag.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

#subset
gridY17.2noTagList <- makeGridY17.2noTag.f(path=path, makeGrid=F,
  sp.val  = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  h.val   = c( "h80"),  #SR steepness
  M.val   = c( "M10", "M08","M06"), #mort
  t.val   = c("t0001"),   #tag weight
  q.val   = c("q0"),           #CPUE q trend % per y
  mix.val = c("x3"),            #tag mixing period (qtrs)
  ind.val = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  sel.val = c("SS","LLd"),           #CPUE standardization method: cluster vs HBF
  ess.val = c("CLRW","ess5"))          #size comp weighting

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2noTag//'
importGrid.f(gridList=gridY17.2noTagList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.2noTagList, mfrowLayout = c(4,2), MSYyLim=c(0,700))








##################################################################################
# gridY17.2Tag - suite of models with tags for OM-ref2
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.2Tag'
gridY17.2TagList <- makeGridY17.2Tag.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.2TagList)
gridY17.2TagList <- makeGridY17.2Tag.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

#subset
gridY17.2TagList <- makeGridY17.2Tag.f(path=path, makeGrid=F,
  sp.val  = c( "R4MvEst"),          #Spatial&Population structure - assessment default single option
  h.val   = c( "h80"),  #SR steepness
  M.val   = c( "M10", "M08","M06"), #mort
  t.val   = c("t0001"),   #tag weight
  q.val   = c("q0"),           #CPUE q trend % per y
  mix.val = c("x3"),            #tag mixing period (qtrs)
  ind.val = c("iH","i10H","iC","i10C"),            #CPUE standardization method: cluster vs HBF  sigma CPUE = 0.3 or 0.1
  sel.val = c("SS","LLd"),           #CPUE standardization method: cluster vs HBF
  ess.val = c("CLRW","ess5"))          #size comp weighting

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2Tag//'
importGrid.f(gridList=gridY17.2TagList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.2TagList, mfrowLayout = c(4,2), MSYyLim=c(0,700))



##################################################################################
# Merge and plot two unbalanced grids...

#gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2noTag//'
#importGrid.f(gridList=gridY17.2noTagList, gridDir=gridDir, covar=F)
#gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2Tag//'
#importGrid.f(gridList=gridY17.2TagList, gridDir=gridDir, covar=F)
OMref17.2grid <- c(gridY17.2TagList, gridY17.2noTagList)

ccDat <- cbind(OMref17.2grid,plotIndices2.f(modList=OMref17.2grid, mfrowLayout = c(4,2), MSYyLim=c(0,700)))
colnames(ccDat)[1] <- 'modID'
# drop failed convergence
cDat <- ccDat[as.numeric(ccDat[,'max.Grad']) < -1, ]

cor(as.numeric(cDat[,'BY']),as.numeric(cDat[,'MSY']))

samp2 <- cDat

titleText <- "\nOMref17.2grid"

N <- nrow(samp2)

#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
layout(matrix(c(1,4,2,4,3,5,6,6), 4, 2, byrow = TRUE))
#par(mar = c(3,4,2,1))
par(mar = c(4,4,4,4))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="SSB(2016)/SSB(MSY) Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="SSB(2016)/SSB0 Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

X <-cbind(as.numeric(samp2[,'MSY']),as.numeric(samp2[,'B_B.MSY']))
print(summary(X))
print(cor(X[,1],X[,2]))
plot(X, xlab="MSY", ylab="SSB/SSBMSY", pch=19, cex=0.5, main="Bivariate Sample Plot")
points(jitter(X[,1], amount=mean(X[,1])/50),jitter(X[,2], amount=mean(X[,2])/50), pch=19, cex=0.01, col='grey')

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID']))/sum(table(samp2[,'modID'])), main="Proportion of each model in OM ensemble", xlab="", ylab='Proportion', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))

title(titleText, outer=TRUE)

plotIndices2.f(modList=as.character(samp2[,1]), mfrowLayout = c(4,2), MSYyLim=c(0,700))











###############################################################################################################################
# OM-ref17.2
# 1 stage bivariate sample on basis of MSY and B/BMSY consistency with assessment

filtered <- gridSamplerBivar(
  gridListPlus = cDat, # model definitions and factors of interest for sampling
  sampleQuant = c('MSY','B_B.MSY'),        # quantities upon which to base the sample
  corr =0., # correlation among characteristics; works but appropriate value not determined

  FSList = list(             # subsample models to achieve the following proportions (specifying too many interactions will have an adverse effect)
  #  FSOpt1 <- c('t0001', 0.5),
  #  FSOpt2 <- c('t10', 0.5)),

  #  FSOpt1 <- c('t0001', 1.)),
  #  FSOpt2 <- c('t10', 0.5)),

  #  FSOpt1 <- c('q0', 0.5),
  #  FSOpt2 <- c('q1', 0.5)),

    FSOpt1 <- c('t0001', 'q0',0.25),
    FSOpt2 <- c('t0001', 'q1',0.25),
    FSOpt3 <- c('t10', 'q0',0.25),
    FSOpt4 <- c('t10', 'q1',0.25)),
  mu    = c(422,0.89),  #distribution mean to attain (roughly)
  #mu    = 0.89,  # B/BMSY mean to attain (roughly)
  sigma = c(3.*refMSYcv,3.*refBYoBMSYcv),  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3, #truncate sample distribution at this many sigma
  nBins = 3, #number of bins covering distribution is a grid of nBins X nBins
  logNorm=T)


OMref17.2List <- filtered[[2]]
details <- filtered[[1]]
details[,1] <- as.factor(details[,1])
samp2 <- details[sample(details[,1], size=2000, replace=T, prob=details[,'wtModList']),]
samp2[,1] <- as.character(samp2[,1])

titleText <- "\n OM-ref17.2"
#titleText <- "\n Alternate OM with sampling CVs 2X SA-ref"

N <- nrow(samp2)

#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
layout(matrix(c(1,4,2,4,3,5,6,6), 4, 2, byrow = TRUE))
#par(mar = c(3,4,2,1))
par(mar = c(4,4,4,4))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="SSB(2016)/SSB(MSY) Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="SSB(2016)/SSB0 Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

X <-cbind(as.numeric(samp2[,'MSY']),as.numeric(samp2[,'B_B.MSY']))
print(summary(X))
print(cor(X[,1],X[,2]))
plot(X, xlab="MSY", ylab="SSB/SSBMSY", pch=19, cex=0.5, main="Bivariate Sample Plot")
points(jitter(X[,1], amount=mean(X[,1])/50),jitter(X[,2], amount=mean(X[,2])/50), pch=19, cex=0.01, col='grey')

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID']))/sum(table(samp2[,'modID'])), main="Proportion of each model in OM ensemble", xlab="", ylab='Proportion', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))

title(titleText, outer=TRUE)

plotIndices2.f(modList=as.character(samp2[,1]), mfrowLayout = c(4,2), MSYyLim=c(0,700))






###############################################################################################################################
# alternate OM,
# 1 stage bivariate sample on basis of MSY and B/BMSY consistency with assessment

filtered <- gridSamplerBivar(
  gridListPlus = cDat, # model definitions and factors of interest for sampling
  sampleQuant = c('MSY','B_B.MSY'),        # quantities upon which to base the sample
  corr =0., # correlation among characteristics; works but appropriate value not determined

  FSList = list(             # subsample models to achieve the following proportions (specifying too many interactions will have an adverse effect)
  #  FSOpt1 <- c('t0001', 0.5),
  #  FSOpt2 <- c('t10', 0.5)),

  #  FSOpt1 <- c('t0001', 1.)),
  #  FSOpt2 <- c('t10', 0.5)),

  #  FSOpt1 <- c('q0', 0.5),
  #  FSOpt2 <- c('q1', 0.5)),

    FSOpt1 <- c('t0001', 'q0',0.25),
    FSOpt2 <- c('t0001', 'q1',0.25),
    FSOpt3 <- c('t10', 'q0',0.25),
    FSOpt4 <- c('t10', 'q1',0.25)),
  mu    = c(422,0.89),  #distribution mean to attain (roughly)
  #mu    = 0.89,  # B/BMSY mean to attain (roughly)
  sigma = c(3.*refMSYcv,3.*refBYoBMSYcv),  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3, #truncate sample distribution at this many sigma
  nBins = 3, #number of bins covering distribution is a grid of nBins X nBins
  logNorm=T)


OMref17.2List <- filtered[[2]]
details <- filtered[[1]]
details[,1] <- as.factor(details[,1])
samp2 <- details[sample(details[,1], size=2000, replace=T, prob=details[,'wtModList']),]
samp2[,1] <- as.character(samp2[,1])

titleText <- "\n OM-ref17.2"
#titleText <- "\n Alternate OM with sampling CVs 2X SA-ref"

N <- nrow(samp2)

#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
layout(matrix(c(1,4,2,4,3,5,6,6), 4, 2, byrow = TRUE))
#par(mar = c(3,4,2,1))
par(mar = c(4,4,4,4))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="SSB(2016)/SSB(MSY) Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="SSB(2016)/SSB0 Distribution \n OM Ensemble; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

X <-cbind(as.numeric(samp2[,'MSY']),as.numeric(samp2[,'B_B.MSY']))
print(summary(X))
print(cor(X[,1],X[,2]))
plot(X, xlab="MSY", ylab="SSB/SSBMSY", pch=19, cex=0.5, main="Bivariate Sample Plot")
points(jitter(X[,1], amount=mean(X[,1])/50),jitter(X[,2], amount=mean(X[,2])/50), pch=19, cex=0.01, col='grey')

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID']))/sum(table(samp2[,'modID'])), main="Proportion of each model in OM ensemble", xlab="", ylab='Proportion', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))

title(titleText, outer=TRUE)

plotIndices2.f(modList=as.character(samp2[,1]), mfrowLayout = c(4,2), MSYyLim=c(0,700))








##################################################################################
#explore some plausibility weighing/sampling options from OM-ref1 & OM-rob-MSY
##################################################################################
#ccDat <- plotIndices2.f(modList=c(gridY17.1MSYList), mfrowLayout = c(4,2), MSYyLim=c(0,700))
#gridList <- gridY17.1List
gridList <- gridY17.2noTagList
gridList <- c(gridY17.2noTagList, gridY17.2TagList)
ccDat <- cbind(gridList,plotIndices2.f(modList=gridList, mfrowLayout = c(4,2), MSYyLim=c(0,700)))
colnames(ccDat)[1] <- 'modID'
# drop failed convergence
cDat <- ccDat[as.numeric(ccDat[,'max.Grad']) < -1, ]

par(mfrow=c(3,1))



# impose prior based on recruitment trend

#reference case slope = -0.0005232924
nRecDevObs <- length(ref2017$recruit$dev[ref2017$recruit$era == "Main"])    #=173
coefficients(lm(ref2017$recruit$dev[ref2017$recruit$era == "Main"] ~ c(1:nRecDevObs)))[2]
#reference case rho = 0.23
cor(ref2017$recruit$dev[ref2017$recruit$era == "Main"][1:(nRecDevObs-1)],ref2017$recruit$dev[ref2017$recruit$era == "Main"][2:(nRecDevObs)])
# ref case rec penalty  = -4.40666e+01
dev <- ref2017$recruit$dev[ref2017$recruit$era == "Main"]
sigR <- 0.6
sum(-log(1/sqrt(2*pi*sigR^2)) + ((dev)^2 - 0.5*sigR^2)/(2*sigR^2))
sum(log(sqrt(2*pi*sigR^2)) + ((dev)^2 - 0.5*sigR^2)/(2*sigR^2))
#72.01249 - difference presumably arises from handling of constant terms and late era constraints and ln bias corrections
sum(log(sqrt(2*pi*sigR^2)) + ((dev)^2 )/(2*sigR^2))


# calculate likelihood of observing recruitment trends for given sigmaR and auto-correlation
par(mfrow=c(2,3), xpd=F)
sigmaR <- 0.6
sd <-NULL
sdSigmaR <- NULL
recPen   <- NULL
#rhoList <- c(0.,0.25,0.5)
rhoList <- c(0.25)
for(rho in rhoList){
  T <- nRecDevObs
  slope <- NULL
  sigmaRSim <- NULL
  for(i in 1:2000){
    omega <- rnorm(T)*0.6
    tau <- omega
    for (t in 2:T){
      tau[t] <- rho*tau[t-1]+omega[t]*sqrt(1-rho^2)
    }
    rSlope <- coefficients(lm(tau ~ c(1:T)))[2]
    slope <- c(rSlope, slope)
    sigmaRSim <- c(sigmaRSim, sd(tau))
    recPen <- c(recPen,sum(-log(2*sigmaR) + tau^2/(2*sigmaR^2)))
  }
  hist(slope,c(-100:1000)*0.0001, xlim=c(-0.005, 0.005), main='Distribution of slopes from 200 timesteps; sigma = 0.6, rho = ' %&% rho, xlab="Trend")
  hist(sigmaRSim, main='Distribution of sigmaR from 200 timesteps; sigma = 0.6, rho = ' %&% rho, xlab="sigmaR")
  hist(recPen, main='Distribution of rec Penalties from 200 timesteps; sigma = 0.6, rho = ' %&% rho, xlab="Recruitment Penalty")

  sd <- c(sd,sd(slope))
  sdSigmaR <- c(sdSigmaR,sd(sigmaRSim))

} #rho loop
sdrTrend25 <- sd

#distribution of MSY, followed by distribution with prior weighting from expected recruitment trend with different levels of rho
summary(cDat[,'rec.Trend'])
hist(sample(cDat[,'MSY'],size=2000,prob=1,replace=T), breaks = c(0:50)*100)
par(mfcol=c(4,2), xpd=F)
pDat <- sample(cDat[,'MSY'],size=2000,replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
#hist(pDat, breaks = c(20:80)*10, main="Full OM-ref1 Grid", xlab="MSY")
hist(pDat, breaks = c(20:80)*10, main="Full OM-rob-MSY Grid", xlab="MSY")

pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(cDat[,'rec.Trend'])^2/(2*sd[1]^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
hist(pDat, breaks = c(20:80)*10, main="Sample - Rec trend prior weight \n rho = " %&% rhoList[1] , xlab="MSY")


pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(cDat[,'rec.Trend'])^2/(2*sd[2]^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
hist(pDat, breaks = c(20:80)*10, main="Sample - Rec trend prior weight \n rho = " %&% rhoList[2] , xlab="MSY")

pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(cDat[,'rec.Trend'])^2/(2*sd[3]^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
hist(pDat, breaks = c(20:80)*10, main="Sample - Rec trend prior weight \n rho = " %&% rhoList[3] , xlab="MSY")



# followed by distribution with prior weighting from expected rec RMSE
pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(cDat[,'rec.RMSE']-mu)^2/(2*sdSigmaR[1]^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
hist(pDat, breaks = c(20:80)*10, main="Sample - Sigma(R) prior weight", xlab="MSY")




# followed by distribution with log-normal prior weighting from MSY +/- error from OM-ref

cvMult = 1
pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(log(cDat[,'MSY']/refMSY))^2/(2*(cvMult*refMSYcv)^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
h <- hist(pDat, breaks = c(20:80)*10, main="Sample - assessment log(MSY) prior weight \n CV Multiplier = " %&% cvMult, xlab="MSY")
lines(1:1000, max(h$counts)*dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)/max(dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)), col=3)

cvMult = 3
pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(log(cDat[,'MSY']/refMSY))^2/(2*(cvMult*refMSYcv)^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
h <- hist(pDat, breaks = c(20:80)*10, main="Sample - assessment log(MSY) prior weight \n CV Multiplier = " %&% cvMult, xlab="MSY")
lines(1:1000, max(h$counts)*dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)/max(dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)), col=3)

cvMult = 5
pDat <- sample(cDat[,'MSY'],size=2000,prob=exp(-(log(cDat[,'MSY']/refMSY))^2/(2*(cvMult*refMSYcv)^2)),replace=T)
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201
h <- hist(pDat, breaks = c(20:80)*10, main="Sample - assessment log(MSY) prior weight \n CV Multiplier = " %&% cvMult, xlab="MSY")
lines(1:1000, max(h$counts)*dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)/max(dnorm(log(1:1000),log(refMSY), sd=cvMult*refMSYcv)), col=3)




# just sample the OM distribution to match the desired MSY distribution
# add sub-sample constraint hierarchy
FSList <- list(
  FSOpt1 <- c('t0001', 'q0',500),
  FSOpt2 <- c('t0001', 'q1',500),
  FSOpt3 <- c('t10', 'q0',500),
  FSOpt4 <- c('t10', 'q1',500)
)

par(mfrow=c(1,1), xpd=T)
cvMult <- 2
sigma <- cvMult*refMSYcv
mu    <- log(refMSY) - 0.5*sigma^2
MSYsamp <- NULL

for(iFS in 1:length(FSList)){
  n <- FSList[[iFS]][length(FSList[[iFS]])]
  modNameList <- strsplit(cDat[,'modID'], split="_")
  keepMods   <- rep(NA,length(modNameList))
  factorList <- matrix(nrow=length(modNameList), ncol=length(modNameList[[1]]),data=unlist(modNameList), byrow=T)
  keepModsList <- FSList[[iFS]][1:(length(FSList[[iFS]])-1)]
  for (i in 1:length(keepMods)){
    keepMods[i] <- prod(keepModsList %in% factorList[i,])
  }
  keepMods <- keepMods==1
  subModList <- cDat[keepMods,]
  prDat <- pnorm(log(as.numeric(subModList[,'MSY'])), mu, sigma)

  for (i in 1:n){
    #rnorm(1,log(refMSY) - 0.5*sigma^2, sd=sigma)
    sdbin <- floor(rnorm(1))
    while (sdbin < -3 | sdbin >2){
      sdbin <- floor(rnorm(1))
    }
    pbin1 <- pnorm(sdbin)
    pbin2 <- pnorm(sdbin+1)

    subDat <- subModList[prDat>pbin1 & prDat <= pbin2,]
    #MSYsamp <- rbind(MSYsamp,sample(as.numeric(cDat[prDat>pbin1 & prDat <= pbin2,c('modID','MSY')]), size=1, prob=prDat[prDat>pbin1 & prDat <= pbin2]))
    MSYsamp <- rbind(MSYsamp,subDat[sample(nrow(subDat), size=1, prob=prDat[prDat>pbin1 & prDat <= pbin2]), ])
    #print(c(sdbin,pbin1,pbin2,MSYsamp))
  }
}
N <- nrow(MSYsamp)

par(mfrow=c(4,1))
h <- hist(as.numeric(MSYsamp[,'MSY']), breaks = c(20:80)*10,
  main="Weighted OM Ensemble MSY Distribution \n OM-ref2; CV Multiplier = " %&% cvMult %&% "; nMods = " %&% length(unique(MSYsamp[,'modID'])),
  xlab="MSY", probability=T)
lines(200:800, dnorm(log(200:800), mu, sd=sigma)/sum(dnorm(log(200:800), mu, sd=sigma)), col=3)
points(exp(c(mu)), c(0), col=2)
summary(as.numeric(MSYsamp[,'MSY']))

h <- hist(as.numeric(MSYsamp[,'B_B.MSY']), breaks = c(0:60)/20,
  main="Weighted OM Ensemble B(2016)/B(MSY) Distribution \n OM-ref2; CV Multiplier = " %&% cvMult %&% "; nMods = " %&% length(unique(MSYsamp[,'modID'])),
  xlab="MSY", probability=T)
points(c(refB_B.MSY), c(0), col=2)
summary(as.numeric(MSYsamp[,'B_B.MSY']))


gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(MSYsamp[,'modID'], "_")))
plot(table(sampledOptions)/N,type='p', col=2, main = "Grid option frequency before and after sub-sampling", ylab="Proportion")
lines(table(gridOptions)/length(gridList))

plot(table(MSYsamp[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')



###############################################################################################################################
# 1 stage sample OM grid on basis of MSY consistency with assessment only

cvMSY <- 1.*refMSYcv
MSYsamp <- gridSampler(
  gridListPlus    = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = refMSY,  #distribution mean to attain (roughly)
  sigma = cvMSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3 #truncate sample distribution at this many sigma
)

samp2 <- MSYsamp
N <- nrow(samp2)

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
par(mar = c(3,4,2,1))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n  OM with MSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="BY/BMSY Distribution \n  OM with MSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID'])), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="BY/B0 Distribution \n  OM with MSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2015)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))
par(def.par)




###############################################################################################################################
#1 stage: sample OM grid on basis of B/BMSY consistency with assessment only

sigmaB_B.MSY <- 1.*refBYoBMSYcv
B_B.MSYsamp <- gridSampler(
  gridListPlus    = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'B_B.MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = refB_B.MSY,  #distribution mean to attain (roughly)
  sigma = sigmaB_B.MSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3 #truncate sample distribution at this many sigma
)

samp2 <- B_B.MSYsamp
N <- nrow(samp2)

#create plots for merged, uhweighted grid...
# par(mfrow=c(4,1))

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
par(mar = c(3,4,2,1))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n  OM with BY/BMSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="BY/BMSY Distribution \n  OM with BY/BMSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID'])), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="BY/B0 Distribution \n  OM with BY/BMSY sampling only;  nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2015)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))
par(def.par)



###############################################################################################################################
# 2 stage - sample OM grid on basis of MSY consistency with assessment, then on basis of B/BMSY

#stage 1 option 1 - sample for an MSY dist

cvMSY <- 3*refMSYcv # reference case SA STD
MSYsamp <- gridSampler(
  gridListPlus    = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',5000),
    FSOpt2 <- c('t0001', 'q1',5000),
    FSOpt3 <- c('t10', 'q0',5000),
    FSOpt4 <- c('t10', 'q1',5000)),
  mu    = refMSY,  #distribution mean to attain (roughly)
  sigma = cvMSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3 #truncate sample distribution at this many sigma
)

#stage 1 option 2 - just truncate the MSY dist
#sigmaTrunc <- 3
#MSYsamp <- cDat[(log(as.numeric(cDat[,"MSY"])) > log(refMSY*(1-sigmaTrunc*cvMSY))) & (log(as.numeric(cDat[,"MSY"])) < log(refMSY*(1+sigmaTrunc*cvMSY))),  ]
samp2 <- MSYsamp

# for B/BMSY CV calculation = 0.59
cor(as.numeric(MSYsamp[,'BY']),as.numeric(MSYsamp[,'BMSY']))



#refB_B.MSY <- 0.89
sigmaB_B.MSY <- 3*refBYoBMSYcv
samp2 <- gridSampler(
  gridListPlus    = MSYsamp, # model definitions and factors of interest for sampling
  sampleQuant = 'B_B.MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),

#  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
#    FSOpt1 <- c('t0001', 'h70',167),
#    FSOpt3 <- c('t10',   'h70',167),
#    FSOpt5 <- c('t0001', 'h80',166),
#    FSOpt7 <- c('t10',   'h80',166),
#    FSOpt9 <- c('t0001', 'h90',167),
#    FSOpt11 <- c('t10',  'h90',167)),

  mu    = refB_B.MSY,  #distribution mean to attain (roughly)
  sigma = sigmaB_B.MSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 2#truncate sample distribution at this many sigma
)

N <- nrow(samp2)

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
par(mar = c(3,4,2,1))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
  main="MSY Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
  main="SSB(2016)/SSB(MSY) Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID'])), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
  main="SSB(2016)/SSB0 Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))
par(def.par)

plotIndices2.f(modList=OMref17.2List, mfrowLayout = c(4,2), MSYyLim=c(0,700))













###############################################################################################################################
# 2 stage - sample OM grid on basis of MSY consistency with assessment, then on basis of B/BMSY

#stage 1 option 1 - sample for an MSY dist

cvMSY <- 2*refMSYcv # reference case SA STD
MSYsamp <- gridSampler(
  gridListPlus    = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',5000),
    FSOpt2 <- c('t0001', 'q1',5000),
    FSOpt3 <- c('t10', 'q0',5000),
    FSOpt4 <- c('t10', 'q1',5000)),
  mu    = refMSY,  #distribution mean to attain (roughly)
  sigma = cvMSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3 #truncate sample distribution at this many sigma
)

#stage 1 option 2 - just truncate the MSY dist
#sigmaTrunc <- 3
#MSYsamp <- cDat[(log(as.numeric(cDat[,"MSY"])) > log(refMSY*(1-sigmaTrunc*cvMSY))) & (log(as.numeric(cDat[,"MSY"])) < log(refMSY*(1+sigmaTrunc*cvMSY))),  ]
samp2 <- MSYsamp

# for B/BMSY CV calculation = 0.59
cor(as.numeric(MSYsamp[,'BY']),as.numeric(MSYsamp[,'BMSY']))



#refB_B.MSY <- 0.89
sigmaB_B.MSY <- 2*refBYoBMSYcv
samp2 <- gridSampler(
  gridListPlus    = MSYsamp, # model definitions and factors of interest for sampling
  sampleQuant = 'B_B.MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = refB_B.MSY,  #distribution mean to attain (roughly)
  sigma = sigmaB_B.MSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 2#truncate sample distribution at this many sigma
)

N <- nrow(samp2)

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
par(mar = c(3,4,2,1))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(20:80)*10,
#  main="MSY Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  main="SSB(2016)/SSB(MSY) Distribution \n Alternative Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:40)/20,
#  main="SSB(2016)/SSB(MSY) Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  main="SSB(2016)/SSB(MSY) Distribution \n Alternative Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID'])), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:40)/40,
#  main="SSB(2016)/SSB0 Distribution \n Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  main="SSB(2016)/SSB(MSY) Distribution \n Alternative Weighted OM Ensemble OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))
par(def.par)





###############################################################################################################################
# 2 stage - sample OM grid on basis of recTrend consistency with simulated recTrend sigR=0.6, rho = 0.25, then on basis of B/BMSY

sdrTrend <- 0.001176811
#sdrTrend <- 3*sdrTrend
recTrendsamp <- gridSampler(
  gridListPlus    = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'rec.Trend',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = 0,  #distribution mean to attain (roughly)
  sigma = sdrTrend,  # sampling sigma (~CV)
  sigmaTrunc = 2, #truncate sample distribution at this many sigma
  logNorm =F
)
#samp2 <- recTrendsamp

refB_B0 <- 0.29
refB_B.MSY <- 0.89
sigmaB_B.MSY <- 0.1
samp2 <- gridSampler(
  gridListPlus    = recTrendsamp, # model definitions and factors of interest for sampling
  sampleQuant = 'B_B.MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = refB_B.MSY,  #distribution mean to attain (roughly)
  sigma = sigmaB_B.MSY,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 2 #truncate sample distribution at this many sigma
)

N <- nrow(samp2)

layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
par(mar = c(3,4,2,1))

pDat <- as.numeric(samp2[,'MSY'])
pDat[pDat>799] <- 799
pDat[pDat<201] <- 201

#h <- hist(pDat, breaks = c(200:800),
#h <- hist(pDat, breaks = c(20:80)*10,
h <- hist(pDat, breaks = c(40:160)*5,
  main="Weighted OM Ensemble MSY Distribution \n OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="MSY", probability=T)
#lines(200:800, dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)/sum(dnorm(log(200:800), refB_B.MSY, sd=sigmaB_B.MSY)), col=3)
points(c(refMSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'MSY']))

#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:300)/100,
#h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:30)/10,
h <- hist(as.numeric(samp2[,'B_B.MSY']), breaks = c(0:60)/20,
  main="Weighted OM Ensemble SSB(2016)/SSB(MSY) Distribution \n OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B(MSY)", probability=T)
points(c(refB_B.MSY), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B.MSY']))

h <- hist(as.numeric(samp2[,'B_B0']), breaks = c(0:60)/20,
  main="Weighted OM Ensemble SSB(2016)/SSB0 Distribution \n OM-ref17.2; nMods = " %&% length(unique(samp2[,'modID'])),
  xlab="B(2016)/B0", probability=T)
points(c(refB_B0), c(0), col=2, pch=15)
summary(as.numeric(samp2[,'B_B0']))

#plot(table(samp2[,'modID']), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
plot(sort(table(samp2[,'modID'])), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')
#plot(cumsum(sort(table(samp2[,'modID']))), main="Frequency of individual models in OM ensemble", xlab="", ylab='Frequency', xaxt='n')


gridOptions <- unlist(c(strsplit(gridList, "_")))
sampledOptions <- unlist(c(strsplit(samp2[,'modID'], "_")))
plot(table(sampledOptions)/sum(table(sampledOptions)),type='p', col=3, pch=15, main = "Grid option frequency before and after filtering", ylab="Proportion")
lines(table(gridOptions)/sum(table(gridOptions)))
par(def.par)











# function to sample an OM grid to achieve particular characteristics in the sampled distribution
gridSampler <- function(
  gridListPlus = cDat, # model definitions and factors of interest for sampling
  sampleQuant = 'MSY',        # quantitiy upon which to base the sample
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 'q0',500),
    FSOpt2 <- c('t0001', 'q1',500),
    FSOpt3 <- c('t10', 'q0',500),
    FSOpt4 <- c('t10', 'q1',500)),
  mu    = 422,  #distribution mean to attain (roughly)
  #mu    = 0.89,  # B/BMSY mean to attain (roughly)
  sigma = 0.11,  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3, #truncate sample distribution at this many sigma
  logNorm=T
){

  samp <- NULL

  for(iFS in 1:length(FSList)){

    nSub <- FSList[[iFS]][length(FSList[[iFS]])]
    modNameList <- strsplit(gridListPlus[,'modID'], split="_")
    keepMods   <- rep(NA,length(modNameList))
    factorList <- matrix(nrow=length(modNameList), ncol=length(modNameList[[1]]),data=unlist(modNameList), byrow=T)
    keepModsList <- FSList[[iFS]][1:(length(FSList[[iFS]])-1)]
    for (i in 1:length(keepMods)){
      keepMods[i] <- prod(keepModsList %in% factorList[i,])
    }
    keepMods <- keepMods==1
    subModList <- gridListPlus[keepMods,]
    if(logNorm){
      prDat <- pnorm(log(as.numeric(subModList[,sampleQuant])), log(mu) - 0.5*sigma^2, sigma)
    } else{
      prDat <- pnorm(as.numeric(subModList[,sampleQuant]), mu, sigma)
    }
    for (i in 1:nSub){
      #rnorm(1,log(refMSY) - 0.5*sigma^2, sd=sigma)
      sdbin <- floor(rnorm(1))
      while (sdbin < -(sigmaTrunc) | sdbin >(sigmaTrunc-1)){
        sdbin <- floor(rnorm(1))
      }
      pbin1 <- pnorm(sdbin)
      pbin2 <- pnorm(sdbin+1)

      subDat <- subModList[prDat>pbin1 & prDat <= pbin2,]
      #MSYsamp <- rbind(MSYsamp,sample(as.numeric(cDat[prDat>pbin1 & prDat <= pbin2,c('modID','MSY')]), size=1, prob=prDat[prDat>pbin1 & prDat <= pbin2]))
      samp <- rbind(samp,subDat[sample(nrow(subDat), size=1, prob=prDat[prDat>pbin1 & prDat <= pbin2]), ])
      #print(c(sdbin,pbin1,pbin2,MSYsamp))
    }
  }
  return(samp)
}


# function to take a bivariate sample of OM models from a grid...
# returns ModList with probability weighitngs...
gridSamplerBivar <- function(
  gridListPlus = cDat, # model definitions and factors of interest for sampling
  sampleQuant = c('MSY','BYoBMSY'),        # quantities upon which to base the sample
  corr =0, # correlation among characteristics (only 0 implemented)
  FSList = list(             # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 0.5),
    FSOpt2 <- c('t10', 0.5)),
  mu    = c(422,0.89),  #distribution mean to attain (roughly)
  #mu    = 0.89,  # B/BMSY mean to attain (roughly)
  sigma = c(0.1,0.1),  #lognormal sampling sigma (~CV)
  sigmaTrunc = 3, #truncate sample distribution at this many sigma
  logNorm=T,
  nBins =3
){

  N <- 10000 # calculated proportions are based on stochastic sampling, so this number should be reasonably big
  samp <- NULL
  wtModList <- NULL

  for(iFS in 1:length(FSList)){ #loop over pre-specified sample frequency requirements

    nSub <- N*as.numeric(FSList[[iFS]][length(FSList[[iFS]])])
    modNameList <- strsplit(gridListPlus[,'modID'], split="_")
    keepMods   <- rep(NA,length(modNameList))
    factorList <- matrix(nrow=length(modNameList), ncol=length(modNameList[[1]]),data=unlist(modNameList), byrow=T)
    keepModsList <- FSList[[iFS]][1:(length(FSList[[iFS]])-1)]
    for (i in 1:length(keepMods)){
      keepMods[i] <- prod(keepModsList %in% factorList[i,])
    }
    keepMods <- keepMods==1
    subModList <- gridListPlus[keepMods,]
    #nBins    <- 3 #each dim
    binWidth <- 2*sigmaTrunc/nBins # in SD units

    sump <- 0
    subSamp <- NULL
    for(x in 1:nBins){
      for(y in 1:nBins){
        #upper and lower bounds of each strata in SD units
        sdx <- c( -(0.5*nBins*binWidth) + (x-1)*binWidth, -(0.5*nBins*binWidth) + x*binWidth )
        sdy <- c( -(0.5*nBins*binWidth) + (y-1)*binWidth, -(0.5*nBins*binWidth) + y*binWidth )

        #probability in strata bivariate normal with correlation
        p <- pmvnorm(lower <- c(sdx[1],sdy[1]), upper <- c(sdx[2],sdy[2]), mean=c(0,0), corr=matrix(c(1,corr,corr,1),2,2))
        sump <- sump+p

        #subset of models that fall in range
        sdMod1 <- (log(as.numeric(subModList[,sampleQuant[1]]))- log(mu[1]) - 0.5*sigma[1]^2 )/sigma[1]
        sdMod2 <- (log(as.numeric(subModList[,sampleQuant[2]]))- log(mu[2]) - 0.5*sigma[2]^2 )/sigma[2]
        subDat <- subModList[(sdMod1 > sdx[1] & sdMod1 < sdx[2]) & (sdMod2 > sdy[1] & sdMod2 < sdy[2]),]

        nSub2 <- nSub*p

        #print(subDat)
        #print(c(p,x,y, sdx, sdy))

        if(!is.null(nrow(subDat))){
          subSamp <- rbind(subSamp, subDat[sample(nrow(subDat), size=nSub2, replace=T), ])
        } else {
          print(c("empty probablity strata: ", sdx, sdy, FSList[[iFS]]))
        }
      }
    }
    samp <- rbind(samp, subSamp)
    print(c("sum P: ", sump))

    #calculate model frequency and convert to proportiosn
    tmp <- table(subSamp[,'modID'])/sum(table(subSamp[,'modID']))
    #reweight the model proportions to account for missing sampling strata, and the expected proportion represented by iFS
    wtModList <- c(wtModList,tmp*nSub/N)

    print(c(FSList[[iFS]], sum(wtModList),dim(tmp)))
    print(names(tmp))

  } #iFS
  #X <-cbind(as.numeric(samp[,sampleQuant[1]]),as.numeric(samp[,sampleQuant[2]]))
  #print(summary(X))
  #print(cor(X[,1],X[,2]))
  #X <- jitter(X)
  #plot(X, xlab="MSY", ylab="SSB/SSBMSY", pch=19, cex=0.01)
  #points(jitter(X[,1], amount=mean(X[,1])/50),jitter(X[,2], amount=mean(X[,2])/50), xlab="MSY", ylab="SSB/SSBMSY", pch=19, cex=0.01, col='grey')
  #z <- kde2d(X[,1], X[,2], n=10)
  #contour(z, drawlabels=FALSE, nlevels=5, add=TRUE, lwd=3, col='white')
  #contour(z, drawlabels=FALSE, nlevels=5, add=TRUE)
  #hist(X[,1],breaks=50)
  #hist(X[,2],breaks=50)
  samp <- unlist(gridListPlus)
  samp <- as.data.frame(samp)

  tmp <- cbind(names(wtModList),wtModList)

  samp <- merge(samp,tmp, by.x=1, by.y=1)
  samp[,2:ncol(samp)] <- apply(samp[,2:ncol(samp)], FUN=as.numeric, MARGIN=c(1,2))
  samp[,1] <- as.character(samp[,1])

  return(list(samp,wtModList))
}

library(MASS)


#Covariance test ... SPB_MSY is not reported in covar matrix; so use sims...
require(mvtnorm)

ref2017SSBCurSD/ref2017SSBCurrent
ref2017SSBMSYSD/ref2017SSBMSY
Bcor <- 0.59

# rmvnorm sigma arguement is covariance matrix
x1 <- rmvnorm(n=10000, mean=c(ref2017SSBCurrent,ref2017SSBMSY), sigma=matrix(c(ref2017SSBCurSD^2,0.0,0.0,ref2017SSBMSYSD^2),2,2))
x2 <- rmvnorm(n=10000, mean=c(ref2017SSBCurrent,ref2017SSBMSY), sigma=matrix(c(ref2017SSBCurSD^2,ref2017SSBCurSD*ref2017SSBMSYSD,ref2017SSBCurSD*ref2017SSBMSYSD,ref2017SSBMSYSD^2),2,2))
#x3 <- rmvnorm(n=10000, mean=c(ref2017SSBCurrent,ref2017SSBMSY), sigma=matrix(c(ref2017SSBCurSD^2,-ref2017SSBCurSD*ref2017SSBMSYSD,-ref2017SSBCurSD*ref2017SSBMSYSD,ref2017SSBMSYSD^2),2,2))
x3 <- rmvnorm(n=10000, mean=c(ref2017SSBCurrent,ref2017SSBMSY), sigma=matrix(c(ref2017SSBCurSD^2,Bcor*ref2017SSBCurSD*ref2017SSBMSYSD,Bcor*ref2017SSBCurSD*ref2017SSBMSYSD,ref2017SSBMSYSD^2),2,2))
cor(x1[,1],x1[,2])
cor(x2[,1],x2[,2])
cor(x3[,1],x3[,2])


ratio <- cbind(x1[,1]/x1[,2],x2[,1]/x2[,2] ,x3[,1]/x3[,2] )
boxplot(ratio)
sd(ratio[,1])/(ref2017SSBCurrent/ref2017SSBMSY)
sd(ratio[,2])/(ref2017SSBCurrent/ref2017SSBMSY)
sd(ratio[,3])/(ref2017SSBCurrent/ref2017SSBMSY)




