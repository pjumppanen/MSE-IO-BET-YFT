#=============================================================================================================================
# Basic R script for making and inspecting the conditioned demonstration case YFT Operating Model OMyft1
#=============================================================================================================================

rootDir <- "H:\\C-offline\\MSE-IO-BET-YFT\\"  #modify for local path

library(TinnRcom)                            #if using TinnR IDE
library (r4ss)                                #R package supporting Stock Synthesis

source(paste(rootDir,"YFT-MSE\\Source\\pasteOperator.R",sep=""))
source(rootDir %&% "OMconditioning\\RStuff\\seasAsYrToDecYr.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\makeGridY3.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\importGrid.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\plotIndices.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\timeSeriesPlots.f.R")


path=rootDir %&% 'OMconditioning\\YFT\\gridY3'
gridY3List <- makeGridY3.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY3List)
#gridY3List <- makeGridY3.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file

# exit this R script and run the SS batch file(s)


# import some summary results & diagnostics from the gridY3List of models
gridDir <- rootDir %&% 'OMconditioning//YFT//gridY3//'
importGrid.f(gridList=gridY3List, gridDir=gridDir, covar=F)

#=============================================================================================================================
# Some results from the Langley et al 2015 assessment for comparison

#Langley 2015 YFT assessment ref case (possibly the pre-release)
ref <- SS_output(dir = rootDir %&% "OMconditioning\\YFT\\LangleyFiles\\ref27Aug", repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=T, forecast=F)
#import and examinde detailed results from a specific SS model configuration
SS_plots(ref, uncertainty=T)

#use boxplots to plot some key model results and diagnostics, marginalized by grid Options
plotIndices.f(modList=gridY3List, mfrowLayout = c(3,2), MSYyLim=c(0,700))

#plot some time series using ref as a reference line
timeSeriesPlots.f(mList = gridY3List, doProj = F, plotRefCase=T, doLegend = F,
     opt = c("R4MvEst",
             "h70", "h80", "h90",
             "M10", "M08","M06",
             "t10","t01","t00",
             "q0","q1"),
     optWt = rep(1,100))




