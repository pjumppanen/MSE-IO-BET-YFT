# Example script demonstrating MSE application for IO- bigeye and yellowfin with Stock Synthesis conditioned OMs

# Note that creating an OM with the new("OMss",...) requires the SS output files
# It is much faster to load() an existing OM than to create it.

#Full list of Performance Measures: (1:15 from user manual)
#c("SBoSB0", "minSBoSB0", "SBoSBMSY", "FoFMSY", "FoFtarg", "GK", "RK","PrSBgt0.2SB0", "PrSBgtSBlim", "Y", "relCPUE", "YoMSY", "APCY", "YcvPct", "PrYlt0.1MSY")

#setwd("C:\\Users\\kol018\\MSE-IO-BET-YFT\\Mseom\\Mseom\\MSE-IO-BET-YFT\\YFT-MSE")  # Set the working directory

rm(list=ls(all=TRUE))                  # Remove all existing objects from environment

#for Tinn-R users only: (may no longer be required)
#.trPaths <- file.path(Sys.getenv("TEMP"), "Tinn-R", c("", "search.txt", "objects.txt", "file.r", "selection.r", "block.r", "lines.r", "reformat-input.r", "reformat-output.r"), fsep="\\")

#########################################################################################################
# ======Prerequisites ===================================================================================

source("Source/seasAsYrToDecYr.f.R")  # functions for converting among date array indexing formats
source("Source/pasteOperator.R")      # replace paste() with %&%
source("Source/MSE_source-OMss.r")    # Reference point calculations
source("Source/Objects-OMss.R")       # Defines and instantiates Operating Models, runs MSE projections
source("Source/MPs.R")                # Defines Management Procedures
source("Source/Diagnostics.r")
source("Source/mseGraphics.R")


#########################################################################################################



#========================================================================================================
# YFT - a more substantive demonstration MSE; 54 equally weighted SS specifications; 108 stochastic replicates
source("Rscripts/Build OMyft2r108.2.R")
print(system.time(OMyft2r108<-new("OMss",OMd, Report=F)))
save(OMyft2r108,file=paste(getwd(),"/Objects/OMyft2r108.RDA",sep=""))
load(file=paste(getwd(),"/Objects/OMyft2r108.RDA",sep=""))

# MSE on a list of candidate MPs
#MPL1 <- c("CC200", "CC400","IT1.50","IT3.50","PT41.tune.9","PT41.100.2","PT41.100.9")
MPL1 <- c("IT3.50","PT41.tune.9")
print(system.time(mseOMyft2r108.MPL1<-new("MSE",OMyft2r108,MPs=MPL1, interval=3, Report=F,UseCluster=0)))
save(mseOMyft2r108.MPL1,file=paste(getwd(),"/Objects/mseOMyft2r108.MPL1.RDA",sep=""))
load(file=paste(getwd(),"/Objects/mseOMyft2r108.MPL1.RDA",sep=""))

# Tables of Performance Measures
tableMSE.f(mseOMyft2r108.MPL1)      # aggregated over fisheries and regions
tableMSEbyRF.f(mseOMyft2r108.MPL1)  # disaggregated by fisheries and regions (not all PMs)

# Compare some specific stats for MPs "CC200" and "PT41.100.2"
# "y10" indicates first 10 year projection period
# 0,5, 0.1, 0.9 =percentiles
c(tableMSE.f(mseOMyft2r108.MPL1)["CC200y10",c("SBoSBMSY0.5","SBoSBMSY0.1","SBoSBMSY0.9")])
c(tableMSE.f(mseOMyft2r108.MPL1)["PT41.100.2y10",c("SBoSBMSY0.5","SBoSBMSY0.1","SBoSBMSY0.9")])

# Bi-variate trade-off plots
plotTO.f(mseObj=list(mseOMyft2r108.MPL1), MPs=MPL1, ylims <- c(2.5,1,1,1))

# Confidence interval plots (Udon-Soba plots)
dev.new(width=6,height=8)
plotUS.f(mseObj=mseOMyft2r108.MPL1,MPs=MPL1, plotStats=c("SBoSB0","RK","Y"))

# Time series plots
dev.new(width=6,height=4)
plotTS.f(mseOMyft2r108.MPL1, plotByRF=F, mwgPlots=T)





gc() # R garbage collection function can free up useful memory
closeAllConnections()  # useful if something crashes while files are open
