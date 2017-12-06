# Example script demonstrating MSE application for IO-BET with Stock Synthesis conditioned OMs


# issues to resolve: check they work on both BET and YFT
# 1) dim length 1 arrays -> karrays
# 2) time mapping for SS year season model vs seasonAsYear model
# 3) data, MP recommendation, and out of date OM time lags
# 4) add vector of observed catches between OM conditioning and first TAC application

#setwd("C:\\MSE-IO-BET-YFT\\")  # Set the working directory
setwd("H:\\C-offline\\MSE-IO-BET-YFT\\gitMirror")  # Set the working directory

rm(list=ls(all=TRUE))                  # Remove all existing objects from environment

#for Tinn-R users only: (may no longer be required)
.trPaths <- file.path(Sys.getenv("TEMP"), "Tinn-R", c("", "search.txt", "objects.txt", "file.r", "selection.r", "block.r", "lines.r", "reformat-input.r", "reformat-output.r"), fsep="\\")

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
# Create OMs



#===== Create & save a test case for the reference case Operating Models OMrefY17.2 ====================================
source('Rscripts\\phase2\\Build OM Model-OMrefY17.2test.R')  # 20 realization test with 3 SS3 configurations

# Create an OM object for the OMd
print(system.time(OMrefY17.2test<-new("OMss",OMd, Report=F)))
OMrefY17.2test@UseCluster <- 0
# Save the OM
save(OMrefY17.2test,file=paste(getwd(),"/Objects/phase2/OMrefY17.2test.RDA",sep=""))



#===== Create & save the reference case Operating Models OMrefY17.2 ====================================
source('Rscripts\\phase2\\Build OM Model-OMrefY17.2.R')  # full 2000 realization Reference case YFT OM proposed Dec 2017.

# Create an OM object for the OMd
print(system.time(OMrefY17.2<-new("OMss",OMd, Report=F)))
OMrefY17.2@UseCluster <- 0
# Save the OM
save(OMrefY17.2,file=paste(getwd(),"/Objects/phase2/OMrefY17.2.RDA",sep=""))










#===== load and possibly modify existing Operating Models and conduct MP tuning ====================================

# Load the previously created test OM
load(file=paste(getwd(),"/Objects/phase2/OMrefY17.2test.RDA",sep=""))

# Run some MPs without tuning
MPList2 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMrefY17.2test.noTune  <- new("MSE",OMrefY17.2test,MPs <- MPList2,interval=3, Report=F,UseCluster=0)))





MPList1 <- c("PT41F.tune.15","PT41.tune.15","IT5.tune.15","CCt")

# tune MPs to YFT tuning objective T1
OMrefY17.2test@UseCluster <- 0
OMrefY17.2test@tunePM     <- "SBoSBMSY0.5"
OMrefY17.2test@tunePMProjPeriod <-1001
OMrefY17.2test@tunePMTarget <- 1.
OMrefY17.2test@CppMethod <- 0
print(system.time(mseOMrefY17.2test.MPT1  <- new("MSE",OMrefY17.2test,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
save(mseOMrefY17.2test.MPT1,file=paste(getwd(),"/Objects/phase2/mseOMrefY17.2test.MPT1.RDA",sep=""))
load(file=paste(getwd(),"/Objects/phase2/mseOMrefY17.2test.MPT1.RDA",sep=""))
plotTS.f(mseOMrefY17.2test.MPT1, plotByRF=F, doWorms=T)  #Time series (worm) plots
plotTS.f(mseOMrefY17.2test.MPT1, plotByRF=F, doWorms=T, param="t10")  #Time series (worm) plots
plotTS.f(mseOMrefY17.2test.MPT1, plotByRF=F, doWorms=T, param="t00")  #Time series (worm) plots





xxxx below here old stuff...






OMyftY17.1.2 <- OMyftY17.1.1
# tune 216 models to YFT tuning objective T2
OMyftY17.1.2@UseCluster <- 0
OMyftY17.1.2@tunePM     <- "SBoSBMSY0.5"
OMyftY17.1.2@tunePMProjPeriod <-1002
OMyftY17.1.2@tunePMTarget <- 1.
OMyftY17.1.2@CppMethod <- 0
print(system.time(mseOMyftY17.1.2  <- new("MSE",OMyftY17.1.2,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1.2,file=paste(getwd(),"/Objects/mseOMyftY17.1.2.RDA",sep=""))
load(file=paste(getwd(),"/Objects/mseOMyftY17.1.2.RDA",sep=""))
plotTS.f(mseOMyftY17.1.2, plotByRF=F, doWorms=T)  #Time series (worm) plots



# tune 216 models to BET tuning objective 2 = YFT objective T3
OMyftY17.1.3 <- OMyftY17.1.1
OMyftY17.1.3@UseCluster <- 0
OMyftY17.1.3@tunePM     <- "GKmean"
OMyftY17.1.3@tunePMProjPeriod <-1001
OMyftY17.1.3@tunePMTarget <- 0.75
OMyftY17.1.3@CppMethod <- 0
print(system.time(mseOMyftY17.1.3  <- new("MSE",OMyftY17.1.3,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1.3,file=paste(getwd(),"/Objects/mseOMyftY17.1.3.RDA",sep=""))
load(file=paste(getwd(),"/Objects/mseOMyftY17.1.3.RDA",sep=""))
plotTS.f(mseOMyftY17.1.3, plotByRF=F, doWorms=T)  #Time series (worm) plots



# tune 216 models to YFT objective T4 (DK made up for illustration purposes)
OMyftY17.1.4 <- OMyftY17.1.1
OMyftY17.1.4@UseCluster <- 0
OMyftY17.1.4@tunePM     <- "SBoSBMSY0.5"
OMyftY17.1.4@tunePMProjPeriod <-1001
OMyftY17.1.4@tunePMTarget <- 1.5
OMyftY17.1.4@CppMethod <- 0
print(system.time(mseOMyftY17.1.4  <- new("MSE",OMyftY17.1.4,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1.4,file=paste(getwd(),"/Objects/mseOMyftY17.1.4.RDA",sep=""))
load(file=paste(getwd(),"/Objects/mseOMyftY17.1.4.RDA",sep=""))
plotTS.f(mseOMyftY17.1.4, plotByRF=F, doWorms=T)  #Time series (worm) plots
plotTS.f(mseOMyftY17.1.4, plotByRF=F, doWorms=T, param="t10")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1.4, plotByRF=F, doWorms=T, param="t01")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1.4, plotByRF=F, doWorms=T, param="t00")  #Time series (worm) plots









#===== load existing Operating Models and conduct MSE with tuned MPs ====================================

#compare MP performance tuned to 216 models against full 2160 realizations...instead of 0.75, got (c(.787,.781,.764)
load(file=paste(getwd(),"/Objects/OMyftY17.1.RDA",sep=""))
MPList2 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1  <- new("MSE",OMyftY17.1,MPs <- MPList2,interval=3, Report=F,UseCluster=0)))
tableMSE.f(mseOMyftY17.1)


#only rerun MSE because of mse Object change for OM subset plotting
load(file=paste(getwd(),"/Objects/OMyftY17.1.1.RDA",sep=""))
MPList2 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1.1t3  <- new("MSE",OMyftY17.1.1,MPs <- MPList2,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1.1t3,file=paste(getwd(),"/Objects/mseOMyftY17.1.1t3.RDA",sep=""))
tableMSE.f(mseOMyftY17.1.1t3)
plotTS.f(mseOMyftY17.1.1t3, plotByRF=F, doWorms=F)  #Time series (worm) plots



#rerun MSE because of mse Object change for OM subset plotting
load(file=paste(getwd(),"/Objects/OMyftY17.1.1.RDA",sep=""))
MPList2 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1.1t1  <- new("MSE",OMyftY17.1.1,MPs <- MPList2,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1.1t1,file=paste(getwd(),"/Objects/mseOMyftY17.1.1t1.RDA",sep=""))
tableMSE.f(mseOMyftY17.1.1t1)
plotTS.f(mseOMyftY17.1.1t1, plotByRF=F, doWorms=T)  #Time series (worm) plots




#test MPs tuned tuned to 216 models against robustness scenario OMyftY17.1selTrend - OM-rob-selTrend
load(file=paste(getwd(),"/Objects/OMyftY17.1selTrend.RDA",sep=""))
MPList.t3 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1selTrend  <- new("MSE",OMyftY17.1selTrend,MPs <- MPList.t3,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1selTrend,file=paste(getwd(),"/Objects/mseOMyftY17.1selTrend.RDA",sep=""))
tableMSE.f(mseOMyftY17.1selTrend)
plotTS.f(mseOMyftY17.1selTrend, plotByRF=F, doWorms=F, param="SS")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1selTrend, plotByRF=F, doWorms=F, param="NS")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1selTrend, plotByRF=F, doWorms=F, param="ST")  #Time series (worm) plots



#test MPs tuned tuned to 216 models against robustness scenario OMyftY17.1tagWt - OM-rob-tagWt
load(file=paste(getwd(),"/Objects/OMyftY17.1tagWt.RDA",sep=""))
MPList.t3 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1tagWt  <- new("MSE",OMyftY17.1tagWt,MPs <- MPList.t3,interval=3, Report=F,UseCluster=0)))
save(mseOMyftY17.1tagWt,file=paste(getwd(),"/Objects/mseOMyftY17.1tagWt.RDA",sep=""))
tableMSE.f(mseOMyftY17.1tagWt)
plotTS.f(mseOMyftY17.1tagWt, plotByRF=F, doWorms=F, param="t15")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1tagWt, plotByRF=F, doWorms=F, param="t10")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1tagWt, plotByRF=F, doWorms=F, param="t01")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1tagWt, plotByRF=F, doWorms=F, param="t001")  #Time series (worm) plots
plotTS.f(mseOMyftY17.1tagWt, plotByRF=F, doWorms=F, param="t00")  #Time series (worm) plots












# Load a previously created OM
load(file=paste(getwd(),"/Objects/OMyftY17.1.test40.RDA",sep=""))
MPList <- c("IT5.15.216.t3","CCt.216.t3")
print(system.time(msePlotTest  <- new("MSE",OMyftY17.1test40,MPs <- MPList,interval=3, Report=F,UseCluster=0)))
tableMSE.f(msePlotTest)
plotTS.f(msePlotTest, plotByRF=F, doWorms=T, param="t00")  #Time series (worm) plots


#test tuned MPs against robustness scenarios
load(file=paste(getwd(),"/Objects/OMyftY17.1.RDA",sep=""))
MPList2 <- c("PT41.15.216.t3","IT5.15.216.t3","CCt.216.t3")
print(system.time(mseOMyftY17.1  <- new("MSE",OMyftY17.1,MPs <- MPList2,interval=3, Report=F,UseCluster=0)))


# Load a previously created OM
load(file=paste(getwd(),"/Objects/OMyftY17.1.RDA",sep=""))

MPList1 <- c("PT41.tune.15","IT5.tune.15","CCt")


OMyftY17.1@UseCluster <- 0
OMyftY17.1@tunePM     <- "SBoSBMSY0.5"
OMyftY17.1@tunePMProjPeriod <-1001
OMyftY17.1@tunePMTarget <- 1.
OMyftY17.1@CppMethod <- 0
print(system.time(mseOMyftY17.1  <- new("MSE",OMyftY17.1,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1, plotByRF=F, doWorms=T)  #Time series (worm) plots


# Run an MSE; "CC001" has fisheries shut off to test unfished equilibrium recovery dynamics
#print(system.time(mseOMyftY17.1 <- new("MSE",OMyftY17.1,MPs<-c("CC001", "CC350", "PT41.100.9", "IT1.00"),interval=3, Report=F,UseCluster=0)))

################################################################################
# test with streamlined OM
# Create a test OM object (20 realizations)
# tunePMProjPeriod =1001 = 2019-2039
# tunePMProjPeriod =1002 = 2024

source('Rscripts\\Build OM Model-yftY17.1test4.R') # 4 OM specifications

# Tune to modified YFT tuning objective 1 (should be 1.0 rather than 1.5)

MPList1 <- c("PT41.tune.15","IT5.tune.15","CCt")

print(system.time(OMyftY17.1test4<-new("OMss",OMd, Report=F)))
OMyftY17.1test4@UseCluster <- 0
OMyftY17.1test4@tunePM     <- "SBoSBMSY0.5"
OMyftY17.1test4@tunePMProjPeriod <-1001
OMyftY17.1test4@tunePMTarget <- 1.

# R-based
OMyftY17.1test4@CppMethod <- 0
print(system.time(mseOMyftY17.1test4.r.t1  <- new("MSE",OMyftY17.1test4,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.r.t1, plotByRF=F, doWorms=T)  #Time series (worm) plots

# Cpp-based
print(system.time(mseOMyftY17.1test4.c.t1  <- new("MSE",OMyftY17.1test4,MPs <- MPList1,interval=3, CppMethod <- 1, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.c.t1, plotByRF=F, doWorms=T)  #Time series (worm) plots


# Tune to modified YFT tuning objective 2

OMyftY17.1test4@UseCluster <- 0
OMyftY17.1test4@tunePM     <- "SBoSBMSY0.5"
OMyftY17.1test4@tunePMProjPeriod <-1002
OMyftY17.1test4@tunePMTarget <- 1.
print(system.time(mseOMyftY17.1test4.r.t2  <- new("MSE",OMyftY17.1test4,MPs <- MPList1[3], interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.r.t2, plotByRF=F, doWorms=T)  #Time series (worm) plots
mseOMyftY17.1test4.r.t2@tune*200

# Tune to BET tuning objective (just to test)

OMyftY17.1test4@tunePM           <- "GKmean"
OMyftY17.1test4@tunePMProjPeriod <- 1001
OMyftY17.1test4@tunePMTarget     <- 0.75
print(system.time(mseOMyftY17.1test4.r.t3 <- new("MSE",OMyftY17.1test4,MPs <- "CCt",interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.r.t3, plotByRF=F, doWorms=T)  #Time series (worm) plots
tableMSE.f(mseOMyftY17.1test4.r.t3)


print(system.time(mseOMyftY17.1test4.r <- new("MSE",OMyftY17.1test4,MPs <- "CCt",interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.r, plotByRF=F, doWorms=T)  #Time series (worm) plots


...reset tuning bounds somewhere...



# Save the OM
save(OMyftY17.1test4,file=paste(getwd(),"/Objects/OMyftY17.1.test4.RDA",sep=""))

MPList0 <- c("PT30.tune.15")
MPList1 <- c("PT41.tune.15")
print(system.time(mseOMyftY17.1test4.c  <- new("MSE",OMyftY17.1test4,MPs <- MPList0,CppMethod=1,interval=3, Report=F,UseCluster=0)))


print(system.time(mseOMyftY17.1test4.r  <- new("MSE",OMyftY17.1test4,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.r, plotByRF=F, doWorms=T)  #Time series (worm) plots

print(system.time(mseOMyftY17.1test4.rI <- new("MSE",OMyftY17.1test4,MPs <- "IT5.tune.15",interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.rI, plotByRF=F, doWorms=T)  #Time series (worm) plots

print(system.time(mseOMyftY17.1test4.rI <- new("MSE",OMyftY17.1test4,MPs <- "CCt",interval=3, Report=F,UseCluster=0)))
plotTS.f(mseOMyftY17.1test4.rI, plotByRF=F, doWorms=T)  #Time series (worm) plots




MPList2 <- c("PT41.tune.15","IT5.tune.15")
print(system.time(mseOMyftY17.1test4 <- new("MSE",OMyftY17.1test4,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))
MPList6 <- c("PT41.tune.9","PT41.tune.15","PT41.tune.05","IT5.tune.15","IT10.tune.15","IT5.tune.05")
print(system.time(mseOMyftY17.1test4 <- new("MSE",OMyftY17.1test4,MPs <- MPList6,interval=3, Report=F,UseCluster=0)))


################################################################################




# Run some tuning MPs
MPList0 <- c("PT41.tune.9")
MPList6 <- c("PT41.tune.9","PT41.tune.15","PT41.tune.05","IT5.tune.15","IT10.tune.15","IT5.tune.05")
print(system.time(mseOMyftY17.1 <- new("MSE",OMyftY17.1,MPs <- MPList1,interval=3, Report=F,UseCluster=0)))





print(system.time(mseOMyftY17.1 <- new("MSE",OMyftY17.1,MPs<-c("CC350","PT41.tune.9"),interval=3, Report=F,UseCluster=0)))



# Plot some key time series, including recruitment and total biomass, which are not part of the MWG standard
plotTS.f(mseOMyftY17.1, plotByRF=F, doWorms=T)  #Time series (worm) plots
plotTS.f(mseOMyftY17.1, plotByRF=T, doWorms=T)  #Time series (worm) plots

# Re-Run the MSE with the C++ projection code this time
print(system.time(mseOMyftY17.1.c <- new("MSE",OMyftY17.1,MPs<-"CC001",interval=3, Report=F,CppMethod=1, UseCluster=0)))
plotTS.f(mseOMyftY17.1.c, plotByRF=F)  #Time series (worm) plots


# Load a previously created OM
load(file=paste(getwd(),"/Objects/mseOMyftY17.1.3.RDA",sep=""))
plotTS.f(mseOMyftY17.1.3, plotByRF=F)  #Time series (worm) plots










stop here....



#load(file="Objects\\OMd.YFTY2r0cpue0n2.RDA")           #loads OMd object named OMd created previously

#source('Rscripts\\Build OM Model-yftY2r6cpue3n20.R')    # 2  SS (7,20) specifications, spanning MSY range; realistic noise; 20 reps total
#load(file="Objects\\OMd.yftY2r6cpue3n20.RDA")          #loads OMd object named OMd created previously

#source('Rscripts\\Build OM Model-yftY2r6cpue3n200.R')    # 2  SS (7,20) specifications, spanning MSY range; realistic noise; 20 reps total
#load(file="Objects\\OMd.yftY2r6cpue3n200.RDA")          #loads OMd object named OMd created previously

#source('Rscripts\\Build OM Model-yftY27r0cpue0n27.R')  # 27 SS specification grid; negligible error; 1 rep per SS file
#load(file="Objects\\OMd.YFT27r0cpue0n27.RDA")         #loads OMd object named OMd created previously

#source('Rscripts\\Build OM Model-yftY27r6cpue3n81.R')  # 27 SS specification grid; realistic noise; 81 reps total
#load(file="Objects\\OMd.yft27r6cpue3n81.RDA")         #loads OMd object named OMd created previously

#str(OMd)  inspect OMd contents



#===== Create, save and/or load an Operating model derived from an OM definition ===========================
OM<-new('OMss',OMd)


#save(OM,file="Objects\\OM.YFTY2r0cpue0n2.RDA")
#save(OM,file="Objects\\OM.YFTY2r6cpue3n20.RDA")
save(OM,file="Objects\\OM.YFTY17.1.RDA")
#save(OM,file="Objects\\OM.YFTY27r6cpue3n81.RDA")

#load(file="Objects\\OM.YFTY2r0cpue0n2.RDA")
#load(file="Objects\\OM.YFTY2r6cpue3n20.RDA")
#load(file="Objects\\OM.YFTY27r0cpue0n27.RDA")
#load(file="Objects\\OM.YFTY27r6cpue3n81.RDA")


gc() #memory check and clean up



#===== Undertake MSE simulations ================================================================

#test new MP return list
system.time(tmse0<-new('MSE',OM,MPs<-c("CC001n", "CE1.0n", "CC100CE1.0n"),interval=3,IE="Umax", loud=F))
system.time(tmse0<-new('MSE',OM,MPs<-c("CE1.0n"),interval=3,IE="Umax", loud=F))


#old MPs below here...


#Check dynamics when fishery shut off (TAC = 1 to collect data)
#system.time(tmse0<-new('MSE',OM,MPs<-c("CC001n"),interval=3,IE="Umax", loud=F))

#(TAC = 249K)
#system.time(tmse0<-new('MSE',OM,MPs<-c("CC249"),interval=3,IE="Umax", loud=F))

#Check MSY from projections consistent with MSY estimate (edit MP CCeq for appropriate TAC; MP allows rebuilding before TAC is applied)
#system.time(tmse0<-new('MSE',OM,MPs<-c("CCeq"),interval=3,IE="Umax", loud=F))

# Compare 3 simple feedback-based MPs
#system.time(tmse0<-new('MSE',OM,MPs<-c("ITarg1.5","ITarg2.5","ITarg3.5"),interval=3,IE="Umax", loud=F))
#system.time(tmse0<-new('MSE',OM,MPs<-c("ITarg1.5","ITarg3.5"),interval=3,IE="Umax", loud=F))

gc() #memory check and clean up



#===== Summarize results ===============================================================================

summary(tmse0)                          # Tabulate results
# Original ABT plotting routines are seemingly not really compatible with Windows
#source("source/plotting.R")
#plot(tmse0)                             # Plot time series and Kobe's

#FLR-based worm plots
plotStockWrapper.f(tmse0)


