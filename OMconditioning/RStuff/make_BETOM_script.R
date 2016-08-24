# Basic R script for making and inspecting the conditioned demonstration case BET Operating Model betOM-1
# based on files from Langley et al 2013 assessment provided by Rishi
# uses R64 3.2.2
#=============================================================================================================================

rootDir <- "H:\\C-offline\\MSE-IO-BET-YFT\\"  #modify for local path

library(TinnRcom)                            #if using TinnR IDE
library (r4ss)                                #R package supporting Stock Synthesis

source(paste(rootDir,"YFT-MSE\\Source\\pasteOperator.R",sep=""))
source(rootDir %&% "OMconditioning\\RStuff\\seasAsYrToDecYr.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\makeGridB1.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\importGrid.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\plotIndices.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\timeSeriesPlots.f.R")


path=rootDir %&% 'OMconditioning\\BET\\gridB1'
gridB1List <- makeGridB1.f(path=path, doHess=F, makeGrid=F)    #make a list of grid elements, but don't create the grid DIR structure
gridB1List
#gridB1List <- makeGridB1.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)


# import some summary results & diagnostics from the gridB1List of models
gridDir <- rootDir %&% 'OMconditioning//BET//gridB1//'
importGrid.f(gridList=gridB1List, gridDir=gridDir, covar=F)


#=============================================================================================================================
# original assessment results
#betOld <- SS_output(dir='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\BET\\betOld', covar=F, ncol=213)
#SS_plots(replist=betOld, uncertainty=F)

# using SS3.23Y files, r4ss Dec 2015
#betNew <- SS_output(dir='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\BET\\betR1h80m10q0', covar=F, ncol=213)
#SS_plots(replist=betNew, uncertainty=F)


#use boxplots to plot some key model results and diagnostics, marginalized by grid Options
plotIndices.f(modList=gridB1List, mfrowLayout = c(2,2),
    SPB_Yr = "SPB_348", refSSB0 = -c(3421, 2492, 4628), refSSBY = -c(786,
    457, 1009), refSSBMSY = -c(1200, 751, 1483), refSSBYoSSBMSY = -c(0.66,
    0.48, 0.72), refSSBYoSSB0 = -c(0.23, 0.17, 0.32), refMSY = -c(402,
    309, 530))



#plot some time series (ignore ref assessment case lines because the grids are essentially the same )
timeSeriesPlots.f(mList = gridB1List, doProj = F, doLegend = F,
     opt = c("R1P1",
             "h70", "h80", "h90",
             "M10", "M08","M06",
             "q0","q1"),
     yrSeasY=200,
     optWt = rep(1,100),
     seasAsYr=99, endSeasAsYr=384, numSeas=4, endYr=2022, endSeas=4,
     xlims=c(1950,2012))


#check every model result
#for(i in 1:length(gridB1List)){
#  rm(tmp)
#  converged <- file.exists(gridDir %&% gridB1List[i] %&% "\\" %&% "ss3_opt.std")
#  print(c(i,gridB1List[i],converged))
#  tmp <- SS_output(dir=gridDir %&% gridB1List[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=converged)
#  par(ask=F)
#  par(mfroc=c(3,2))
#  SS_plots(tmp, uncertainty=converged, plot=c(3,4,6,8,11,12),png=F)
#}


gridB1FullList <- makeGridB1.f (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridD', makeGrid=F, splitMasterBatch=T, doHess=T,
  sp.val  = c( "R1P1"), #1 regions, 1 population, no mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M10", "M08","M06"),#mort
  q.val   = c( "q0", "q1")) # longline catchability



# these functions are incomplete for BET...

k2sm.f(mList = gridB1FullList, doProj = F, doLegend = F,
     opt = c("R1P1", #1 regions, 1 population, no mixing
             "h70","h80","h90", #SR steepness
             "M10", "M08","M06",#mort
             "q0", "q1"), # longline catchability
     optWt = rep(1,9))

plotIndices.f(modList=gridB1FullList, mFile = F, mfrowLayout = c(2,2),
    SPB_Yr = "SPB_348", refSSB0 = c(0, 0, 0), refSSBY = c(0,0,0), refSSBMSY = c(0, 0, 0), refSSBYoSSBMSY = c(1.92,0.87,02.22), refSSBYoSSB0 = c(0., 0., 0.), refMSY = c(188,98,207))











stop - below here copied from YFT




gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridA//'

#look at every gridA model result
for(i in 1:length(gridAList)){
  rm(tmp)
  converged <- file.exists(gridDir %&% gridAList[i] %&% "\\" %&% "ss3_opt.std")
  print(c(i,gridAList[i],converged))
  tmp <- SS_output.f(dir=gridDir %&% gridAList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=converged)
  par(ask=T)
  SS_plots(tmp, uncertainty=converged, plot=c(3,4,6,8,11,12))
}




#extract the relevant results of gridA
gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridA//'
for(i in 1:length(gridAList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% gridAList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F)
  converged <- file.exists(gridDir %&% gridAList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(gridAList[i],tmp)

}

for(i in 1:length(gridAList)){
  print(c(i,gridAList[i],converged))
}
                                                                               


makdGridDlist(splitMasterBatch=T)  #set up the grid
gridDListFull <- makeGridD.f(makeGrid=F)          #just make a list of the grid elements
gridDList <- gridDListFull[73:144]

#extract the relevant results of gridD

gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridD//'
for(i in 1:length(gridDList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% gridDList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
  converged <- file.exists(gridDir %&% gridDList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(gridDList[i],tmp)

}

for(i in 1:length(gridDList)){
  print(c(i,gridDList[i],converged))
}

tmp <- SS_output.f(dir=gridDir %&% gridDList[1], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
SS_plots(tmp, uncertainty=F, plot=c(3,4,6,8,11,12))

                                                                               
plotIndices.f(modList=gridDListFull)



  
plotIndices.f(modList=gridDfList)


# import missing elements of gridD
tmpList1 <- makeGridD.f (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridD', makeGrid=F, splitMasterBatch=T, doHess=T,
  sp.val  = c( "R2P1Hi","R2P1Lo"), #2 regions, 1 population, Lo or Hi mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M12","M10", "M08"),      #mort              
  Ucv.val =c("Ucv15","Ucv30"),  #cpue weight
  CL.val = c("CL5","CL1"), # CL weighting"CL20",
  Rcv.val=c("Rcv3")) # rec CV

tmpList2 <- makeGridD.f (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridD', makeGrid=F, splitMasterBatch=T, doHess=T,
  sp.val  = c( "R2P1Hi","R2P1Lo"), #2 regions, 1 population, Lo or Hi mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M06"),      #mort              
  Ucv.val =c("Ucv15","Ucv30"),  #cpue weight
  CL.val = c("CL5","CL1"), # CL weighting"CL20",
  Rcv.val=c("Rcv6","Rcv3")) # rec CV

tmpList <- c(tmpList1,tmpList2)

gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridD//'
for(i in 1:length(tmpList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% tmpList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
  converged <- file.exists(gridDir %&% tmpList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(tmpList[i],tmp)
}



gridDFullList <- makeGridD.f (path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridD', makeGrid=F, splitMasterBatch=T, doHess=T,
  sp.val  = c( "R2P1Hi","R2P1Lo"), #2 regions, 1 population, Lo or Hi mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M12","M10", "M08","M06"),      #mort              
  Ucv.val = c("Ucv15","Ucv30"),  #cpue weight
  CL.val  = c("CL5","CL1"), # CL weighting"CL20",
  Rcv.val = c("Rcv6","Rcv3")) # rec CV


  
k2sm.f(mList = gridDFullList, doProj = F, doLegend = F, 
     opt = c("R2P1Hi","R2P1Lo",
             "h70", "h80", "h90", 
             "M10", "M08", 
             "Ucv30", 
             "Rcv6", 
             "CL1"),
     optWt = rep(1,10))  

plotIndices.f(modList=gridDFullList)


gridDfList <- makeGridD.f(path='H:\\C-offline\\MSE-IO-BET-YFT\\OMconditioning\\YFT\\gridD', makeGrid=F, splitMasterBatch=F, doHess=F,
  sp.val  = c( "R2P1Hi","R2P1Lo"), #2 regions, 1 population, Lo or Hi mixing
  h.val   = c( "h70","h80","h90"), #SR steepness
  M.val   = c( "M10", "M08","M06"),      #mort              
  Ucv.val =c("Ucv30"),  #cpue weight
  CL.val = c("CL1"), # CL weighting"CL20",
  Rcv.val=c("Rcv6")) # rec CV

  
k2sm.f(mList = gridDfList, doProj = F, doLegend = F, 
     opt = c("R2P1Hi","R2P1Lo",
             "h70", "h80", "h90", 
             "M10", "M08","M06","M12", 
             "Ucv30","Ucv15", 
             "Rcv6","Rcv3", 
             "CL1","CL5"),
     optWt = rep(1,15))  

plotIndices.f(modList=gridDfList)

#grid derived from the assessment: includes tags and environmental variables
gridLList <- makeGridL.f(doHess=F,makeGrid=F)  

tmpList <- gridLList
gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridL//'
for(i in 1:length(tmpList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% tmpList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
  converged <- file.exists(gridDir %&% tmpList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(tmpList[i],tmp)
}

plotIndices.f(modList=gridLList)
k2sm.f(mList = gridLList, doProj = F, doLegend = F, 
     opt = c("R2P1Es",
             "h70", "h80", "h90", 
             "M10", "M08","M06", 
             "t10","t01","t00"), 
     optWt = rep(1,9))  


#grid derived from the assessment: environmental variables are omitted (but tags included)
gridMList <- makeGridM.f(doHess=F)  

tmpList <- gridMList
gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridM//'
for(i in 1:length(tmpList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% tmpList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
  converged <- file.exists(gridDir %&% tmpList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(tmpList[i],tmp)
}

plotIndices.f(modList=gridMList)

k2sm.f(mList = gridMList, doProj = F, doLegend = F, 
     opt = c("R2P1Lo","R2P1Hi",
             "h70", "h80", "h90", 
             "M10", "M08","M06", 
             "t10","t01","t00"), 
     optWt = rep(1,11))  

#pool grid L and M models
gridLMList <- c(gridLList,gridMList) 
plotIndices.f(modList=gridLMList)

k2sm.f(mList = gridLMList, doProj = F, doLegend = F, 
     opt = c("R2P1Lo","R2P1Hi","R2P1Es",
             "h70", "h80", "h90", 
             "M10", "M08","M06", 
             "t10","t01","t00"), 
     optWt = rep(1,100))  

#No environment, estimated movement models
gridNList <- makeGridN.f(doHess=F, makeGrid=F)

tmpList <- gridNList
gridDir <- 'H://C-offline//MSE-IO-BET-YFT//OMconditioning//YFT//gridN//'
for(i in 1:length(tmpList)){
  rm(tmp)
  tmp <- SS_output.f(dir=gridDir %&% tmpList[i], repfile = "Report.sso", compfile = "CompReport.sso", ncol=213, covar=F,forecast=F)
  converged <- file.exists(gridDir %&% tmpList[i] %&% "\\" %&% "ss3_opt.std") 
  M <- tmp$endgrowth$M
  tmp <- list(tmp$derived_quants, tmp$recruit, tmp$maximum_gradient_component,tmp$likelihoods_used,tmp$cpue, tmp$Length_comp_Eff_N_tuning_check, tmp$timeseries[,1:8], tmp$parameters[,1:9],M,
         tmp$equil_yield, converged )
  names(tmp) <- c('derived_quants', 'recruit', 'maximum_gradient_component','likelihoods_used','cpue', 'Length_comp_Eff_N_tuning_check', 'timeseries', 'parameters','M',
          'equil_yield', 'converged')

  assign(tmpList[i],tmp)
}

plotIndices.f(modList=gridNList)

k2sm.f(mList = gridNList, doProj = F, doLegend = F, 
     opt = c("R2P1NE",
             "h70", "h80", "h90", 
             "M10", "M08","M06", 
             "t10","t01","t00"), 
     optWt = rep(1,100))  

