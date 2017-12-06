#=============================================================================================================================
# Basic R script for making and inspecting the conditioned YFT Operating Model OMyft17.2 and associated robustness scenarios
#=============================================================================================================================
#online
rootDir <- "H:\\C-offline\\MSE-IO-BET-YFT\\"  #modify for local path
#offline
#rootDir <- "C:\\MSE-IO-BET-YFT\\"  #modify for local path

library(mvtnorm)
library(MASS)
library(TinnRcom)                            #if using TinnR IDE
library (r4ss)                                #R package supporting Stock Synthesis - if something's broken it might be becuase this changed
library(PerformanceAnalytics)                 # for chart.Correlation
source(paste(rootDir,"OMconditioning\\RStuff\\phase2\\pasteOperator.R",sep=""))

#source(paste(rootDir,"Source\\pasteOperator.R",sep=""))
source(rootDir %&% "OMconditioning\\RStuff\\seasAsYrToDecYr.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.2Tag.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\makeGridY17.2noTag.f.R")
source(rootDir %&% "OMconditioning\\RStuff\\phase2\\gridSamplerBivar.f.R")
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

#SS_plots(ref2017, uncertainty=F)
ref$derived_quants[ref$derived_quants$LABEL == 'TotYield_MSY',]$Value*4/1000
max(ref$equil_yield$Catch)


# OM17.2 reference case
# tag
# noTag
# combine and weight
#



##################################################################################
# gridY17.2noTag - suite of models without tags for OM-refY17.2
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.2noTag'
gridY17.2noTagList <- makeGridY17.2noTag.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.2noTagList)

#gridY17.2noTagList <- makeGridY17.2noTag.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)


gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2noTag//'
importGrid.f(gridList=gridY17.2noTagList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.2noTagList, mfrowLayout = c(4,2), MSYyLim=c(0,700))




##################################################################################
# gridY17.2Tag - suite of models with tags for OM-refY17.2
path=rootDir %&% 'OMconditioning\\YFT\\gridY17.2Tag'
gridY17.2TagList <- makeGridY17.2Tag.f(path=path, doHess=F, makeGrid=F)    #check the list of grid elements, but don't create the grid DIR structure
print(gridY17.2TagList)

#gridY17.2TagList <- makeGridY17.2Tag.f(path=path, doHess=F, makeGrid=T)    #create the grid structure for running a batch file
# exit this R script and run the SS batch file(s)

gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2Tag//'
importGrid.f(gridList=gridY17.2TagList, gridDir=gridDir, covar=F)

plotIndices2.f(modList=gridY17.2TagList, mfrowLayout = c(4,2), MSYyLim=c(0,700))



##################################################################################
# Merge and plot the two grids to make one unbalanced aggregate...

#gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2noTag//'
#importGrid.f(gridList=gridY17.2noTagList, gridDir=gridDir, covar=F)
#gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2Tag//'
#importGrid.f(gridList=gridY17.2TagList, gridDir=gridDir, covar=F)
OMref17.2grid <- c(gridY17.2TagList, gridY17.2noTagList)
gridList <- OMref17.2grid
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
# Bivariate sample of OMref17.2grid on basis of MSY and B/BMSY consistency with assessment and assumed variance

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
  mu    = c(422,0.89),  #distribution means to aim for
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


#need to ensure that the "correct", i.e. assessment-preferred CPUE is the one used for the MP
#this is achieved by making sure this series is the first model in the list... (iH and q0)
i <- 1
while(!prod((c("iH","q0") %in% unlist(strsplit(names(OMref17.2List[i]), split="_")))) {
  i <- i+1
  if(i>length(OMref17.2List)){print("oops - better check your list for CPUE series..."); break}
}
if(i>1){
  tmp <- OMref17.2List[1]
  OMref17.2List[1] <- OMref17.2List[i]
  OMref17.2List[i] <- tmp
  names(OMref17.2List) <- names(OMref17.2List[c(i,1:(i-1),(i+1):length(OMref17.2List))])
}

#save just the list of  weights (and model names as names)
save(OMref17.2List,file="OMref17.2List.RDA")

# save the full SS3 files, so they don't have to be reloaded from ASCII to recreate OMs ...too big
#gridDir <- rootDir %&% 'OMconditioning//YFT//gridY17.2//'
#importGrid.f(gridList=names(OMref17.2List), gridDir=gridDir, covar=F, keepEverything=T)
#save(list=names(OMref17.2List),file="OMref17.2ListMods.RDA")





