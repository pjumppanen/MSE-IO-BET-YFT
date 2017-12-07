
# function to take a bivariate sample of OM models from a grid...
# returns a list object - [[1]] List of model names and [[2]] probability weightings for sampling in OM formation

gridSamplerBivar <- function(
  gridListPlus = cDat,              # model definitions and factors of interest for sampling
  sampleQuant = c('MSY','BYoBMSY'), # quantities upon which to base the sample
  corr =0,                          # correlation among characteristics
  FSList = list(                    # subsample models to achieve the following proportions (specifying too many will fail somehow)
    FSOpt1 <- c('t0001', 0.5),
    FSOpt2 <- c('t10', 0.5)),
  mu    = c(422,0.89),              # distribution mean to attain (roughly)
  sigma = c(0.1,0.1),               # lognormal sampling sigma (~CV)
  sigmaTrunc = 3,                   # truncate sample distribution at this many sigma
  logNorm=  T,                      # lognormal is the only option still supported
  nBins  =  3)                      # bivariate sampling on a grid of dim nBins X nBins (going too fine is pointless because model grid is coarse)
{

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





