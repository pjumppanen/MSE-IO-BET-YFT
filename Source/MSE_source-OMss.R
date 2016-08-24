#Modifications to original ABT MSE code include:
#  -several initialization functions were not relevant because the conditioned SS inputs provide values
#  -NSN and SSN partitioning was replaced by N (total) and SSN (ABT formulation did not really do SSN-NSN split properly)
#  -spatial targeting dropped; Effort distribution based on 'terminal F' in assessments is used instead
#  -Recruitment occurs in every (nsubyrs) season (as does SSB calculations)
#  -Annual Rec devs are auto-correlated and distributed according to input params (not SSB)
#  -growth, maturity, M, etc also seasonal
#  -Order of events altered to be consistent with SS: Movement, Recruitment, M+F, Graduation
#  -Many karray operations were altered (in various ways) to prevent unwanted dim dropping for npop=1
#        -these should probably be replaced with pkg karray to be more elegant and robust to future dim issues
#
# ===================================================================================================================
# ==== MSE source code ==============================================================================================
# ===================================================================================================================
cat("Installing and loading needed libraries")
cat("\n")
# Not actually used but C++ interface files require Oarray in the case of ADT arrays
# with starting indices other than 1
if(!require("Oarray"))install.packages("Oarray",repos="https://cloud.r-project.org")
if(!require("keep"))install.packages("keep",repos="https://cloud.r-project.org")
if(!require("stringr"))install.packages("stringr",repos="https://cloud.r-project.org")
if(!require("r4ss"))install.packages("r4ss",repos="https://cloud.r-project.org")
if(!require("ggplot2"))install.packages("r4ss",repos="https://cloud.r-project.org")
if(!require("reshape2"))install.packages("r4ss",repos="https://cloud.r-project.org")


library(parallel)
library(abind)
library(stringr)


#load Mseom module and R interface code
# Load the library
if (version$os == "mingw32")
{
  # Running in Windows
  LibName      <- "Mseom"
  LibExtension <- ".dll"

} else
{
  # Running in Linux
  LibName      <- "libMseom"
  LibExtension <- ".so"
}

# Only load library if not already loaded. Loading more than once results in
# R mis-behaving and crashing
if (is.na(match(LibName,  attr(getLoadedDLLs(), "names"))))
{
  LibPath <- paste("./lib/", LibName, LibExtension, sep="")
  dyn.load(LibPath)
}


source("./lib/OmB_R_interface.r")
source("./lib/Om_R_interface.r")


popdyn.domov <- function(Ntemp, movtemp)
{
  # P A R  x  P A R R
  nareas <- dim(movtemp)[4]

  apply(karray(Ntemp, c(dim(Ntemp),nareas)) * movtemp, c(1,2,4), sum)
}


Cpp_MSYrefs <- function(SimNums, .Object, nyears=40)
{
  nsims <- length(SimNums)

  Obj <- Om.create(nsims,
                   .Object@npop,
                   .Object@nages,
                   .Object@nsubyears,
                   .Object@nareas,
                   .Object@nfleets,
                   .Object@Recsubyr)

  Report        <- as.integer(1)
  maxit         <- as.integer(1000)
  ntargets      <- as.integer(length(.Object@targpop))
  MinPar        <- as.double(rep(NA, nsims))
  MSY           <- as.double(rep(NA, nsims))
  BMSY          <- as.double(rep(NA, nsims))
  SSBMSY        <- as.double(rep(NA, nsims))
  SSBMSY_B0     <- as.double(rep(NA, nsims))
  N             <- karray(as.double(NA),c(nsims, .Object@npop, .Object@nages, .Object@nsubyears + 1, .Object@nareas))
  NBefore       <- karray(as.double(NA),c(nsims, .Object@npop, .Object@nages, .Object@nsubyears + 1, .Object@nareas))
  SSN           <- karray(as.double(NA),c(nsims, .Object@npop, .Object@nages, .Object@nsubyears, .Object@nareas))
  C             <- karray(as.double(NA),c(nsims, .Object@npop, .Object@nages, .Object@nsubyears, .Object@nareas, .Object@nfleets))
  SSBA          <- karray(as.double(NA),c(nsims, .Object@npop))

  Om.nt.findMSYrefs(Obj,
                    Report,
                    .Object@ECurrent[SimNums,,,],
                    .Object@q[SimNums,],
                    .Object@R0[SimNums,],
                    .Object@M[SimNums,,,1],
                    .Object@mat[SimNums,,,1],
                    .Object@Idist[SimNums,,,],
                    .Object@Len_age[SimNums,,,1],
                    .Object@Wt_age[SimNums,,,1],
                    .Object@sel[SimNums,,],
                    .Object@mov[SimNums,,,,,],
                    .Object@h[SimNums,],
                    .Object@Recdist[SimNums,,],
                    .Object@SRrel,
                    N,
                    NBefore,
                    SSN,
                    C,
                    SSBA,
                    ntargets,
                    as.integer(.Object@targpop),
                    as.integer(nyears),
                    MinPar,
                    MSY,
                    BMSY,
                    SSBMSY,
                    SSBMSY_B0,
                    maxit)

  B0   <- OmB.get.nt.B0(Obj)
  SSB0 <- OmB.get.nt.SSB0(Obj)

  refs <- karray(NA, dim=c(9, nsims))

  print("Effort multiplier at MSY:")
  print(exp(MinPar))

  for (cn in 1:nsims)
  {
    SimNum    <- SimNums[cn]
    refs[,cn] <- MSYreferencePoints(.Object@ECurrent[SimNum,,,],
                                    .Object@sel[SimNum,,],
                                    .Object@M[SimNum,,,1],
                                    C[cn,,,,,],
                                    SSB0[cn,],
                                    B0[cn,],
                                    N[cn,,,,],
                                    NBefore[cn,,,,],
                                    SSN[cn,,,,],
                                    .Object@Wt_age[SimNum,,,1],
                                    .Object@targpop,
                                    .Object@nsubyears,
                                    .Object@npop,
                                    .Object@nages,
                                    .Object@nareas,
                                    .Object@nfleets)
  }

  # delete the object
  Om.destroy(Obj)

  return (refs)
}


R_MSYrefs <- function(SimNums, UseCluster, .Object, nyears=40)
{
  nCores    <- detectCores()
  nOMfiles  <- length(SimNums)

  if (nOMfiles < nCores)
  {
    nCores <- nOMfiles
  }

  if (UseCluster && (nCores > 1))
  {
    cl <- makeCluster(nCores)
    clusterExport(cl, c("popdyn","popdyn.domov","MSYreferencePoints"))
    MSYrefs <-parSapply(cl, SimNums, FUN=getMSYrefs, .Object, AsCluster=UseCluster, nyears=70)
    stopCluster(cl)
    printClusterOutput(SimNums)
  }
  else
  {
    MSYrefs <-sapply(SimNums, FUN=getMSYrefs, .Object, nyears=70)
  }

  return (MSYrefs)
}


printClusterOutput <- function(SimNums)
{
  for (sim in SimNums)
  {
    StdOutFileName <- paste("StdOutFile", sim, ".txt", sep="")

    if (file.exists(StdOutFileName))
    {
      Con <- file(StdOutFileName, "rt")
      writeLines(readLines(Con))
      close(Con)
      unlink(StdOutFileName)
    }
  }
}


getMSYrefs <- function(sim, .Object, nyears=40, AsCluster=FALSE, toly=1e-1)
{
  if (AsCluster)
  {
    require(keep)

    # In a cluster context we re-direct output to file so we can then play
    # it back in the host
    StdOutFileName <- paste("StdOutFile", sim, ".txt", sep="")
    sink(file=StdOutFileName, append=FALSE, type=c("output", "message"))
  }

  test<-optimize(popdyn,
                 interval=log(c(0.001,10.)),
                 reportIndicators=FALSE,
                 npop=.Object@npop,
                 nages=.Object@nages,
                 nyears=nyears,
                 nsubyears=.Object@nsubyears,
                 nareas=.Object@nareas,
                 nfleets=.Object@nfleets,
                 R0=.Object@R0[sim,],
                 M=.Object@M[sim,,,],
                 mat=.Object@mat[sim,,,],
                 Idist=.Object@Idist[sim,,,],
                 Len_ageVec=.Object@Len_age[sim,1,,1],
                 Wt_ageVec=.Object@Wt_age[sim,1,,1],
                 #Wt_age_SBVec=.Object@Wt_age_SB[sim,1,,1],
                 q=.Object@q[sim,],
                 nSpawnPerYr = length(.Object@Recsubyr),
                 a=.Object@a,
                 b=.Object@b,
                 Recdist=.Object@Recdist[sim,,],
                 ECurrent=.Object@ECurrent[sim,,,],
                 sel=.Object@sel[sim,,],
                 E=.Object@E[sim,,],
                 mov=.Object@mov[sim,,,,,],
                 Recsubyr=.Object@Recsubyr, h=.Object@h[sim,],
                 Recdevs=.Object@Recdevs[sim,,],
                 SRrel=.Object@SRrel,
                 targpop=.Object@targpop,
                 MSYyear=.Object@nyears,
                 tol=toly)

  best <- popdyn(test$minimum,
                 reportIndicators=TRUE,
                 npop=.Object@npop,
                 nages=.Object@nages,
                 nyears=nyears,
                 nsubyears=.Object@nsubyears,
                 nareas=.Object@nareas,
                 nfleets=.Object@nfleets,
                 R0=.Object@R0[sim,],
                 M=.Object@M[sim,,,],
                 mat=.Object@mat[sim,,,],
                 Idist=.Object@Idist[sim,,,],
                 Len_ageVec=.Object@Len_age[1,1,,1],
                 Wt_ageVec=.Object@Wt_age[1,1,,1],
                 #Wt_age_SBVec=.Object@Wt_age_SB[1,1,,1],
                 q=.Object@q[sim,],
                 nSpawnPerYr = length(.Object@Recsubyr),
                 a=.Object@a,
                 b=.Object@b,
                 Recdist=.Object@Recdist[sim,,],
                 ECurrent=.Object@ECurrent[sim,,,],
                 sel=.Object@sel[sim,,],
                 E=.Object@E[sim,,],
                 mov=.Object@mov[sim,,,,,],
                 Recsubyr=.Object@Recsubyr,
                 h=.Object@h[sim,],
                 Recdevs=.Object@Recdevs[sim,,],
                 SRrel=.Object@SRrel,
                 targpop=.Object@targpop,
                 MSYyear=.Object@nyears)

  print("Effort multiplier at MSY:")
  print(exp(test$minimum))

  if (AsCluster)
  {
    sink()
  }

  return(best)
}

#--------------- MSYreferencePoints function------------------------------------------

MSYreferencePoints <- function(ECurrent,
                               sel,
                               M,
                               C,
                               SSB0,
                               B0,
                               N,
                               NBefore,
                               SSN,
                               Wt_age,
                               #Wt_age_SB,
                               targpop,
                               nsubyears,
                               npop,
                               nages,
                               nareas,
                               nfleets)
{
  # Assume M,C,SSN,Wt_age are single year slices. In addition assume N and NBefore
  # have nsubyear + 1 entries for the subyear index with the +1 th entry
  # corresponding to the first time step of the following year. In all cases the
  # simulation index is removed / not present.
  MSY    <- sum(karray(C[targpop,,,,], c(length(targpop), nages, nsubyears, nareas, nfleets)) * karray(Wt_age[targpop,], c(length(targpop),nages,nsubyears,nareas,nfleets)))
  BMSY   <- sum(karray((NBefore[targpop,,1,]), c(length(targpop),nages,nareas)) * karray(Wt_age[targpop,], c(length(targpop),nages,nareas)))
  softmp <- karray(NA, dim=c(nages,nsubyears,nareas,nfleets))

  for(im in 1:nsubyears)
  {
    for(ir in 1:nareas)
    {
      for(ifleets in 1:nfleets)
      {
        softmp[,im,ir,ifleets] <- ECurrent[im,ir,ifleets] * sel[ifleets,]
      }
    }
  }

  sof <- apply(softmp, sum, MARGIN=c(1:3))

  rm(softmp)

  # DK: prod(N, sel, ECurrent for MSY year)   TC uses season 1 only; take average over seasons instead
  # this is probably wrong and not a good option anyway
  VBMSY <- sum(karray(rep(sof,each=length(targpop)), c(length(targpop),nages,nsubyears,nareas)) *
               karray((N[targpop,,1:nsubyears,]), c(length(targpop),nages,nsubyears,nareas)) *
               karray(rep(Wt_age[targpop,],nsubyears), c(length(targpop),nages,nsubyears,nareas)))

  # All this stuff used to calculate FMSY (summed over regions, mean over age and season)
  NsoRbyPAM <- apply(NBefore[,,,], FUN=sum, MARGIN=c(1:3))

  # Note that as we are in steady state numbers at y + 1 subyear 1 should be the same as
  # numbers at y subyear 1 . The values at nsubyears + 1 are incomplete because they need
  # to wrap around to the next year first. As such, we base the nsubyears + 1 case on
  # the first subyear.
  NsoRbyPAM[,,nsubyears + 1] <- apply(NBefore[,,1,], FUN=sum, MARGIN=c(1:2))

  ZsoRbyPAM <- karray(NA, dim=c(npop,nages-2,nsubyears))

  ZsoRbyPAM[,,1] <- -log(NsoRbyPAM[,2:(nages-1),2] / NsoRbyPAM[,1:(nages-2),1])
  ZsoRbyPAM[,,2] <- -log(NsoRbyPAM[,2:(nages-1),3] / NsoRbyPAM[,1:(nages-2),2])
  ZsoRbyPAM[,,3] <- -log(NsoRbyPAM[,2:(nages-1),4] / NsoRbyPAM[,1:(nages-2),3])
  ZsoRbyPAM[,,4] <- -log(NsoRbyPAM[,2:(nages-1),5] / NsoRbyPAM[,1:(nages-2),4])

  MbyPAM    <- karray(rep(M[,1:(nages-2)],times=nsubyears), dim=c(npop,nages-2,nsubyears))
  FsoRbyPAM <- ZsoRbyPAM - MbyPAM / nsubyears

  FMSY1     <- mean(FsoRbyPAM[,2:27,])  # 2:27 = true ages 1:26 (1:26 used by SS)

  #potential change to integrated biomass calculation
  SSBMSY    <- sum(karray(SSN[targpop,,1,],c(length(targpop),nages,nareas)) * karray(Wt_age[targpop,], c(length(targpop),nages,nareas)))
  #SSBMSY     <- sum(karray(NBefore[targpop,,1,],c(length(targpop),nages,nareas)) * karray(Wt_age_SB[targpop,], c(length(targpop),nages,nareas)))

  UMSY      <- MSY / VBMSY
  SSBMSY_B0 <- SSBMSY / sum(SSB0[targpop])

  return(c(MSY,BMSY,VBMSY,SSBMSY,UMSY,FMSY1,SSBMSY_B0,SSB0,B0))
}

#--------------- popdyn function------------------------------------------

popdyn <- function(par,
                   reportIndicators,
                   npop,
                   nages,
                   nyears,
                   nsubyears,
                   nareas,
                   nfleets,
                   R0,
                   M,
                   mat,
                   Idist,
                   Len_ageVec,
                   Wt_ageVec,
                   #Wt_age_SBVec,
                   q,
                   nSpawnPerYr=4,
                   a,
                   b,
                   Recdist,
                   ECurrent,
                   sel,
                   E,
                   mov,
                   Recsubyr,
                   h,
                   Recdevs,
                   SRrel,
                   targpop=NA,
                   MSYyear=1,
                   loud=F)
{
  #Effort multiplier for ECurrent
  totF <- exp(par)

  # Need to get initial equilibrium numbers right because SSB0 and rec depend on it
  # NBefore recorded before M, hence the M=0 for a=1; integral for lastAge added (important for low M scenatios)
  Madvanced           <- karray(as.double(NA),c(npop,nages))
  Madvanced[1:npop,1] <- 0.0
  Madvanced[,2:nages] <- M[,1:(nages-1),1]
  surv                <- t(exp(-apply(Madvanced[,1:nages], c(1), cumsum) / nsubyears))

  #infinite sum for plus group
  surv[,nages] <- surv[,nages-1]*exp(-Madvanced[,nages]/nsubyears)/(1-exp(-Madvanced[,nages]/nsubyears))

  Wt_age    <- karray(NA, c(npop,nages,nyears))
  #Wt_age_SB <- karray(NA, c(npop,nages,nyears))

  ind    <- as.matrix(expand.grid(1:npop, 1:nages, 1:nyears))
  indo   <- as.matrix(expand.grid(1:npop, 1:nages, 1:nyears))

  indo[,3] <- rep(MSYyear, npop * nages * nyears)
  M        <- karray(M[indo], c(npop,nages,nyears))
  mat      <- karray(mat[indo], c(npop,nages,nyears))

  indAPY      <- as.matrix(expand.grid(1:nages, 1:npop, 1:nyears))
  Wt_age[ind]    <- karray(rep(Wt_ageVec, times=npop * nyears), dim=c(nages,npop,nyears))[indAPY]
  #Wt_age_SB[ind] <- karray(rep(Wt_age_SBVec, times=npop * nyears), dim=c(nages,npop,nyears))[indAPY]

  #NBefore probably not required in equilibrium calcs, but added for ease of reading
  N_Y       <- karray(NA,c(npop,nages,nsubyears + 1,nareas))      # only need aggregated catch for these purposes
  NBefore_Y <- karray(NA,c(npop,nages,nsubyears + 1,nareas))      # only need aggregated catch for these purposes
  SSN_Y     <- karray(NA,c(npop,nages,nsubyears,nareas))          # only need aggregated catch for these purposes
  SSB_Y     <- karray(NA,c(npop,nages,nsubyears,nareas))          # only need aggregated catch for these purposes
  B_Y       <- karray(NA,c(npop,nages,nsubyears,nareas))          # only need aggregated catch for these purposes
  C_Y       <- karray(NA,c(npop,nages,nsubyears,nareas,nfleets))  # Catch

  Z    <- karray(NA,c(npop,nages,nareas))
  FD   <- karray(NA,c(nfleets,nareas))            # F distribution
  FM   <- karray(NA,c(npop,nages,nareas,nfleets)) # Fishing Mortality

  y <- 1
  m <- 1

  PAYMR <- as.matrix(expand.grid(1:npop,1:nages,y,m,1:nareas))    # Set up some karray indexes
  PAMR  <- PAYMR[,c(1,2,3,5)]
  PA    <- PAYMR[,c(1,2)]
  P     <- PAYMR[,c(1)]
  PAR   <- PAYMR[,c(1,2,5)]
  PAY   <- PAYMR[,c(1,2,3)]

  if (nareas > 1)
  {
    # Calculate virgin Numbers
    NBefore_Y[PAMR] <- surv[PA] * R0[P] * Idist[PAR]
  }
  else
  {
    NBefore_Y[PAMR] <- surv[PA] * R0[P]
  }

  SSN_Y[PAMR] <- NBefore_Y[PAMR] * mat[PAY]

  # potential change to combined maturity and mass-at-age
  # Calculate spawning stock biomass
  SSB_Y[PAMR] <- SSN_Y[PAMR] * Wt_age[PAY]
  #SSB_Y[PAMR] <- NBefore_Y[PAMR] * Wt_age_SB[PAY]

  # Calculate total biomass
  B_Y[PAMR] <- NBefore_Y[PAMR] * Wt_age[PAY]

  SSB0 <- apply(SSB_Y[,,m,], 1, sum)
  B0   <- apply(B_Y[,,m,], 1, sum)

  # Calculate spawning stock biomass per recruit
  SSBpR <- SSB0 / R0

  # Ricker SR params
  bR <- log(5.0 * h) / (0.8 * SSB0)

  # Ricker SR params
  aR <- exp(bR * SSB0) / SSBpR

  #minimizer quite robust to this value
  NBefore_Y[PAMR] <- NBefore_Y[PAMR] * 0.3
  SSB_Y[PAMR]     <- SSB_Y[PAMR] * 0.3

  PAYMRF2 <- as.matrix(expand.grid(1:npop, 1:nages, y, m, 1:nareas, 1:nfleets))
  PAMRF2  <- PAYMRF2[,c(1,2,4,5,6)]
  PAMR2   <- PAYMRF2[,c(1,2,4,5)]
  PARF2   <- PAYMRF2[,c(1,2,5,6)]
  PAR2    <- PAYMRF2[,c(1,2,5)]
  PAY2    <- PAYMRF2[,c(1,2,3)]
  FA2     <- PAYMRF2[,c(6,2)]
  F2      <- PAYMRF2[,c(6)]
  MRF2    <- PAYMRF2[,c(4,5,6)]

  for (y in 1:nyears)
  {
    PAYMRF2[,3] <- y
    PAY2[,3]    <- y

    for (m in 1:nsubyears)
    {
      PAYMRF2[,4] <- m
      PAMRF2[,3]  <- m
      PAMR2[,3]   <- m
      MRF2[,1]    <- m

      # do recruitment
      if (m %in% Recsubyr)
      {
        # ie every qtr for YFT
        SSN_Y[,,m,] <- NBefore_Y[,,m,] * karray(rep(mat[,,y],times=nareas), dim=c(npop,nages,nareas))

        #potential change to combined spawning and Wt_age
        SSBA_Y      <- apply(SSN_Y[,,m,] * karray(Wt_age[,,y], dim=c(npop,nages,nareas)), 1, sum, na.rm=T)
        #SSBA_Y      <- apply(NBefore_Y[,,m,] * karray(Wt_age_SB[,,y], dim=c(npop,nages,nareas)), 1, sum, na.rm=T)

        # recruit fish
        for(pp in 1:npop)
        {
          if (SRrel[pp] == 1)
          {
            # Beverton-Holt recruitment
            rec <- ((0.8 * R0[pp] * h[pp] * SSBA_Y[pp]) / (0.2 * SSBpR[pp] * R0[pp] * (1.0 - h[pp]) + (h[pp] - 0.2) * SSBA_Y[pp]))
          }
          else
          {
            #readline("warning: Ricker not tested")
            # Most transparent form of the Ricker uses alpha and beta params
            rec <- aR[pp] * SSBA_Y[pp] * exp(-bR[pp] * SSBA_Y[pp])
          }

          NBefore_Y[pp,1,m,] <- rec * Recdist[pp,]
        }
      }

      # move fish (moved to before F+M)
      if (nareas > 1)
      {
        N_Y[,,m,] <- popdyn.domov(Ntemp=NBefore_Y[,,m,], movtemp=mov[,,m,,])
      }
      else
      {
        N_Y[,,m,] <- NBefore_Y[,,m,]
      }

      FM[PARF2]     <- totF * ECurrent[MRF2] * sel[FA2] * q[F2]
      FM[is.na(FM)] <- 0
      Ftot          <- apply(FM, c(1,2,3), sum)
      Z[PAR2]       <- Ftot[PAR2] + M[PAY2] / nsubyears

      C_Y[PAMRF2]  <- N_Y[PAMR2] * (1.0 - exp(-Z[PAR2])) * (FM[PARF2] / Z[PAR2])

      N_Y[,,m,] <- N_Y[,,m,] * exp(-Z)

      # Age fish
      NBefore_Y[,nages,m + 1,]          <- N_Y[,nages - 1,m,] + N_Y[,nages,m,]
      NBefore_Y[,2:(nages - 1),m + 1,]  <- N_Y[,1:(nages - 2),m,]
      NBefore_Y[,1,m + 1,]              <- 0
      N_Y[,,m + 1,]                     <- NBefore_Y[,,m + 1,]

    } # season loop

    if (y != nyears)
    {
      NBefore_Y[,,1,] <- NBefore_Y[,,nsubyears + 1,]
    }
  } # year loop

  if (reportIndicators)
  {
    return (MSYreferencePoints(ECurrent,
                               sel,
                               M[,,1],
                               C_Y,
                               SSB0,
                               B0,
                               N_Y,
                               NBefore_Y,
                               SSN_Y,
                               Wt_age[,,1],
                               #Wt_age_SB[,,1],
                               targpop,
                               nsubyears,
                               npop,
                               nages,
                               nareas,
                               nfleets))
  }
  else
  {
    BMSY <- sum(karray(C_Y[targpop,,,,], c(length(targpop),nages,nsubyears,nareas,nfleets)) *
                karray(Wt_age[targpop,,nyears], c(length(targpop),nages,nsubyears,nareas,nfleets)))

    return(-log(BMSY))
  }
}
