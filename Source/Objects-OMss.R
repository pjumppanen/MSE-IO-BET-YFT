# =============================================================================
# ==== YFT MSE object classes modified to use conditioned SS OMs===============
# =============================================================================



# ------------------ Operating model definition object ------------------------
# This object has been merged with the OM definition class and OM class, so that
# key parameters can be extracted from the SS output files. It includes minimal
# set of inputs for dynamics and observation models most parameter specifications
# are drawn from the Stock Synthesis results files

setClass("OMd",representation(

               # Description
               Name             = "character",
               Label             = "character",
               Date             = "character",
               Author           = "character",
               Notes            = "character",
               PrimarySource    = "character",
               npop             = "integer",
               proyears         = "integer",
               targpop          = "integer",
               seed             = "integer",
               recentPerFirst   = "integer",
               recentPerLast    = "integer",
               seasonCEDist     = "integer",
               SSRootDir        = "character",
               SBlim            = "numeric",
               Flim             = "numeric",
               OMList           = "list",
               nsimPerOMFile    = "karray",
               CppMethod        = "numeric",
               UseCluster       = "numeric",

               nfleets          = "integer",
               ReccvTin         = "karray",
               RecACTin         = "karray",
               ReccvRin         = "numeric",
               NInitCV          = "numeric",   # additional noise on initial N(a) (CV on age 1)
               NInitCVdecay     = "numeric",   # exponential decay on CV on initial N(a) = exp(NinitCV*(a-1))
               selExpRange      = "numeric",   # sel temporal variability exponentoscillates with a sin wave rangeing between these values
               selAgeRange      = "numeric",   # 0=no age shift, 2 means (discretized) sine wave shift of sel vector between - 2 and + 2 age class
               selWLRange       = "karray",    # sel temporal variability wavelength range (0.0625 = quarter wavelength in 25 years, 0.5=2 full cycles in 25 years

               nsubyears        = "integer",
               firstCalendarYr  = "integer",
               lastCalendarYr   = "integer",
               lastSeas         = "integer",
               firstSeas        = "integer",
               firstSSYr        = "integer",
               firstMPYr        = "integer",   # MP kicks in here (projections in the intervening period may require an observed catch input to handle data lags)
               MPDataLag        = "integer",   # MP kicks in here (projections in the intervening period may require an observed catch input to handle data lags)
               catchBridge      = "karray",    # known catch history between last assessment year in OM and firstMPYr (length of min 0 to max firstMPYr-lastCalendarYr-1)
               catchBridgeCV    = "numeric",   # error to add onto catch for bridge years with unknown catch

               indexFisheries   = "integer",   # vector of fisheries to aggregate for index calculations (e.g. LL CPUEfleets)

                                               #Observation errors (original ABT code had a separate observation class)
               TACEcv           = "karray",    # fleet-specific lognormal errors on TAC/TAE (independent among fleets and seasons)
               Ccv              = "numeric",   # Observation error
               Cbcv             = "numeric",   # bias in total annual catches
               nCAAobs          = "numeric",   # Number of annual catch at age (CAA) observations
               nCALobs          = "numeric",   # Number of annual catch-at-length (CAL) observations
               Lcv              = "numeric",
               Ibeta            = "numeric",   # Hyperstability parameter I^beta
               Icv              = "numeric",   # Observation error in relative abundance indices
               IACin            = "numeric",   # cpue index autorrelation error
               ITrendin         = "numeric",   # cpue index catchability trend (for observation error only, this parameter does not cause increases in F)
               Mbcv             = "numeric",   # Bias in observation of natural mortality rate
               Kbcv             = "numeric",   # Bias in estimation of growth parameters
               t0bcv            = "numeric",
               Linfbcv          = "numeric",
               LFCbcv           = "numeric",   # Bias in observation of length at first capture (LFC)
               LFSbcv           = "numeric",   # Bias in observation of length at full selection (LFS)
               FMSYbcv          = "numeric",   # Bias in observaton of FMSY
               FMSY_Mbcv        = "numeric",   # Bias in ratio of FMSY/M
               BMSY_B0bcv       = "numeric",   # Bias ratio of BMSY/B0
               ageMbcv          = "numeric",   # Bias in observation of age at 50% maturity and
               Dbcv             = "numeric",   # Bias in observation of current stock depletion
               Dcv              = "numeric",   # Imprecision in observation of current stock depletion
               Btbcv            = "numeric",   # Bias in observation of current stock biomass
               Btcv             = "numeric",   # Imprecision in observation of current stock biomass
               Ftbcv            = "numeric",   # Bias in observation of current fishing mortality rate
               Ftcv             = "numeric",   # Imprecision in observation of current fishing mortality rate
               hbcv             = "numeric",   # Bias in observation of steepness
               Recbcv           = "numeric",   # Bias in observation of recent recrutiment
               IMSYbcv          = "numeric",   # Bias in observation of target CPUE (CPUE @ MSY)
               MSYbcv           = "numeric",   # Bias in observation of target catch (MSY)
               BMSYbcv          = "numeric",   # Bias in observation of target biomass (BMSY)

               # management procedure tuning
               tunePMProjPeriod = "numeric",
               tunePM           = "character", # name of performance measure to tune to
               tunePMTarget     = "numeric",   # level of performance measure that tuning should achieve
               tuneTol          = "numeric",   # relative level of precision required in tuning
               tuneLogDomain    = "numeric"    # log base 10 of search domain for tuning solution. It is a log domain to improve search dynamic range
             ),
             prototype = list(
               tuneTol       = 0.01,
               tuneLogDomain = c(-3,0.5)
             )
)





#------------------- Operating model object -------------------------------------------------------------------------------------

# An OMss object is created from a minimal OMd (Operating Model definition) object
# with most specifications imported from the list of SS3 assessment models in the OMd

setClass("OMss",representation(

              # Description
              Name            = "character",
              Label           = "character",
              Date            = "character",
              Author          = "character",
              Notes           = "character",
              PrimarySource   = "character",
              npop            = "integer",
              proyears        = "integer",
              targpop         = "integer",
              seed            = "integer",
              recentPerFirst  = "integer",
              recentPerLast   = "integer",
              seasonCEDist    = "integer",
              SSRootDir       = "character",
              SBlim           = "numeric",
              Flim            = "numeric",
              OMList          = "list",
              OMid            = "karray",
              nsimPerOMFile   = "karray",
              CppMethod       = "numeric",
              UseCluster      = "numeric",
              UseMSYss        = "numeric",  # not 0 means use the SS-based reference points rather than the R or Cpp verions
              indexFisheries  = "integer",  # vector of fisheries to aggregate for index calculations (e.g. LL CPUEfleets)

              #ReccvTin        = "karray",   # dimmed below
              #RecACTin        = "karray",   # dimmed below
              ReccvRin        = "numeric",  # dimmed below

              # MSE Dimensions
              nsim            = "integer",
              nages           = "integer",
              nyears          = "integer",
              nsubyears       = "integer",
              nareas          = "integer",

              #Calendar time labels for time series
              lastCalendarYr   = "integer",
              yrLabels        = "karray",
              seasonLabels    = "karray",
              yrSeasLabels    = "karray",
              firstMPYr       = "numeric",
              catchBridge     = "karray",  # known catch history between last assessment year in OM and firstMPYr (length of min 0 to max firstMPYr-lastCalendarYr-1)
              catchBridgeCV   = "numeric", # error to add onto catch for bridge years with unknown catch
              MPDataLag       = "numeric", # The lag in number of years between assessment data availability and the timing of the assessment/HCR calculation; +1 means data is one year behind
              NInitCV         = "numeric", # additional noise on initial N(a) (CV on age 1)
              NInitCVdecay    = "numeric", # exponential decay on CV on initial N(a) = exp(NinitCV*(a-1))
              selExpRange     = "numeric", # sel temporal variability exponentoscillates with a sin wave rangeing between these values
              selAgeRange     = "numeric", # 0=no age shift, 2 means (discretized) sine wave shift of sel vector between - 2 and + 2 age class
              selWLRange      = "karray",  # sel temporal variability wavelength range (0.0625 = quarter wavelength in 25 years, 0.5=2 full cycles in 25 years
              selTSSign       = "karray",  # first exp then age; sim and fishery-specific values for initial direction of sel sin wave functions
              selTSWaveLen    = "karray",  # first exp then age; sim and fishery-specific values for sel TS WaveLen

              # OM parameters extracted from assessment outputs
              SRrel           = "integer", # Stock-recruitment relationship type
              h               = "karray",  # steepness
              recgrad         = "karray",  # underlying gradient % yr-1
              ReccvT          = "karray",  # CV of aggregate recruitment deviations (i.e. Rec(t))
              RecAC           = "karray",  # recruitment auto-correlation
              ReccvR          = "karray",  # CV of recruitment spatial deviations
              #Recsubyr        = "integer", # Sub-year in which recruitment occurs, area in which recruitment occurs
              Linfmu          = "karray",  # Mean growth parameters
              Kmu             = "karray",
              t0              = "numeric",
              Ksd             = "karray",
              Kgrad           = "karray",
              Linfsd          = "karray",  # Interannual variability in growth
              Linfgrad        = "karray",  # growth mean trajectory % yr-1
              a               = "numeric", # Weight - Length conversion W=aL^b
              b               = "numeric",
              ageMmu          = "karray",
              ageM            = "karray",  # Age-at-maturity
              ageMsd          = "karray",  # Interannual variability
              ageMgrad        = "karray",  # gradient % yr-1
              D               = "karray",  # Current stock depletion
              SSBcur          = "karray",  # Current stock abundance
              Size_area       = "karray",  # Size of regions
              mov             = "karray",  # Markov movement matrix for all fish
              Mmov            = "karray",  # Markov movement matrix for mature fish
              movvar          = "karray",  # Inter-simulation variability in movement
              movsd           = "karray",  # Interannual-variability in movement
              movgrad         = "karray",  # Gradient changes in area gravity weights
              Mmovvar         = "karray",  # Inter-simulation variability in mature movement
              Mmovsd          = "karray",  # Interannual-variability in mature movement
              Mmovgrad        = "karray",  # Gradient changes in area gravity weights
              excl            = "karray",  # Exclusion matrix [0,1] depending on whether the stock can go there
              nfleets         = "integer", # Number of fleets,
              nCPUE           = "integer", # Number of CPUE series,
              age05           = "karray",  # age at 5% vulnerability
              Vmaxage         = "karray",  # vulnerability of oldest age class
              AFS             = "karray",  # age at full selection
              Fsd             = "karray",  # Interannual variability in F
              Fgrad           = "karray",  # Final gradient in F yr-1
              Area_names      = "character", # Area definitions (polygons)
              Area_defs       = "list",

              # Simulation data -------
              E               = "karray",
              dFfinal         = "karray",
              q               = "karray",
              sel             = "karray",
              CPUEsel         = "karray",
              mat             = "karray",
              Recdevs         = "karray",
              R0              = "karray",  # Recruitment deviations and unfished recruitment
              M               = "karray",  # Natural instantaneous mortality rate.
              Z               = "karray",  # Total instantaneous mortality rate.
              B               = "karray",  # Biomass
              SSB             = "karray",  # Spawning stock biomass
              NBeforeInit     = "karray",
              C               = "karray",  # Catches taken
              CAA             = "karray",  # Catch at age taken
              MSY             = "numeric", # Maximum sustainable yield
              BMSY            = "numeric", # Biomass at maximum sustainable yield
              VBMSY           = "numeric", # Vulnerable biomass at maximum sustainable yield
              SSBMSY          = "numeric",
              UMSY            = "numeric",
              FMSY1           = "numeric", # Fishing mortality rate at maximum sustainable yield
              SSB0            = "numeric",
              B0              = "numeric",
              SSBMSY_SSB0     = "numeric",
              IMSY            = "numeric", # Relative abundance index at maximum sustainable yield

              Idist           = "karray",  # Initial distribution (mean unfished N by region given recruitment by region and movement)

              # DK changes for YFT
              ECurrent        = "karray",  #recent,e.g. 2 year, period for which relative effort among fleets is assumed to remain constant for MSY and projections
              CMCurrent       = "karray",  #relative catch among fleets in recent period for TAC allocations
              EByQtrLastYr    = "karray",  #last year of effort, for projecting from N[last assessment qtr] -> N[first projection qtr]
              Recsubyr        = "integer", # Sub-years in which recruitment occurs, area in which recruitment occurs?
              Recdist         = "karray",  # mean spatial distribution for recruitment...
              initIDev        = "karray",  # initial deviation for autocorrelated aggregate annual CPUE errors

              # Reference points from SS for comparison with MSE model (note potential definition differences from MSE model, particularly due to seasonality)
              SSB0ss          = "numeric", # SSB0 extract from SS
              B0ss            = "numeric", # B0 extract from SS
              MSYss           = "numeric", # MSY from SS output
              BMSYss          = "numeric", # total biomass at MSY from SS output
              SSBMSYss        = "numeric", # SSBMSY from SS
              FMSYss          = "numeric", # FMSY from SS
              F_FMSYss        = "karray",  # F/FMSY from SS
              Frepss          = "karray",  # Frep from SS

              # DK additions for SS3 YFT situation
              Len_age         = "karray",  #extract from SS
              Len_age_mid     = "karray",  #extract from SS
              Wt_age          = "karray",  #extract from SS
              Wt_age_SB       = "karray",  #extract from SS
              Wt_age_mid      = "karray",  #extract from SS
              Css             = "karray",  #Catch SPYF from SS
              CAAFss          = "karray",  #Catch SAYF from SS
              SSBAss          = "karray",  #age aggregated SSB
              CBss            = "karray",  #Catch biomass SPY from SS
              Bss             = "karray",  #biomass from SS
              Recss           = "karray",  #recruitment from SS
              RecYrQtrss      = "karray",  #quarterly recruitment from SS
              NLLss           = "karray",  #Longline-selected Numbers over all ages,  populations by area and regions,for aggregate abundance index
              NLLIss          = "karray",  #Longline-selected Numbers over all ages, areas, populations for aggregate abundance index
              CPUEobsMR       = "karray",  #observed CPUE series in SS model by area by season
              CPUEobsY        = "karray",  #observed CPUE series in SS model summed over area and season
              CPUEFleetNums   = "karray",  #fishery corresponding to each CPUE series (for selectivity)
              CPUEFleetAreas  = "karray",  #areas corresponding to each CPUE series
              Edist           = "karray",  #Distribution of relative F in final year by seasn, region and fishery (redundant with region)

#Observation errors (original ABT code had a separate observation class)
#xxx Remove the redundant ones
               TACEcv         = "karray",  # fleet-specific lognormal errors on TAC/TAE (independent among fleets and seasons)
               Ccv            = "numeric", # Observation error in total annual catches
               Cbcv           = "numeric", # Bias in total annual catches
               nCAAobs        = "numeric", # Number of annual catch at age (CAA) observations
               nCALobs        = "numeric", # Number of annual catch-at-length (CAL) observations
               Lcv            = "numeric",
               Ibeta          = "numeric", # Hyperstability parameter I^beta
               Icv            = "numeric", # Observation error in relative abundance indices
               IAC            = "karray",  # cpue index autorrelation error
               ITrend         = "karray",  # cpue catchability trend % per annum compounded, -1 means infer the value from the assessment model
               Mbcv           = "numeric", # Bias in observation of natural mortality rate
               Kbcv           = "numeric", # Bias in estimation of growth parameters
               t0bcv          = "numeric",
               Linfbcv        = "numeric",
               LFCbcv         = "numeric", # Bias in observation of length at first capture (LFC)
               LFSbcv         = "numeric", # Bias in observation of length at full selection (LFS)
               FMSYbcv        = "numeric", # Bias in observaton of FMSY
               FMSY_Mbcv      = "numeric", # Bias in ratio of FMSY/M
               BMSY_B0bcv     = "numeric", # Bias in ratio of BMSY/B0
               ageMbcv        = "numeric", # Bias in observation of age at 50% maturity and
               Dbcv           = "numeric", # Bias and imprecision in observation of current stock depletion
               Dcv            = "numeric",
               Btbcv          = "numeric", # Bias and imprecision in observation of current stock biomass
               Btcv           = "numeric",
               Ftbcv          = "numeric", # Bias and imprecision in observation of current fishing mortality rate
               Ftcv           = "numeric",
               hbcv           = "numeric", # Bias in observation of steepness
               Recbcv         = "numeric", # Bias in observation of recent recrutiment
               IMSYbcv        = "numeric", # Bias in observation of target CPUE (CPUE @ MSY)
               MSYbcv         = "numeric", # Bias in observation of target catch (MSY)
               BMSYbcv        = "numeric", # Bias in observation of target biomass (BMSY)

               # management procedure tuning
               tunePMProjPeriod = "numeric",
               tunePM           = "character", # name of performance measure to tune to
               tunePMTarget     = "numeric",   # level of performance measure that tuning should achieve
               tuneTol          = "numeric",   # relative level of precision required in tuning
               tuneLogDomain    = "numeric"    # log base 10 of search domain for tuning solution. It is a log domain to improve search dynamic range
              ))

setMethod("initialize", "OMss", function(.Object,OMd, Report=F, UseMSYss=0)
{
  print("start OMss definition")

  #SS reports in 1000 individuals and tonnes), which should be consistent with MSE
  #Confirm rate translations from seasonal to annual
  #ABT MSE inputs annual M definition, assumes constant over seasons, and applies M/nSubyears
  #ABT uses year-specific age attributes in seasonal model (i.e annual age structure: growth assumed static within years)
  #YFT uses season-specific age attributes, i.e. (quarterly age structure; quarter-specific biology)
  #Both models input annual mortality rate definitions, but YFTs are quarterly age-specific i.e. multiply seasonAsYear SS input by 4 )

  if (class(OMd) != 'OMd')
  {
    print(paste('Could not build operating model:',deparse(substitute(OMd)),'not of class OMd'))
    stop()
  }

  if((OMd@firstMPYr - OMd@lastCalendarYr) < 2)
  {
    print("Aborted: firstMPYr - OMd@firstCalendarYr must be >= 2 ")
    stop()
  }


  .Object@Name          <- OMd@Name
  .Object@Label         <- OMd@Label
  .Object@Date          <- OMd@Date
  .Object@Author        <- OMd@Author
  .Object@Notes         <- OMd@Notes
  .Object@PrimarySource <- OMd@PrimarySource
  .Object@npop          <- OMd@npop
  .Object@proyears      <- OMd@proyears
  .Object@targpop       <- OMd@targpop
  .Object@seed          <- OMd@seed
  .Object@recentPerFirst<- OMd@recentPerFirst
  .Object@recentPerLast <- OMd@recentPerLast
  .Object@seasonCEDist  <- OMd@seasonCEDist
  .Object@SSRootDir     <- OMd@SSRootDir
  .Object@SBlim         <- OMd@SBlim
  .Object@Flim          <- OMd@Flim
  .Object@OMList        <- OMd@OMList
  .Object@nsimPerOMFile <- OMd@nsimPerOMFile
  .Object@OMid          <- karray(rep(NA, sum(OMd@nsimPerOMFile), dim=1))
  .Object@CppMethod     <- OMd@CppMethod
  .Object@UseCluster    <- OMd@UseCluster
  .Object@UseMSYss      <- UseMSYss

  .Object@NInitCV      <- OMd@NInitCV
  .Object@NInitCVdecay <- OMd@NInitCVdecay
  .Object@selExpRange  <- OMd@selExpRange
  .Object@selAgeRange  <- OMd@selAgeRange
  .Object@selWLRange   <- OMd@selWLRange

  .Object@lastCalendarYr <- OMd@lastCalendarYr
  .Object@MPDataLag      <- OMd@MPDataLag
  .Object@catchBridge    <- OMd@catchBridge
  .Object@catchBridgeCV  <- OMd@catchBridgeCV
  .Object@Ccv            <- OMd@Ccv
  .Object@Cbcv           <- OMd@Cbcv
  .Object@Icv            <- OMd@Icv
  .Object@Ibeta          <- OMd@Ibeta
  .Object@Btcv           <- OMd@Btcv
  .Object@Btbcv          <- OMd@Btbcv
  .Object@Mbcv           <- OMd@Mbcv
  .Object@Linfbcv        <- OMd@Linfbcv
  .Object@MSYbcv         <- OMd@MSYbcv
  .Object@BMSYbcv        <- OMd@BMSYbcv
  .Object@IMSYbcv        <- OMd@IMSYbcv
  .Object@FMSYbcv        <- OMd@FMSYbcv
  .Object@FMSY_Mbcv      <- OMd@FMSY_Mbcv
  .Object@nCAAobs        <- OMd@nCAAobs
  .Object@ageMbcv        <- OMd@ageMbcv
  .Object@indexFisheries <- OMd@indexFisheries

  .Object@tunePMProjPeriod  <- OMd@tunePMProjPeriod
  .Object@tunePM            <- OMd@tunePM
  .Object@tunePMTarget      <- OMd@tunePMTarget
  .Object@tuneTol           <- OMd@tuneTol
  .Object@tuneLogDomain     <- OMd@tuneLogDomain

  set.seed(.Object@seed)

  cat("Loading OM files from SS outputs")

  nOMfiles <- length(.Object@OMList)

  .Object@nsim       <- as.integer(sum(.Object@nsimPerOMFile))
  .Object@nsubyears  <- as.integer(OMd@nsubyears)
  .Object@Area_names <- c("NW","SW","SE","NE")             # Not yet used
  .Object@TACEcv     <- OMd@TACEcv

  .Object@MSY         <- rep(as.double(NA), .Object@nsim)
  .Object@BMSY        <- rep(as.double(NA), .Object@nsim)
  .Object@VBMSY       <- rep(as.double(NA), .Object@nsim)
  .Object@SSBMSY      <- rep(as.double(NA), .Object@nsim)
  .Object@UMSY        <- rep(as.double(NA), .Object@nsim)
  .Object@FMSY1       <- rep(as.double(NA), .Object@nsim)
  .Object@SSB0        <- rep(as.double(NA), .Object@nsim)
  .Object@B0          <- rep(as.double(NA), .Object@nsim)
  .Object@SSBMSY_SSB0 <- rep(as.double(NA), .Object@nsim)

  MSYss       <- NULL
  BMSYss      <- NULL
  SSBMSYss    <- NULL
  SSB0ss      <- NULL
  B0ss        <- NULL
  F_FMYSss    <- NULL
  FMYSss      <- NULL

  MSY_projection_sims <- karray(NA,dim=c(nOMfiles,2))

  par(mfrow=c(3,3)) #Plot some SS results while importing

  for(iom in 1:nOMfiles)
  {
    if (exists("ssMod"))
    {
      rm(ssMod)
    }

    #indices to fill karrays with OM-specific number of simulations
    if (iom == 1)
    {
      firstSimNum <- 1

    }
    else
    {
      firstSimNum <- sum(.Object@nsimPerOMFile[1:(iom-1)]) + 1
    }

    lastSimNum <- sum(.Object@nsimPerOMFile[1:iom])

    print(c("importing ", .Object@OMList[iom]))

    ssMod <- SS_output(dir=OMd@SSRootDir %&% .Object@OMList[iom], covar=F, ncols=213,forecast=F) # ssoutput.f no longer req'd, r4ss fixed

    if (iom == 1)
    {
      #define karrays based on these inputs; dimensions must remain the same for other files

      .Object@nfleets   <- as.integer(sum(ssMod$IsFishFleet)) # P nfleets excluding surveys (and some fleets are misleading - should be re-configured as time blocks in selectivity)

      if (.Object@nfleets != OMd@nfleets)
      {
        print("nfleet inconsistency error")
        stop()
      }

      .Object@nCPUE     <- ssMod$nfleets - ssMod$nfishfleets
      .Object@nareas    <- as.integer(ssMod$nareas)

      #time translations:
      #YFT 101 =1972 Q1 used to define MSE startyr (for plotting purposes)
      #YFT 273 =2015 Q1 is first projection yr
      #BET 101 = 1951 Q1
      #BET 345 =2013 Q1 is first proj yr

      .Object@nyears    <- as.integer((ssMod$endyr - OMd@firstSSYr + 1) / .Object@nsubyears)
      nagesss           <- length(ssMod$endgrowth$Age) #note MSE age 1 = SS age 0
      .Object@nages     <- as.integer(nagesss)                   #instead of plus group
      allyears          <- .Object@nyears + .Object@proyears

      .Object@yrLabels       <- karray(OMd@firstCalendarYr - 1 + 1:allyears)
      .Object@seasonLabels   <- karray((1:.Object@nsubyears) / .Object@nsubyears - 0.5 * (1 / .Object@nsubyears))
      .Object@yrSeasLabels   <- karray(rep((OMd@firstCalendarYr - 1 + 1:allyears), each=.Object@nsubyears) + rep(.Object@seasonLabels, times=allyears))
      .Object@firstMPYr      <- OMd@firstMPYr
      .Object@MPDataLag      <- OMd@MPDataLag
      .Object@catchBridge    <- OMd@catchBridge
      .Object@catchBridgeCV  <- OMd@catchBridgeCV

      .Object@M         <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))
      .Object@h         <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop))
      .Object@ReccvT    <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop))
      .Object@ReccvR    <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nareas))

      .Object@Recsubyr  <- c(1:4)  #recruit every quarter

      .Object@Recdist   <- karray(NA, dim=c(.Object@nsim, .Object@npop, ssMod$nareas)) #recruitment distribution by areas

      .Object@Recdevs      <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears * length(.Object@Recsubyr)))
      .Object@Len_age      <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))      #cm
      .Object@Len_age_mid  <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))      #cm
      .Object@Wt_age       <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))      #kg
      .Object@Wt_age_SB    <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))      #kg
      .Object@Wt_age_mid   <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))      #kg
      .Object@mat          <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, allyears))
      .Object@sel          <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nfleets, .Object@nages))
      .Object@CPUEsel      <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nages, .Object@nCPUE))
      .Object@selTSSign    <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nfleets, 2))
      .Object@selTSWaveLen <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nfleets, 2))
      .Object@mov          <- karray(as.double(1.0), dim=c(.Object@nsim, .Object@npop, .Object@nages, .Object@nsubyears, .Object@nareas, .Object@nareas)) # Initialise to 1 because C++ model requires valid mov
      .Object@Mmov         <- karray(as.double(0.0), dim=c(.Object@nsim, .Object@npop, .Object@nages, .Object@nsubyears, .Object@nareas, .Object@nareas))
      .Object@Idist        <- karray(as.double(1.0), dim=c(.Object@nsim, .Object@npop, .Object@nages, .Object@nareas)) # Initialise to 1 because C++ model requires valid Idist
      .Object@R0           <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop))
      .Object@Css          <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears,.Object@nfleets))      # Catch numbers
      .Object@CAAFss       <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nages, allyears,.Object@nfleets))     # Catch at age by fleets
      .Object@SSBAss       <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears))                      # age aggregated SSB
      .Object@CBss         <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears))                      # Catch biomass
      .Object@Bss          <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears))                      # biomass
      .Object@Recss        <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears))                      # recruitment
      .Object@RecYrQtrss   <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, allyears*.Object@nsubyears))    # quarterly recruitment
      .Object@NLLss        <- karray(as.double(NA), dim=c(.Object@nsim, allyears,.Object@nsubyears, .Object@nareas))  # Longline-selected Numbers over all ages, areas, populations for aggregate abundance index
      .Object@NLLIss       <- karray(as.double(NA), dim=c(.Object@nsim, allyears))                                    # Longline-selected Numbers over all ages, areas, populations for aggregate abundance index
      .Object@NBeforeInit  <- karray(as.double(NA), dim=c(.Object@nsim, .Object@npop, .Object@nages, .Object@nareas)) # Initial population at the beginning of projection before movement and mortality for reporting
      .Object@q            <- karray(as.double(0.0), dim=c(.Object@nsim, .Object@nfleets))
      .Object@E            <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nfleets, allyears))

      #should probably force nCPUE = nareas (?)
      .Object@CPUEFleetNums  <- karray(ssMod$fleet_ID[!ssMod$IsFishFleet], dim=c(.Object@nCPUE))
      .Object@CPUEFleetAreas <- karray(ssMod$fleet_area[!ssMod$IsFishFleet], dim=c(.Object@nCPUE))

      .Object@CPUEobsMR     <- karray(as.double(NA), dim=c(.Object@nsim, allyears, .Object@nsubyears, .Object@nCPUE))
      .Object@CPUEobsY      <- karray(as.double(NA), dim=c(.Object@nsim, allyears))
      .Object@initIDev      <- karray(as.double(NA), dim=c(.Object@nsim))
      .Object@ECurrent      <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nsubyears, .Object@nareas, .Object@nfleets))
      .Object@CMCurrent     <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nsubyears, .Object@nareas, .Object@nfleets))
      .Object@EByQtrLastYr  <- karray(as.double(NA), dim=c(.Object@nsim, .Object@nsubyears, .Object@nareas, .Object@nfleets))

      .Object@F_FMSYss <- karray(as.double(NA), dim=c(.Object@nsim, allyears))
      .Object@Frepss   <- karray(as.double(NA), dim=c(.Object@nsim, allyears))

      .Object@IAC      <- karray(as.double(NA), dim=c(.Object@nsim))
      .Object@ITrend   <- karray(as.double(1.0), dim=c(.Object@nsim, allyears))

    }

    #index each sim with an OM specification that can be parsed later
    .Object@OMid[firstSimNum:lastSimNum] <- .Object@OMList[[iom]]

    # --- set up M -------
    #assumes the last year vector is the one to use in projections, ignores sex, morphs, temporal varying growth)
    #for some reason SS M matrix exceeds the final year and begins to deviate...must be used for something else
    # not sure why the original YFT specification seems unduly complicated, check its still ok
    nColOffset <- 3

    #SS does not report M in final time step for some reason (i.e. should be ssMod$endyr - ssMod$startyr +1)
    .Object@M[firstSimNum:lastSimNum,,,] <- .Object@nsubyears * karray(rep(as.numeric(ssMod$M_at_age[ssMod$endyr - ssMod$startyr ,nColOffset+1:(nagesss)]),
                                                                           each=.Object@nsimPerOMFile[iom] * .Object@npop),
                                                                       dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))

    #SS reported M(last age)=NA for some reason
    .Object@M[firstSimNum:lastSimNum,, .Object@nages,] <- .Object@M[keep(firstSimNum:lastSimNum),, .Object@nages - 1, ]

    # ---- Stock-recruit relationships -------
    .Object@SRrel <- as.integer(1) #BH: steepness ssMod$parameters$Label="SR_BH_steep"
    .Object@h[firstSimNum:lastSimNum,] <- karray(ssMod$parameters[ssMod$parameters$Label == "SR_BH_steep",]$Value,
                                                 dim=c(.Object@nsimPerOMFile[iom], .Object@npop))

    #assumes all sims and pops identical for given OMFile
    .Object@h[firstSimNum:lastSimNum,] <- karray(ssMod$parameters[ssMod$parameters$Label == "SR_BH_steep",]$Value,
                                                 dim=c(.Object@nsimPerOMFile[iom], .Object@npop))

    #assumes all sims and pops identical for given OMFile
    #.Object@ReccvT[((iom - 1) * nsimPerOMFile + 1):(iom * nsimPerOMFile),] <- karray(ssMod$parameters[ssMod$parameters$Label == "SR_sigmaR",]$Value, dim=c(nsimPerOMFile, .Object@npop))

    # Rec sigma - if ReccvTin is negative use the SS file, if non-positive use ReccvTin
    if (OMd@ReccvTin[1] >= 0)
    {
      .Object@ReccvT[firstSimNum:lastSimNum,] <- karray(OMd@ReccvTin,
                                                        dim=c(.Object@nsimPerOMFile[iom], .Object@npop))
    }
    else
    {
      .Object@ReccvT[firstSimNum:lastSimNum,] <- karray(ssMod$parameters[ssMod$parameters$Label == "SR_sigmaR",]$Value,
                                                        dim=c(.Object@nsimPerOMFile[iom], .Object@npop))
    }

    #assumes all sims and pops identical for given OMFile
    if (.Object@nareas > 1)
    {
      tmp <- NULL

      for(r in 1:.Object@nareas)
      {
        tmp <- c(tmp, ssMod$parameters$Value[ssMod$parameters$Label == paste("RecrDist_Area_", r, sep="")])
      }

      # note that lower bound parameter value is input for areas which have zero rec by definition in SS (~0.1%) - keep this to prevent NA issues
      tmp <- exp(tmp)
      tmp <- tmp / sum(tmp)
      .Object@Recdist[firstSimNum:lastSimNum,,] <- karray(rep(tmp, each=.Object@nsimPerOMFile[iom]),
                                                          dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nareas))
    }
    else
    {
      .Object@Recdist[] <- 1
    }

    RecACT <- karray(OMd@RecACTin, dim=c(.Object@nsim, .Object@npop))  #input value; could extract from each OM file

    .Object@ReccvRin                          <- OMd@ReccvRin
    .Object@ReccvR[firstSimNum:lastSimNum,,]  <- karray(OMd@ReccvRin, dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nareas))
    .Object@IAC[firstSimNum:lastSimNum]       <- karray(OMd@IACin, dim=lastSimNum - firstSimNum + 1)  #input value; could extract from each OM file

    # Trend in CPUE observation error (Multiplier)
    if (OMd@ITrendin < 0)
    {
      #the catchability trend is determined from the q value in each file name independently
      print("Warning: file-specific index trends assume the single digit following a single q in the filename defines the trend")
      qVal <- as.numeric(substr(OMd@OMList[iom], str_locate(OMd@OMList[iom], "q") + 1, str_locate(OMd@OMList[iom], "q") + 1)) #Should probably revise to 2 digits on next grid iteration
      .Object@ITrend[firstSimNum:lastSimNum, (.Object@nyears + 1):(.Object@nyears + .Object@proyears)] <- rep(cumprod(rep(1 + 0.01 * qVal, .Object@proyears)), times=length(firstSimNum:lastSimNum))
    }
    else
    {
      #use the input value for all sims
      .Object@ITrend[firstSimNum:lastSimNum,(.Object@nyears + 1):(.Object@nyears + .Object@proyears)] <- rep(cumprod(rep(1 + 0.01 * OMd@ITrendin, .Object@proyears)), times=length(firstSimNum:lastSimNum))
    }

    #autocorrelated rec devs for projection years
    rndDevs <- karray(rep(.Object@ReccvT[firstSimNum:lastSimNum],
                          times=.Object@npop * allyears * length(.Object@Recsubyr)) * rnorm(.Object@nsimPerOMFile[iom] * .Object@npop * allyears * length(.Object@Recsubyr)),
                      dim=c(.Object@nsimPerOMFile[iom], .Object@npop, allyears * length(.Object@Recsubyr)))

    for (t in (.Object@nyears * length(.Object@Recsubyr) + 1):(allyears * length(.Object@Recsubyr)))
    {
      rndDevs[,,t] <- RecACT[keep(firstSimNum:lastSimNum),] * rndDevs[,,t-1] + rndDevs[,,t] * sqrt(1 - RecACT[keep(firstSimNum:lastSimNum),] ^ 2)
    }

    .Object@Recdevs[firstSimNum:lastSimNum,,] <- exp(rndDevs[,,] - 0.5 * .Object@ReccvT[keep(firstSimNum:lastSimNum)] ^ 2)

    #Len_age relationships
    #nColOffset <- 4
    #.Object@Len_age[firstSimNum:lastSimNum,,,] <- karray(rep(as.numeric(ssMod$growthseries[nColOffset + 1:nagesss]), each=.Object@nsimPerOMFile[iom] * .Object@npop),
    #                                                     dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))
    .Object@Len_age[firstSimNum:lastSimNum,,,]     <- karray(rep(ssMod$endgrowth$Len_Beg, each=.Object@nsimPerOMFile[iom] * .Object@npop),
         dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))
    .Object@Len_age_mid[firstSimNum:lastSimNum,,,] <- karray(rep(ssMod$endgrowth$Len_Mid, each=.Object@nsimPerOMFile[iom] * .Object@npop),
         dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))

    #weight-length parameters
    # appears that SS does not apply them to the mean Len(a), but integrates to account for the non-linearity - use SS mass inputs instead
    .Object@a <- ssMod$parameters[ssMod$parameters$Label == "Wtlen_1_Fem",]$Value
    .Object@b <- ssMod$parameters[ssMod$parameters$Label == "Wtlen_2_Fem",]$Value

    #This mass-at-age conversoin appears to be overly simplified:
    #.Object@Wt_age[firstSimNum:lastSimNum,,,]     <- .Object@a * .Object@Len_age[keep(firstSimNum:lastSimNum),,,] ^ .Object@b
    #.Object@Wt_age_mid[firstSimNum:lastSimNum,,,] <- .Object@a * .Object@Len_age_mid[keep(firstSimNum:lastSimNum),,,] ^ .Object@b

    # Wt_age_SB includes Wt_age * maturity (useful for dispoportional fecundity)
    .Object@Wt_age_SB[firstSimNum:lastSimNum,,,]  <- karray(rep(as.numeric(ssMod$wtatage[1,7:(7+.Object@nages-1)]), each=.Object@nsimPerOMFile[iom] * .Object@npop),
         dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))
    .Object@Wt_age[firstSimNum:lastSimNum,,,]     <- karray(rep(as.numeric(ssMod$wtatage[2,7:(7+.Object@nages-1)]), each=.Object@nsimPerOMFile[iom] * .Object@npop),
         dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))
    .Object@Wt_age_mid[firstSimNum:lastSimNum,,,] <- karray(rep(as.numeric(ssMod$wtatage[3,7:(7+.Object@nages-1)]), each=.Object@nsimPerOMFile[iom] * .Object@npop),
         dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))

    #Maturity - product of age and length-based maturity seems to cover both options
    #But does not add up to the calculation below; presumably SS does an integration for maturity as well)
    #.Object@mat[firstSimNum:lastSimNum,,,] <- karray(rep(ssMod$endgrowth$Age_Mat * ssMod$endgrowth$Len_Mat, each=.Object@nsimPerOMFile[iom] * .Object@npop),
    #                                                 dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, allyears))

    .Object@mat[firstSimNum:lastSimNum,,,] <- .Object@Wt_age_SB[firstSimNum:lastSimNum,,,]/.Object@Wt_age[firstSimNum:lastSimNum,,,]

    #.Object@Wt_age_SB[firstSimNum:lastSimNum,,,]ssMod$endgrowth$Age_Mat * ssMod$endgrowth$Len_Mat

    #extract CPUE Observations
    print("Warning: missing CPUE values (after time series start) are substituted with preceding value(i.e. not a good idea if there are many missing obs)")
    print("Warning: Observed CPUE series are adopted from the first SS model in the OMList")
    #substitute missing CPUE values for expedience (assumes SS output in temporal order)

    # all OMs use the CPUE data defined in the first SS model list
    # This is required because catchability trends in the SS assessment models are implemented with fake CPUE series
    # but the MP only has one set of data to use

    if (iom == 1)
    {
      totalSims <- sum(.Object@nsimPerOMFile)

      for (i in 1:nrow(ssMod$cpue))
      {
        yrSeas <- seasAsYrToMSEYrSeas.f(seasAsYr=ssMod$cpue$Yr[i], endSeasAsYr = ssMod$endyr, numSeas = OMd@nsubyears, endYr = OMd@lastCalendarYr, endSeas = OMd@lastSeas, firstSeasAsYr = OMd@firstSSYr, firstSeas = OMd@firstSeas) #generic
        #.Object@CPUEobsMR[firstSimNum:lastSimNum,yrSeas[1]:.Object@nyears,yrSeas[2]:.Object@nsubyears, as.numeric(substr(ssMod$cpue$Name[i], 7, 7))] <- ssMod$cpue$Obs[i]
        .Object@CPUEobsMR[,yrSeas[1]:.Object@nyears,yrSeas[2]:.Object@nsubyears, as.numeric(substr(ssMod$cpue$Name[i], 7, 7))] <- rep(ssMod$cpue$Obs[i],times=totalSims)
      }

      #annual
      .Object@CPUEobsY[] <- apply(.Object@CPUEobsMR[,,,], FUN=sum, MARGIN=c(1,2))
    }

    #calculate initial aggregate annual CPUE index deviate for aggregate autocorrelation of index
    lastYrIndices <- (max(ssMod$cpue$Yr) - 3):max(ssMod$cpue$Yr)
    .Object@initIDev[firstSimNum:lastSimNum] <- rep(log(sum(ssMod$cpue$Exp[(ssMod$cpue$Yr %in% lastYrIndices)]) /
                                                    (sum(ssMod$cpue$Obs[(ssMod$cpue$Yr %in% lastYrIndices)]))),
                                                    times=.Object@nsimPerOMFile[iom])

    sel     <- karray(NA, c(.Object@nsimPerOMFile[iom], .Object@nfleets, .Object@nages))
    CPUEsel <- karray(NA, c(.Object@nsimPerOMFile[iom], .Object@nages, .Object@nCPUE))

    #fishery selectivity
    for (ff in 1:.Object@nfleets)
    {
      nColOffset <- 7
      sel[,ff,]  <- karray(rep(as.numeric(ssMod$ageselex[ssMod$ageselex$label == ssMod$endyr %&% "_"  %&% ff %&% "Asel", c(nColOffset + 1:nagesss)]),
                               each=.Object@nsimPerOMFile[iom]),
                           dim=c(.Object@nsimPerOMFile[iom], .Object@nages))
      if (Report)
      {
        if (ff == 1)
        {
          plot(sel[1,ff,], type='l', main="Fishery Sel, iom=" %&% iom)
        }
        else
        {
          lines(sel[1,ff,], col=ff)
        }
      }
    }

   .Object@sel[firstSimNum:lastSimNum,,] <- sel

    #CPUE selectivity
    for(iCPUE in 1:.Object@nCPUE)
    {
      nColOffset        <- 7
      CPUEsel[,,iCPUE]  <- karray(rep(as.numeric(ssMod$ageselex[ssMod$ageselex$label == ssMod$endyr %&% "_"  %&% .Object@CPUEFleetNums[iCPUE] %&% "Asel",c(nColOffset + 1:nagesss)]),
                                      each=.Object@nsimPerOMFile[iom]),
                                  dim=c(.Object@nsimPerOMFile[iom],.Object@nages))

      if (Report)
      {
        if (iCPUE == 1)
        {
          plot(CPUEsel[1,,iCPUE], type='l', main="CPUE Sel, iom=" %&% iom)
        }
        else
        {
          lines(CPUEsel[1,,iCPUE], col=iCPUE)
        }
      }
    }

    .Object@CPUEsel[firstSimNum:lastSimNum,,] <- CPUEsel

    mov <- karray(0, dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nsubyears, .Object@nareas, .Object@nareas))

    if (.Object@nareas > 1)
    {
      #import SS movement parms
      ssmov <- ssMod$movement

      for (i in 1:nrow(ssmov))
      {
        from        <- ssmov$Source_area[i]
        to          <- ssmov$Dest_area[i]
        nColOffSet  <- 6
        moveVec     <- ssmov[i, (nColOffSet + 1):(nColOffSet + nagesss)]
        moveVec     <- as.numeric(moveVec) # MSE maxage older than SS

        for (is in 1:.Object@nsimPerOMFile[iom])
        {
          for (ip in 1:.Object@npop)
          {
            #Note that OM is configured for seasonal movement, but SS uses mean annual migration; hence the repetition for each OM season
            for(im in 1:.Object@nsubyears)
            {
              mov[is,ip,,im,from,to] <- moveVec
            }
          }
        }
      }

      .Object@mov[firstSimNum:lastSimNum,,,,,] <- mov

  #readline(prompt="Warning: remove hi movvement rate test: ")
  #.Object@mov[firstSimNum:lastSimNum,,,,,] <- mov*0+0.25 #uniform mixing test

      #Initial distribution of N -not really relevant when adopting SS initial N
      Idist <- karray(rep(0.0, times=(.Object@nsimPerOMFile[iom] * .Object@npop * .Object@nages * .Object@nareas)),
                      c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nareas))

      for (i in 1:30)
      {
        #at least nages/nsubyears
        for (mm in 1:.Object@nsubyears)
        {
          #mean unfished equilibrium distribution (extra age calculaitons are redundant)
          if (i == 1 & mm == 1)
          {
            Idist[,ip,1,] <- .Object@Recdist[keep(firstSimNum:lastSimNum),ip,]
          }

          Idist <- projection.domov(Ntemp=karray(Idist, dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nareas)),
                                    movtemp=karray(.Object@mov[firstSimNum:lastSimNum,,,mm,,], dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nareas, .Object@nareas)))

          for(ip in 1:.Object@npop)
          {
            Idist[,ip,2:(.Object@nages),] <- Idist[,ip,1:(.Object@nages - 1),]
            Idist[,ip,1,]                 <- as.double(.Object@Recdist[keep(firstSimNum:lastSimNum),ip,])
          }
        }
      }

      if (Report)
      {
        plot(Idist[1,1,,1], col=1, type='l', ylim=c(0, max(Idist[1,1,,])), xlab="age", ylab="Prop", main="Init N Dist by region")
        lines(Idist[1,1,,2], col=2)
        lines(Idist[1,1,,3], col=3)
        lines(Idist[1,1,,4], col=4)
      }

      .Object@Idist[firstSimNum:lastSimNum,,,] <- Idist
    }  # end if (nareas > 1)

    # R0 extracted from SS
    R0ss                                <- ssMod$parameters[ssMod$parameters$Label == "SR_LN(R0)",]$Value  # in thousands
    .Object@R0[firstSimNum:lastSimNum,] <- karray(rep(exp(R0ss), each=.Object@nsimPerOMFile[iom] * .Object@npop),
                                                  dim=c(.Object@nsimPerOMFile[iom], .Object@npop))

    #total numbers and spawning stock numbers
    SSN     <- karray(NA, c(.Object@npop, .Object@nages, allyears + 1, .Object@nsubyears, .Object@nareas))
    NBefore <- karray(NA, c(.Object@npop, .Object@nages, allyears + 1, .Object@nsubyears, .Object@nareas))

    # Catch (numbers, then mass in MSE framework)
    CssTmp <- karray(NA, c(.Object@npop, .Object@nages, allyears, .Object@nsubyears, .Object@nareas, .Object@nfleets))

    #YFT input is in mass, except LL fleets in numbers (3,7,10,11,13) (though not #25 fresh tuna LL apparently?)
    #check if catage in numbers or mass or mixed: appears to be in numbers:
    #sum(ssMod$catage[ssMod$catage$Fleet==1 & ssMod$catage$Yr==272, 11:39])  = 936.9    # YFT.dat = 13598 = 15kg per gillnet fish
    #sum(ssMod$catage[ssMod$catage$Fleet==3 & ssMod$catage$Yr==272, 11:39])  = 1.355    # YFT.dat = 1.355

    # SS Catch
    CssTmp[,,1:.Object@nyears,,,] <- 0.0

    # SS Catch-at-age
    CAss <- ssMod$catage[as.numeric(ssMod$catage$"Yr") >= OMd@firstSSYr,] # "Area"   "Fleet"  "Gender" "XX"     "XX"     "Morph"  "Yr"     "Seas"   "XX"     "Era"    "0"

    for (i in 1:nrow(CAss))
    {
      #convert SS season (defined as year) to MSE year, season
      yrSeas <- seasAsYrToMSEYrSeas.f(seasAsYr = CAss$Yr[i], endSeasAsYr = ssMod$endyr, numSeas = OMd@nsubyears, endYr = OMd@lastCalendarYr, endSeas = OMd@lastSeas, firstSeasAsYr = OMd@firstSSYr, firstSeas = OMd@firstSeas) #generic
      nColOffset <- 10

      #PAYMRF
      CssTmp[1,1:nagesss,yrSeas[1],yrSeas[2],CAss$Area[i],CAss$Fleet[i]] <- as.numeric(CAss[i,nColOffset+1:nagesss])
    }

    .Object@Css[keep(firstSimNum:lastSimNum),,,]    <- rep(apply(CssTmp, c(1,3,6), sum), each=.Object@nsimPerOMFile[iom])
    .Object@CAAFss[keep(firstSimNum:lastSimNum),,,] <- rep(apply(CssTmp, c(2,3,6), sum), each=.Object@nsimPerOMFile[iom])
    .Object@CBss[keep(firstSimNum:lastSimNum),,]    <- rep(apply(CssTmp * karray(.Object@Wt_age_mid[1,,,], c(.Object@npop,.Object@nages, allyears,.Object@nsubyears,.Object@nareas,.Object@nfleets)), c(1,3), sum), each=.Object@nsimPerOMFile[iom])

    #SS N(age)
    # YFT Nss <- ssMod$natage[as.numeric(ssMod$natage$"Yr")>=101 & ssMod$natage$"Beg/Mid"=="B",]  #in thousands

    Nss <- ssMod$natage[as.numeric(ssMod$natage$"Yr") >= OMd@firstSSYr & ssMod$natage$"Beg/Mid"=="B",]  #in thousands

    for (i in 1:nrow(Nss))
    {
      #convert SS season (defined as year) to MSE year, season
      yrSeas <- seasAsYrToMSEYrSeas.f(seasAsYr=Nss$Yr[i], endSeasAsYr = ssMod$endyr, numSeas = OMd@nsubyears, endYr = OMd@lastCalendarYr, endSeas = OMd@lastSeas, firstSeasAsYr = OMd@firstSSYr, firstSeas = OMd@firstSeas) #generic

      nColOffset <- 11

      NBefore[,1:nagesss,yrSeas[1],yrSeas[2],Nss$Area[i]] <- as.numeric(Nss[i,nColOffset + 1:nagesss])
      SSN[,,yrSeas[1],yrSeas[2],Nss$Area[i]] <- NBefore[,,yrSeas[1],yrSeas[2],Nss$Area[i]] * .Object@mat[firstSimNum,,,yrSeas[1]]
    }

    .Object@NBeforeInit[firstSimNum:lastSimNum,,,] <- rep(NBefore[,,.Object@nyears,1,], each=.Object@nsimPerOMFile[iom])

    #add noise to NBeforeInit
    CVMatBySPA  <- rep(.Object@NInitCV * exp(-.Object@NInitCVdecay * (0:(.Object@nages - 1))), each=.Object@nsimPerOMFile[iom] * .Object@npop)
    CVMatBySPAR <- karray(rep(CVMatBySPA, times=.Object@nareas), dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nareas))
    Ndevs       <- exp(CVMatBySPAR * rnorm(prod(dim(.Object@NBeforeInit[keep(firstSimNum:lastSimNum),,,]))) - 0.5 * CVMatBySPAR*CVMatBySPAR)

    .Object@NBeforeInit[firstSimNum:lastSimNum,,,] <- karray(.Object@NBeforeInit[keep(firstSimNum:lastSimNum),,,], dim=c(.Object@nsimPerOMFile[iom], .Object@npop, .Object@nages, .Object@nareas)) * Ndevs

    .Object@SSBAss[firstSimNum:lastSimNum,,]       <- rep(apply(karray(SSN[,,keep(1:allyears),.Object@nsubyears,], c(.Object@npop,.Object@nages,allyears,.Object@nareas)) * karray(.Object@Wt_age[firstSimNum,,,1:allyears], c(.Object@npop,.Object@nages,allyears,.Object@nareas)), c(1,3), sum), each=.Object@nsimPerOMFile[iom])

    # Bss mean over seasons
    #.Object@Bss[firstSimNum:lastSimNum,,]          <- rep(apply(NBefore[,,keep(1:allyears),,] * karray(.Object@Wt_age[firstSimNum,,,keep(1:allyears)], c(.Object@npop,.Object@nages,allyears,.Object@nsubyears,.Object@nareas)), c(1,3), sum), each=.Object@nsimPerOMFile[iom])
    BbyM                                           <- apply(NBefore[,,keep(1:allyears),,] * karray(.Object@Wt_age[firstSimNum,,,keep(1:allyears)], c(.Object@npop,.Object@nages,allyears,.Object@nsubyears,.Object@nareas)), c(1,3,4), sum)
    BbyY                                           <- apply(BbyM, c(1,2), mean)
    .Object@Bss[firstSimNum:lastSimNum,,]          <- rep(BbyY, each=.Object@nsimPerOMFile[iom])
    .Object@Recss[firstSimNum:lastSimNum,,]        <- rep(apply(NBefore[,1,keep(1:allyears),,], c(2), sum), each=.Object@nsimPerOMFile[iom])

    for (yyy in 1:allyears)
    {
      FirstIdx <- (yyy - 1) * .Object@nsubyears + 1
      LastIdx  <- yyy * .Object@nsubyears

      .Object@RecYrQtrss[firstSimNum:lastSimNum,,FirstIdx:LastIdx] <- rep(apply(NBefore[,1,yyy,keep(1:.Object@nsubyears),], c(1,2), sum), each=.Object@nsimPerOMFile[iom])
    }

    #All this stuff used to calculate Frep (summed over regions, mean over age and season)
    NsoRbyPAYM   <- karray(NA, dim=c(.Object@npop,.Object@nages,allyears+1,.Object@nsubyears))
    NsoRbyPAYM[] <- apply(NBefore, FUN=sum, MARGIN=c(1:4))

    ZsoRbyPAYM       <- karray(NA, dim=c(.Object@npop,.Object@nages-2,allyears,.Object@nsubyears))
    ZsoRbyPAYM[,,,1] <- -log(NsoRbyPAYM[,2:(.Object@nages - 1),1:allyears,2]       / NsoRbyPAYM[,1:(.Object@nages - 2),1:allyears,1])
    ZsoRbyPAYM[,,,2] <- -log(NsoRbyPAYM[,2:(.Object@nages - 1),1:allyears,3]       / NsoRbyPAYM[,1:(.Object@nages - 2),1:allyears,2])
    ZsoRbyPAYM[,,,3] <- -log(NsoRbyPAYM[,2:(.Object@nages - 1),1:allyears,4]       / NsoRbyPAYM[,1:(.Object@nages - 2),1:allyears,3])
    ZsoRbyPAYM[,,,4] <- -log(NsoRbyPAYM[,2:(.Object@nages - 1),2:(allyears + 1),1] / NsoRbyPAYM[,1:(.Object@nages - 2),1:allyears,4])

    FsoRbyPAYM <- ZsoRbyPAYM - karray(.Object@M[firstSimNum,,1:(.Object@nages-2),] / .Object@nsubyears, dim=c(.Object@npop,.Object@nages-2,allyears,.Object@nsubyears))

    # 2:27 = true ages 1:26 (1:26 used by SS)
    .Object@Frepss[firstSimNum:lastSimNum,] <- rep(apply(FsoRbyPAYM[,2:27,,], mean, MARGIN=c(3)), each=.Object@nsimPerOMFile[iom])

    rm(NsoRbyPAYM)
    rm(ZsoRbyPAYM)

    #includes two abundance index karrays:
    #Iobs is LL selected annual numbers aggregated over seasons and regions
    #CPUEobs is partitioned by season and fishery (one fishery per region in initial YFT implementation)
    NLLssByAYMR <- apply(NBefore[,keep(1:.Object@nages),keep(1:allyears),,], MARGIN=c(2,3,4,5), FUN=sum,na.rm=T)

    # Longline-selected Numbers by area and season for partitioned abundance index
    for(isubyears in 1:.Object@nsubyears)
    {
      for(iCPUE in 1:.Object@nCPUE)
      {
        .Object@NLLss[firstSimNum:lastSimNum,,isubyears,iCPUE] <- rep(apply(NLLssByAYMR[,,isubyears,.Object@CPUEFleetAreas[iCPUE]] * karray(rep(.Object@CPUEsel[firstSimNum,,iCPUE], times=allyears), dim=c(.Object@nages,allyears)), MARGIN=c(2), FUN=sum), each=.Object@nsimPerOMFile[iom])
      }
    }

    .Object@NLLIss[firstSimNum:lastSimNum,] <- apply(.Object@NLLss[keep(firstSimNum:lastSimNum),,,], MARGIN=c(1,2), FUN=sum)

    rm(NLLssByAYMR)

    # some of these should be retained from SS, some could be recalculated with popdyn, but should check consistency
    SSBss <- (ssMod$recruit$spawn_bio)  #ss yr=13-272 in tonnes
    Recss <- (ssMod$recruit$exp_recr)   # in thousands

    #SS F/FMSY - xxx need to resolve inconsistencies if series to be usefully merged
    for (ymse in 1:.Object@nyears)
    {
      .Object@F_FMSYss[firstSimNum:lastSimNum,ymse] <- rep(mean(ssMod$Kobe$F.Fmsy[seasAsYrToMSEYr.f(seasAsYr=ssMod$Kobe$Year, endSeasAsYr = ssMod$endyr, numSeas = OMd@nsubyears, endYr = OMd@lastCalendarYr, endSeas = OMd@lastSeas, firstSeasAsYr = OMd@firstSSYr, firstSeas = OMd@firstSeas)==ymse]),.Object@nsimPerOMFile[iom])
    }

    #MSE works in individuals and kg (or 1000 individuals and tonnes)
    #SS MSY, SSBMSY, etc (no obvious BMSY in SS outputs)
    MSYss           <- 4.*ssMod$derived_quants$Value[ssMod$derived_quants$LABEL == "TotYield_MSY"]  #SS MSY is quarterly
    SSBMSYss        <- ssMod$derived_quants$Value[ssMod$derived_quants$LABEL == "SSB_MSY"]
    .Object@SSB0ss  <- c(.Object@SSB0ss, rep(ssMod$derived_quants$Value[ssMod$derived_quants$LABEL == "SSB_Unfished"], .Object@nsimPerOMFile[iom]))
    .Object@B0ss    <- c(.Object@B0ss, rep(ssMod$derived_quants$Value[ssMod$derived_quants$LABEL == "TotBio_Unfished"], .Object@nsimPerOMFile[iom]))
    .Object@FMSYss  <- c(.Object@FMSYss, rep(ssMod$derived_quants$Value[ssMod$derived_quants$LABEL == "Fstd_MSY"], .Object@nsimPerOMFile[iom]))

    # ABT MSE assumes TAC split envenly among quarters
    # YFT q defined as 1 (or zero for fisheries that are time blocks that do not operate in 2014)
    # YFT assumes a default annual distribution for Effort with a constant proportion by season, fishery and area strata
    # This is true for aggregate annual TAC or TAE-managed fisheries, but not if there is a mix (or fishery-specific controls)
    q <- karray(1.0, dim=c(.Object@nsimPerOMFile[iom], .Object@nfleets))

    #YFT Effort and Catch-in-Mass distribution (by season region and fishery for one year)  for future effort/catch allocations
    EByQtrLastYr  <- karray(0.0, c(.Object@nsubyears, .Object@nareas, .Object@nfleets))
    ECurrent      <- karray(0.0, c(.Object@nsubyears, .Object@nareas, .Object@nfleets))
    CMCurrent     <- karray(0.0, c(.Object@nsubyears, .Object@nareas, .Object@nfleets))

    # Note that YFT fisheries are always defined on the basis of a specific region, such that
    # fishery and area distinctions are redundant for many purposes. The F definitions below
    # appear to be scaled relative to the season duration, such that:
    # if duration is 3 months, ssMod$seasDuration = 0.25
    # Fs reported as "F:_fleetnumber" and plotted are annuallized, but only applied for seasDuration
    # meaning that the F value in absolute practical terms is F * seasduration!
    # Ms and growth appear to be re-scaled separate to this, and the seasduration introduces an
    # error in terms of mid-year sizes:
    # i.e. midyearSize is t=t+0.125 rather than t+0.5

    lastSSyrSeas   <- ssMod$endyr
    numRecentSeas   <- .Object@recentPerLast - .Object@recentPerFirst + 1 # last timestep indicated by recentPerFirst=0

    # calculate mean catch and effort by season and fleet
    if ((.Object@seasonCEDist) & (numRecentSeas %% OMd@nsubyears ==0))
    {
      recentPeriodYrs    <- numRecentSeas/.Object@nsubyears

      for (ifleets in 1:.Object@nfleets)
      {
        for (iy in 1:recentPeriodYrs)
        {
          for (im in 1:.Object@nsubyears)
          {
            r <- ssMod$"fleet_area"[ifleets]

            # one row for each area, but fleet only operates in one area so can sum over areas
            YrRows                      <- lastSSyrSeas - .Object@recentPerFirst - (recentPeriodYrs * OMd@nsubyears) + ((iy - 1) * OMd@nsubyears) + im
            EByQtrLastYr[im,r,ifleets]  <- sum(ssMod$timeseries[ssMod$timeseries$Yr == YrRows, "F:_" %&% ifleets == names(ssMod$timeseries)])
            ECurrent[im,r,ifleets]      <- ECurrent[im,r,ifleets] + (1.0 / recentPeriodYrs) * EByQtrLastYr[im,r,ifleets]

            yrSeas <- seasAsYrToMSEYrSeas.f(seasAsYr      = YrRows,
                                            endSeasAsYr   = ssMod$endyr,
                                            numSeas       = OMd@nsubyears,
                                            endYr         = OMd@lastCalendarYr,
                                            endSeas       = OMd@lastSeas,
                                            firstSeasAsYr = OMd@firstSSYr,
                                            firstSeas     = OMd@firstSeas) #generic

            CMCurrent[im,r,ifleets] <- CMCurrent[im,r,ifleets] + sum((1.0 / recentPeriodYrs) *
                                                                 colSums(CssTmp[,,yrSeas[1], yrSeas[2], r, ifleets]) *
                                                                 colSums(.Object@Wt_age_mid[1,,,seasAsYrToMSEYrSeas.f(YrRows)[1]]))
          }
        }
      }
    }
    else #calculate mean catch and effort that remains constant over seasons
    {
      if ((.Object@seasonCEDist) & (numRecentSeas %% OMd@nsubyears !=0))
      {
        readline("Warning: recent C&E dists are not seasonal because recentPeriod is not a multiple of nsubyears. Press ENTER to continue")
      }

      for (ifleets in 1:.Object@nfleets)
      {
        for (iseason in .Object@recentPerLast:.Object@recentPerFirst)
        {
          r <- ssMod$"fleet_area"[ifleets]
          # one row for each area, but fleet only operates in one area so can sum over areas
          YrRows <- lastSSyrSeas - iseason
          EThisSeason <- ssMod$timeseries[ssMod$timeseries$Yr == YrRows, "F:_" %&% ifleets == names(ssMod$timeseries)]

          ECurrent[1:OMd@nsubyears,r,ifleets] <- ECurrent[1:OMd@nsubyears,r,ifleets] + (1.0 /numRecentSeas) * EThisSeason

          yrSeas <- seasAsYrToMSEYrSeas.f(seasAsYr      = YrRows,
                                          endSeasAsYr   = ssMod$endyr,
                                          numSeas       = OMd@nsubyears,
                                          endYr         = OMd@lastCalendarYr,
                                          endSeas       = OMd@lastSeas,
                                          firstSeasAsYr = OMd@firstSSYr,
                                          firstSeas     = OMd@firstSeas) #generic

          CMCurrent[1:OMd@nsubyears,r,ifleets] <- CMCurrent[1:OMd@nsubyears,r,ifleets] + sum(1.0 / (numRecentSeas) *
                                                                                         colSums(CssTmp[,,yrSeas[1], yrSeas[2], r, ifleets]) *
                                                                                         colSums(.Object@Wt_age_mid[1,,,seasAsYrToMSEYrSeas.f(YrRows)[1]]))
        }
      }
    }

    rm(CssTmp)
    rm(SSN)
    rm(NBefore)

    EByQtrLastYr  <- EByQtrLastYr * ssMod$seasdurations      #see note aboive
    ECurrent      <- ECurrent * ssMod$seasdurations

    #plot the F by fleet time series to see if the SS output is sensible
    par(mfrow=c(1,1), ask=F)

    for (ifleets in c(1:.Object@nfleets))
    {
      f <- "F:_" %&% ifleets
      y <- "Yr"
      d <- ssMod$timeseries[, c(f,y) ]
      d <- d[d[,1] > 0,] * ssMod$seasdurations

      if (Report)
      {
        if (ifleets == 1)
        {
          plot( d[,2],d[,1], type='l', ylim=c(0,2.5))
        }
        else
        {
          lines( d[,2],d[,1],col=ifleets)
        }
      }
    }

    .Object@q[firstSimNum:lastSimNum,] <- q

    for (isim in firstSimNum:lastSimNum)
    {
      .Object@ECurrent[isim,,,]     <- ECurrent
      .Object@CMCurrent[isim,,,]    <- CMCurrent
      .Object@EByQtrLastYr[isim,,,] <- EByQtrLastYr
    }

    MSY_projection_sims[iom,] <- c(firstSimNum, lastSimNum)

    .Object@MSYss    <- c(.Object@MSYss,   rep(MSYss, .Object@nsimPerOMFile[iom]))
    .Object@SSBMSYss <- c(.Object@SSBMSYss,rep(SSBMSYss, .Object@nsimPerOMFile[iom]))

  } #iom loop

  print("Doing MSY projections")

  # do MSYref projections in parallel
  if (.Object@CppMethod != 0)
  {
    MSYrefs <- Cpp_MSYrefs(MSY_projection_sims[,1], .Object, nyears=70)
  }
  else
  {
    MSYrefs <- R_MSYrefs(MSY_projection_sims[,1], .Object@UseCluster != 0, .Object, nyears=70)
  }

  for(iom in 1:nOMfiles)
  {
    indices <- seq(MSY_projection_sims[iom,1], MSY_projection_sims[iom,2])

    if (.Object@UseMSYss==0)
    {
      .Object@MSY[indices]         <- MSYrefs[1,iom]
      .Object@BMSY[indices]        <- MSYrefs[2,iom]
      .Object@VBMSY[indices]       <- MSYrefs[3,iom]
      .Object@SSBMSY[indices]      <- MSYrefs[4,iom]
      .Object@UMSY[indices]        <- MSYrefs[5,iom]
      .Object@FMSY1[indices]       <- MSYrefs[6,iom]
      .Object@SSBMSY_SSB0[indices] <- MSYrefs[7,iom]
      .Object@SSB0[indices]        <- MSYrefs[8,iom]
      .Object@B0[indices]          <- MSYrefs[9,iom]

    }
    else
    {
      .Object@MSY[indices]         <- .Object@MSYss[iom]
      .Object@BMSY[indices]        <- NA #.Object@BMSYss[iom]
      .Object@VBMSY[indices]       <- NA
      .Object@SSBMSY[indices]      <- .Object@SSBMSYss[iom]
      .Object@UMSY[indices]        <- NA
      .Object@FMSY1[indices]       <- .Object@FMSYss[iom]
      .Object@SSBMSY_SSB0[indices] <- .Object@SSBMSYss[iom]/.Object@SSB0ss[iom]
      .Object@SSB0[indices]        <- .Object@SSB0ss[iom]
      .Object@B0[indices]          <- .Object@B0ss[iom]
    }
  }

  # selectvity temporal variability
  .Object@selTSSign[]                       <- round(runif(prod(dim(.Object@selTSSign))))
  .Object@selTSSign[.Object@selTSSign < 1]  <- -1   # -1 or 1 indicate initial sin wave trend
  .Object@selTSWaveLen[]                    <- runif(prod(dim(.Object@selTSWaveLen))) * (.Object@selWLRange[2] - .Object@selWLRange[1]) + .Object@selWLRange[1]

  print("End OMss definition")

  gc()

  return(.Object)
})
#### End OMss Class




#------------------------------------Class MSE----------------------------------------------------------------------
#modifification of Class MSE that takes operating model of class OMss as input

setClass("MSE",representation(

               # Description
               OMName           = "character",
               Label            = "character",
               Date             = "character",
               Author           = "character",
               Notes            = "character",
               PrimarySource    = "character",
               seed             = "integer",
               CppMethod        = "numeric",
               UseCluster       = "numeric",
               OMid             = "karray",
               SBlim            = "numeric",
               Flim             = "numeric",
               TACTime          = "numeric",
               rULim            = "numeric",

                # Dimensions
               nsim             = "integer",
               npop             = "integer",
               nages            = "integer",

               # MSE dimensions
               nyears           = "integer",
               nsubyears        = "integer",
               nareas           = "integer",
               nfleets          = "integer",
               proyears         = "integer", # Proyears
               nMPs             = "integer", # number of management procedures
               targpop          = "integer",

               # tuning vector for MP tuning
               tune             = "numeric",

               # projection period for tuning diagnostic
               tunePMProjPeriod = "numeric",

               # Observation model
               Cimp             = "numeric",
               Cb               = "numeric",
               Cerr             = "karray",
               Iimp             = "numeric",
               Ibeta            = "numeric",
               Ierr             = "karray",
               IAC              = "karray",
               ITrend           = "karray",
               CPUEerr          = "karray",
               nCAAobs          = "numeric",
               nCALobs          = "numeric",
               Lcv              = "numeric",
               Mb               = "numeric",
               Kb               = "numeric",
               t0b              = "numeric",
               Linfb            = "numeric",
               LFCb             = "numeric",
               LFSb             = "numeric",
               FMSYb            = "numeric",
               FMSY_Mb          = "numeric",
               BMSY_B0b         = "numeric",
               ageMb            = "numeric",
               Dimp             = "numeric",
               Db               = "numeric",
               Derr             = "karray",
               Btimp            = "numeric",
               Btb              = "numeric",
               Bterr            = "karray",
               Ftimp            = "numeric",
               Ftb              = "numeric",
               Fterr            = "karray",

               hb               = "numeric",
               Recbcv           = "numeric",
               IMSYb            = "numeric",
               MSYb             = "numeric",
               BMSYb            = "numeric",

               # Management quantities
               CM               = "karray",
               CMbyF            = "karray",
               C_MSY            = "karray",
               B_B0             = "karray",
               SSB_SSB0         = "karray",
               B_BMSY           = "karray",
               SSB_SSBMSY       = "karray",
               F_FMSY           = "karray",
               B                = "karray",
               SSB              = "karray",
               TAC              = "karray",
               MPs              = "character",

               #Calendar time labels for time series
               yrLabels         = "karray",
               seasonLabels     = "karray",
               yrSeasLabels     = "karray",
               firstMPYr        = "numeric",
               lastCalendarYr   = "numeric",

               catchBridge      = "karray",
               catchBridgeCV    = "numeric",

               MPDataLag        = "numeric",
               Rec              = "karray",
               RecYrQtr         = "karray",
               Css              = "karray",
               CBss             = "karray",
               CAAFss           = "karray",
               SSBAss           = "karray",
               Bss              = "karray",
               Recss            = "karray",
               RecYrQtrss       = "karray",
               NLLIss           = "karray",
               F_FMSYss         = "karray", #import from SS
               Frepss           = "karray",
               CPUEobsR         = "karray", #CPUE indices by R
               IobsArchive      = "karray", #CPUE indices
               IobsRArchive     = "karray"  #CPUE indices by R
 ))



runProjection <- function(.Object, OM, projSims, CPUEobsR, TACEErrorAll, MPs, interval, Report, CppMethod, UseCluster, EffortCeiling, TACTime, rULim, tune)
{
  nsim      <- length(projSims)
  npop      <- OM@npop
  nyears    <- OM@nyears
  proyears  <- OM@proyears
  nages     <- OM@nages
  nsubyears <- OM@nsubyears
  nareas    <- OM@nareas
  nfleets   <- OM@nfleets
  nCPUE     <- OM@nCPUE
  targpop   <- as.integer(OM@targpop)
  allyears  <- nyears + proyears
  nMPs      <- length(MPs)

  NMass     <- karray(as.double(NA),c(nsim,npop,nages,allyears+1,nsubyears,nareas))  # N in mass
  Z         <- karray(as.double(NA),c(nsim,npop,nages,allyears+1,nsubyears,nareas))  # Z aggregated over fleets
  SSBA      <- karray(as.double(NA),c(nsim,npop,allyears))                           # SSB aggregated over ages archive
  NLLI      <- karray(as.double(NA),c(nsim,allyears))                                # Longline-selected Numbers over all ages, areas, populations for aggregate abundance index
  NLL       <- karray(as.double(NA),c(nsim,allyears,nsubyears,nareas))               # Longline-selected Numbers over all ages, areas, populations for aggregate abundance index
  Fdist     <- karray(as.double(NA),c(nsim,npop,nfleets,nareas))                     # current F dist by fleet
  CAAF      <- karray(as.double(NA),c(nsim,nages,allyears,nfleets))                  # Catch at age by fleets
  B         <- karray(as.double(NA),c(nsim,npop,allyears))                           # Biomass annual aggregate series
  Rec       <- karray(as.double(NA),c(nsim,npop,allyears))                           # Recruitment series
  RecYrQtr  <- karray(as.double(NA),c(nsim,npop,allyears * nsubyears))               # Quaterly series

  initYear <- nyears  # last assessment year

  # subset all model parameters needed by the model run
  Len_age       <- OM@Len_age[keep(projSims),,,]
  Len_age_mid   <- OM@Len_age_mid[keep(projSims),,,]
  Wt_age        <- OM@Wt_age[keep(projSims),,,]
  #Wt_age_SB     <- OM@Wt_age_SB[keep(projSims),,,]
  Wt_age_mid    <- OM@Wt_age_mid[keep(projSims),,,]
  h             <- OM@h[keep(projSims),]
  R0            <- OM@R0[keep(projSims),]
  Idist         <- OM@Idist[keep(projSims),,,]
  M             <- OM@M[keep(projSims),,,]
  mat           <- OM@mat[keep(projSims),,,]
  Recdevs       <- OM@Recdevs[keep(projSims),,]
  Recdist       <- OM@Recdist[keep(projSims),,]
#    Linf          <- OM@Linf[keep(projSims)]
#    K             <- OM@K[keep(projSims)]
  q             <- OM@q[keep(projSims),]
  EByQtrLastYr  <- OM@EByQtrLastYr[keep(projSims),,,]
  ECurrent      <- OM@ECurrent[keep(projSims),,,]
  CMCurrent     <- OM@CMCurrent[keep(projSims),,,]
  sel           <- OM@sel[keep(projSims),,]
  mov           <- OM@mov[keep(projSims),,,,,]
  CPUEsel       <- OM@CPUEsel[keep(projSims),,]
  CPUEobsY      <- OM@CPUEobsY[keep(projSims),]
  TACEError     <- TACEErrorAll[keep(projSims),]

  CAAF[,,1:initYear,]                 <- OM@CAAFss[keep(projSims),,1:initYear,]
  SSBA[,,1:initYear]                  <- OM@SSBAss[keep(projSims),,1:initYear]
  B[,,1:initYear]                     <- OM@Bss[keep(projSims),,1:initYear]
  Rec[,,1:initYear]                   <- OM@Recss[keep(projSims),,1:initYear]
  RecYrQtr[,,1:(initYear*nsubyears)]  <- OM@RecYrQtrss[keep(projSims),,1:(initYear*nsubyears)]
  NLLI[,1:initYear]                   <- OM@NLLIss[keep(projSims),1:initYear]
  NLL[,1:initYear,,]                  <- OM@NLLss[keep(projSims),1:initYear,,]

  # q to keep CPUE on original scale; should move to OM initialization
  qCPUE <- rep(NA, nsim)

  for(isim in 1:nsim)
  {
    qCPUE[isim] <- sum(CPUEobsY[isim,1:nyears][!is.na(CPUEobsY[isim,1:nyears])]) /
                   sum(NLLI[isim,1:nyears][!is.na(CPUEobsY[isim,1:nyears])])
  }

  SSB0 <- karray(OM@SSB0[keep(projSims)],dim=c(nsim,npop))

  # Calculate spawning stock biomass per recruit
  SSBpR   <- SSB0 / R0

  # Ricker SR params
  bR      <- log(5.0 * h) / (0.8 * SSB0)
  aR      <- exp(bR * SSB0) / SSBpR

  # historical simulation of one year to re-create N for first projection year
  # ---------------------------------------------------------------------------
  cat("Re-running final year of OM")
  cat("\n")

  y           <- initYear
  mm          <- nsubyears
  nSpawnPerYr <- length(OM@Recsubyr)  # rec Index for more than 1 rec per year

  NullRecSpatialDevs <- karray(as.double(1.0),c(nsim,npop,nareas))
  N_Y                <- karray(NA,c(nsim,npop,nages,nsubyears + 1,nareas))
  NBefore_Y          <- karray(NA,c(nsim,npop,nages,nsubyears + 1,nareas))
  SSN_Y              <- karray(as.double(NA),c(nsim,npop,nages,nsubyears,nareas))
  NLL_Y              <- karray(as.double(NA),c(nsim,nsubyears,OM@nCPUE))
  NLLI_Y             <- karray(as.double(NA),c(nsim))
  C_Y                <- karray(as.double(NA),c(nsim,npop,nages,nsubyears,nareas,nfleets))
  SSBA_Y             <- karray(as.double(NA),c(nsim,npop))
  M_Y                <- M[,,,y]
  mat_Y              <- mat[,,,y]
  Len_age_Y          <- Len_age[,,,y]
  Len_age_mid_Y      <- Len_age_mid[,,,y]
  Wt_age_Y           <- Wt_age[,,,y]
  #Wt_age_SB_Y        <- Wt_age_SB[,,,y]
  Wt_age_mid_Y       <- Wt_age_mid[,,,y]
  RecdevInd_Y        <- (y - 1) * nSpawnPerYr + 1
  Recdevs_Y          <- Recdevs[,,keep(RecdevInd_Y:(RecdevInd_Y + nSpawnPerYr - 1))]
  #All this stuff used to calculate Frep (summed over regions, mean over age and season)
  NsoRbySPAMInit     <- karray(NA, dim=c(nsim,npop,nages,nsubyears + 1))

  # Initialise starting population
  N_Y[,,,1,]       <- OM@NBeforeInit[keep(projSims),,,]
  NBefore_Y[,,,1,] <- OM@NBeforeInit[keep(projSims),,,]

  # Matrix index arrays used in historic and projection R code
  SPAYMRF <- as.matrix(expand.grid(1:nsim, 1:npop, 1:nages, initYear, 1, 1:nareas, 1:nfleets))
  SPAMRF  <- SPAYMRF[,c(1,2,3,5,6,7)]
  SMRF    <- SPAYMRF[,c(1,5,6,7)]
  SPARF   <- SPAYMRF[,c(1,2,3,6,7)]
  SFA     <- SPAYMRF[,c(1,7,3)]
  SF      <- SPAYMRF[,c(1,7)]
  SARF    <- SPAYMRF[,c(1,3,6,7)]
  SAR     <- SPAYMRF[,c(1,3,6)]
  SRF     <- SPAYMRF[,c(1,6,7)]

  SPAYMR  <- as.matrix(expand.grid(1:nsim, 1:npop, 1:nages, initYear, 1, 1:nareas))
  SPAMR   <- SPAYMR[,c(1,2,3,5,6)]
  SPAY    <- SPAYMR[,c(1,2,3,4)]
  SPAR    <- SPAYMR[,c(1,2,3,6)]
  SA      <- SPAYMR[,c(1,3)]

  if (.Object@CppMethod != 0) # use C++ Baranov sub-routine
  {
    Obj <- Om.create(nsim,
                     npop,
                     nages,
                     nsubyears,
                     nareas,
                     nfleets,
                     OM@Recsubyr)

    OmB.nt.initialiseParameters(Obj,
                                M_Y,
                                R0,
                                mat_Y,
                                Idist,
                                Wt_age_Y,
                                h)

    OmB.nt.runHistoric(Obj,
                       as.double(1.0),
                       q,
                       EByQtrLastYr,
                       R0,
                       M_Y,
                       mat_Y,
                       Idist,
                       Len_age_Y,
                       Wt_age_Y,
                       sel,
                       mov,
                       h,
                       Recdist,
                       Recdevs_Y,
                       NullRecSpatialDevs,
                       OM@SRrel,
                       N_Y,
                       NBefore_Y,
                       SSN_Y,
                       C_Y,
                       SSBA_Y)
  }
  else # Use R-based projection code
  {
    FM <- karray(NA, dim=c(nsim,npop,nages,nareas,nfleets))
    Z  <- karray(NA, dim=c(nsim,npop,nages,nareas))

    for (mm in 1:nsubyears)
    {
      SPAMRF[,4] <- mm
      SPAMR[,4]  <- mm
      SMRF[,2]   <- mm

      SSN_Y[,,,mm,] <- NBefore_Y[,,,mm,] * karray(rep(mat[,,,y],times=nareas), c(nsim,npop,nages,nareas))
      #potential change to _SB
      SSBA_Y        <- apply(SSN_Y[,,,mm,] * karray(Wt_age[,,,y], dim=c(nsim,npop,nages,nareas)), c(1,2), sum, na.rm=T)
      #SSBA_Y        <- apply(NBefore_Y[,,,mm,] * karray(Wt_age_SB[,,,y], dim=c(nsim,npop,nages,nareas)), c(1,2), sum, na.rm=T)

      if (mm %in% OM@Recsubyr)
      {
        for(pp in 1:npop)
        {
          # ie every qtr for YFT
          RecdevInd <- (y - 1) * nSpawnPerYr + mm

          # recruit fish
          if (OM@SRrel[pp] == 1)
          {
            # Beverton-Holt recruitment
            rec <- Recdevs[,pp,RecdevInd] * ((0.8 * R0[,pp] * h[,pp] * SSBA_Y[,pp]) / (0.2 * SSBpR[,pp] * R0[,pp] * (1 - h[,pp]) + (h[,pp] - 0.2) * SSBA_Y[,pp]))
          }
          else
          {
            # Most transparent form of the Ricker uses alpha and beta params
            rec <- Recdevs[,pp,RecdevInd] * aR[,pp] * SSBA_Y[,pp] * exp(-bR[,pp] * SSBA_Y[,pp])
          }

          NBefore_Y[,pp,1,mm,] <- rec * Recdist[,pp,]
        }
      }

      # move fish
      if (nareas > 1)
      {
        N_Y[,,,mm,] <- projection.domov(Ntemp=karray(NBefore_Y[,,,mm,], dim=c(nsim,npop,nages,nareas)),
                                        movtemp=karray(mov[,,,mm,,], dim=c(nsim,npop,nages,nareas,nareas)))
      }
      else
      {
        N_Y[,,,mm,] <- NBefore_Y[,,,mm,]
      }

      # Apply M and F
      FM[SPARF] <- EByQtrLastYr[SMRF] * sel[SFA] * q[SF]
      Ftot      <- apply(FM, c(1,2,3,4), sum)
      Z[SPAR]   <- Ftot[SPAR] + M[SPAY] / nsubyears

      C_Y[SPAMRF] <- N_Y[SPAMR] * (1.0 - exp(-Z[SPAR])) * (FM[SPARF] / Z[SPAR])

      N_Y[,,,mm,] <- N_Y[,,,mm,] * exp(-Z[,,,])

      # Age fish
      NBefore_Y[,,nages,mm + 1,]          <- N_Y[,,nages - 1,mm,] + N_Y[,,nages,mm,]
      NBefore_Y[,,2:(nages - 1),mm + 1,]  <- N_Y[,,1:(nages - 2),mm,]
      NBefore_Y[,,1,mm + 1,]              <- 0
      N_Y[,,,mm + 1,]                     <- NBefore_Y[,,,mm + 1,]
    }
  }

  # calculate LL selected numbers for the year for the abundance index
  NLLbySAMR <- apply(N_Y[,,keep(1:nages),1:nsubyears,], MARGIN=c(1,3,4,5), FUN=sum, na.rm=T)

  for (isubyears in 1:nsubyears)
  {
    for (iCPUE in 1:OM@nCPUE)
    {
      NLL_Y[1:nsim,isubyears,iCPUE] <- apply(NLLbySAMR[,,isubyears,OM@CPUEFleetAreas[iCPUE]] * CPUEsel[,,iCPUE], FUN=sum, MARGIN=1, na.rm=T)
    }
  }

  NLLI_Y <- apply(NLL_Y, FUN=sum, MARGIN=1, na.rm=T)

  # copy results back into historic data
  SSBA[,,y] <- SSBA_Y
  NLLI[,y]  <- NLLI_Y
  NLL[,y,,] <- NLL_Y

  E <- apply(ECurrent, sum, MARGIN=c(1,2,4))

  #first N needs to be re-initialized for each MP
  NInit_Y <- N_Y[,,,nsubyears + 1,]

  # Store results ...
  .Object@CM[1:nMPs,projSims,,y]     <- rep(apply(C_Y * karray(Wt_age_mid[,,,y], c(nsim,npop,nages,nsubyears,nareas,nfleets)), c(1,2), sum),each=nMPs)
  .Object@CMbyF[1:nMPs,projSims,,y,] <- rep(apply(C_Y * karray(Wt_age_mid[,,,y], c(nsim,npop,nages,nsubyears,nareas,nfleets)), c(1,2,6), sum),each=nMPs)

  CAAF[,,y,]  <- apply(C_Y, c(1,3,6), sum)
  #B[,,y]      <- apply(NBefore_Y[,,,keep(1:nsubyears),] * karray(Wt_age[,,,y], c(nsim,npop,nages,nsubyears,nareas)), c(1,2), sum)
  BbyS        <- apply(NBefore_Y[,,,keep(1:nsubyears),] * karray(Wt_age[,,,y], c(nsim,npop,nages,nsubyears,nareas)), c(1,2,4), sum)
  B[,,y]      <- apply(BbyS, c(1,2), mean)
  Rec[,,y]    <- apply(NBefore_Y[,,1,keep(1:nsubyears),], c(1,2), sum)

  FirstIdx <- (y - 1) * nsubyears + 1
  LastIdx  <- y * nsubyears

  RecYrQtr[,,FirstIdx:LastIdx] <- apply(NBefore_Y[,,1,keep(1:nsubyears),], c(1,2,3), sum)

  # Save data for Frep calculation
  NsoRbySPAMInit[,,,1:nsubyears] <- apply(NBefore_Y[,,,keep(1:nsubyears),], FUN=sum, MARGIN=c(1:4))

  # Run projections
  # ---------------------------------------------------------------------------

  cat("\n Running projections \n")

  firstMPy <- nyears + .Object@firstMPYr - .Object@lastCalendarYr
  upyrs    <- firstMPy + (0:(floor(OM@proyears / interval))) * interval  # the years in which there are MP updates (e.g. every three years)

  #xxx zzz need to load historical Catch and CL (depending on MP), bypassed for now

  CAAInit     <- karray(NA, dim=c(nsim, nages, nyears))
  CAAIndInit  <- karray(NA, dim=c(nsim, nages, nyears)) #Catch-at-age for the abudnance index fleets

  if (min(apply(CAAF, c(1,2,3), sum) < 0, na.rm=T))
  {
    print("CA<0 a")
    browser()
  }

  CAAInit[,,1:nyears]    <- sampCatch(apply(CAAF[,,keep(1:nyears),], c(1,2,3),sum), .Object@nCAAobs)                    # need time series for MP
  CAAIndInit[,,1:nyears] <- sampCatch(apply(CAAF[,,keep(1:nyears),OM@indexFisheries], c(1,2,3), sum), .Object@nCAAobs) # need time series for MP

  # xxx zzz need to load historical from OM eventually
  nCALbins <- 30

  #DK: quick and dirty fix
  #CAL_bins   <- seq(0, max(Linf), length.out=nCALbins)
  CAL_bins   <- seq(0, 200, length.out=nCALbins - 1)
  CAL_bins   <- c(CAL_bins, CAL_bins[nCALbins - 1] * 2)     #add big +bin
  CALInit    <- karray(NA, dim=c(nsim, nCALbins, nyears))
  CALIndInit <- karray(NA, dim=c(nsim, nCALbins, nyears))

  #DK change to use Len_age rather than recalc
  #CAL                     <- makeCAL(CAA, Linf=Linf[,1,1:nyears], K=K[,1,1:nyears], t0=t0[1], CAL_bins)
  #CAL[,,1:nyears]         <- makeCAL(CAA[,,1:nyears,drop=F], Len_age[1,1,,1:nyears], CAL_bins)
  CALInit[,,1:nyears]     <- makeCAL(CAAInit[,,1:nyears], Len_age_mid[1,1,,1:nyears], CAL_bins)
  #CALLL10[,,1:nyears]     <- makeCAL(CAALL10[,,1:nyears,drop=F], Len_age[1,1,,1:nyears], CAL_bins)
  CALIndInit[,,1:nyears]  <- makeCAL(CAAIndInit[,,1:nyears], Len_age_mid[1,1,,1:nyears], CAL_bins)

  for (MP in 1:nMPs)
  {
    #MP loop
    #Need to re-initialize some arrays
    if (.Object@CppMethod != 0)
    {
      N_Y[,,,nsubyears + 1,]        <- NInit_Y
      NBefore_Y[,,,nsubyears + 1,]  <- NInit_Y
    }
    else
    {
      N_Y[,,,1,]        <- NInit_Y
      NBefore_Y[,,,1,]  <- NInit_Y
    }

    NsoRbySPAM <- NsoRbySPAMInit

    .Object@CM[MP,keep(projSims),,]     <- OM@CBss[keep(projSims),,]
    .Object@F_FMSY[MP,keep(projSims),] <- OM@Frepss[keep(projSims),] / OM@FMSY1[keep(projSims)]

    CAA     <- CAAInit
    CAAInd  <- CAAIndInit
    CAL     <- CALInit
    CALInd  <- CALIndInit

    set.seed(OM@seed) #repeat the same stochastic history for each MP to keep the comparison fair

    cat(paste(MP,"/",nMPs," Running MSE for ",MPs[MP],"    lastSim ", projSims[length(projSims)], sep=""))  # print a progress report
    cat("\n")
    flush.console()                                             # update the console

    # initial TAC and TAE values required for first MP application
    TAC   <- karray(apply(CMCurrent, sum, MARGIN=c(1)))      # refresh the MP store of TAC among simulations
    TAE   <- karray(rep(0, nfleets * nsim), dim=c(nsim, nfleets))
    TACE  <- cbind(TAE, TAC)

    if (.Object@CppMethod != 0)
    {
      Om.nt.beginProjection(Obj, as.double(rep(log(0.001), nfleets)))
    }

    for (y in (nyears + 1):(nyears + proyears))
    {
      SPAYMRF[,4] <- y
      SPAY[,4]    <- y

      #projection years loop includes repeat of first year for rec and initializtion
      if (.Object@CppMethod == 0)
      {
        cat(".")
      }

      selTS <- sel

      if (OM@selAgeRange >= 1)
      {
        #add an age shift to the selectivity noise
        selAge <- karray(rep(round(OM@selTSSign[keep(1:nsim),,2] * sin((y - nyears - 1) * OM@selTSWaveLen[keep(1:nsim),,2]) * OM@selAgeRange), nages), dim=dim(sel)[1:2])

        for (si in 1:nsim)
        {
          for (fi in 1:nfleets)
          {
            if (selAge[si,fi] > 0)
            {
              selTS[si,fi,(selAge[si,fi] + 1):nages]  <- sel[si,fi,1:(nages - selAge[si,fi])]
              selTS[si,fi,1:selAge[si,fi]]            <- sel[si,fi,1]
            }
            else
            {
              if (selAge[si,fi] < 0)
              {
                selTS[si,fi,1:(nages + selAge[si,fi])]          <- sel[si,fi,(1 - selAge[si,fi]):nages]
                selTS[si,fi,(nages + selAge[si,fi] + 1):nages]  <- sel[si,fi,nages]
              }
            }
          }
        }
      }

      selExp    <- exp(OM@selTSSign[keep(1:nsim),,1] * sin((y-nyears-1) * OM@selTSWaveLen[keep(1:nsim),,1]) * (OM@selExpRange))
      selTS     <- selTS[keep(1:nsim),,] ^ rep(selExp[keep(1:nsim),],nages)
      CPUEselTS <- CPUEsel

      for (ai in 1:nareas)
      {
        CPUEselTS[1:nsim,,ai] <- selTS[keep(1:nsim),OM@indexFisheries[ai],]
      }

      #demo plot of selectivity temporal variability
      # browser()
      # fii <- 2 # 2= LL region 1, 5=PSLS region 1
      # plot(sel[1,fii,], type='l', main='LL R1; y = ' %&% (y-61), ylab='Age', xlab='Sel',lwd=1)
      # #plot(sel[1,fii,], type='l', main='PSLS R1; y = ' %&% (y-61), ylab='Age', xlab='Sel')
      # lines(selTS[1,fii,], col=fii, lty=2)

      if ((y %in% c(nyears + proyears)) && (MP == 1))
      {
        #plot some original and  modified selectivities
        par(mfrow=c(4,4))

        if(Report){
          for (fi in 1:nfleets)
          {
            plot(sel[1,fi,], type='l', main='sel')
            lines(selTS[1,fi,], col=2)
          }

          for (ai in 1:min(16,nCPUE))
          {
            plot(CPUEsel[1,,ai], type='l', main="CPUE sel")
            lines(CPUEselTS[1,,ai], col=3)
          }
        }
      }

      # set the bridging Catches between the last OM year and the first MP year
      firstMPy <- nyears + .Object@firstMPYr - .Object@lastCalendarYr

      if (y < firstMPy)
      {
        yBridge <- y - nyears

        if ((yBridge <= length(.Object@catchBridge)) && (sum(.Object@catchBridge) > 0))
        {
          # use the known aggregate catch imported from the OMd
          TAC   <- .Object@catchBridge[yBridge]
        }
        else
        {
          #use the previous aggregate catch+error
          TAC   <- TAC * exp(rnorm(length(TAC)) * .Object@catchBridgeCV - 0.5 * .Object@catchBridgeCV ^ 2)
        }

        TAE   <- karray(rep(0, nfleets * nsim), dim=c(nsim, nfleets))
        TACE  <- cbind(TAE, TAC)
      }

      if (y %in% upyrs)
      {
        #CPUEobsY based on merger of observed and simulated CPUE (proportional to IObs in projections)
        #calculate CPUE up to but not including the current year (MP does not necessarily have access to this value depending on data lag)
        CPUEobsY[,nyears:(y - 1)] <- qCPUE * (NLLI[,nyears:(y - 1)] ^ .Object@Ibeta[keep(projSims)]) * .Object@Ierr[keep(projSims),nyears:(y - 1)]
        Iobs                      <- CPUEobsY[,1:(y - 1)]

        # Operate MP S P A Y M R          #MP applied in these years only

        # Simulate imperfect information --------------------------------------

        if (y >=firstMPy) #DK change ... xxx why +1 ?
        {
          # DK change in relation to redefined upyrs
          #update data
          if (y == firstMPy)
          {
            nuy <- (nyears + 1):(firstMPy - 1)
          }
          else
          {
            nuy <- (upyrs[match(y, upyrs) - 1]):(y - 1)
          }

          nCAA    <- sampCatch(apply(CAAF[,,keep(nuy),], c(1,2,3), sum), .Object@nCAAobs)

          set.seed(OM@seed + y * nsim)

          nCAAInd <- sampCatch(apply(CAAF[,,keep(nuy),OM@indexFisheries], c(1,2,3), sum), .Object@nCAAobs)

          set.seed(OM@seed + y * nsim + 1)

          CAA     <- abind(CAA, nCAA, along=3)
          CAAInd  <- abind(CAAInd, nCAAInd, along=3)
          #CAL     <- abind(CAL, makeCAL(nCAA, Linf=OM@Linf[,1,nuy], K=K[,1,nuy], t0=t0[1], CAL_bins), along=3)
          CAL     <- abind(CAL, makeCAL(nCAA, Len_age_mid[1,1,,nuy], CAL_bins), along=3)
          CALInd  <- abind(CALInd, makeCAL(nCAAInd, Len_age_mid[1,1,,nuy], CAL_bins), along=3)
        }

        Cobs <- apply(.Object@CM[MP,keep(projSims),,1:(y - 1)], c(1,3), sum) * .Object@Cerr[keep(projSims),1:(y - 1)]

        # xxx zzz MP rate parameters presumably need to be harmonized for YFT as well
        # xxx zzz need to finish growth curve substitution flagged ###
        # missing dimension problem temporarily commented out...

        pset<-list("Cobs"     = Cobs[,1:(y - OM@MPDataLag)],
                   #"K"        = K[,1,y - 1]*.Object@Kb[keep(projSims)],
                   #"Linf"     = Linf[,1,y - 1]*.Object@Kb[keep(projSims)],
                   "t0"       = rep(OM@t0[1],nsim),
                   "M"        = M[,1,,(y - 1)]*.Object@Mb[keep(projSims)],
                   "MSY"      = OM@MSY[keep(projSims)] * .Object@MSYb[keep(projSims)],
                   "BMSY"     = OM@BMSY[keep(projSims)] * .Object@BMSYb[keep(projSims)],
                   "UMSY"     = OM@UMSY[keep(projSims)] * .Object@FMSYb[keep(projSims)],
                   "a"        = rep(OM@a, nsim),
                   "b"        = rep(OM@b, nsim),
                   "nages"    = nages,
                   "Mat"      = mat[,1,,1:(y - 1)],

                   #need to lag data available for HCR by appropriate amount (OM@MPDataLag)
                   "CMCsum"   = apply(CMCurrent, sum, MARGIN=1),
                   "UMSY_PI"  = OM@UMSY[keep(projSims)],
                   "Iobs"     = Iobs[,keep(1:(y - OM@MPDataLag))],
                   "CAA"      = CAA[,,keep(1:(y - OM@MPDataLag))],
                   "CAL"      = CAL[,,keep(1:(y - OM@MPDataLag))],
                   "CALInd"   = CALInd[,,keep(1:(y - OM@MPDataLag))],
                   "CAL_bins" = CAL_bins,
                   "prevTACE" = TACE,
                   "y"        = y - OM@MPDataLag,
                   "tune"     = tune[MP])

        simset <- karray(1:nsim, c(nsim))

        if (.Object@UseCluster != 0)
        {
          TACE <- t(parSapply(cluster, simset, FUN=runMP, get(MPs[MP]), pset))
          printClusterOutput(simset)
        }
        else
        {
          TACE <- t(sapply(simset, FUN=get(MPs[MP]), pset))
        }
      }

      # Unpack MP TACs and TAEs
      TAC    <- TACE[,ncol(TACE)]  # TAC is final entry
      TAEbyF <- karray(TACE[,1:(ncol(TACE) - 1)], dim=c(nsim,nfleets))
      #if the fleet has a TAE, this vector is used to exclude the fleet from the TAC extractions

      # Start of annual projection
      #---------------------------------------------------------------------
      #Spatial devs in rec (affect spatial distribution but not total; streamlined implementationsame for all sims, pops, area)
      recSpatialDevs <- karray(exp(OM@ReccvR[keep(projSims)] * rnorm(length(OM@Recdist[keep(projSims),,]))), dim=dim(OM@Recdist[keep(projSims),,]))
      recSpatialDevs <- recSpatialDevs / karray(rep(apply(recSpatialDevs, FUN=mean, MARGIN=c(1,2)), nareas),dim=dim(OM@Recdist[keep(projSims),,]))

      if (.Object@CppMethod != 0)
      {
        M_Y           <- M[,,,y]
        mat_Y         <- mat[,,,y]
        Len_age_Y     <- Len_age[,,,y]
        Len_age_mid_Y <- Len_age_mid[,,,y]
        Wt_age_Y      <- Wt_age[,,,y]
        #Wt_age_SB_Y   <- Wt_age_SB_Y[,,,y]
        Wt_age_mid_Y  <- Wt_age_mid[,,,y]
        RecdevInd_Y   <- (y - 1) * nSpawnPerYr + 1
        Recdevs_Y     <- Recdevs[,,keep(RecdevInd_Y:(RecdevInd_Y + nSpawnPerYr - 1))]

        Om.nt.projection(Obj,
                         y,
                         as.integer(if (Report) 1 else 0),
                         EffortCeiling,
                         TAC,
                         TAEbyF,
                         TACEError,
                         ECurrent,
                         CMCurrent,
                         q,
                         R0,
                         M_Y,
                         mat_Y,
                         Idist,
                         Len_age_Y,
                         Wt_age_Y,
                         Wt_age_mid_Y,
                         selTS,
                         mov,
                         h,
                         Recdist,
                         Recdevs_Y,
                         recSpatialDevs,
                         OM@SRrel,
                         N_Y,
                         NBefore_Y,
                         SSN_Y,
                         C_Y,
                         SSBA_Y,
                         as.integer(100));
      }
      else
      {
        # distribute TAC by season and fleet for all sims, for all fishries that do not have TAEs
        isTACFleet          <- karray(NA, dim=dim(CMCurrent))
        SMRFim              <- as.matrix(expand.grid(1:nsim, 1:nsubyears, 1:nareas, 1:nfleets))
        SFim                <- SMRFim[,c(1,4)]
        isTACFleet[SMRFim]  <- karray(rep((!TAEbyF[1,]) * 1.0, each=nsim), dim=c(nsim, nfleets))[SFim]  # exclude TAC fleets
        TACbySMRF           <- TAC * CMCurrent * isTACFleet / apply(CMCurrent * isTACFleet, sum, MARGIN=c(1))

        for (mm in 1:nsubyears)
        {
          SPAYMRF[,5]  <- mm
          SPAMRF[,4]   <- mm
          SPAMR[,4]    <- mm
          SMRF[,2]     <- mm

          N_Y[,,,mm,]   <- NBefore_Y[,,,mm,]
          SSN_Y[,,,mm,] <- NBefore_Y[,,,mm,] * karray(rep(mat[,,,y],times=nareas), c(nsim,npop,nages,nareas))
          #potential change
          SSBA_Y        <- apply(SSN_Y[,,,mm,] * karray(rep(Wt_age[,,,y],times=nareas), dim=c(nsim,npop,nages,nareas)), c(1,2), sum, na.rm=T)
          #SSBA_Y        <- apply(NBefore_Y[,,,mm,] * karray(rep(Wt_age_SB[,,,y],times=nareas), dim=c(nsim,npop,nages,nareas)), c(1,2), sum, na.rm=T)

          # do recruitment
          if(mm %in% OM@Recsubyr)
          {
            for(pp in 1:npop)
            {
              # ie every qtr for YFT
              RecdevInd <- (y - 1) * nSpawnPerYr + mm

              # recruit fish
              if (OM@SRrel[pp] == 1)
              {
                # Beverton-Holt recruitment
                rec <- Recdevs[,pp,RecdevInd] * ((0.8 * R0[,pp] * h[,pp] * SSBA_Y[,pp]) / (0.2 * SSBpR[,pp] * R0[,pp] * (1 - h[,pp]) + (h[,pp] - 0.2) * SSBA_Y[,pp]))
              }
              else
              {
                # Most transparent form of the Ricker uses alpha and beta params
                rec <- Recdevs[,pp,RecdevInd] * aR[,pp] * SSBA_Y[,pp] * exp(-bR[,pp] * SSBA_Y[,pp])
              }

              N_Y[,pp,1,mm,]        <- rec * Recdist[,pp,] * recSpatialDevs[,pp,]
              NBefore_Y[,pp,1,mm,]  <- N_Y[,pp,1,mm,]
            }
          }

          # move fish (order of events altered from original)
          if (nareas > 1)
          {
            N_Y[,,,mm,] <- projection.domov(Ntemp=karray(N_Y[,,,mm,], dim=c(nsim,npop,nages,nareas)),
                                            movtemp=karray(mov[,,,mm,,], dim=c(nsim,npop,nages,nareas,nareas)))
          }

          #---------------------------------------------------------------------
          # Use the new population dynamics for mix of TACs and TAEs
          # Pope's approximation; values up to ~0.6 may be substantially closer
          # to Baranov depending on TAC, TAE and M

          #TACTime <- 0.5

          CNTACbySPARF <- karray(0, dim=c(nsim,npop,nages,nareas,nfleets))
          CNTAEbySPARF <- karray(0, dim=c(nsim,npop,nages,nareas,nfleets))

          # Fishing mort for TAE-managed fleets
          FTAE    <- karray(0, dim=c(nsim,npop,nages,nareas,nfleets)) #by fleet
          FTAEAgg <- karray(0, dim=c(nsim,npop,nages,nareas))         #aggregated over fleets

          # Define some index matrices
          SPARFim <- SPAYMRF[,c(1,2,3,6,7)]
          SFim    <- SPAYMRF[,c(1,7)]
          SMRFim  <- SPAYMRF[,c(1,5,6,7)]

          # TAE Fs by fleet
          FTAE[SPARFim] <- ECurrent[SMRFim] * TAEbyF[SFim] * selTS[SFA] * q[SF] * TACEError[SFim]

          #sum TAE Fs over fleets
          FTAEAgg <- apply(FTAE, MARGIN=c(1:4), sum)

          #---------------------------------------------------------------------
          # first half timestep natural mort and F for TAEs (before TAC)

          CNTAEbySPARF[SPARFim] <- FTAE[SPARFim] / (FTAEAgg[SPAR] + M[SPAY] / nsubyears) * (1.0 - exp(-TACTime * (M[SPAY] / nsubyears + FTAEAgg[SPAR]))) * N_Y[SPAMR]
          N_Y[SPAMR]            <- N_Y[SPAMR] * exp(-TACTime * (M[SPAY] / nsubyears + FTAEAgg[SPAR]))

          #---------------------------------------------------------------------
          # mid-year TAC extraction

          CNTACbySPAR <- karray(0.0, dim=c(nsim,npop,nages,nareas))

          # skip if all TACs = 0
          if ((sum(isTACFleet) > 0) && (TAC > 0.0))
          {
            # xxx still need to add TAC implementation error
            # xxx some of this is probably incorrect for multi-population case

            # Vulnerable biomass and numbers for each fishery by pop
            VBbySPARF    <- karray((N_Y[SPAMR] * Wt_age_mid[SPAY] * selTS[SFA]), dim=c(nsim,npop,nages,nareas,nfleets))

            # VB summed over ages (and populations)
            VBbySRF      <- apply(VBbySPARF,sum,MARGIN=c(1,4,5))

            # U for each fishery and region independently (TAC could be unachievable , i.e. U>1)
            UbySRFtest                 <- TACbySMRF[,mm,,] / VBbySRF[]
            UbySRFtest[UbySRFtest > 9] <- 9 #bound ridiculous U to prevent exp(U) -> inf

            # U for each age by region and fishery  (could be unachievable, i.e. U>1)
            UbySARFtest       <- karray(NA, dim=c(nsim,nages,nareas,nfleets))
            UbySARF           <- karray(NA, dim=c(nsim,nages,nareas,nfleets))
            UbySARFtest[SARF] <- UbySRFtest[SRF] * selTS[SFA] * q[SF] * TACEError[SFim]

            # U for each age by region aggregated over fisheries
            UbySARtest <- apply(UbySARFtest, sum, MARGIN=c(1:3))

            # ad hoc limit on U
            UbySAR                     <- UbySARtest
            UbySAR[UbySARtest > rULim] <- rULim+(1-rULim-0.3)*(1-exp(-UbySARtest[UbySARtest > rULim]+rULim))

            # rULim=50% original: proportional to U=0.5, then approaches an asymptote of 0.88;
            #if(rULim == 0.5) UbySAR[UbySARtest > 0.5] <- exp(UbySARtest[UbySARtest > 0.5]) / (1.0 + exp(UbySARtest[UbySARtest > 0.5])) - 0.122
            # rULim=30% - seemingly best (coupled with the TACtime=0.01); proportional to U=0.3, then approaches an asymptote of 0.28
            #if(rULim == 0.3) UbySAR[UbySARtest > 0.3] <- exp(UbySARtest[UbySARtest > 0.3]) / (1.0 + exp(UbySARtest[UbySARtest > 0.3])) - 0.28
            #test  H65
            #UbySAR[UbySARtest > 0.65] <- exp(UbySARtest[UbySARtest > 0.65]) / (1.0 + exp(UbySARtest[UbySARtest > 0.65])) - 0.01
            #test  H10
            #UbySAR[UbySARtest > 0.1] <- exp(UbySARtest[UbySARtest > 0.1]) / (1.0 + exp(UbySARtest[UbySARtest > 0.1])) - 0.43
            #test  H99
            #UbySAR[UbySARtest > 0.99] <- 0.99

            # rescale U for each fishery as a proportion of the fishery-aggregate U (only relevant for those that exceed U limit)
            UbySARF[SARF] <- UbySAR[SAR] * UbySARFtest[SARF] / UbySARtest[SAR]

            # TAC catch in numbers (U identical for all pops, xxx probably not correct for multiple pops)
            CNTACbySPARF[SPARFim] <- UbySARF[SARF] * N_Y[SPAMR]

            # aggregate catch over fleets
            CNTACbySPAR <- karray(apply(CNTACbySPARF, sum, MARGIN=c(1:4)), dim=c(nsim,npop,nages,nareas))

            # Update N post-TAC extraction
            N_Y[SPAMR] <- N_Y[SPAMR] - CNTACbySPAR[SPAR]
          }

          #---------------------------------------------------------------------
          # Second TAE extraction (and M), following TAC extraction

          # Catch from first TAE + second TAE extraction
          CNTAEbySPARF[SPARFim] <- CNTAEbySPARF[SPARFim] + FTAE[SPARFim] / (FTAEAgg[SPAR] + M[SPAY] / nsubyears) * (1.0 - exp(-(1-TACTime) * (M[SPAY] / nsubyears + FTAEAgg[SPAR]))) * N_Y[SPAMR]
          N_Y[SPAMR]            <- N_Y[SPAMR] * exp(-(1 - TACTime) * (M[SPAY] / nsubyears + FTAEAgg[SPAR]))
          CNTAEbySPAR           <- karray(apply(CNTAEbySPARF, sum, MARGIN=c(1:4)), dim=c(nsim,npop,nages,nareas))

          # Total Catch in numbers from TAE and TAC
          C_Y[SPAMRF] <- CNTAEbySPARF[SPARF] + CNTACbySPARF[SPARFim]

          if (Report)
          {
            # aggregate catch in mass for rep 1
            CMTACsum <- apply(CNTACbySPAR[,,,] * karray(rep(Wt_age_mid[,,,y],nareas),dim=c(nsim,npop, nages, nareas)),sum,MARGIN=1)
            CMTAEsum <- apply(CNTAEbySPAR[,,,] * karray(rep(Wt_age_mid[,,,y],nareas),dim=c(nsim,npop, nages, nareas)),sum,MARGIN=1)

            print("CMass TAC, TAE, combined:")
            print(CMTACsum)
            print(CMTAEsum)
            print(CMTACsum + CMTAEsum)
          }

          # end harvest calculations
          #---------------------------------------------------------------------

          #  age individuals
          NBefore_Y[,,nages,mm + 1,]          <- N_Y[,,nages - 1,mm,] + N_Y[,,nages,mm,]
          NBefore_Y[,,2:(nages - 1),mm + 1,]  <- N_Y[,,1:(nages - 2),mm,]
          NBefore_Y[,,1,mm + 1,]              <- 0
        } # season loop
      }

      # calculate LL selected numbers for the year for the abundance index
      NLLbySAMR <- apply(N_Y[,,keep(1:nages),1:nsubyears,], MARGIN=c(1,3,4,5), FUN=sum, na.rm=T)

      for (isubyears in 1:nsubyears)
      {
        for (iCPUE in 1:OM@nCPUE)
        {
          NLL_Y[1:nsim,isubyears,iCPUE] <- apply(NLLbySAMR[,,isubyears,OM@CPUEFleetAreas[iCPUE]] * CPUEselTS[,,iCPUE], FUN=sum, MARGIN=1, na.rm=T)
        }
      }

      NLLI_Y <- apply(NLL_Y, FUN=sum, MARGIN=1, na.rm=T)

      # copy results back into historic data
      SSBA[,,y] <- SSBA_Y
      NLLI[,y]  <- NLLI_Y
      NLL[,y,,] <- NLL_Y

      # Store results ...
      .Object@CM[MP,projSims,,y]     <- apply(C_Y * karray(Wt_age_mid[,,,y], c(nsim,npop,nages,nsubyears,nareas,nfleets)), c(1,2), sum)
      .Object@CMbyF[MP,projSims,,y,] <- rep(apply(C_Y * karray(Wt_age_mid[,,,y], c(nsim,npop,nages,nsubyears,nareas,nfleets)), c(1,2,6), sum))

      CAAF[,,y,]  <- apply(C_Y, c(1,3,6), sum)


      #B[,,y]      <- apply(NBefore_Y[,,,keep(1:nsubyears),] * karray(Wt_age_mid[,,,y], c(nsim,npop,nages,nsubyears,nareas)), c(1,2), sum)
      BbyS        <- apply(NBefore_Y[,,,keep(1:nsubyears),] * karray(Wt_age[,,,y], c(nsim,npop,nages,nsubyears,nareas)), c(1,2,4), sum)
      B[,,y]      <- apply(BbyS, c(1,2), mean)
      Rec[,,y]    <- apply(NBefore_Y[,,1,keep(1:nsubyears),], c(1,2), sum)

      FirstIdx <- (y - 1) * nsubyears + 1
      LastIdx  <- y * nsubyears

      RecYrQtr[,,FirstIdx:LastIdx] <- apply(NBefore_Y[,,1,keep(1:nsubyears),], c(1,2,3), sum)

      # Calculate Frep
      NsoRbySPAM[,,,nsubyears + 1]      <- apply(NBefore_Y[,,,1,], FUN=sum, MARGIN=c(1:3))
      Frep                              <- findFrep(NsoRbySPAM, M[,,,y - 1], nsim, npop, nages, nsubyears)
      .Object@F_FMSY[MP,projSims,y - 1] <- Frep / OM@FMSY1[keep(projSims)]

      # Save data for next Frep calculation
      NsoRbySPAM[,,,1:nsubyears] <- apply(NBefore_Y[,,,keep(1:nsubyears),], FUN=sum, MARGIN=c(1:4))

      # set up next year starting point
      N_Y[,,,1,]       <- N_Y[,,,nsubyears + 1,]
      NBefore_Y[,,,1,] <- NBefore_Y[,,,nsubyears + 1,]

      #---------------------------------------------------------------------
      # End of annual projection

    } # projection year loop

    # Calculate Frep
    NsoRbySPAM[,,,nsubyears + 1]  <- apply(NBefore_Y[,,,1,], FUN=sum, MARGIN=c(1:3))
    Frep                          <- findFrep(NsoRbySPAM, M[,,,y], nsim, npop, nages, nsubyears)
    .Object@F_FMSY[MP,projSims,y] <- Frep / OM@FMSY1[keep(projSims)]

    # Store results ...
    # archive timing may not be entirely consistent with SS for all quantitities,
    # but should be internally consistent
    # recalculate CPUE so last years can be reported even if they are outside of MP years

    OM@CPUEobsY[projSims,nyears:y]    <- qCPUE * (NLLI[,nyears:y] ^ .Object@Ibeta[keep(projSims)]) * .Object@Ierr[keep(projSims),nyears:y]
    Iobs                              <- OM@CPUEobsY[keep(projSims),1:y]
    .Object@IobsArchive[MP,projSims,] <- Iobs

    NLLR                                   <- apply(NLL[keep(1:nsim),1:y,,], sum, MARGIN=c(1,2,4))
    CPUEobsR[1:nsim,1:y,]                  <- qCPUE * (NLLR ^ .Object@Ibeta[keep(1:nsim)]) * rep(.Object@Ierr[keep(projSims),], nCPUE)
    IobsR                                  <- CPUEobsR[,1:y,]
    .Object@IobsRArchive[MP,projSims,1:y,] <- IobsR[keep(1:nsim),1:y,]

    # note that not everything is summarized with respect to sub-populations

    .Object@SSB_SSB0[MP,projSims,,]     <- SSBA / apply(SSB0, 1, sum)
    .Object@B_B0[MP,projSims,,]         <- apply(karray(B[,targpop,], dim=c(nsim,length(targpop),allyears)), c(1,3), sum) / OM@B0[keep(projSims)]
    .Object@C_MSY[MP,keep(projSims),,]  <- karray(.Object@CM[MP,projSims,,], dim=c(nsim,length(targpop),allyears)) / OM@MSY[keep(projSims)]
    .Object@B_BMSY[MP,projSims,]        <- apply(karray(B[,targpop,], dim=c(nsim,length(targpop),allyears)), c(1,3), sum) / OM@BMSY[keep(projSims)]
    .Object@SSB_SSBMSY[MP,projSims,]    <- apply(karray(SSBA[,targpop,], dim=c(nsim,length(targpop),allyears)), c(1,3), sum) / OM@SSBMSY[keep(projSims)]
    .Object@Rec[MP,projSims,]           <- apply(karray(Rec[,targpop,], dim=c(nsim,length(targpop),allyears)), c(1,3), sum)
    .Object@RecYrQtr[MP,projSims,]      <- apply(karray(RecYrQtr[,targpop,], dim=c(nsim,length(targpop),allyears * nsubyears)), c(1,3), sum)

    cat("\n")
  } # end of MP loop

  print("end MP")
  if (.Object@CppMethod != 0)
  {
    Om.destroy(Obj)
  }

  .Object@MPs <- MPs

  return (.Object)
}


setMethod("initialize", "MSE", function(.Object, OM, MPs, interval=3, Report=F, CppMethod=NA, UseCluster=NA, EffortCeiling = as.double(20.0), TACTime = 0.5, rULim=0.5)
{
  if (class(OM) != 'OMss')
  {
    print(paste('Could not run MSE:',deparse(substitute(OM)),'not of class OMss'))
    stop()
  }

  fixed_MPs     <- c()
  fixed_MPs_Idx <- c()
  tuned_MPs     <- c()
  tuned_MPs_Idx <- c()
  cn            <- 1

  for (MP in MPs)
  {
    MP_class <- class(get(MP))

    if (MP_class == 'IO_MP')
    {
      # Normal MP
      fixed_MPs     <- c(fixed_MPs, MP)
      fixed_MPs_Idx <- c(fixed_MPs_Idx, cn)

    } else if (MP_class == 'IO_MP_tune')
    {
      # Tunable MP
      tuned_MPs     <- c(tuned_MPs, MP)
      tuned_MPs_Idx <- c(tuned_MPs_Idx, cn)

    } else
    {
      print(paste('Could not run MSE:',deparse(substitute(MP)),'not of class IO_MP or IO_MP_tune'))
      stop()
    }

    cn <- cn + 1
  }

  .Object@OMName        <- OM@Name
  .Object@Label         <- OM@Label
  .Object@SBlim         <- OM@SBlim
  .Object@Flim          <- OM@Flim
  .Object@TACTime       <- TACTime
  .Object@rULim         <- rULim
  .Object@Date          <- OM@Date
  .Object@Author        <- OM@Author
  .Object@Notes         <- OM@Notes
  .Object@PrimarySource <- OM@PrimarySource
  .Object@F_FMSYss      <- OM@F_FMSYss
  .Object@UseCluster    <- OM@UseCluster
  .Object@OMid          <- OM@OMid

  if(!is.na(UseCluster)) OM@UseCluster <- UseCluster

  if (is.na(CppMethod))
  {
    .Object@CppMethod     <- OM@CppMethod
  }
  else
  {
    if (CppMethod != OM@CppMethod)
    {
      print("Warning: Projection sub-routine differs from that used for MSY calculations")
      .Object@CppMethod <- CppMethod
    }
  }

  # -------------------------------------------------------------------------

  cat("Constructing karrays")
  cat("\n")

  flush.console()

  .Object@yrLabels       <- OM@yrLabels
  .Object@seasonLabels   <- OM@seasonLabels
  .Object@yrSeasLabels   <- OM@yrSeasLabels
  .Object@lastCalendarYr <- OM@lastCalendarYr
  .Object@firstMPYr      <- OM@firstMPYr
  .Object@MPDataLag      <- OM@MPDataLag
  .Object@catchBridge    <- OM@catchBridge
  .Object@catchBridgeCV  <- OM@catchBridgeCV

  # Dimensions  S P A Y M R
  nallsims  <- OM@nsim
  npop      <- OM@npop
  nyears    <- OM@nyears
  proyears  <- OM@proyears
  nages     <- OM@nages
  nsubyears <- OM@nsubyears
  nareas    <- OM@nareas
  nfleets   <- OM@nfleets
  nCPUE     <- OM@nCPUE
  targpop   <- as.integer(OM@targpop)
  allyears  <- nyears + proyears
  nMPs      <- length(MPs)

  .Object@nsim      <- nallsims
  .Object@npop      <- npop
  .Object@nyears    <- nyears
  .Object@proyears  <- proyears
  .Object@nages     <- nages
  .Object@nsubyears <- nsubyears
  .Object@nareas    <- nareas
  .Object@nfleets   <- nfleets
  .Object@targpop   <- targpop
  .Object@nMPs      <- nMPs
  .Object@tune      <- rep(1, times=nMPs)

  .Object@tunePMProjPeriod <- OM@tunePMProjPeriod

  if (identical(OM@tunePMProjPeriod, numeric(0)))
  {
    .Object@tunePMProjPeriod = .Object@proyears - (.Object@firstMPYr - .Object@lastCalendarYr)
  }

  allowed_years <- c(1,3,5,10,20,.Object@proyears - (.Object@firstMPYr - .Object@lastCalendarYr),1001,1002)

  if (!any(allowed_years == .Object@tunePMProjPeriod))
  {
    stop(paste("tunePMProjPeriod must be one of ", paste(allowed_years, collapse=","), sep=""))
  }

  # Define karrays
  # ---------------------------------------------------------------------------

  # Management variables
  .Object@CM          <- karray(NA,c(nMPs,nallsims,npop,allyears))
  .Object@CMbyF       <- karray(NA,c(nMPs,nallsims,npop,allyears,nfleets))
  .Object@C_MSY       <- karray(NA,c(nMPs,nallsims,npop,allyears))
  .Object@SSB_SSB0    <- karray(NA,c(nMPs,nallsims,npop,allyears))
  .Object@B_B0        <- karray(NA,c(nMPs,nallsims,npop,allyears))
  .Object@TAC         <- karray(NA,c(nMPs,nallsims,allyears))
  .Object@F_FMSY      <- karray(NA,c(nMPs,nallsims,allyears))
  .Object@B_BMSY      <- karray(NA,c(nMPs,nallsims,allyears))
  .Object@SSB_SSBMSY  <- karray(NA,c(nMPs,nallsims,allyears))
  .Object@Rec         <- karray(NA,c(nMPs,nallsims,allyears))
  .Object@IobsArchive <- karray(NA,c(nMPs,nallsims,allyears))           #observed CPUE
  .Object@IobsRArchive<- karray(NA,c(nMPs,nallsims,allyears,nareas))    #observed CPUE by R
  .Object@RecYrQtr    <- karray(NA,c(nMPs,nallsims,allyears*nsubyears))

  # Generate observation errors
  # ---------------------------------------------------------------------------

  set.seed(OM@seed)

  .Object@Cimp  <- runif(nallsims, OM@Ccv[1], OM@Ccv[2])
  .Object@Cb    <- trlnorm(nallsims, 1, OM@Cbcv)
  .Object@Cerr  <- karray(trlnorm(nallsims * allyears, rep(.Object@Cb, allyears), rep(.Object@Cimp, allyears)), c(nallsims, allyears))
  .Object@Iimp  <- runif(nallsims, OM@Icv[1], OM@Icv[2])
  .Object@Css   <- OM@Css

  IrndDevs          <- karray(rnorm(nallsims * allyears, 0, rep(.Object@Iimp, allyears)), c(nallsims, allyears))
  IrndDevs[,nyears] <- OM@initIDev[] #initial dev based on historical obs

  for (t in (nyears + 1):allyears)
  {
    IrndDevs[,t] <- OM@IAC * IrndDevs[,t - 1] + IrndDevs[,t] * sqrt(1.0 - OM@IAC ^ 2)
  }

  #season X area CPUE indices - autocorrelated error structure only properly implemented for aggregate
  CPUErndDevs <- karray(rnorm(nallsims * allyears * nsubyears, 0, rep(.Object@Iimp, allyears * nsubyears)), c(nallsims, allyears, nsubyears))

  #xxx Initial deviate not correlated to historical obs
  for (iYr in (nyears + 1):allyears)
  {
    for (iSeas in 1:nsubyears)
    {
      if (iSeas == 1)
      {
        iLastYrIndex    <- iYr - 1
        iLastSeasIndex  <- nsubyears
      }
      else
      {
        iLastYrIndex    <- iYr
        iLastSeasIndex  <- iSeas - 1
      }

      CPUErndDevs[,iYr,iSeas] <- OM@IAC * CPUErndDevs[,iLastYrIndex, iLastSeasIndex] + CPUErndDevs[,iYr, iSeas] * sqrt(1.0 - OM@IAC ^ 2)
    }
  }

  .Object@Ierr    <- exp(IrndDevs[,]) * OM@ITrend[,]
  .Object@CPUEerr <- exp(CPUErndDevs[,,])
  .Object@Ibeta   <- exp(runif(nallsims, log(OM@Ibeta[1]), log(OM@Ibeta[2])))
  CPUEobsR        <- as.karray(apply(OM@CPUEobsMR,MARGIN=c(1,2,4),sum))

  .Object@Btimp   <- runif(nallsims, OM@Btcv[1], OM@Btcv[2])
  .Object@Btb     <- trlnorm(nallsims, 1, OM@Btbcv)
  .Object@Bterr   <- karray(trlnorm(nallsims* allyears, rep(.Object@Btb, allyears), rep(.Object@Btimp, allyears)), c(nallsims, allyears))

  .Object@Mb      <- trlnorm(nallsims, 1, OM@Mbcv)
  .Object@Kb      <- trlnorm(nallsims, 1, OM@Kbcv)
  .Object@Linfb   <- trlnorm(nallsims, 1, OM@Linfbcv)
  .Object@t0b     <- rep(1, nallsims)

  .Object@MSYb    <- trlnorm(nallsims, 1, OM@MSYbcv)
  .Object@BMSYb   <- trlnorm(nallsims, 1, OM@BMSYbcv)
  .Object@IMSYb   <- trlnorm(nallsims, 1, OM@IMSYbcv)
  .Object@FMSYb   <- trlnorm(nallsims, 1, OM@FMSYbcv)
  .Object@FMSY_Mb <- trlnorm(nallsims, 1, OM@FMSY_Mbcv)

  .Object@nCAAobs <- ceiling(runif(nallsims, OM@nCAAobs[1], OM@nCAAobs[2]))

  .Object@ageMb   <- trlnorm(nallsims, 1, OM@ageMbcv)

  # Fleet-specific Error multiplier for TAC and TAE applications
  TACEErrorAll    <- karray(exp(rnorm(nallsims * OM@nfleets) * rep(OM@TACEcv - 0.5 * OM@TACEcv ^ 2, each=nallsims)), c(nallsims,nfleets))

  # break the processing up into a smaller numbers of sims to avoid excessive
  # memory usage
  ncores    <- detectCores()
  nsimsleft <- nallsims
  lastSim   <- 0

  if (.Object@UseCluster != 0)
  {
    cluster <- makeCluster(ncores)

    if (length(MP_FunctionExports) > 0)
    {
      clusterExport(cluster, MP_FunctionExports)
    }
  }

  RNG_state <- .GlobalEnv$.Random.seed

  HasTuning <- !identical(OM@tunePM, character(0))

  if (HasTuning && (length(tuned_MPs) > 0))
  {
    for (idx in tuned_MPs_Idx)
    {
      print(paste("tuning ", MPs[idx]))

      # define an optimisation function for obtaining desired MP tuning
      # and use a log10 based tuning argument to provide a wide dynamic
      # range for the minimisation domain
      opt_fn <- function(tuneLog10)
      {
        .GlobalEnv$.Random.seed <- RNG_state

        tune      <- 10 ^ tuneLog10
        nsimsleft <- nallsims
        lastSim   <- 0

        while (nsimsleft > 0)
        {
          nsim      <- if (nsimsleft > ncores) ncores else nsimsleft
          firstSim  <- lastSim + 1
          lastSim   <- firstSim + nsim - 1
          projSims  <- firstSim:lastSim

          .Object   <- runProjection(.Object, OM, projSims, CPUEobsR, TACEErrorAll, c(MPs[idx]), interval, Report, CppMethod, UseCluster, EffortCeiling, TACTime, rULim, c(tune))

          nsimsleft <- nsimsleft - nsim
        }

        tuneValue <- tableMSE.f(.Object, MPsSub=c(1))[OM@tunePM][paste(MPs[idx], "y", .Object@tunePMProjPeriod, sep=""),]
        #tuneError <- abs(OM@tunePMTarget - tuneValue)
        tuneError <- (OM@tunePMTarget - tuneValue)^2

        print(paste("target:", OM@tunePMTarget, ", value:", tuneValue, ", tuning error:", tuneError, ", tune:", tune))

        return (tuneError)
      }

      res <- optimise(opt_fn, interval=OM@tuneLogDomain, tol=OMd@tuneTol)

      .Object@tune[idx] <- 10 ^ res$minimum
    }
  }

  .GlobalEnv$.Random.seed <- RNG_state

  while (nsimsleft > 0)
  {
    nsim      <- if (nsimsleft > ncores) ncores else nsimsleft
    firstSim  <- lastSim + 1
    lastSim   <- firstSim + nsim - 1
    projSims  <- firstSim:lastSim
    tune      <- rep(1.0, times=length(MPs))

    .Object   <- runProjection(.Object, OM, projSims, CPUEobsR, TACEErrorAll, MPs, interval, Report, CppMethod, UseCluster, EffortCeiling, TACTime, rULim, .Object@tune)

    nsimsleft <- nsimsleft - nsim
  }

  if (exists("cluster"))
  {
    stopCluster(cluster)
  }

  .Object
})


runMP <- function(sim, MP, pset)
{
  require(keep)

  # In a cluster context we re-direct output to file so we can then play
  # it back in the host
  StdOutFileName <- paste("StdOutFile", sim, ".txt", sep="")
  sink(file=StdOutFileName, append=FALSE, type=c("output", "message"))

  Result <- MP(sim, pset)

  sink()

  return (Result)
}


findFrep <- function(NsoRbySPAM, M_Y, nsim, npop, nages, nsubyears)
{
  ZsoRbySPAM                 <- karray(NA, dim=c(nsim,npop,nages-2,nsubyears))
  ZsoRbySPAM[,,,1:nsubyears] <- -log(NsoRbySPAM[,,2:(nages - 1),2:(nsubyears + 1)] / NsoRbySPAM[,,1:(nages - 2),1:nsubyears])

  FsoRbySPAM <- ZsoRbySPAM - karray(M_Y[,,1:(nages-2)] / nsubyears, dim=c(nsim,npop,nages-2,nsubyears))
  Frep       <- apply(FsoRbySPAM[,,2:27,], mean, MARGIN=c(1))  # 2:27 = true ages 1:26 (1:26 used by SS)

  return (Frep)
}


cv <- function(x)
{
  sd(x) / mean(x)
}


sdconv <- function(m, sd)
{
  # get log normal standard deviation from transformed space mean and standard deviation
  (log(1.0 + ((sd ^ 2) / (m ^ 2)))) ^ 0.5
}


mconv <- function(m, sd)
{
  # get log normal mean from transformed space mean and standard deviation
  log(m) - 0.5 * log(1.0 + ((sd ^ 2) / (m ^ 2)))
}


alphaconv <- function(m, sd)
{
  m * (((m * (1.0 - m)) / (sd ^ 2)) - 1.0)
}


betaconv <- function(m, sd)
{
  (1.0 - m) * (((m * (1.0 - m)) / (sd ^ 2)) - 1.0)
}


trlnorm <- function(reps, mu, cv)
{
  return(rlnorm(reps, mconv(mu, mu * cv), sdconv(mu, mu * cv)))
}


sampCatch <- function(Csamp, nSamp)
{
  out    <- karray(NA, dim(Csamp))
  nsim   <- dim(Csamp)[1]
  nages  <- dim(Csamp)[2]
  nyears <- dim(Csamp)[3]

  for (ss in 1:nsim)
  {
    for (yy in 1:nyears)
    {
      Csampo <- Csamp[ss,,yy]

      if (sum(Csampo < 0.0) > 0.0) browser()
      if (sum(Csampo) == 0) Csampo <- rep(1.0 / nages, nages)

      out[ss,,yy] <- ceiling(rmultinom(1, size=nSamp[ss], Csampo)) #small value added for reproducibility in rnd number table position)
    }
  }

  return(out)
}

makeCAL <- function(CAA, LenAtAge, CAL_bins, CALsd=0.05)
{
  ny       <- dim(CAA)[3]
  na       <- dim(CAA)[2]
  ns       <- dim(CAA)[1]
  CALmu    <- -0.5 * CALsd ^ 2
  nCALbins <- length(CAL_bins)
  CAL      <- karray(NA, dim=c(ns,nCALbins,ny))

  for (i in 1:ns)
  {
    for (j in 1:ny)
    {
      ages      <- rep(1:na, CAA[i,,j]) + runif(sum(CAA[i,,j]), -0.5, 0.5)  # vector of the age of the samples (with error)
      lengths   <- LenAtAge[round(ages)] * exp(rnorm(sum(CAA[i,,j]), CALmu, CALsd))
      CAL[i,,j] <- c(0, hist(lengths, breaks=CAL_bins, plot=F)$counts)
    }
  }

  return(CAL)
}

projection.domov <- function(Ntemp, movtemp)
{
  # S P A R  x  S P A R R
  nareas <- dim(movtemp)[length(dim(movtemp))]

  #return dim = SPAR
  apply(karray(Ntemp, c(dim(Ntemp),nareas)) * movtemp, MARGIN=c(1,2,3,5), sum)
}
