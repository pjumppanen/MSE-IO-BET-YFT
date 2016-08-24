# ============================================================================================================================================
# === Operating Model definition object 'OMd' ============================================================================
# ============================================================================================================================================

# --- Create a blank OM definition object --------------------------------------
OMd<-new('OMd')

# --- Description --------------------------------------------------------------
OMd@Name          <- "OMbet0c"
OMd@Label         <- "OMbet0c" #useful for changing graphics labels
OMd@Date          <- "Jun 2016"
OMd@Author        <- "D.Kolody"
OMd@Notes         <- "Low error, no stochastic variability 18 SS model ensemble specification for testing"
OMd@PrimarySource <- "test model using C++ projections"
OMd@CppMethod     <- 1 # 1 = Use C++ Baranov solution, 0 = Use R based Baranov approximation
OMd@UseCluster    <- 0 # 1 = Use cluster of R processes, 0 = Use single R process


# --- Specifications -----------------------------------------------------------
OMd@npop          <- as.integer(1) #multiple stocks not supported at this time
OMd@nfleets       <- as.integer(12) #number of fisheries; needs to be consistent with SS input files
OMd@SSRootDir     <- "./OMconditioning/BET/gridB1/"          # Root dir for SS results outputs
OMd@SBlim         <- 0.5
OMd@Flim          <- 1.3

source("./OMconditioning/RStuff/makeGridB1.f.R")
gridB1List <- makeGridB1.f(makeGrid=F)
OMd@OMList        <- as.list(gridB1List) #[c(1,2,17,18)]) #("R1P1_h70_M06_q0Mhigh") #gridB1List[c(1,18)])
OMd@nsimPerOMFile <- array(rep(1,length(OMd@OMList)),dim=length(OMd@OMList))   # Number of simulations per each SS specification file (vector of length OMList allows differental weighting, i.e. c(10,50,25)
OMd@proyears      <- as.integer(50)              # Number projection years
OMd@targpop       <- as.integer(1)               # summary stats by population; irrelevant for single stock case
OMd@seed          <- as.integer(1)               # rnd seed
OMd@recentPerFirst<- as.integer(4)               # number of most recent seasons to include in "recent" C and E definition counting backward (0 means use last season of assessment)
OMd@recentPerLast <- as.integer(19)               # number of most recent seasons to include in "recent" C and E definition counting backward
OMd@seasonCEDist  <- as.integer(0)               # 0/1 - 1=use seasonal pattern of C/E (recentPeriod must be multiple of recentPeriod) 0 = all seasons equal

#some assessment specific time mapping requirements to align SS years as quarters with OM year-seasons and real-time
OMd@nsubyears      <- as.integer(4)
OMd@lastSeas       <- as.integer(4)        # just in case whole years not used
OMd@firstSeas      <- as.integer(1)        # just in case whole years not used
OMd@firstSSYr      <- as.integer(101)      # SS equivalent of firstYrToPlot
OMd@firstCalendarYr<- as.integer(1951)     # used to convert to actual years
OMd@lastCalendarYr <- as.integer(2012)     # used to convert to actual years
OMd@firstMPYr      <- as.integer(2018)     # MP management kicks in here, but projections using an observed catch or assumed effort may be required between lastCalendarYr and firstMPYr
OMd@MPDataLag      <- as.integer(3)        # The lag in number of years between assessment data availability and the timing of the assessment/HCR calculation; +1 means data is one year behind
OMd@catchBridge    <- as.karray(c(1,1,1,1,1))    # no catch
#OMd@catchBridge    <- as.karray(c(10000,200000,10000,200000,10000))    # known catch history between last assessment year in OM and firstMPYr (length of min 0 to max firstMPYr-lastCalendarYr-1)
#OMd@catchBridge    <- as.karray(c(-999))    # known catch history between last assessment year in OM and firstMPYr (length of min 0 to max firstMPYr-lastCalendarYr-1)
OMd@catchBridgeCV  <- 0.1                #error to add onto catch for bridge years with unknown catch

# fleet numbers to include in size comp sampling for MPs, i.e. probably the same as the survey (CPUE) fleets
OMd@indexFisheries   <- as.integer(2,3,4)       # LL CPUE fleets

#These are dimmed in the OM
OMd@ReccvTin <-array(rep(0.,OMd@npop), dim=OMd@npop)       # Temporal variability in recruitment aggregate over regions
OMd@ReccvRin <- 0.                                        # Spatial variability in recruitment for all pops, areas, sims (i.e. streamlined implemented, because its probably a low priority)
OMd@RecACTin <-array(rep(0.5,OMd@npop), dim=OMd@npop)       # Recruitment autocorrelation (for regional aggregate)
OMd@NInitCV      <- 0.
OMd@NInitCVdecay <- 0.
OMd@selExpRange  <- 0.             # sel temporal variability exponent - oscillates with a sin wave rangeing between +/- exp(selExpRange)
OMd@selAgeRange  <- 0.               # 0=no age shift, 2 means (discretized) sine wave shift of sel vector between - 2 and + 2 age class
OMd@selWLRange   <- array(c(0.0625,0.5))   # sel temporal variability wavelength range (0.0625 = quarter wavelength in 25 years, 0.5=2 full cycles in 25 years

#Observation errors (original ABT code had a separate observation class)
OMd@TACEcv <-array(rep(0.,OMd@nfleets)) # fleet-specific lognormal errors on TAC/TAE (independent among fleets and seasons)
OMd@Ccv    <-c(0., 0.0001)
OMd@Icv    <-c(0., 0.0001)
OMd@Dcv    <-c(0., 0.0001)
OMd@Btcv   <-c(0., 0.0001)
OMd@Ftcv   <-c(0., 0.0001)
OMd@Cbcv   <- 0.
OMd@Mbcv   <- 0.
OMd@LFCbcv <- 0.
OMd@LFSbcv <- 0.
OMd@ageMbcv<- 0.
OMd@Ftbcv  <- 0.
OMd@Recbcv <- 0.
OMd@IMSYbcv<- 0.
OMd@MSYbcv <- 0.
OMd@BMSYbcv<- 0.
OMd@hbcv   <- 0.
OMd@Btbcv  <- 0.
OMd@Dbcv   <- 0.
OMd@Kbcv       <- 0.0001
OMd@t0bcv      <- 0.0001
OMd@Linfbcv    <- 0.0001
OMd@FMSYbcv    <- 0.
OMd@FMSY_Mbcv  <- 0.
OMd@BMSY_B0bcv <- 0.
OMd@nCAAobs    <- c(100,101)
OMd@nCALobs    <- 1000
OMd@Lcv        <- c(0., 0.)
OMd@Ibeta      <- c(0.999, 1.0001) # hyperstability parm, cv # c(0.66,1.5) #exp(runif(nsim,log(0.75),log(1.25))) #check definition
OMd@IACin      <- 0.5             # cpue autocorrelation
OMd@ITrendin   <- -1               # cpue trend % per annum compounded, negative means the trend is extracted from the assessment model filname, i.e. q1 = 1%






# --- save object --------------------------------------------------------------
save(OMd,file=paste(getwd(),"/Objects/OMd.OMbet0c.RDA",sep=""))


# ==========================================================================================================================
# End of build sript ========================================================================================================
# ==========================================================================================================================
