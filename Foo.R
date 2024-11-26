######################
# coMet_ConfigFile.R #
######################

# This file is a config file for the coMet_Pipeline.
# It as to be modified using any text editor to fulfill the user and will be sourced in the main script. 

# The parameters with no tag can be modified at choice. 
# The parameters with the tag /FIXED/ are computed from the previous parameters and SHOULD NOT be modified. It still possible but can lead to bugs or spurious behavior. 
# The parameters with the tag /CONST/ MUST NOT be modified and are placed in this file for the sake of debugging and further upgrading of the script. 

# ---------------------------- #
##### SCENARIO NAME CHOICE #####
# ---------------------------- #

# What name will be given to each of the outputs of the pipeline. It is important to choose it carefully and it will be useful
# if the user wants to compare multiple scenarios. This name can be either a integer or a character string but lust not contains any special caracters as "/" ...

Scenario <- "Foo"

# How many cores of your computer should the script use ? 
Ncore <- 10

#--------------------------------#
##### COMMUNITIES PARAMETERS #####
#--------------------------------#

#####------------------ Replication ----------------------#####

Ngen <- 2000  # Number of generation (To be determine in function of the computation time). 
Nrep <- 3    # Number of replicates for each set of initial parameters.
Ntree <- 5  # Number of Null Model trees replicate 

# Type of randomization to perform
Randomization <- 2  # If "1", tip-shuffling is performed across the whole phylogenie
                    # If "2", tip-shuffling is performed across the species at least one time across the sample (The default)

#####---------------- Grid definition --------------------#####

Nmin <- 2      # Minimum number of species per community (Cannot be less than 2)

Nhabitat = 2   # /CONST/ Number of habitat

Nip <- 50   	 # Number of individuals per plot (A plot is a square on the simulation map).
 
sq_x = 3       # Number of fragmented habitat "squares" along x axis.
sq_y = sq_x    # |FIXED| Number of fragmented habitat "squares" along y axis.
    # For now, the sq_x and sq_y are set equal to ensure a square GRID.

sq_size = 14 # Size of fragmented habitat "squares". 
    # CAUTION: The minimum is 10 or the cores and ecotones will be overlapping and overwrite themselves.
    # CAUTION: To ensure a symmetrical repartition of the plots, sq-size must be an EVEN number. Otherwise, a warning will be thrown but the script will continue. 
    # Look carefully the GRID plotting to verify that no plots overlap. 

nbe_sq_r = 1:sq_x      # |CONST| Number of fragmented habitat "squares" along x axis (as a vector)
nbe_sq_c = 1:sq_y      # |CONST| Number of fragmented habitat "squares" along y axis (as a vector)

L <- ((sq_x*2)+1)*sq_size  # |FIXED| Length of the grid.
W <- ((sq_y*2)+1)*sq_size  # |FIXED| Width of the grid.
  # L and W are already splited for a future possibility of non-square matrices.
  # sq_x * 2 (because we have fragmented and continuous habitat) + 1 (because we have one more cotinuous than fragmented habitat) * sq_size (number of plots by "squares")

NP <- L*W              # |CONST| Number of plots.
NG <- Ngen*L*W*Nip     # |CONST| Number of steps during the simulation:
# Because we want one generation to be the replacement of a number of individuals equal to the number of individuals
# in the community, we multiply the number of generations by the number of plots and the number of individuals by plots. 

HabitatSize = sq_size       # |CONST| Size of the fragmented habitat.       
nbe_sample = sq_x*sq_y*16   # |CONST| Number of samples for the parameters entered. Used to check the validity and make some tests on further results for robustness.
  # 16 = 2 Habitats * 2 SampleTypes * 4 Samples.

#####---------------- Assembly rules --------------------#####

S <- 0		        # Intensity of environmental filtering.          
Competition <- 0  # Intensity of stabilizing mechanisms.  
equalizing <- 1   # Intensity of equalizing mechanisms.   

stoc <- 1-Competition-equalizing-S   # |CONST| Intensity of the stochastic processes (No trait structure / Must be comprised between 0 and 1).

# REMINDER:
  # Stabilizing mechanisms: Rare species, functionally different -> Competition and trait overdispersion.
  # Equalizing mechanisms: Common species, functionally identical -> Niche conservatism and trait aggregation.

#####---------------- Dispersion mode ------------------#####

deltadisp <- 1     # Mean value of inter distance plot dispersal/migration. The peak of the curve :Probability of recruitment VS Distance
betadisp <- 1.0    # Determines the shape of the distribution of the dispersion kernel. 
                        # The tail is longer for lower betadisp. 
                        # beta=1 is equivalent of a exponential distribution and beta=2 is equivalent of a lognormal distribution
                        # Here, with beta = 1, The probability of close recruitment is far more superior than the probability of far recruitment. 
Ms <- 0.0001   	   # Migration rate from 'source' meta community (see after)
Mno <- 0.8         # Probability of recruitment on-site (from the same plot) --> 1 - Mno is the proportion of recruitment from outside the plot.

#####--------- Meta Community Diversity ----------------#####

  # The meta community is a community OUTSIDE the grid of plots we are observing. It represents a pool of potential migrating species. 

Nsp <- 1000   # Expected species count S after generation of a log-series.
NSs <- 800    # The maximum number of species abundance classes to consider.
MetaComMode <- 1 # 0 = Relative abundance of species from the source meta community respect the fisher log distribution
                 # 1 = Relative abundance of species are manually modified to be identical (initial species abundance set to 1 for all species)  

#####--------- Trait(s) Simulation ---------------------#####

MODE_DIff_Traits <- 1   # Are different traits are involved in different mechanisms ? 
                           # == 1 : Traits are involved in different mechanisms.
                           # != 1 : Traits are NOT involved in different mechanisms. In that case, only the Equalizing mechanisms are taken into account. 

TraitsEqualizing <- 3    # Number of trait involved in equalizing mechanisms.
TraitsStablizing <- 3    # Number of trait involved in stabilizing mechanisms.
TraitsNeutral <- 3       # Number of traits that are neutral according to assembly rules (not involved in any assembly mechanisms).

Traits <-                # |CONST| Number of total traits for each species.
  ifelse(
    MODE_DIff_Traits != 1,
    TraitsEqualizing,
    TraitsEqualizing + TraitsStablizing + TraitsNeutral)  

  # Sigma is given for each species and characterize the level of ecological specialization: there is a trade-off between the ecological amplitude and competitive fitness.
    SigmaValues <- 2
      
    # If 0, Sigma of each species is determined by a trait evolved under a Brownian motion.
    # A small sigma gives specialist species, a big sigma gives generalist species. 

    # propvarerr must be a vector of size equal to the number of trait
propvarerr <- rep(0,Traits) # Measurement error [0,1] to intrinsic species trait value due to within species variation. 
                            # == 0 : No measurement error (Intra-specific variation) of the species traits.
KlimSup <- 0                # If defining a inferior limit of K (i.e. phylogenetic signal).

# When processes are 100% equalizing, the fitness for one of the two habitat off all the species can be set to 1 to simulate stochastic processes 
Stochastic_Fit <- 0  # 1: The continuous habitat is set to "stochastic" / 2 
                    # Anything different from 1 or 2 will be ignored.

#####--------- Habitat structure ---------------------#####

HabOpt     <- c(-0.5,0.5)    # Input of the optimal values for each different habitat (Here, two habitats).
HabOptRand <- 0           # == 0: In each habitat, all the traits receive the same value of "HabOpt".
                          # != 0: Traits optimums are chosen randomly in "HabOpt", so they can vary inside each habitat. 

# ------------------------------------------------------------ #
##### -------- Metrics to compute in coMet_Metrics ------- #####
# ------------------------------------------------------------ #

# Set metrics to TRUE to compute it, Set it to FALSE to not compute it.

PD <- TRUE
FJ <- TRUE 
PSV <- TRUE
PDab <- TRUE 
MPD <- TRUE 
MNTD <- TRUE 
ED <- TRUE 
PDb <- TRUE 
UniPhylo <- TRUE 
S <- TRUE 
PIst <- TRUE 
PCD <- FALSE

# ---- #
Simpson <- TRUE 

# -------------------------------------------------------------# 
##### -------- Tests on the parameters ------------------- #####
# -------------------------------------------------------------# 

# Is it a square matrix ?  
if (L != W) {stop("STOP: Non-square grid not yet implemented. Please enter an identical Length (L) and Width (W) of Grid or number of squares in x and y (sq_x / sq_y).")}

# Is the square size odd ? 
if (sq_size %% 2 != 0){warning("WARNING: You entered an ODD sq-size, therefore the plots could not be symmetrical. Have a carefull look at the GRID plotting.")}

# Is the square size < 9? 
if (sq_size < 10){stop("STOP: You entered a sq_size < 10, therefore the plots cores and ecotones will overlap and mess the results. Please enter a sq_size >= 10.")}

# Is the number of habitat optimum equal to the number of habitat ? 
if (length(HabOpt) != 2 ){stop("STOP: Number of habitat optimums entered different from the number of habitat.")}

# If Stochastic fit is set up to 1 or 2, are the equalizing processes set to 100% ? 
if ((Stochastic_Fit == 2 || Stochastic_Fit == 1) && equalizing !=1 ){stop("STOP: Stochastic Fit is set up but Equilizing processes are not set up to 100% ")}


# --------------------------------------------------------- # 
##### -------- DEBUGGING PARAMETERS ------------------- #####
# --------------------------------------------------------- # 