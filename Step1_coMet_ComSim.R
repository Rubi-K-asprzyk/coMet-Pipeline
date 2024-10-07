#!/usr/bin/env Rscript

# ----------------------------------------#
# Simulation Script: G.Dauby / K.Thibault #
# --------------------------------------- #

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_COMSIM.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 
cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))

# ------------------------------ #
##### STEP 0: INITIALISATION #####
# ------------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

#### . Packages . ####

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})
  
# Install/load tons of packages.
p_load("doParallel", # Allow parallel computation
       "fields",     # Allow the display of the grid of plots
       "ape",        # Multiple tools of phylogenetic analyses
       "untb",       # Allow the creation of the fisherian ecosystem
       "TreeSim",    # Allow the creation of the phylogenetic trees.
       "dplyr",      # Allow the use of the plyr syntax
       "tidyr",      # Allow the use of function "pivot_longer"
       "stringr",    # Allow the use of function "str_remove"
       "ggplot2",    # Graphical representation tools
       "magrittr",   # Allow the use of function "set_colnames"
       "argparser",   # Add a command_line parser 
       "progressr"
      )

# Set a global option to disable dplyr::summarise() warnings when grouping. 
options(dplyr.summarise.inform = FALSE)
options(scipen = 999)

# Call the function to generate the habitat
source("Utilities/Generate_Hab.R") # The function is called "generate_hab"

#### . Argument Parser . ####

# Create the parser
arg_parser <- arg_parser("Simulate ecological communities using the coMet pipeline.", name = "coMet_ComSim.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = 1, help = "Configuration file for coMet_ComSim.R")

# Add the parameter to have the Array_Task_ID of the cluster.
arg_parser <- add_argument(arg_parser, arg = "--Rep", short = "-R", nargs = Inf, help = "$SLURM_ARRAY_TASK_ID")

# Add the debug mode to use as an optional argument, with a default. 
arg_parser <- add_argument(arg_parser, arg = "--debug", short = "-d", nargs = 1 , help = "Use the debug mode: Save the initial data (before the C code) to test eventual bugs.", default = FALSE)

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Load them into the environnement
source(Parameters$Config)

# Print a message to inform that the right configuration file was correctly loaded. 
cli::cli_alert_info(paste0("Working on configuration file: ",Parameters$Config," - TASK ID: ",Parameters$Rep))

# ---------------------------------- #
##### STEP 1: INITIAL PARAMETERS #####
# ---------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

##### ------------ 1.A: PLOT GRID -------------------- #####

  #### . Habitat Structure . ####

hab_structure <- 
  Generate_Hab (                # Use of the function "generate_hab" to create the GRID of habitat.
    L = L,                      # Length
    W = W,                      # Width
    sq_x = sq_x,                # Number of fragmented habitat "squares" along x axis
    sq_y = sq_y,                # Number of fragmented habitat "squares" along y axis
    sq_size = sq_size,          # Size of fragmented habitat "squares"
    Nhabitat = Nhabitat         # Number of different habitat
)     

# Modify the values of L and W to be the new ones, created by the function generate_hab to respect the tore shape of the simulation
L_Tore <- as.integer(ncol(hab_structure))
W_Tore <- as.integer(nrow(hab_structure))

# Re computation of NP and NG
NP <- L_Tore*W_Tore            # |CONST| Number of plots.
NG <- Ngen*L_Tore*W_Tore*Nip   # |CONST| Number of steps during the simulation:

Nhabitat <- length(table(hab_structure))   # Number of different habitat (different colors "Squares") : Here, we have two habitats. 
HabitatP <- as.vector(hab_structure)       # Values of the habitat for each plot. 

    # ----- # 
  
  #### . Habitat Optimums . ####

HabTrait <- matrix(NA, nrow = Nhabitat, ncol = Traits)  # Create a matrix to stack the trait optimum for each habitat and each trait
rownames(HabTrait) <- paste0("Habitat-",1:Nhabitat)     # Change the row names
colnames(HabTrait) <- paste0("Trait-",1 :Traits)        # Change the column names
  
if (HabOptRand == 0){      # Set optimum habitat to have the same optimum values for all traits in one habitat.
  # For each habitat
  for (i in 1:Nhabitat){HabTrait[i,] <- rep(HabOpt[i],Traits)}  # Fill the habitat optimum for each trait
} else { # Random optimum for each habitat. 
  # For each trait (each column)
  for (j in 1:Traits) {HabTrait[,j] <- sample(HabOpt, Nhabitat, replace = FALSE)} # Randomly choose a trait optimum in the vector of possible optimal habitat
}
    # ----- # 

  #### . Plot Optimums . ####

HabitatOpTrait <- matrix(NA, nrow = NP, ncol = Traits) %>%  # Create a matrix to stack the trait optimum for each plot and each trait
  set_colnames(paste0("Trait-",1:Traits))

for (i in 1:Nhabitat) {
  for (j in 1:Traits) {                                  
    HabitatOpTrait[which(HabitatP == i), j] <- HabTrait[i, j] # Give the optimal trait value for the habitat in each plot.
  } 
}
  
    # ----- # 

##### ------------ 1.B: DISPERSION FUNCTION -------------------- #####

{
alphadisp <- deltadisp*gamma(2./betadisp)/gamma(3./betadisp) # Quelle est cette fonction ?
dmax <- as.integer(sqrt(L_Tore*L_Tore+W_Tore*W_Tore)/2)      # Max distance is divided by two to respect the tore shape of the Grid. 
df <- seq(0.1,dmax,0.1)                                      # Multiplication of the number of distance unit by 10 ? 
fEcum <- matrix(0, dmax*10 , 3)                              # Creation of the matrix of distances
fEcum[,3] <- df                  # Fill the last column with the extended units of distances
fcum <- 0                        # Initialize fcum, the cumulated frequency)
    
for (i in 1:length(df)) {                            
  fcum <- fcum + exp(-1*((df[i])/alphadisp)^betadisp)  # At each distance, sum the new computed frequency
  fEcum[i,1] <- fcum                                   # Add the computed summed frequency
  fEcum[i,2] <- exp(-1*((df[i])/alphadisp)^betadisp)   # Add the computed frequency 
}
  
fEcum[,1] <- fEcum[,1]/fcum           # Report the summed frequencies between 0 and 1
fEcum[,2] <- fEcum[,2]/sum(fEcum[,2]) # Report the absolute frequencies as they sum to 1

# Rename the data
colnames(fEcum) <- c("Freq_Cum","Freq_Abs","Dist_Units")

# Que fait tout ceci ? ---
fEcumR <- as.vector(fEcum[,1])       
fEr <- 0        # ?
fEcumPar <- c() # ? 
  for (i in 5:length(fEcumR)) {
  fEr=fEr+fEcum[i,2]  
  fEcumPar=c(fEcumPar,fEr)
}
fEcumPar=fEcumPar/max(fEcumPar)
# ---------------------- #

# Save the values of the function of dispersion
write.csv(fEcum, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step0_Dispersion_Values.csv")) 

}

    # ----- # 

##### ------------ 1.C: SAVING -------------------- #####

# --- Save the initial parameters --- #

options(scipen = 999) # Remove the scientific notation. To revert it : options(scipen = 0) 
  
# Here, the initial parameters are saved to verify and keep track of what values were used in the simulation.
Parameters_file <- as.matrix(c("Scenario" = Scenario,"Ntree" = Ntree, "Nip" = Nip, "Ngen" = Ngen, "L" = L,"L_Tore" = L_Tore,"W" = W,"W_Tore" = W_Tore, "NP" = NP, "NG" = NG, "Nrep" = Nrep,
                  "S" = S, "Competition" = Competition, "equalizing" = equalizing, "stoc" = stoc,
                  "deltadisp" = deltadisp, "betadisp" = betadisp, "Ms" = Ms, "Mno" = Mno,
                  "Nsp" = Nsp, "NSs" = NSs,
                  "MODE_DIff_Traits" = MODE_DIff_Traits,"MetaComMode"= MetaComMode,"Stochastic_Fit"= Stochastic_Fit, "TraitsEqualizing" = TraitsEqualizing, "TraitsStablizing" = TraitsStablizing, "TraitsNeutral" = TraitsNeutral, "Traits" = Traits,
                  "SigmaValues" = SigmaValues, "KlimSup" = KlimSup, "Nhabitat" = Nhabitat, "HabOpt" = HabOpt))

# Change the column name
colnames(Parameters_file) <- "Parameters"
  
# Save the Scenario Parameters
write.csv(Parameters_file, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Scenario",Scenario,"_Parameters.csv")) 
  
# Display a message
cat(rule(left = "Script initialized", line_col = "green", line = "~", col = "br_green")) 

# ------------------------------------- #
##### STEP 2: COMMUNITY SIMULATIONS #####
# ------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Display a message
cat(rule(left = paste0("BEGINNING OF THE SIMULATIONS: ",Sys.time()), line_col = "red", line = "-", col = "br_red"))
    
# Keep in memory the time already used since the opening of the script to later only have the time of the computation. 
ptm <- proc.time() 

# We will not use any loop, each replicate will be done on a different core in an array manner. 
Rep <- Parameters$Rep

# Print a Message
cat(rule(left = paste0("Starting simulation ",Rep," of Scenario ",Scenario), line_col = "white", line = "-", col = "orange"))

# Initiate the list "Full Results" to keep the results.
FullResults <- list()

    # ----- #

  #### . Metacommunity . ####

# Generation of metacommunity following log-series abundance distribution and Fisher's alpha = FISHER
MetaCom <- fisher.ecosystem(1000000, Nsp, NSs)  # The first argument is the number of wanted individuals in the metacommunity

NS <- length(MetaCom)                           # Real number of species we have at the end of the generation of the metacommunity

# If MetaComMode == 0, let the values respect the fisher log distribution
# If MetaComMode == 1, change the values of species abundance to be equal
if(MetaComMode == 1){
  MetaCom[] <- 1
}

# Change format for saving
MetaComData <- as.data.frame(MetaCom) %>%
  setNames(c("sp","ab"))

# Change the species notation to be consistent with further analysis
MetaComData$sp <- str_replace_all(MetaComData$sp,"sp.","t")

# Compute the relative abundance of each species
MetaComData <- MetaComData %>%
  select(sp,ab) %>%
  mutate(freq = ab/sum(ab))

# Save the meta community
write.csv(MetaComData, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step0_SourceMetaCom.csv"))

    # ----- #

  #### . Empty Matrices . ####

# Matrix of selection coefficient of each species for each habitat
SelectionEsp <- matrix(NA, nrow = NS, ncol = Nhabitat)

# Matrix of fitness coefficient of each specie in each habitat based on traits involved in equalizing mechanisms.
FitnessEsp <- matrix(NA, nrow = NS, ncol = Nhabitat)

# Matrix of traits optimum for each species
TraitsEsp <- matrix(NA, nrow = NS, ncol = Traits)

# Matrix of traits sigma for each species, simulating ecological amplitude of the species around its optimum
SigmaEsp <- matrix(NA, nrow = NS, ncol = Traits)

    # ----- #

  #### . Phylogenetic Tree . ####

rtree <- sim.bd.taxa(
      n = length(MetaCom),    # Number of tips (from the metacommunity)
      numbsim = 1,            # Number of tree we wanted
      lambda = 0.2,           # Speciation rate
      mu = 0,                 # Extinction rate
      complete = FALSE        # Extinct lineages are suppressed
    )[[1]]

    # ----- #

  #### . Trait Evolution . ####

# If we don't want any inferior limit to K (Klimsup == 0) : rTrueTrait == rTrait

for (i in 1:Traits) {
  if (KlimSup == 0)
    rtruetrait <- rTraitCont(rtree, model = "BM", sigma = 0.1 )  # - Intrinsic species trait value
                                # Model Brownian
                                # Sigma = 0.1 : Random deviation of the trait for each branch.

  if (KlimSup != 0 && propvarerr[i] == 0){
    rtruetrait <- blomtest(rtree, 1.5, 0.005)
    rtruetrait <- rtruetrait / sd(rtruetrait)
  }

  # Give the value of rtruetrait to rtrait
  rtrait <- rtruetrait

  # add optional measurement error / or intraspecific trait variation
  if (propvarerr[i] > 0) {
    Nmeasures = 4
    rn <- matrix(ncol = Nmeasures, nrow = NS)
    for (meas in 1:Nmeasures)
      rn[, meas] <-
      (1 - propvarerr) * rtruetrait + propvarerr * rnorm(NS, mean = 0, sd = sqrt(Nmeasures))
    rtrait <-
      apply(rn, 1, mean) # estimated species trait value based on Nmeasures
    sern <-
      apply(rn, 1, sd) / sqrt(Nmeasures) #sd on the estimated species trait value
  }

  # Give the trait to each species
  TraitsEsp[, i] <- rtrait
  rownames(TraitsEsp) <-  names(rtruetrait)

}

# Rename the data
colnames(TraitsEsp) <- paste0("Trait_",c(1:ncol(TraitsEsp)))

# Are the traits involved in different mechanisms ? --- Trait values have to be vectorized to be imported in C
  # Yes the are
if (MODE_DIff_Traits == 1) {
   MatTraitsEsp <- as.vector(round(TraitsEsp[, 1:TraitsStablizing], 3)) # Only traits values involved in stabilizing mechanisms are imported in C
  # No they are not
} else {
   MatTraitsEsp <- as.vector(round(TraitsEsp, 3))
}

    # ----- #

  #### . Species Specialization . ####

if (SigmaValues > 0 & KlimSup == 0) { # If sigma is constant across the tree (Create a matrix filled with the value of sigma)
  rownames(SigmaEsp) <-  names(rtruetrait)
  SigmaEsp[,] = SigmaValues

} else {   # If sigma value of each species with an evolved trait in the phylogeny is variable
  for (j in 1:Traits)  {
    SigmaEsp2 <- rTraitCont(rtree, model = "BM")     # intrinsic species trait value
    SigmaEsp2 <- as.matrix(SigmaEsp2/sd(SigmaEsp2))
    SigmaEsp[,j] <- exp(-SigmaEsp2)
  }
}

    # ----- #

  #### . Species Frequencies . ####

names(MetaCom) <- rtree$tip.label
Fs <- as.numeric(MetaCom / sum(MetaCom))
names(Fs) <- names(MetaCom)

    # ----- #

##### ------------ 2.A.a : Simulation Coefficients  ------------ #####

    # ----- #

#### . Selection Coefficient . ####

# This selection coefficient is based on the stabilizing traits
# SelectionEsp is initially empty
rownames(SelectionEsp) <- names(MetaCom)   # Nrow = NSpecies (NS) / Ncol = NHabitat

# Computation of the Selection Coefficient
for (i in 1:Nhabitat) {
  for (j in 1:NS) {
    # Initialization of a keeper for the selection coefficient
    TR <- 1
    # Loop across all the traits if MODE_DIff_Traits != 1 or only across the stabilizing traits
    for (k in 1:ifelse(MODE_DIff_Traits != 1, Traits, TraitsStablizing)) {
      # Computation of the selection coefficient
      SEL <- 1 / (SigmaEsp[j, k]) * exp(-0.5 * (((HabTrait[i,k] - TraitsEsp[j, k]) ^ 2) / SigmaEsp[j, k] ^ 2))
      TR <- TR * SEL
    } # End of the loop on Traits
    SelectionEsp[j, i] <- TR
  } # End of the loops on the Species
} # End of the loop on habitats

# Standardization of the Selection Coefficient
for (i in 1:length(SelectionEsp[1, ])) {
  SelectionEsp[, i] <- round(SelectionEsp[, i] / max(SelectionEsp[, i]), 2) }

# Rename the dataframe
colnames(SelectionEsp) <- paste0("habitat_", seq(1, Nhabitat, 1))


# Selection coefficient for each species within each habitat transformed into vector for exporting into C
MatSelEsp <- as.vector(round(SelectionEsp, 2))

  # ----- #

#### . Fitness Coefficient . ####

# FitnessEsp is initially empty
rownames(FitnessEsp) <- names(MetaCom)

# Computation of the Specific Fitness
if (MODE_DIff_Traits == 1) {
  for (i in 1:Nhabitat) {
    for (j in 1:NS) {
      # Initialization of a keeper for the selection coefficient
      TR <- 1
      for (k in (TraitsStablizing + 1):(TraitsEqualizing + TraitsStablizing)) {
        # Computation of the selection coefficient
        SEL <- 1 / (SigmaEsp[j, k]) * exp(-0.5 * ((( HabTrait[i,k] - TraitsEsp[j, k]) ^ 2) / SigmaEsp[j, k] ^ 2))
        TR = TR * SEL
      } # End of the loop on Traits
      FitnessEsp[j, i] <- TR
    } # End of the loops on the Species
  } # End of the loop on habitats

  # Standardization
  for (i in 1:length(FitnessEsp[1, ]))
    FitnessEsp[, i] <- round(FitnessEsp[, i] / max(FitnessEsp[, i]), 4)

# If the traits are not involved in different mechanisms
} else { FitnessEsp <- SelectionEsp }

# Rename the dataframe
colnames(FitnessEsp) <- paste0("habitat_", seq(1, Nhabitat, 1))

#### .. Stochastic_Fit call .. ####

  # Stochastic Fit is used to modify the fitness of species to 1 for one of the habitat to simulate a stochastic selection for this habitat,
  # because in that case, all species will have the same fitness and therefore chance to set up on this habitat.

if(Stochastic_Fit == 1){FitnessEsp[,1] <- 1}
if(Stochastic_Fit == 2){FitnessEsp[,2] <- 1}

  # ----- #

#### . Competition Coefficient . ####

# Maximum functional distance between two Species
MAXFunc <- max(dist(TraitsEsp[, 1:ifelse(MODE_DIff_Traits != 1, Traits, TraitsStablizing)], method = "euclidean"))

# Computation of the Matrix
InteractionsMat_vec <- foreach(k = 1:Nhabitat, .combine = c) %do% {

  # Compute interaction matrix
  Fit_Mat <- outer(FitnessEsp[,k],FitnessEsp[,k],FUN="/") # Matrix of fitness ratio
  # Modify values > 1 to 1
  Fit_Mat[Fit_Mat>1] <- 1
  # Vectorisation to be input in C
  Fit_Mat <- as.vector(round(Fit_Mat, 2))
  # return the wanted object
  return(Fit_Mat)

}

##### ------------ 2.B: C-SCRIPT CALL  -------------------- #####

##### ------------ 2.B.a : Variable Initialization  ------------ #####

  #### . Input Variable . ####

# Realized dispersal distribution
disp <- vector(mode = "integer", length = trunc(sqrt(L_Tore * L_Tore + W_Tore * W_Tore) / 2))
# Realized dispersal distribution after recruitment
dispRec <- vector(mode = "integer", length = trunc(sqrt(L_Tore * L_Tore + W_Tore * W_Tore) / 2))
# Number of dispersal item which were made between plot of same habitat (relation between habitat aggregation and dispersal abilities)
samehab <- 0
# Number of dispersal item which were made between plot of different habitat (relation between habitat aggregation and dispersal abilities)
diffhab <- 0
# Species Occurrence after simulation
X <- vector(mode = "integer", length = (NS * W_Tore * L_Tore))
# Species abundances per habitat across generations during simulation
Xhsg <- vector(mode = "integer", length = (Nhabitat * NS * (Ngen+1)))
# Habitat of each plot
habitatp <- vector(mode = "integer", length = W_Tore * L_Tore)
# Vector exporting coefficient selection that has been entered by R. To verify that the matrix "SelectionEsp" was used
SelectExport <- vector(mode = "double", length = NS)
# Number of failed recruitment
failedRec <- 0
# Number of successful recruitment
succesRec <- 0
# Total sum of the probability of being recruited
sumProbRecAll <- 0.
# Sum of the probability of being Successfully recruited
sumProbRecSucc <- 0.
# Sum of the probability of NOT being recruited
sumProbaRecFail <- 0.
# Matrix of probability of being recruited
MatProbRec <- vector(mode = "double", length = (10 * Ngen))
# Matrix of count of recruitment
MatCountRec <- vector(mode = "integer", length = (10 * Ngen))

  # ----- #

#### . Variable Download . ####

# If the Debug mode mode is activated, all the variables used for the C code to run will be saved in separate folders by replicates.
# Therefore, if some replicates show problems, the same initial values could be used and replicated to check if the problem is from the C code or these initials values.

if (Parameters$debug == T){
  write.csv(NS, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_NS.csv"))
  write.csv(Fs, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_Fs.csv"))
  write.csv(HabitatSize, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_HabitatSize.csv"))
  write.csv(MAXFunc, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_MAXFunc.csv"))
  write.csv(SelectExport, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_SelectExport.csv"))
  write.csv(SelectionEsp, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_SelectionEsp.csv"))
  write.csv(MatSelEsp, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_MatSelEsp.csv"))
  write.csv(MatTraitsEsp, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_MatTraitsEsp.csv"))
  write.csv(fEcumR, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_fEcumR.csv"))
  write.csv(fEcumPar, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_fEcumPar.csv"))
  write.csv(InteractionsMat_vec, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_InteractionsMat_vec.csv"))
  write.csv(HabitatP, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_HabitatP.csv"))
  write.csv(MetaCom, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_MetaCom.csv"))
  write.csv(TraitsEsp, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_TraitsEsp.csv"))
  write.csv(FitnessEsp, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_FitnessEsp.csv"))
  write.csv(hab_structure, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_hab_structure.csv"))
  write.tree(rtree, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Debug/Rep_",Rep,"/Scenario",Scenario,"_rtree.tree"), )
}

  # ----- #

#### . Variable Storage . ####

storage.mode(Nip) = "integer"
storage.mode(W_Tore) = "integer"
storage.mode(L_Tore) = "integer"
# storage.mode(W) = "integer"
# storage.mode(L) = "integer"
storage.mode(NG) = "double"
storage.mode(disp) = "integer"
storage.mode(dispRec) = "integer"
storage.mode(samehab) = "integer"
storage.mode(diffhab) = "integer"
storage.mode(NS) = "integer"
storage.mode(Fs) = "double"
storage.mode(X) = "integer"
storage.mode(Xhsg) = "integer"

storage.mode(habitatp) <- "integer"
storage.mode(SelectExport) <- "double"
storage.mode(HabitatSize) = "integer"

storage.mode(HabitatP) = "integer"
storage.mode(MatSelEsp) = "double"
storage.mode(Nhabitat) = "integer"
storage.mode(Ms) = "double"

storage.mode(MatTraitsEsp) = "double"
storage.mode(Traits) = "integer"
storage.mode(TraitsStablizing) = "integer"
storage.mode(Competition) = "double"
storage.mode(MAXFunc) = "double"
storage.mode(failedRec) = "integer"
storage.mode(succesRec) = "integer"
storage.mode(fEcumR) = "double"
storage.mode(Mno) = "double"
storage.mode(fEcumPar) = "double"

storage.mode(sumProbRecAll) = "double"
storage.mode(sumProbRecSucc) = "double"
storage.mode(sumProbaRecFail) = "double"
storage.mode(Ngen) = "integer"
storage.mode(MatProbRec) = "double"
storage.mode(MatCountRec) = "integer"

storage.mode(InteractionsMat_vec) = "double"
storage.mode(equalizing) = "double"
storage.mode(stoc) = "double"

test = 0
test2 = 0
rerror = 0
storage.mode(test) = "double"
storage.mode(test2) = "integer"
storage.mode(rerror) = "integer"

  # ----- #

##### ------------ 2.B.a : Simulations ------------ #####

  #### . C loading . ####

# Load the C-Model
# dyn.load("Utilities/modele14.dll")   # // On Windows
dyn.load("Utilities/coMet_C_Model.so") # // On Linux

# Compute the Simulation and return a Large list with the results
result = .C(
  "main",
  rNS = NS,
  rNip = Nip,
  rNG = NG,
  rW = W_Tore,
  rL = L_Tore,
  X = X,
  habp = habitatp,
  rS = S,
  rFs = Fs,
  rdisp = disp,
  rsamehab = samehab,
  rdiffhab = diffhab,
  rHabitatP = HabitatP,
  rMatSelEsp = MatSelEsp,
  rNhabitat = Nhabitat,
  rSelecExport = SelectExport,
  rMs = Ms,
  rMatTraitsEsp = MatTraitsEsp,
  rTraits = ifelse(MODE_DIff_Traits == 1, TraitsStablizing, Traits),
  rCompetition = Competition,
  rMAXFunc = MAXFunc,
  rfailedRec = failedRec,
  rsuccesRec = succesRec,
  RfEcumR = fEcumR,
  rMno = Mno,
  rfEcumPar = fEcumPar,
  RdispRec = dispRec,
  rsumProbRecAll = sumProbRecAll,
  rsumProbRecSucc = sumProbRecSucc,
  rsumProbaRecFail = sumProbaRecFail,
  rNgen = Ngen,
  rMatProbRec = MatProbRec,
  rMatCountRec = MatCountRec,
  rInteractionsMat = InteractionsMat_vec,
  requalizing = equalizing,
  Xhsg = Xhsg,
  rstoc = stoc,
  rtest = test,
  rtest2 = test2,
	rerror = rerror)

  #### . C unloading . ####

# Unload the C-Model
dyn.unload("Utilities/coMet_C_Model.so")

    # ----- #

##### ------------ 2.B.c : Replicate results ------------ #####

# -- Get some results -- #

# Number of failed recruitment
NbeFailed <- result$rfailedRec
# Number of successful recruitment
NbeSuccess <- result$rsuccesRec
# Number of successful recruitment from the same habitat
samehab <- result$rsamehab
# Number of successful recruitment from a different habitat
diffhab <- result$rdiffhab

# Number of dispersal events in function of the distance
disp <- round(result$rdisp / sum(result$rdisp), 4)
# Number of successfull recruitements based on the distance
dispRec <- round(result$RdispRec / sum(result$RdispRec), 4)
# Create a matrix combining the two values
Rec_by_Dist <- data.frame("Prop_Rec" = disp, "Prop_success_Rec" = dispRec) %>% mutate("Dist" = rownames(.), .before = 1)

# -- Create a matrix (TabFull) containing the number of individuals of each species, in each plot

# Create the dataframe from the vector of results X
TabFull <- as.data.frame(matrix(result$X,nrow=W_Tore*L_Tore, byrow=F))
# Change the colnames to be the vector of species.
colnames(TabFull) <- names(MetaCom)

# -- Create a matrix (TabHSG) containing the number of individuals of each species, in each habitat for each generation

# # Create an array to contains Nhabitats*Nspecies*Ngenerations values
# TabHSG <- array(data = integer(), dim=c(Nhabitat, NS, Ngen+1), dimnames=list(hab=1:Nhabitat, sp=1:NS, gen=0:Ngen))
#
# # Transfer data from Xhsg to fill TabHsg
# j=1
# for(h in 1:Nhabitat) {
#   for(s in 1:NS) {
#     for(g in 1:(1+Ngen)) {
#       TabHSG[h,s,g] = result$Xhsg[j]
#       j=j+1
#     }
#   }
# }
#
# # Transform the array into a dataframe
# TabHSG <- as.data.frame.table(TabHSG)

#-- Create "You" and "CountRec" that contains Mean recruitment probability (You) and Number of recruitment processes (CountRec) at each generation (row) for

# (1) Failed recruitment from MetaCom,
# (2) Successful recruitment from MetaCom,
# (3) Failed recruitment when no migration (within-plot recruitment)
# (4) successful recruitment when no migration (within-plot recruitment),
# (5) Failed recruitment when migration,
# (6) Successful recruitment when migration
# (7) Failed recruitment when the recruit come from the same habitat,
# (8) Successful recruitment when the recruit come from the same habitat,
# (9) Failed recruitment when the recruit come from a different habitat,
# (10) Successful recruitment when the recruit come from a different habitat,

if (Ngen > 1) {

  # Create "You" -> A probability Matrix
  Y <- result$rMatProbRec
  You <- vector(length = 10)
  for (j in 1:Ngen)
  {
    c = Y[(((j - 1) * (10)) + 1):(j * (10))]
    You = rbind(You, c)
  }
  You <- You[-1, ]
  You = as.matrix(You)

  # Create "CountRec" -> A "Count" Matrix
  Y <- result$rMatCountRec
  CountRec = vector(length = 10)

  for (j in 1:Ngen){
    c = Y[(((j - 1) * (10)) + 1):(j * (10))]
    CountRec = rbind(CountRec, c)
  }
  CountRec <- CountRec[-1, ]
  CountRec <- as.matrix(CountRec)

} # End of the creation of the two matrices

# Removal of the old data
rm(result, InteractionsMat_vec, Y)

# --- Processing and plotting statistics of recruitment across generation --- #

# Change the colnames of the recruitment results
colnames(CountRec) <-
c("nbe_failed_metacom", "nbe_success_metacom",
  "nbe_failed_within_plot", "nbe_success_within_plot",
  "nbe_failed_between_plot", "nbe_success_between_plot",
  "nbe_failed_same_hab", "nbe_success_same_hab",
  "nbe_failed_diff_hab", "nbe_success_diff_hab")

# Recruitment results
summary_stats <- CountRec %>%
  as_tibble() %>%
  mutate(gen = 1:nrow(.), .before = 1)

# Species traits results
TraitsEsp_tb <- TraitsEsp %>%
  as_tibble() %>%
  mutate(sp = row.names(TraitsEsp), .before = 1)

# Rename it
names(TraitsEsp_tb)[2:ncol(TraitsEsp_tb)] <- paste0("Trait_", seq(1, ncol(TraitsEsp), 1))

# Species fitness results
FitnessEsp_tb <- FitnessEsp %>%
  as_tibble() %>%
  mutate(sp = row.names(FitnessEsp))

# Rename it
names(FitnessEsp_tb)[1:Nhabitat] <- paste0("habitat_", seq(1, Nhabitat, 1))

# -- Save all the individuals results of each Replicates -- #
write.csv(summary_stats, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/summary_stats/summary_stats",Rep,".csv"))
write.csv(TabFull, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TabFull/TabFull",Rep,".csv"))
write.csv(TraitsEsp_tb, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TraitsEsp/TraitsEsp",Rep,".csv"))
write.csv(FitnessEsp_tb, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/FitnessEsp/FitnessEsp",Rep,".csv"))
write.csv(Rec_by_Dist, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/DispRec/DispRec",Rep,".csv"))
write.tree(rtree, file = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/PhyloTree/PhyloTree",Rep,".csv"))

# Print the time (hour) along the replicate that is finished
cat(rule(left = paste0("-- END: Simulation ",Rep," / ",Sys.time()," --"), line_col = "white", line = " ", col = "grey"))

  #### . /!\ End of the loop /!\ . ####

cat(rule(left = paste0("END OF THE SIMULATIONS: ",Sys.time()), line_col = "red", line = "-", col = "br_red"))  
# Show the elapsed time of the simulation
print(proc.time() - ptm)
cat(rule(line_col = "red", line = "-", col = "br_red"))  

# End of the Simulations

cat(rule(left = "SCRIPT COMET_COMSIM_1.R ENDING", line_col = "red", line = "-", col = "br_red")) 

# --------------- #
