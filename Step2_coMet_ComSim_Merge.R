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

# Call the function to generate the habitat
source("Utilities/Generate_Hab.R") # The function is called "generate_hab"

#### . Local Mode . ####

Parameters <- list()
Parameters$Config <- "coMet_ConfigFiles/Equalizing100_Sig5.R"

#### . Argument Parser . ####

# Create the parser
arg_parser <- arg_parser("Simulate ecological communities using the coMet pipeline.", name = "coMet_ComSim.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = 1, help = "Configuration file for coMet_ComSim.R")

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Load them into the environnement
source(Parameters$Config)

# Print a message to inform that the right configuration file was correctly loaded. 
cli::cli_alert_info(paste0("Working on configuration file: ",Parameters$Config))

# # ---------------------------------- #
# ##### STEP 1: INITIAL PARAMETERS #####
# # ---------------------------------- #
# 
# ##### ______________________________________________________________________________________________________________________________________ #####
# 
# ##### ------------ 1.A: PLOT GRID -------------------- #####
# 
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
NP <- L_Tore*W_Tore              # |CONST| Number of plots.
NG <- Ngen*L_Tore*W_Tore*Nip     # |CONST| Number of steps during the simulation:

Nhabitat <- length(table(hab_structure))   # Number of different habitat (different colors "Squares") : Here, we have two habitats.
HabitatP <- as.vector(hab_structure)       # Values of the habitat for each plot.

##### ------------ 2.C: RESULT MERGING  -------------------- #####

# Display a message
cat(rule(left = "Merging the simulations outputs", line_col = "green", line = "~", col = "br_green"))

   # ----- #

#### . Loading and Merging . ####

  # Reloading of each of the intermediate results.

# -- Recruitment  -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/summary_stats"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(Recruitment_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {
  
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.csv(File[j],row.names = 1) %>% # Load the file
    mutate(rep = rep,.before = 1) # Append the replicate
  # remove the file as it goes
  # file.remove(File[j]) ; Raw

})

# -- Community abundance data   -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TabFull"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(Community_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.csv(File[j],row.names = 1) %>% # Load the file
    mutate(rep = rep,.before = 1) # Append the replicate
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
})

# -- Species Traits -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TraitsEsp"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(TraitsEsp_Total  <- foreach(j = 1:length(File), .combine = full_join) %do% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.csv(File[j],row.names = 1) %>% # Load the file
    mutate(rep = rep,.before = 1) # Append the replicate
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
})

# -- Species Fitness -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/FitnessEsp"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(FitnessEsp_Total  <- foreach(j = 1:length(File), .combine = full_join) %do% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.csv(File[j],row.names = 1) %>% # Load the file
    mutate(rep = rep,.before = 1) # Append the replicate
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
})

# -- Species Fitness -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/DispRec"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(Rec_by_Dist_Total  <- foreach(j = 1:length(File), .combine = full_join) %do% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.csv(File[j],row.names = 1) %>% # Load the file
    mutate(rep = rep,.before = 1) # Append the replicate
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
})

# -- PhyloTrees -- #
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/PhyloTree"),recursive = T, full.names = T)
# Loop across the found files to load them
suppressMessages(PhyloTrees_Total  <- foreach(j = 1:length(File)) %do% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  Raw <- read.tree(File[j])
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
})

# Change the class to multiphylo
class(PhyloTrees_Total) <- "multiPhylo"

# Change the names of each tree to be the replicate number. 
for(j in 1:length(File)) {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(File[j], '\\d+')
  # Load the file
  names(PhyloTrees_Total)[j] <- paste0("rep",rep)
  # remove the file as it goes
  # file.remove(File[j]) ; Raw
}


#### . Downloading . ####

# Saving the Raw Data.
write.csv(Recruitment_Total, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Result_Count.csv"))
write.csv(TraitsEsp_Total, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Traits.csv"))
write.csv(Community_Total, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Community_Abundance.csv"))
write.csv(FitnessEsp_Total, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Fitness.csv"))
write.csv(Rec_by_Dist_Total, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Rec_by_Dist.csv"))
write.tree(PhyloTrees_Total,
           file = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_PhyloTrees.tree"),
           append = FALSE,
           digits = 10,
           tree.names = TRUE)

# ---------------------------- #
##### STEP 3: SUB-SAMPLING #####
# ---------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####
  
# Display a message
cat(rule(left = "Creating the sub_samples", line_col = "green", line = "~",col = "br_green"))
  
#### . Grid . ####

# Create a table with characteristics of each site (X, Y coordinates, habitat, ecotone/core, sample n°)
sites <- data.frame(cbind(X = rep(1:(L_Tore), times = L_Tore), Y = rep(1:(W_Tore), each = W_Tore)))

# Gives the position of each site in table "sites" according to X,Y coordinates: pos[X, Y]
pos <- matrix(1:((L_Tore)*(W_Tore)), ncol=(L_Tore)) 

# Gives the value of the habitat corresponding to each position
for(i in 1:Nhabitat){sites$hab[which(hab_structure == i)] <- i}

# Change the str() of sites$hab to be a character
sites$hab <- as.character(sites$hab)

# Create a variable to determine if the site is core or ecotone 
sites$SampleType = NA # Define whether a site is core (1) or ecotone (2)
sites$sample = NA     # Define n° of sample to which a site contributes (total = 144 samples)

# Create a counter for number sample.
sam = 1 

# Create a vector to find each of the possible center of cores as a compt token.
Center <-seq(1,max(nbe_sq_c)*2, by = 2) 

# Define samples for core hab 2 (Fragmented)  
for(x in (Center*(sq_size)+1)){
  for(y in (Center*(sq_size)+1)){   # For now, OK ! Find the left bottom of the right up corner.
    for( x2 in x:(x+2))     for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam       # Right up corner
    for( x2 in (x-3):(x-1)) for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam + 1   # Left up corner
    for( x2 in x:(x+2))     for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 2   # Bottom right corner
    for( x2 in (x-3):(x-1)) for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 3   # Bottom left corner
    sam = sam + 4
  }
}

# Define samples for core hab 1 (Continuous)
for(x in (Center+1)*sq_size+1){
  for(y in Center*(sq_size)+1){       # For now, OK ! Find the left bottom of the right up corner.
    #! In all the "inner" cases, when we are not at the edge of the tore.
    if (x != max((Center+1)*sq_size+1)){  
      for( x2 in x:(x+2))     for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam       # Right up corner
      for( x2 in (x-3):(x-1)) for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam + 1   # Left up corner
      for( x2 in x:(x+2))     for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 2   # Bottom right corner
      for( x2 in (x-3):(x-1)) for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 3   # Bottom left corner
      sam = sam + 4
      #When we are not at the edge of the tore.
    } else {
      for( x2 in 0:3)     for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam       # Right up corner
      for( x2 in (x-3):(x-1)) for( y2 in y:(y+2))     sites$sample[pos[x2,y2]] = sam + 1   # Left up corner
      for( x2 in 0:3)     for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 2   # Bottom right corner
      for( x2 in (x-3):(x-1)) for( y2 in (y-1):(y-3)) sites$sample[pos[x2,y2]] = sam + 3   # Bottom left corner
      sam = sam + 4
    }
  }
}

# Give the value to the plots
sites$SampleType[sites$sample < sam] = 1 # 1 for core habitat, 2 for ecotones

# Save this sam for after
sam_eco <- sam

# Define samples for Ecotone hab 2 (Fragmented)
for(x in Center*(sq_size)+1) for(y in Center*(sq_size)+1){   # Left side
  for( y2 in (y:(y+8)-(sq_size/2))) sites$sample[pos[x-(sq_size/2),y2]] = sam
  sam = sam + 1
}
for(x in ((Center+1)*sq_size)-(sq_size/2)) for(y in ((Center+1)*sq_size)-(sq_size/2)){ # Right side
  for( y2 in y:(y-8)) sites$sample[pos[x,y2]] = sam
  sam = sam + 1
}
for(y in (Center*sq_size)-(sq_size/2)+1) for(x in ((Center+1)*sq_size)-(sq_size/2)){ # Bottom side 
  for( x2 in x:(x-8)) sites$sample[pos[x2,y]] = sam
  sam = sam + 1
}
for(y in (Center+1)*sq_size-(sq_size/2)) for(x in Center*(sq_size)-(sq_size/2)+1){ # Up side
  for( x2 in x:(x+8)) sites$sample[pos[x2,y]] = sam
  sam = sam + 1
}

# Define samples for Ecotone hab 1 (Continuous)
for(x in Center*(sq_size)) for(y in Center*(sq_size)+1){   # Left side
  for( y2 in (y:(y+8)-(sq_size/2))) sites$sample[pos[x-(sq_size/2),y2]] = sam
  sam = sam + 1
}
for(x in ((Center+1)*sq_size)-(sq_size/2)+1) for(y in ((Center+1)*sq_size)-(sq_size/2)){ # Right side
  for( y2 in y:(y-8)) sites$sample[pos[x,y2]] = sam
  sam = sam + 1
}
for(y in (Center*sq_size)-(sq_size/2)) for(x in ((Center+1)*sq_size)-(sq_size/2)){ # Bottom side 
  for( x2 in x:(x-8)) sites$sample[pos[x2,y]] = sam
  sam = sam + 1
}
for(y in (Center+1)*sq_size-(sq_size/2)+1) for(x in Center*(sq_size)-(sq_size/2)+1){ # Up side
  for( x2 in x:(x+8)) sites$sample[pos[x2,y]] = sam
  sam = sam + 1
}

# Give the value to the plots
sites$SampleType[sites$sample < sam & sites$sample >= sam_eco] = 2

# Plot it to verify that everithing is fine
# plot(sites$X,sites$Y,col =sites$hab)
# points(sites$X,sites$Y,col =sites$sample, pch = 16)

# Save the global Site_Data
write.csv(sites,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step2_Sites_Data.csv"))

    # ----- #

#### . Communities . ####

# "Tabfull" represents our simulated communities. It is an Abundance dataset with species in columns and plots in rows
# We have "Nrep" communities of L x W plots to transform into pre_abs data or sampled data.
# "Sites" represents the coordinates and data about the sampled plots we want for the analysis

# Initialization of the number of sampled communities
N_samp <- max(nbe_sq_r)^2 * 16   # nb_sq_r^2 group * 4 ecotones * 4 cores 
  
# Create a loop to work on each community created for each replicates. 
# We wanted to create the sampled communities as well as the presence absence version.
Sampled_Data <- foreach(i = unique(Community_Total$rep)) %do% { library("dplyr")
    
  # Extract the data for one community
  TabFull <- Community_Total[Community_Total$rep == i,] %>%  # Select the wanted lines
  subset(select = (-c(rep)))
    
  # Create the matrix to collect data of sampled plots <- Abundance data
  Abund_Data <- matrix(
    nrow=N_samp,
    ncol=ncol(TabFull))       # Ncol is the number of species 
    
  # Give back to "SampleComm" the names of the species
  colnames(Abund_Data) <- colnames(TabFull)
    
# Create a matrix storing characteristics of each sample
  # Copy the first nbe_sq_r*16 rows because we only need this number
  samples <- sites[1:N_samp,] 
  
  # Find the mean "X" position of the samples (constituted by 9 plots)
  samples$X <- tapply(sites$X, INDEX = sites$sample, FUN="mean")
    
  # Find the mean "Y" position of the samples (constituted by 9 plots)
  samples$Y <- tapply(sites$Y, INDEX = sites$sample, FUN="mean")
    
  # Find the habitat (2 or 1) of the samples (constituted by 9 plots)
  samples$hab <- tapply(as.numeric(sites$hab), INDEX = sites$sample, FUN="mean")     # It's not a "proper" mean, it is just a manner of having the habitat
    
  # Find the sample type (core or ecotone) of the samples (constituted by 9 plots)
  samples$SampleType <- tapply(sites$SampleType, INDEX = sites$sample, FUN="mean", na.rm=T)
  
  # Give a number to all the samples, from 1 to N_Samp
  samples$sample <- tapply(sites$sample, INDEX = sites$sample, FUN="mean", na.rm=T)
  
  # Add the rep ID. 
  samples <-  mutate(samples, rep = i, .before = X)
    
  # -- Create the abundance data for the composite samples -- #
  for(sp in 1:ncol(TabFull)) {  # For each species                                             
    Abund_Data[,sp]  <- tapply(TabFull[,sp], INDEX = sites$sample, FUN="sum") # Sum the abundance of each species in every 144 sampled composed by 9 plots
    }
  
  # Create a copy
  Pre_Abs_Data <- Abund_Data
  
  # Add the rep ID. 
  Abund_Data <- as.data.frame(Abund_Data) %>% mutate(rep = i, sample = samples$sample, .before = 1) # Abundance
  
  # Transform it into Presence_Absence Data
  Pre_Abs_Data[Pre_Abs_Data > 0] <- 1
  
  # Add the rep ID. 
  Pre_Abs_Data <- as.data.frame(Pre_Abs_Data) %>% mutate(rep = i, sample = samples$sample, .before = 1) # Pre_Abs data
  
  # Compute the species richness
  SR_Data <- apply(dplyr::select(Pre_Abs_Data,-rep,-sample),1,sum, na.rm = T)
  # Add the rep ID. 
  SR_Data <- as.data.frame(SR_Data) %>% mutate(rep = i, sample = samples$sample, .before = 1)  # Abundance
  
  # Create the output list
  Sampled_Data <- list("Abund_Data" = Abund_Data,
                       "Pre_Abs_Data" = Pre_Abs_Data,
                       "Sites_Data" = samples,
                       "SR_Data" = SR_Data)
  
} # End of the foreach loop on sampled data.

    # ----- #

#### . Merging and Saving . ####
  
# Merging all the data from each replicate and saving them. 
Abundance_Total  <- lapply(Sampled_Data, function(x) x[["Abund_Data"]]) %>% bind_rows() %>%
  write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Abundance.csv"))

Occurrence_Total <- lapply(Sampled_Data, function(x) x[["Pre_Abs_Data"]]) %>% bind_rows() %>% 
  write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Occurence.csv"))

Sites_Total  <- lapply(Sampled_Data, function(x) x[["Sites_Data"]]) %>% bind_rows() %>%   # Sites_Data are the same for each replicates but still saved with the replicate index for continuity. 
  write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Sites.csv"))
  
SR_Total  <- lapply(Sampled_Data, function(x) x[["SR_Data"]]) %>% bind_rows() %>% 
  write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_SR.csv"))
  
# --------------- #

cat(rule(left = "SCRIPT COMET_COMSIM.R ENDING", line_col = "red", line = "-", col = "br_red")) 

# --------------- #
