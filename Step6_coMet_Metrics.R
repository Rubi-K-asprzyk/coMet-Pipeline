#!/usr/bin/env Rscript

##### ______________________________________________________________________________________________________________________________________ #####

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_METRICS.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

# ---------------------------------------------------------------------- #
# coMet_Metrics : Computes a variety of pĥylogenetic diversity metrics.  #
# ---------------------------------------------------------------------- #

# This script is the main part of the coMet_pipeline and cannot be used without the outputs of the coMet_ComSim script.
# It is used to compute various diversity metrics from ecological community matrices created under coMet scenarios. 
# By default, it comes with 4 different computation possibilities to test for different sub-scenarios:
  # A = Intra-Core: Metrics are computed for each core of habitats.
  # B = Inter-Core: Metrics are computed between cores of habitats.
  # C = Inter-Ecotone: Metrics are computed between each ecotones of habitats.
  # D = Intra_Ecotone: Metrics are computed for each ecotones of habitats.

# As input, this script takes: 
  # 1: MANDATORY - The coMet_ConfigFile(s) of the scenario we want to analyze. 
  # This script will automatically navigate through the repository path to find all the needed files. 

# If we want to compute metrics for an empirical data set, please check for the individual code for each metric. 

# ------------------------------ #
##### STEP 0: INITIALISATION #####
# ------------------------------ #

cat(rule(left = "SCRIPT BEGINNING", line_col = "red", line = "-", col = "br_red"))

##### ______________________________________________________________________________________________________________________________________ #####

  #### . Packages . ####

cat(rule(left = "Loading packages and functions ...", line_col = "green", line = "-", col = "br_green"))

# Install/load pacman 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");require(pacman)})
# Install/load tons of packages
p_load("dplyr",
       "doParallel",
       "ape",
       "janitor",
       "data.table",
       "tictoc",
       "purrr",
       "rlist",
       "geiger",
       "picante",
       "gdata",
       "phylobase",
       "adiv",
       "spacodiR",
       "ade4",
       "vegan",
       "argparser")

cli_alert_success("Packages correctly loaded !")

# ----- #

  #### . Local Mode . ####

#Parameters <- list()
#Parameters$Config <- "coMet_ConfigFiles/Stochastic100_Sig5.R"
#Parameters$Rep <- 3

# ----- #

  #### . Argument Parser . ####

# Create the parser
arg_parser <- arg_parser("Check the communities computed by coMet_ComSim and create an output pdf.", name = "coMet_ComSim_Check.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = 1, help = "Configuration file for coMet_ComSim.R")

# Add the parameter to have the Array_Task_ID of the cluster.
arg_parser <- add_argument(arg_parser, arg = "--Rep", short = "-R", nargs = Inf, help = "$SLURM_ARRAY_TASK_ID")

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Extraction of the replicate number #
Rep <- Parameters$Rep

# Print a message to inform that the right configuration file was correctly loaded. 
cli::cli_alert_info(paste0("Working on configuration file: ",Parameters$Config," - TASK ID: ",Parameters$Rep))

# Load them into the environnement
source(Parameters$Config)

# Display some information about the config-file to be sure that all is OK.
cat(col_red("Scenario: "),Scenario,"\n")
cat(col_yellow("Number of replicates: "),Nrep,"\n")
cat(col_yellow("Number of samples: "),nbe_sample,"\n")
cat(col_yellow("Number of individuals per plot: "),Nip,"\n")
cat(col_yellow("Number of fragmented habitats: "),sq_x,"*",sq_y,"\n")
cat(col_yellow("Square size: "),sq_size,"\n")

# Print a message
cli_alert_success("coMet_ConfigFile.R correctly read !")

  #### . Creation of functions . ####

cli_alert_info("Creation of the needed functions ... ") ; cli_ul() # Opening a list container

# Call the functions needed to compute each of the metrics
source(paste0(getwd(),"/Utilities/PW_to_Vector.R"))   # Function to transform a PW square matrix into a vector with coordinates

# -- Create a function to create to columns based on an INDEX -- #   
    
Index <- function(Data){

    # Extract the index
    Index <- Data[,"Sample"]
    # Split the characters
    Sample_A <- as.integer(gsub("-[0-9]+$","",Index))
    Sample_B <- as.integer(gsub("^[0-9]+-","",Index))
    # Add those vectors to the intial data frame
  Data <- cbind("Sample" = Data[,"Sample"],Sample_A,Sample_B,Data[,- which(colnames(Data) == "Sample")])

  # Return Data
  return(Data)
  }

  cli_li("DONE: Index")
  
# -- Create a function to transform a pairwise matric into a vector shape-- #   
  
PW_to_Vector <- function(PW_matrix, Colname, Diag = T, Type = "Intra"){
     
     if (Type == "Intra"){  
       
       if (Diag == T) {Compound_new <- upper.tri(PW_matrix , diag = TRUE)}
       if (Diag == F) {Compound_new <- upper.tri(PW_matrix , diag = FALSE)}
       i2 <- which(Compound_new, arr.ind=TRUE)
       Compound_final <- data.frame(Sample = paste0(rownames(PW_matrix)[i2[,1]],"-",colnames(PW_matrix)[i2[,2]]),
                                    Sample_A = rownames(PW_matrix)[i2[,1]], 
                                    Sample_B = colnames(PW_matrix)[i2[,2]], 
                                    Colname = PW_matrix[Compound_new]) %>%
       stats::setNames(nm = c("Sample","Sample_A","Sample_B",Colname))
       
       return(Compound_final) } # End of type "Intra"
     
     if (Type == "Inter"){
       Compound_new <- matrix(data = TRUE, nrow = nrow(PW_matrix), ncol = ncol(PW_matrix))
       i2 <- which(Compound_new, arr.ind=TRUE)
       Compound_final <- data.frame(Sample = paste0(rownames(PW_matrix)[i2[,1]],"-",colnames(PW_matrix)[i2[,2]]),
                                    Sample_A = rownames(PW_matrix)[i2[,1]], 
                                    Sample_B = colnames(PW_matrix)[i2[,2]], 
                                    Name = PW_matrix[Compound_new]) %>%
       stats::setNames(nm = c("Sample","Sample_A","Sample_B",Colname))
     
       return(Compound_final)  } # End of type "Inter"
     
   } 
   
   cli_li("DONE: PW_to_Vector")
     
# -- Create of a function to add the MetaData for the metrics each time they are computed -- #
  
  # Metric is the Metric Computed. It could be either a vector of alpha values or a matrix of pairwise Beta Values
  # Type: Is the metric an Alpha or Beta metric ? 
  # Name is a character string used to name the output file (Example: PD)
  # Rep is the replicate we are currently working with
  # Diag is a boolean to determine if we want or not the intra plot values for Beta metrics
  
MetaData <- function(Metric, Type, Name, Rep, Diag = F){
    
  if (Type == "Alpha"){
      
    # -- Transformation / MetaData -- # 
      
    # Create a lite version of Sites_Data 
    Sites <- Sites_Data %>% 
      dplyr::filter(rep == Rep) %>% # Select the replicate
      dplyr::select(!c(X,Y)) # Remove the coordinates of the plot  
      
    Metric <- Metric %>%
      # Add the sample number
      tibble::rownames_to_column("Sample") %>%
      # Get the MetaData
      left_join(y = Sites, join_by("Sample" == "sample")) %>% dplyr::rename(Rep = rep, Habitat = hab, Spatial = SampleType) %>% # Get the MetaData for the Sample
      # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
      mutate(Spatial = case_when(Spatial == 1 ~ "C", Spatial == 2 ~ "B")) %>%
      mutate("Group" = paste0(Spatial,Habitat)) %>%
      # Transform what needs to be transformed into factors
      mutate(across(c(Rep,Sample,Group,Habitat,Spatial),as.factor)) %>% 
      # Reorder the columns in a more friendly way
      dplyr::select(Rep,Sample,Group,Habitat,Spatial,all_of(Name))
  
    } # End of Type : Alpha 
    
  if (Type == "Beta"){
      
    # -- Transformation / MetaData -- # 
      
    # Create a lite version of Sites_Data 
    Sites <- Sites_Data %>% 
      dplyr::filter(rep == Rep) %>% # Select the replicate
      dplyr::select(!c(X,Y)) # Remove the coordinates of the plot
      
    Metric <- PW_to_Vector(Metric, Colname = Name, Type = "Intra", Diag = Diag) %>% # Transform the pairwise matrix into a vector
        # Get the MetaData
        left_join(y = Sites, join_by("Sample_A" == "sample")) %>% dplyr::rename(Rep = rep, Habitat_A = hab, Spatial_A = SampleType) %>% # Get the MetaData for Sample_A
        left_join(y = Sites, join_by("Rep" == "rep","Sample_B" == "sample")) %>% dplyr::rename(Habitat_B = hab, Spatial_B = SampleType) %>% # Get the MetaData for Sample_B
        # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
        mutate(Spatial_A = case_when(Spatial_A == 1 ~ "C", Spatial_A == 2 ~ "B")) %>%
        mutate(Spatial_B = case_when(Spatial_B == 1 ~ "C", Spatial_B == 2 ~ "B")) %>% 
        # Create a last column "Group" that is a combination of Habitat and Spatial and "Type" for the intra comparisons (Same Habitat and Spatial)
        mutate("Group" = paste0(Spatial_A,Habitat_A,"_",Spatial_B,Habitat_B)) %>%
        # Transform what needs to be stransformed into factors
        mutate(across(c(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B),as.factor)) %>% 
        # Reorder the columns in a more friendly way
        dplyr::select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,all_of(Name))
      
    } # End of Type : Alpha 
 
  # Return the Results   
  return(Metric)
} 
  
  cli_li("DONE: MetaData ")

# Print a message
cli_end() ; cli_alert_success("Functions correctly created !")

# ------------------------------------------- #
##### STEP 1: LOADING THE COMMUNITY FILES #####
# ------------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "LOADING OF THE COMMUNITY FILES", line_col = "green", line = "-", col = "br_green"))

# ----- #

#### . Loading files . ####

cat(rule(left = "Loading the community files ...", line_col = "green", line = "-", col = "br_green"))

# -- COMMUNITIES -- #

# Community abundance data
Abundance_Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Abundance.csv"),row.names = 1) %>%
  mutate(across(c(rep,sample),as.factor))

# Community occurrence data
Pre_Abs_Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Occurence.csv"),row.names = 1) %>%
  mutate(across(c(rep,sample),as.factor))

# Sites_data
Sites_Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Sites.csv"),row.names = 1) %>%
  mutate_all(as.factor)

# Specific richness. 
Richness_Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_SR.csv"),row.names = 1) %>% 
  mutate(across(c(rep,sample),as.factor))

# Remove samples with species richness < Nmin
  # Find the samples
  To_remove <- which(Richness_Data$SR_Data < Nmin)
  
  # Remove the samples from Abundance, Presence-absence data and Richness Data if necessary
  if (is_empty(To_remove) == FALSE){
  Pre_Abs_Data <- Pre_Abs_Data[-To_remove,]
  Abundance_Data <- Abundance_Data[-To_remove,]
  Richness_Data <-  Richness_Data[-To_remove,] 
  # Print a message to inform how many communities were eventually removed. 
  cat(col_red("WARNING: "),length(To_remove)," communities of less than ",Nmin," species were removed out of ",nbe_sample*Nrep," total communities accross all replicates.","\n")
  } else {cat(col_green("OK: "),"No communities were less than ",Nmin," species.","\n")}
  
# -- TREES -- #

# Load the multiphylo file of observed trees
List_Phylo <-read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_PhyloTrees.tree"))
  
# Load the multiphylo files of null-model trees
List_Phylo_NULL <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TreeDistances/Scenario",Scenario,"_Step3_MultiPhylo_NullModel",Rep,".tree"))

# If they exist (i.e Randomization == 2), load the Observed Trimmed trees.
List_Phylo_Trimmed <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TreeDistances/Scenario",Scenario,"_Step3_ObservedTrimmed_NullModel",Rep,".tree"))

# Detect the size of the Null Model
NNM <- length(List_Phylo_NULL)

# if(length(List_Phylo) != Nrep){
#   stop("Number of Replicates of observed tree files different from the number of Replicates ; NREP = ",Nrep," - FILES = ",length(List_Phylo))
# } else if (length(List_Phylo_NULL) != Nrep){
#   stop("Number of Replicates of null model tree files different from the number of Replicates ; NREP = ",Nrep," - FILES = ",length(List_Phylo_NULL))
# } else {cli_li(cat(col_yellow("DONE"),": Tree data","\n"))}

cat(rule(left = "Loading the community files : DONE", line_col = "green", line = "-", col = "br_green"))

# --- Abundance Data ---- #

Sampled_Abund_Data_Inter <- Abundance_Data[Abundance_Data$rep == Rep,] %>%
    remove_empty(.,which = "cols")

# --- Occurrence Data ---- #

Sampled_Pre_Abs_Data_Inter <- Pre_Abs_Data[Pre_Abs_Data$rep == Rep,] %>%
    remove_empty(.,which = "cols") 

# --- Richness Data --- #

Sampled_Richness_Data <- Richness_Data[Richness_Data$rep == Rep,] %>%
    remove_empty(.,which = "cols") 

  # ---- #

# If Randomization == 2 (i.e tip-shuffling is only done on the species present in the communities), 
  # modify the Sampled_Data sets to be the one that are trimmed from the absent species. 

# Keep only the species present in the replicate / Trimmed tree if it exist
if(Randomization == 2){
  
  # -- Occurrence Data -- #  
  
  Sampled_Pre_Abs_Data_Inter <- Sampled_Pre_Abs_Data_Inter %>%
      # Keep only the species that are present in the trimmed tree
      dplyr::select("rep","sample",List_Phylo_Trimmed$tip.label)
  
  # -- Abundance Data -- # 
  
  Sampled_Abund_Data_Inter <- Sampled_Abund_Data_Inter %>%
      # Keep only the species that are present in the trimmed tree
      dplyr::select("rep","sample",List_Phylo_Trimmed$tip.label) 

  # -- Phylo Tree -- # 
  List_Phylo <- List_Phylo_Trimmed
  
}

# Change the rownames of each dataframes to be the sample numbers
rownames(Sampled_Abund_Data_Inter) <- Sampled_Abund_Data_Inter$sample
rownames(Sampled_Pre_Abs_Data_Inter) <- Sampled_Pre_Abs_Data_Inter$sample

# Display a message
cat(rule(left = "Script initialized", line_col = "green", line = "~", col = "br_green"))

# ----------------------------------- #
##### STEP 2: METRICS COMPUTATION #####
# ----------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "BEGINNING OF THE METRICS COMPUTATION", line_col = "green", line = "-", col = "br_green"))

  #### . Local Mode . ####
# 
# PD <- TRUE # Done
# FJ <- TRUE # Done
# PSV <- TRUE # Done
# PDab <- TRUE # Done
# MPD <- TRUE # Done
# MNTD <- TRUE # Done
# ED <- TRUE # Done
# PDb <- TRUE # Done
# UniPhylo <- TRUE # Done
# COMDIST <- FALSE # Done
# S <- TRUE # Done
# Dab <- TRUE # Done
# PIst <- TRUE # Done
# PCD <- TRUE # Done
# DISC <- FALSE # Done
# #
# Simpson <- TRUE

# -- PD (Faith 1992) and avPD (Clarke 2001) ----

# Do we want to compute PD and avPD ? 
if (PD == FALSE) {
  cli_alert_warning("PD and avPD not computed")
} else {

# "PD"   = "PHYLOGENETIC DIVERSITY"         = Sum of total branch lengths connecting species together
# "avPD" = "AVERAGE PHYLOGENETIC DIVERSITY" = Sum of total branch lengths connecting species together DIVIDED by species richness
  
# WARNING: 
# If the root is to be included in all calculations (include.root=TRUE), the PD of all samples will include the branch length connecting taxa in those samples and the root node of the supplied tree. The root of the supplied tree may not be spanned by any taxa in the sample.
# Single-species samples will be assigned a PD value equal to the distance from the root to the present.

cat(rule(left = "Computation of PD and avPD", line_col = "yellow", line = "-", col = "br_yellow"))
 
# --- OBSERVED --------------------------------------------------------------- #
tic() 

# OBSERVED: PD and avPD 
PD_Obs <- picante::pd(Sampled_Pre_Abs_Data_Inter[,-(1:2)], tree = List_Phylo, include.root = T) %>%
    # Add the avPD values
    mutate("avPD"= PD / SR) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = c("PD","avPD","SR"), Rep = Rep)
  
# OBSERVED avPD
avPD_Obs <- select(PD_Obs,-c(PD,SR))

# OBSERVED avPD
PD_Obs <- select(PD_Obs,-c(avPD,SR))

# Message
cat(col_yellow("DONE:"),"Observed PD and avPD", "\n") ; toc()

# ------------------------------- #

# --- NULL-MODEL ------------------------------------------------------------- #

# NULL: PD and avPD 
PD_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees
      
    # -- Metric Computation -- #  
    Metric <- picante::pd(Sampled_Pre_Abs_Data_Inter[,-(1:2)],tree = List_Phylo_NULL[[j]], include.root = T) %>%
      # Add the avPD values
      mutate("avPD"= PD / SR) %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("PD_NM",j) := PD) %>%
      dplyr::rename(!!paste0("avPD_NM",j) := avPD) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = c(paste0("PD_NM",j),paste0("avPD_NM",j),"SR"), Rep = Rep)
    }

# NULL: avPD
avPD_Null <- select(PD_Null,-c(starts_with("PD"),SR))

# NULL: PD
PD_Null <- select(PD_Null,-c(starts_with("avPD"),SR))

# Message
cat(col_yellow("DONE:"),"Null_Model PD and avPD","\n")
toc()

# --- COMBINE AND SAVE -------------------------------------------------------------- #

PD_Total <- full_join(PD_Obs, PD_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
  write.csv(PD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PD/Scenario",Scenario,"_PD_Faith_",Rep,".csv"))

avPD_Total <- full_join(avPD_Obs, avPD_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
  write.csv(avPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PD/Scenario",Scenario,"_avPD_Faith_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: PD and avPD", line_col = "yellow", line = "-", col = "br_yellow")) 
# --- #

} # End of PD and avPD

# --------------------------------------------------------------------------------------------------------------------------------------------- #
          
# -- F (Izsak 2000) and J (Izsak 2000) ----

# Do we want to compute F and J ? 
if (FJ == FALSE) {
  cli_alert_warning("F and J not computed")
} else {
 
# Load the needed function
source(paste0(getwd(),"/Utilities/Diversity_Metrics/F_extensive_quadratic_entropy.R"))  # Function to compute the F quadratic entropy

  # "F" = "EXTENSIVE QUADRATIC ENTROPY"   = Sum of pairwise distance
  # "J" = "INTENSIVE QUADRATIC ENTROPY"   = Average distance between two randomly chosen species from a community
    # "Q" is not computed because it is very similar to MPD
    
cat(rule(left = "F and J", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic() 

# Load the observed pairwise distances
Obs_PW <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step4_Intra_PW_Distances_Observed",Rep,".csv"), row.names = 1)


 # Create the NULL_PW by using the Obs_PW_Trimmed and the tips_label from the list_Phylo_NULL trees. 
Null_PW <- 
    foreach(NM = 1:length(List_Phylo_NULL)) %do% {
    # Make a copy
    Null <- Obs_PW
    # Extract the species found in the tip-labels of the Trimmed Trees (That are, therefore species that are at least present one time in the dataframe)
    Null_Sp <- List_Phylo$tip.label
    # Un-triangle the matrix to make the function work
    upperTriangle(Null[-(1:2)] ) <- lowerTriangle(Null[-(1:2)] , byrow=TRUE)  
    # Filter and Select NULL for the species present in the tip.labels (that can be low if randomization == 2)
    Null <- Null[Null_Sp,c("Sp","tree_NM",Null_Sp)]
    # Randomize the species with the random trees tips.labels
    Null$Sp <-  List_Phylo_NULL[[NM]]$tip.label
    colnames(Null)[-(1:2)] <- List_Phylo_NULL[[NM]]$tip.label
    rownames(Null) <- List_Phylo_NULL[[NM]]$tip.label
    return(Null)
  }

# OBSERVED: F and J
F_Obs <- Ex_quad_entropy(Occurence_data = Sampled_Pre_Abs_Data_Inter[,-(1:2)], 
                            Phylo_tree = List_Phylo, 
                            PW = Obs_PW[,-(1:2)], FUN = sum ) %>% 
    as.data.frame() %>%
    set_names("F_") %>%
  
  # Add the J_ values
  mutate("J_" = F_ / Sampled_Richness_Data$SR_Data^2) %>%
  # Add the MetaData
  MetaData(Type = "Alpha", Name = c("F_","J_"), Rep = Rep)
  

# OBSERVED J_
J_Obs <- select(F_Obs,-F_)

# OBSERVED F_
F_Obs <- select(F_Obs,-J_)

# Message
cat(col_yellow("DONE:"),"Observed F and J", "\n") ; toc()      
 
# --- NULL-MODEL ------------------------------------------------------------- #
tic()  

# NULL: F and J 
F_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees
    
    # -- Metric Computation -- #  
    Metric <- Ex_quad_entropy(Occurence_data = Sampled_Pre_Abs_Data_Inter[,-(1:2)], 
                              Phylo_tree = List_Phylo_NULL[[j]], 
                              PW = Null_PW[[j]][,-(1:2)], FUN = sum ) %>% as.data.frame() %>% set_names("F_") %>%
      # Add the J_ values
      mutate("J_" = F_ / Sampled_Richness_Data$SR_Data^2) %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("F_NM",j) := F_) %>%
      dplyr::rename(!!paste0("J_NM",j) := J_) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = c(paste0("F_NM",j),paste0("J_NM",j)), Rep = Rep)
  }

# NULL: J
J_Null <- select(F_Null,-c(starts_with("F_")))

# NULL: F
F_Null <- select(F_Null,-c(starts_with("J_")))

# Message
cat(col_yellow("DONE:"),"Null_Model F and J","\n")
toc()

# --- COMBINE AND SAVE -------------------------------------------------------------- #

F_Total <- full_join(F_Obs, F_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(F_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/F/Scenario",Scenario,"_F_Izsak_",Rep,".csv"))

J_Total <- full_join(J_Obs, J_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(J_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/F/Scenario",Scenario,"_J_Izsak_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: F and J", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of F and J

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# -- PSV, PSR, and PSE (Helmus, 2007) ----

# Do we want to compute PSV, PSR and PSE ? 
if (PSV == FALSE) {
  cli_alert_warning("PSV, PSR and PSE not computed")
} else {
 
  # "PSV" = "PHYLOGENETIC SPECIES VARIABILITY"  = Variability in an unmeasured neutral trait or the relative amount of unshared branch length.
  # "PSR" = "PHYLOGENETIC SPECIES RICHNESS"     = Variability in an unmeasured neutral trait multiplied by species richness
  # "PSE" = "PHYLOGENETIC SPECIES EVENESS"      = Abundance weighted PSV
          
  cat(rule(left = "PSV, PSR and PSE", line_col = "yellow", line = "-", col = "br_yellow"))
    
# --- OBSERVED --------------------------------------------------------------- #
tic()     

# OBSERVED: PSV
PSV_Obs <- picante::psv(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo, compute.var = F) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "PSVs", Rep = Rep)
  
# OBSERVED: PSR
PSR_Obs <- picante::psr(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo, compute.var = F) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "PSR", Rep = Rep)
  
# OBSERVED: PSE
PSE_Obs <- picante::pse(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "PSEs", Rep = Rep)
  
# Message
cat(col_yellow("DONE:"),"Observed PSV, PSR and PSE", "\n") ; toc()   
      
# --- NULL-MODEL ------------------------------------------------------------- #

tic()  

# PSV: Null_Model
PSV_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees

    # -- Metric Computation -- #  
    Metric <- picante::psv(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo_NULL[[j]], compute.var = F) %>%
      # Set the name of the Metric
      rename(!!paste0("PSV_NM",j) := PSVs) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("PSV_NM",j), Rep = Rep)
}

# PSR: Null_Model
PSR_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees
    
    # -- Metric Computation -- #  
    Metric <- picante::psr(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo_NULL[[j]], compute.var = F) %>%
      # Set the name of the Metric
      rename(!!paste0("PSR_NM",j) := PSR) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("PSR_NM",j), Rep = Rep)
} 

# PSE: Null_Model
PSE_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees
    
    # -- Metric Computation -- #  
    Metric <- picante::pse(Sampled_Abund_Data_Inter[,-(1:2)], tree = List_Phylo_NULL[[j]]) %>%
      # Set the name of the Metric
      rename(!!paste0("PSE_NM",j) := PSEs) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("PSE_NM",j), Rep = Rep)
} 

# Message
cat(col_yellow("DONE:"),"Null_Model PSV, PSR and PSE","\n")
toc()

# --- COMBINE AND SAVE ------------------------------------------------------- #

PSV_Total <- full_join(PSV_Obs, PSV_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
  write.csv(PSV_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV/Scenario",Scenario,"_PSV_Helmus_",Rep,".csv"))
  
PSR_Total <- full_join(PSR_Obs, PSR_Null, by = join_by(Rep, Sample, Group,  Habitat, Spatial)) 
  write.csv(PSR_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV/Scenario",Scenario,"_PSR_Helmus_",Rep,".csv"))
  
PSE_Total <- full_join(PSE_Obs, PSE_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
  write.csv(PSE_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV/Scenario",Scenario,"_PSE_Helmus_",Rep,".csv"))

 # --- #
cat(rule(left = "DONE: PSV, PSR and PSE", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of PSV, PSR and PSE

# --------------------------------------------------------------------------------------------------------------------------------------------- #
  
# -- PDab (Vellend 2010), ΔnPD (Barker 2002) and avPDab (Tucker 2017) ----
 
if (PDab == FALSE) {
  cli_alert_warning("PDab, DeltaNPD and avPDab not computed")
} else {
  
# "PDab"   = "ABUNDANCE WEIGHTED PHYLOGENETIC DIVERSITY ss VELLEND" = sum of total branch length connecting all the species together where branches are scaled by proportional abundances of subtending species
# "ΔnPD"   = "ABUNDANCE WEIGHTED PHYLOGENETIC DIVERSITY ss BARKER" = sum of total branch length connecting all the species together where branches are scaled by abundances of subtending species
# "avPDab"   = " AVERAGE ABUNDANCE WEIGHTED PHYLOGENETIC DIVERSITY" = sum of total branch length connecting all the species together where branches are scaled by proportional abundances of subtending species, divided by the number of species 
  
cat(rule(left = "PDab, Delta_n_PD and avPDab", line_col = "yellow", line = "-", col = "br_yellow"))

# Load the needed function  
source(paste0(getwd(),"/Utilities/Diversity_Metrics/Weighted_Faith-Swenson.R")) # Function to compute PDab values ( Nathan Swenson )
  
# --- OBSERVED --------------------------------------------------------------- #
tic()     
  

# OBSERVED: PDab
PDab_Obs <- weighted.faith(my.phylo = List_Phylo, my.sample = Sampled_Abund_Data_Inter[,-(1:2)], type = "V" ) %>%
    # Add the avPDab values
    mutate("avPDab" = PDab / Sampled_Richness_Data$SR_Data^2) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = c("PDab","avPDab"), Rep = Rep) 
  
# OBSERVED avPDab
avPDab_Obs <- select(PDab_Obs,-PDab)

# OBSERVED PDab
PDab_Obs <- select(PDab_Obs,-avPDab)

# OBSERVED: DnPD
DnPD_Obs <- weighted.faith(my.phylo = List_Phylo, my.sample = Sampled_Abund_Data_Inter[,-(1:2)], type = "B" ) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "DnPD", Rep = Rep)
  
# Message
cat(col_yellow("DONE:"),"Observed PDab, DnPD and avPDab", "\n") ; toc()
  
# --- NULL-MODEL ------------------------------------------------------------- #
tic()   

# NULL: PDab and avPDab
PDab_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees
    
    # -- Metric Computation -- #  
    Metric <- weighted.faith(Sampled_Abund_Data_Inter[,-(1:2)], my.phylo = List_Phylo_NULL[[j]], type = "V" ) %>%
      # Add the avPDab values
      mutate("avPDab" = PDab / Sampled_Richness_Data$SR_Data^2) %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("PDab_NM",j) := PDab) %>%
      dplyr::rename(!!paste0("avPDab_NM",j) := avPDab) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = c(paste0("PDab_NM",j),paste0("avPDab_NM",j)), Rep = Rep)
} 

# NULL: avPDab
avPDab_Null <- select(PDab_Null,-c(starts_with("PDab")))

# NULL: PDab
PDab_Null <- select(PDab_Null,-c(starts_with("avPDab")))

# NULL: DnPD
DnPD_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- weighted.faith(Sampled_Abund_Data_Inter[,-(1:2)], my.phylo = List_Phylo_NULL[[j]], type = "B" ) %>%
      # Rename 
      dplyr::rename(!!paste0("DnPD_NM",j) := DnPD) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("DnPD_NM",j), Rep = Rep)
  }
  
# Message
cat(col_yellow("DONE:"),"Null_Model PDab, DnPD and avPDab","\n")
toc()

# --- COMBINE AND SAVE ------------------------------------------------------- #

PDab_Total <- full_join(PDab_Obs, PDab_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(PDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab/Scenario",Scenario,"_PDab_Vellend_",Rep,".csv"))

DnPD_Total <- full_join(DnPD_Obs, DnPD_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(DnPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab/Scenario",Scenario,"_DnPD_Barker_",Rep,".csv"))

avPDab_Total <- full_join(avPDab_Obs, avPDab_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(avPDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab/Scenario",Scenario,"_avPDab_Tucker_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: PDab, DnPD and avPDab", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of PDab, DnPD and avPDab 

# --------------------------------------------------------------------------------------------------------------------------------------------- #  

# -- MPD, MPDab and InterMPDab (Clarke & Warwick 1998) ----

if (MPD == FALSE) {
  cli_alert_warning("MPD, MPDab and InterMPDab not computed")
} else {
  
# Load the needed function
source(paste0(getwd(),"/Utilities/Diversity_Metrics/MPD_picante.R"))  # Function to compute the F quadratic entropy

  # "Δp"   =  Mean phylogenetic distance between distinct species
  # "MPD         = "MEAN PAIRWISE DISTANCES" = Mean phylogenetic distances between all pair of species occurring within a site
    # For now, only MPD is set up because ( problem with the weighted mean thingy in the custom fonction)
  # "MPDab"      = "ABUNDANCE WEIGHTED MEAN PAIRWISE DISTANCES" = Abundance weighted mean phylogenetic distances between all pair of species occurring within a site / Mean phylogenetic distance between INDIVIDUALS
  # "InterMPDab" = "INTERSPECIFIC ABUNDANCE WEIGHTED MEAN PAIRWISE DISTANCES" = Abundance weighted mean phylogenetic distances between all pair of species occurring within a site / Mean phylogenetic distance between INDIVIDUALS from different species

cat(rule(left = "MPD, MPDab and InterMPDab", line_col = "yellow", line = "-", col = "br_yellow"))

# --- OBSERVED --------------------------------------------------------------- #
tic()

# Load the observed pairwise distances
Obs_PW <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TreeDistances/Scenario",Scenario,"_Step4_Intra_PW_Distances_Observed",Rep,".csv"), row.names = 1)

# Create the NULL_PW by using the Obs_PW_Trimmed and the tips_label from the list_Phylo_NULL trees. 
Null_PW <-   
  foreach(NM = 1:length(List_Phylo_NULL)) %do% {
      # Make a copy
      Null <- Obs_PW
      # Extract the species found in the tip-labels of the Trimmed Trees (That are, therefore species that are at least present one time in the dataframe)
      Null_Sp <- List_Phylo$tip.label
      # Un-triangle the matrix to make the function work
      upperTriangle(Null[-(1:2)] ) <- lowerTriangle(Null[-(1:2)] , byrow=TRUE)  
      # Filter and Select NULL for the species present in the tip.labels (that can be low if randomization == 2)
      Null <- Null[Null_Sp,c("Sp","tree_NM",Null_Sp)]
      # Randomize the species with the random trees tips.labels
      Null$Sp <-  List_Phylo_NULL[[NM]]$tip.label
      colnames(Null)[-(1:2)] <- List_Phylo_NULL[[NM]]$tip.label
      rownames(Null) <- List_Phylo_NULL[[NM]]$tip.label
      return(Null)
    }

# OBSERVED: MPD

# WARNING ! MPD = MPD intra
MPD_Obs <- mpd_V2(Sample = Sampled_Abund_Data_Inter[,-c(1:2)], PW_distances =  Obs_PW[,-(1:2)], abundance.weighted = F, intra = T) %>% 
    as.data.frame() %>%
    # Set the name of the Metric
    dplyr::rename(MPD = ".") %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "MPD", Rep = Rep)
  
# OBSERVED: DeltaP

# WARNING ! DP = MPD 
Dp_Obs <- mpd_V2(Sample = Sampled_Abund_Data_Inter[,-c(1:2)], PW_distances =  Obs_PW[,-(1:2)], abundance.weighted = F, intra = F) %>% 
    as.data.frame() %>%
    # Set the name of the Metric
    dplyr::rename(Dp = ".") %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "Dp", Rep = Rep)
  
# Message
cat(col_yellow("DONE:"),"MPD and Dp", "\n") ; toc()
  
# --- NULL-MODEL ------------------------------------------------------------- #
tic()  

# NULL: MPD 
MPD_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- mpd_V2(Sample = Sampled_Abund_Data_Inter[,-(1:2)],
                     PW_distances = Null_PW[[j]][,-(1:2)],abundance.weighted = F, intra = T)  %>% 
      as.data.frame() %>% 
      # Set the name of the Metric
      dplyr::rename(!!paste0("MPD_NM",j) := .) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("MPD_NM",j), Rep = Rep)
  } 

# NULL: MPD 
Dp_Null <- 
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- mpd_V2(Sample = Sampled_Abund_Data_Inter[,-(1:2)],
                     PW_distances = Null_PW[[j]][,-(1:2)],abundance.weighted = F, intra = F)  %>% 
      as.data.frame() %>% 
      # Set the name of the Metric
      dplyr::rename(!!paste0("Dp_NM",j) := .) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("Dp_NM",j), Rep = Rep)
  } 

# Message
cat(col_yellow("DONE:"),"Null_Model MPD and Dp","\n")
toc()
  
# --- COMBINE AND SAVE ------------------------------------------------------- #

MPD_Total <- full_join(MPD_Obs, MPD_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
# write.csv(MPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MPD/Scenario",Scenario,"_MPD_Clarke_",Rep,".csv"))

Dp_Total <- full_join(Dp_Obs, Dp_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
# write.csv(Dp_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MPD/Scenario",Scenario,"_Dp_Hardy_",Rep,".csv"))
  
# TEST - TEST - TEST #

# Remove any unwanted columns
  # Raw_Data <- Dp_Total %>%
  #   select(-any_of(c("Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial")))
  # 
  # # Extract the Metric Name
  # Name <- colnames(Raw_Data)[4]
  # 
  # # Transform the data to add a column for the null vs obs data
  # Raw_Trans <- Raw_Data %>%
  #   pivot_longer(
  #     cols = contains(Name, ignore.case = TRUE, vars = NULL),
  #     names_to = "OBSvsNM", values_to = "Value")
  # 
  # # Change the values contained in the new column
  # Raw_Trans$OBSvsNM[grep("NM",Raw_Trans$OBSvsNM)] <- paste0(Name,"_NM")
  # 
  # # Transform all the stuff into factor
  # Raw_Trans <-  Raw_Trans %>%
  #   mutate_at(vars(Rep,Sample,Group,), factor) %>%
  #   as.data.frame()
  # 
  # Raw_Data <- Raw_Data %>%
  #   mutate_at(vars(Rep,Sample,Group,), factor) %>%
  #   as.data.frame()
  # 
  # # ----- #
  # 
  # #### . Observed VS Null Model . ###
  # 
  # 
  # # Compute summary statistics
  # Sum_Raw <- Raw_Trans %>%
  #   group_by(OBSvsNM) %>%  # Group_by to have two groups in a tibble for faceting
  #   get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups
  #   # Draw the summary statistics
  #   ggsummarytable(x = "OBSvsNM",           # Sub scenario
  #                  y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
  #                  digits = 4,                        # Number of digits
  #                  size = 4,                          # Size of the text
  #                  ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot()
  # 
  # # Create the density plot
  # Dens_Raw <- Raw_Trans %>%
  #   # Plotting
  #   ggplot(aes(x = Value, group = OBSvsNM, color = OBSvsNM, fill = OBSvsNM)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
  #   geom_density(alpha = 0.3) +
  #   theme(legend.position= "right",
  #         axis.title.x=element_blank(),
  #         legend.title =element_blank())


# --- #
cat(rule(left = "DONE: MPD and Dp", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of MPD

# --------------------------------------------------------------------------------------------------------------------------------------------- #  
  
# -- MNTD and MNTDab (Webb et al, 2002,2008) ----

if (MNTD == FALSE) {
  cli_alert_warning("MNTD and MNTDab not computed")
} else {
  
#   # "MNTD"    =  MEAN NEAREST TAXON DISTANCE = Mean phylogenetic distance to the nearest relative for all species occurring within a site
#   # "MNTDab"  =  ABUNDANCE WEIGHTED MEAN NEAREST TAXON DISTANCE = Mean phylogenetic distance to the nearest relative for all species occurring within a site weighted by abundance  
  
cat(rule(left = "MNTD and MNTDab", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()
  
# Load the observed pairwise distances
Obs_PW <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step4_Intra_PW_Distances_Observed",Rep,".csv"), row.names = 1)

# Create the NULL_PW by using the Obs_PW_Trimmed and the tips_label from the list_Phylo_NULL trees. 
Null_PW <-   
  foreach(NM = 1:length(List_Phylo_NULL)) %do% {
    # Make a copy
    Null <- Obs_PW
    # Extract the species found in the tip-labels of the Trimmed Trees (That are, therefore species that are at least present one time in the dataframe)
    Null_Sp <- List_Phylo$tip.label
    # Un-triangle the matrix to make the function work
    upperTriangle(Null[-(1:2)] ) <- lowerTriangle(Null[-(1:2)] , byrow=TRUE)  
    # Filter and Select NULL for the species present in the tip.labels (that can be low if randomization == 2)
    Null <- Null[Null_Sp,c("Sp","tree_NM",Null_Sp)]
    # Randomize the species with the random trees tips.labels
    Null$Sp <-  List_Phylo_NULL[[NM]]$tip.label
    colnames(Null)[-(1:2)] <- List_Phylo_NULL[[NM]]$tip.label
    rownames(Null) <- List_Phylo_NULL[[NM]]$tip.label
    return(Null)
  }

# Make the pairwise matrix an inter-specific pairwise distance matrix.
Meta <- Obs_PW[,1:2]
PW <- Obs_PW[,-c(1:2)]
PW[upper.tri(PW,diag = T)] <- NA
# Un-triangle the matrix to make the function work
upperTriangle(PW) <- lowerTriangle(PW, byrow=TRUE)
Obs_PW_Inter  <- cbind(Meta,PW)

Null_PW_Inter <- lapply(seq_along(Null_PW), function(y){
    Meta <- Null_PW[[y]][,1:2]
    PW <- Null_PW[[y]][,-c(1:2)]
    PW[upper.tri(PW,diag = T)] <- NA
    # Un-triangle the matrix to make the function work
    upperTriangle(PW) <- lowerTriangle(PW, byrow=TRUE)
    PW <- cbind(Meta,PW)
  })

# OBSERVED: MNTD 
MNTD_Obs <-  mntd(samp = Sampled_Abund_Data_Inter[,-c(1:2)], dis = Obs_PW_Inter[,-(1:2)], abundance.weighted = F)  %>%
    as.data.frame() %>%
    # Set the name of the Metric
    dplyr::rename(MNTD = ".") %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "MNTD", Rep = Rep)
  
# OBSERVED: MNTDab 
MNTDab_Obs <- mntd(samp = Sampled_Abund_Data_Inter[,-c(1:2)], dis = Obs_PW_Inter[,-(1:2)], abundance.weighted = T)  %>%
    as.data.frame() %>%
    # Set the name of the Metric
    dplyr::rename(MNTDab = ".") %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = "MNTDab", Rep = Rep)

# Message
cat(col_yellow("DONE:"),"MNTD and MNTDab observed", "\n") ; toc() 

# --- NULL-MODEL ------------------------------------------------------------- #
tic()  

# NULL: MNTD
MNTD_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- mntd(samp = Sampled_Abund_Data_Inter[,-c(1:2)], dis = Null_PW_Inter[[j]][,-(1:2)], abundance.weighted = F) %>%
      as.data.frame() %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("MNTD_NM",j) := .) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("MNTD_NM",j), Rep = Rep)
  }

# NULL: MNTDab
MNTDab_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- mntd(samp = Sampled_Abund_Data_Inter[,-c(1:2)], dis = Null_PW_Inter[[j]][,-(1:2)], abundance.weighted = T) %>%
      as.data.frame() %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("MNTDab_NM",j) := .) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = paste0("MNTDab_NM",j), Rep = Rep)
  }

# Message
cat(col_yellow("DONE:"),"Null_Model MNTD and MNTDab","\n")
toc()

# --- COMBINE AND SAVE ------------------------------------------------------- #

MNTD_Total <- full_join(MNTD_Obs, MNTD_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(MNTD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MNTD/Scenario",Scenario,"_MNTD_Webb_",Rep,".csv"))

MNTDab_Total <- full_join(MNTDab_Obs, MNTDab_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(MNTDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MNTD/Scenario",Scenario,"_MNTDab_Webb_",Rep,".csv"))

# ---------------------------- #

# --- #
cat(rule(left = "DONE: MNTD and MNTDab", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of MNTD and MNTDab

# --------------------------------------------------------------------------------------------------------------------------------------------- #    
  
# -- ED, AED and MED (Tucker 2017) ----

if (ED == FALSE) {
  cli_alert_warning("ED, AED and MED not computed")
} else {
  
  # "ED"  = "SUMMED EVOLUTIONARY DISTINCTIVENESS"   = Sum of species evolutionary distinctiveness
  # "AED" = "SUMMED ABUNDANCE WEIGHTED EVOLUTIONARY DISTINCTIVENESS"   = Sum of species evolutionary distinctiveness weighted by species abundances    
  # "MED" = "MEAN EVOLUTIONNARY DISTINCTIVENESS" = Mean species evolutionary distinctiveness divided by species richness
  
cat(rule(left = "ED, AED and MED", line_col = "yellow", line = "-", col = "br_yellow"))
 
# library(ecoPD) is not supported anymore, so I copy-paste the source code in a script
source(paste0(getwd(),"/Utilities/Diversity_Metrics/AED_ecoPD.R"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()

# Compute the Evolutionary distinctiveness and Abundance weighted evolutionary distinctiveness. 

# Total_ED: The Evolutionary distinctiveness for all the species. 
Total_ED <- picante::evol.distinct(tree = List_Phylo, type = "fair.proportion")

# Here; we have the ED computed for all species, we need to retrieve the data for the species present in each samples
# Sample_ED: The Evolutionary distinctiveness for the species present in each sample. 
Sample_ED <- apply(Sampled_Abund_Data_Inter[,-c(1:2)],1,function(z){
      # Find the species present in the plot
      Sp.names <- names(z[which(z > 0)])
      # Find the ED
      ED <- Total_ED[which(Total_ED$Species %in% Sp.names),2]
    })
  
# OBSERVED: ED and MED
ED_Obs <-  sapply(Sample_ED,sum) %>%
    as.data.frame() %>%
    # Set the name of the Metric
    dplyr::rename(ED = ".") %>%
    # Mutate the MED
    mutate("MED" = sapply(Sample_ED,mean)) %>%
    # Add the MetaData
    MetaData(Type = "Alpha", Name = c("ED","MED"), Rep = Rep)
  

# OBSERVED MED
MED_Obs <- select(ED_Obs,-ED)

# OBSERVED ED
ED_Obs <- select(ED_Obs,-MED)

# OBSERVED: AED

# Add the abundance community data
AED_format <- phylo4d(x = List_Phylo,tip.data = t(Sampled_Abund_Data_Inter[,-c(1:2)])) # Species needed to be in row
  
# -- Metric Computation -- #  
Metric <- aed(AED_format) %>%
  sapply(sum) %>%
  as.data.frame() %>%
  # Set the name of the Metric
  rename(AED = ".")
  
  # Change the rownames
  rownames(Metric) <- Sampled_Abund_Data_Inter$sample 
    
# Add the MetaData
AED_Obs <- Metric %>% MetaData(Type = "Alpha", Name = "AED", Rep = Rep)
  
# Message
cat(col_yellow("DONE:"),"ED, AED and MED observed", "\n") ; toc() 
  
# --- NULL-MODEL ------------------------------------------------------------- #
tic()  
  
# ED : Null_Model values
  # Total_ED: The Evolutionary distinctiveness for all the species. 

Total_ED_NM <- foreach(NM = 1:length(List_Phylo_NULL)) %do% {
  # Change the colnames and rownames of Obs_PW[[x]] with the randomized rownames and colnames of List_Phylo_NULL
  Null <- Total_ED
  Null$Species <- List_Phylo_NULL[[NM]]$tip.label
  rename_with(Null, .fn = ~ paste0("ED_NM",NM), .cols = w)
  return(Null)
}

# Sample_ED_NM: The Evolutionary distinctiveness for the species present in each sample. 
Sample_ED_NM <- foreach(j = 1:NNM) %do% { # Loop over the Null trees 
    # Apply a function on each of the rows (each sample)
    ED <- apply(Sampled_Abund_Data_Inter[,-c(1:2)],1,function(z){
    # Find the species present in the plot
    Sp.names <- names(z[which(z > 0)])
    # Find the ED
    ED <- Total_ED_NM[[j]][which(Total_ED_NM[[j]]$Species %in% Sp.names),2]
  })
}    


# NULL: ED and MED
ED_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # -- Metric Computation -- #  
    Metric <- sapply(Sample_ED_NM[[j]],sum) %>%
      as.data.frame() %>%
      # Add the J_ values
      mutate("MED" = sapply(Sample_ED_NM[[j]],mean)) %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("ED_NM",j) := .) %>%
      dplyr::rename(!!paste0("MED_NM",j) := MED) %>%
      # Add the MetaData
      MetaData(Type = "Alpha", Name = c(paste0("ED_NM",j),paste0("MED_NM",j)), Rep = Rep)
}

# NULL: MED
MED_Null <- select(ED_Null,-c(starts_with("ED_")))

# NULL: ED
ED_Null <- select(ED_Null,-c(starts_with("MED_")))

# NULL: AED
AED_Null <-
  foreach(j = 1:NNM, .combine = merge) %do% { # Loop over the Null trees 
    
    # The function uses the rownames as an index number, therefore, change the column "sample" to the rownames. 
    AED_format <- phylo4d(x = List_Phylo_NULL[[j]], t(Sampled_Abund_Data_Inter[,-c(1:2)])) # Species needed to be in rows
    
    # -- Metric Computation -- #  
    Metric <- aed(AED_format) %>%
      sapply(sum) %>%
      as.data.frame() %>%
      # Set the name of the Metric
      dplyr::rename(!!paste0("AED_NM",j) := .)
    
    # Change the rownames
    rownames(Metric) <- Sampled_Abund_Data_Inter$sample
    
    # Add the MetaData
    Metric <- Metric %>% MetaData(Type = "Alpha", Name = paste0("AED_NM",j), Rep = Rep)
  }

# --- COMBINE AND SAVE ------------------------------------------------------- #

ED_Total <- full_join(ED_Obs, ED_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(ED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED/Scenario",Scenario,"_ED_Tucker_",Rep,".csv"))

MED_Total <- full_join(MED_Obs, MED_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(MED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED/Scenario",Scenario,"_MED_Tucker_",Rep,".csv"))

AED_Total <- full_join(AED_Obs, AED_Null, by = join_by(Rep, Sample, Group, Habitat, Spatial)) 
write.csv(AED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED/Scenario",Scenario,"_AED_Tucker_",Rep,".csv"))

# ---------------------------- #

# --- #
cat(rule(left = "DONE: ED, AED and MED", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of ED, AED and MED

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# -- Faith's PDbeta (Nipperess, 2010) ----

# Do we want to compute the PDbeta (as well as avPD) ? 
if (PDb == FALSE) {
  cli_alert_warning("PD_Beta not computed")
} else {
  
  # "PDbeta" = BETA DIVERSITY OF FAITH = Proportional shared branch length between two communities (Or percentage of the tree that is covered when we link all the species in the assemblage together)   
  
# Source the needed function
source(paste0(getwd(),"/Utilities/Diversity_Metrics/Phylosim_Nipperess.R"))
  
cat(rule(left = "Computation of PDBeta", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()
  
# OBSERVED: PDbeta
PDBeta_Obs <- as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo, incidence = T, method = "faith"))) %>%
    MetaData(Type = "Beta", Name = "PDb", Rep = Rep, Diag = F)
  
# Message
cat(col_yellow("DONE:"),"PD_Beta observed", "\n") ; toc() 
  
# --- RANDOM --------------------------------------------------------------------------- #  
tic()
  
# PDBeta_Null: Observed values
PDBeta_Null <- foreach(j = 1:NNM,.combine = merge) %do% { # Loop over the Null trees( "Merge" is used instead of "full join" because it is silent )

  # -- Metric Computation -- #  
    Metric <- as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo_NULL[[j]], incidence = T, method = "faith"))) %>%
    MetaData(Type = "Beta", Name = paste0("PDb_NM",j), Rep = Rep, Diag = F)  
  
  }
  
# Message
cat(col_yellow("DONE:"),"Null_Model PD_Beta","\n")
toc()  

# --- COMBINE AND SAVE ------------------------------------------------------- #
  
# Join the data
PDBeta_Total <- full_join(PDBeta_Obs,PDBeta_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B)) 
  write.csv(PDBeta_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDb/Scenario",Scenario,"_PDbeta_Nipperess_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: PD_Beta", line_col = "yellow", line = "-", col = "br_yellow")) 
  
} # End of PDbeta

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# -- PHYLOSOR and UNIFRAC (Leprieur 2012) ----

# Do we want to compute the PD (as well as avPD) ? 
if (UniPhylo == FALSE) {
  cli_alert_warning("PhyloSor and UniFrac not computed")
} else {
  
  # "PHYLOSOR"     = PHYLOGENETIC SORENSEN INDEX                       = Proportional shared branch length between two communities 
  # "PHYLOSORab"   = PHYLOGENETIC SORENSEN INDEX WEIGHTED BY ABUNDANCE = Proportional shared branch length between two communities weighted by abundance
  # "PHYLOSORturn" = TURNOVER COMPONENT OF PHYLOSOR                    = Proportional shared branch length between two communities explained by turnover
  
  # /!\ WARNING /!\ = We remove the Intra-Plot values of PHYLOSORab and UNIFRACab
  
# Source the needed function
source(paste0(getwd(),"/Utilities/Diversity_Metrics/Phylosim_Nipperess.R"))
source(paste0(getwd(),"/Utilities/Diversity_Metrics/Unifrac_Phylosor_Leprieur_2012.R")) # Function to compute Unifrac and Phylosor
  
cat(rule(left = "Phylosor and Unifrac (and related)", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()

# Create a lite version of Sites_Data 
Sites <- Sites_Data %>% 
  filter(rep == Rep) %>% # Select the replicate
  select(!c(X,Y)) # Remove the coordinates of the plot

# PhyloFroc_Obs: Observed values of Phylosor and Unifrac
PhyloFroc_Obs <- as.data.frame(beta.pd.decompo(com = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], tree = List_Phylo)) %>%
      tibble::rownames_to_column("Sample") %>% 
      Index(.) %>% 
      mutate(across(starts_with("Sample"), as.factor)) %>%
      # Get the MetaData
      left_join(y = Sites, join_by("Sample_A" == "sample")) %>% rename(Rep = rep, Habitat_A = hab, Spatial_A = SampleType) %>% # Get the MetaData for Sample_A
      left_join(y = Sites, join_by("Rep" == "rep","Sample_B" == "sample")) %>% rename(Habitat_B = hab, Spatial_B = SampleType) %>% # Get the MetaData for Sample_B
      # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
      mutate(Spatial_A = case_when(Spatial_A == 1 ~ "C", Spatial_A == 2 ~ "B")) %>%
      mutate(Spatial_B = case_when(Spatial_B == 1 ~ "C", Spatial_B == 2 ~ "B")) %>% 
      # Create a last column "Group" that is a combination of Habitat and Spatial and "Type" for the intra comparisons (Same Habitat and Spatial)
      mutate("Group" = paste0(Spatial_A,Habitat_A,"_",Spatial_B,Habitat_B)) %>%
      # Transform what needs to be stransformed into factors
      mutate(across(c(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B),as.factor))
    
# Phylosor
PhyloSor_Obs <- PhyloFroc_Obs %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,"PhyloSor" =  betadiv.PhyloSor)

# Phylosor_Turn
PhyloSor_Turn_Obs <- PhyloFroc_Obs %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,"PhyloSor_Turn" =  betadiv.PhyloSor_turn)

# Phylosor_Ab
PhyloSor_Ab_Obs <-  as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo, incidence = F, method = "sorensen"))) %>%
    MetaData(Type = "Beta", Name = "PhyloSor_Ab", Rep = Rep, Diag = F)   

# UniFrac
UniFrac_Obs <- PhyloFroc_Obs %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,"UniFrac" =  betadiv.UniFrac)

# Unifrac_Turn
UniFrac_Turn_Obs <- PhyloFroc_Obs %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,"UniFrac_Turn" =  betadiv.UniFrac_turn)

# Unifrac_Ab
UniFrac_Ab_Obs <- as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo, incidence = F, method = "jaccard"))) %>%
  MetaData(Type = "Beta", Name = "UniFrac_Ab", Rep = Rep, Diag = F)

# Message
cat(col_yellow("DONE:"),"PhyloSor and UniFrac observed", "\n") ; toc() 
    
# --- NULL-MODEL ------------------------------------------------------------- #
tic()
    
# PhyloFroc_Null
PhyloFroc_Null <- 
  foreach(j = 1:NNM, .verbose = F) %do% {
      
    # Create a lite version of Sites_Data 
    Sites <- Sites_Data %>% 
      filter(rep == Rep) %>% # Select the replicate
      select(!c(X,Y)) # Remove the coordinates of the plot
      
  Metric <- as.data.frame(beta.pd.decompo(com = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], tree = List_Phylo_NULL[[j]]))  %>%
    tibble::rownames_to_column("Sample") %>% Index(.) %>%
    mutate(across(starts_with("Sample"), as.factor)) %>%
    # Get the MetaData
    left_join(y = Sites, join_by("Sample_A" == "sample")) %>% rename(Rep = rep, Habitat_A = hab, Spatial_A = SampleType) %>% # Get the MetaData for Sample_A
    left_join(y = Sites, join_by("Rep" == "rep","Sample_B" == "sample")) %>% rename(Habitat_B = hab, Spatial_B = SampleType) %>% # Get the MetaData for Sample_B
    # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2) 
    mutate(Spatial_A = case_when(Spatial_A == 1 ~ "C", Spatial_A == 2 ~ "B")) %>%
    mutate(Spatial_B = case_when(Spatial_B == 1 ~ "C", Spatial_B == 2 ~ "B")) %>% 
    # Create a last column "Group" that is a combination of Habitat and Spatial and "Type" for the intra comparisons (Same Habitat and Spatial)
    mutate("Group" = paste0(Spatial_A,Habitat_A,"_",Spatial_B,Habitat_B)) %>%
    # Transform what needs to be stransformed into factors
    mutate(across(c(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B),as.factor))
  
  } 

PhyloSor_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    metric <- PhyloFroc_Null[[j]] %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,betadiv.PhyloSor) %>%
      rename_with(.fn = ~ paste0("PhyloSor_NM",j), .cols = betadiv.PhyloSor)
  }

PhyloSor_Turn_Null <-   
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    metric <- PhyloFroc_Null[[j]] %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,betadiv.PhyloSor_turn) %>% # Discard the column SR
      rename_with(.fn = ~ paste0("PhyloSor_Turn_NM",j), .cols = betadiv.PhyloSor_turn)
  }

PhyloSor_Ab_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo_NULL[[j]], incidence = F, method = "sorensen"))) %>%
    MetaData(Type = "Beta", Name = paste0("PhyloSor_Ab_NM",j), Rep = Rep, Diag = F)
  }
        
UniFrac_Null <-  
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    metric <- PhyloFroc_Null[[j]] %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,betadiv.UniFrac)  %>% # Discard the column SR
      rename_with(.fn = ~ paste0("UniFrac_NM",j), .cols = betadiv.UniFrac)
  } 

UniFrac_Turn_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    metric <- PhyloFroc_Null[[j]] %>% select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,betadiv.UniFrac_turn)  %>% # Discard the column SR
      rename_with(.fn = ~ paste0("UniFrac_Turn_NM",j), .cols = betadiv.UniFrac_turn)
  }

UniFrac_Ab_Null <- 
   foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(phylosim(x = Sampled_Pre_Abs_Data_Inter[,-c(1:2)], phy = List_Phylo_NULL[[j]], incidence = F, method = "jaccard"))) %>%
    MetaData(Type = "Beta", Name = paste0("UniFrac_Ab_NM",j), Rep = Rep, Diag = F)
    
  }
        
# Message
cat(col_yellow("DONE:"),"Null_Model PhyloSor and UniFrac","\n")
toc()  

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Join the data
PhyloSor_Total <- full_join(PhyloSor_Obs,PhyloSor_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B)) 
  write.csv(PhyloSor_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_PhyloSor_Bryant_",Rep,".csv"))

PhyloSor_Turn_Total <- full_join(PhyloSor_Turn_Obs,PhyloSor_Turn_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B)) 
  write.csv(PhyloSor_Turn_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_PhyloSor_Turn_Leprieur_",Rep,".csv"))

PhyloSor_Ab_Total <- full_join(PhyloSor_Ab_Obs,PhyloSor_Ab_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
  write.csv(PhyloSor_Ab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_PhyloSor_Ab_Nipperess_",Rep,".csv"))

UniFrac_Total <- full_join(UniFrac_Obs,UniFrac_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
  write.csv(UniFrac_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_UniFrac_Lozupone_",Rep,".csv"))

UniFrac_Turn_Total <- full_join(UniFrac_Turn_Obs,UniFrac_Turn_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
  write.csv(UniFrac_Turn_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_UniFrac_Turn_Leprieur_",Rep,".csv"))

UniFrac_Ab_Total <- full_join(UniFrac_Ab_Obs,UniFrac_Ab_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
  write.csv(UniFrac_Ab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc/Scenario",Scenario,"_UniFrac_Ab_Leprieur_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: PhyloSor and UniFrac", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Phylosor and Unifrac

# --------------------------------------------------------------------------------------------------------------------------------------------- #        

# -- S_METRICS (Pavoine & Ricotta, 2014) ----

# Do we want to compute the PD (as well as avPD) ? 
if (S == FALSE) {
  cli_alert_warning("S_Metrics not computed")
} else {
  
  #  "S_Jaccard"      = Generalization of Jaccard measure of species turnover
  #  "S_Ochial"       = Generalization of Ochial measure of species turnover
  #  "S_Sorensen"     = Generalization of Sorensen measure of species turnover
  #  "S_Sokal-Sneath" = Generalization of Sokal & Sneath measure of species turnover
  #  "S_Beta"         = Generalization of Rao's Dab measure of species turnover
  
cat(rule(left = "S_Metrics ", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()     

# Compute the dissimilarity values. 

# We need, as the function require that the dissimilarity between species as rows and columns in the same order as in comm,
# Therefore, we will need to reorder all the columns and rows. 
sp_order <- colnames(Sampled_Pre_Abs_Data_Inter) 
sp_order <- sp_order [! sp_order  %in% c("rep","sample")]
  
# Sokal-sneath
SS <- dsimTree(phyl = List_Phylo, rootedge = NULL, method = 1, type = "similarity")
SS <- SS[sp_order,sp_order]
# Jaccard
J <- dsimTree(phyl = List_Phylo, rootedge = NULL, method = 2, type = "similarity")
J <- J[sp_order,sp_order]
# Sorensen
S <- dsimTree(phyl = List_Phylo, rootedge = NULL, method = 3, type = "similarity")
S <- S[sp_order,sp_order]
# Ochiai  
O <- dsimTree(phyl = List_Phylo, rootedge = NULL, method = 4, type = "similarity")
O <- O[sp_order,sp_order]
# S_Beta
SB <- dsimTree(phyl = List_Phylo, rootedge = NULL, method = 5, type = "similarity")
SB <- SB[sp_order,sp_order]
  
# Create a list 
Dissim_Obs<- list(SS,J,S,O,SB) %>% setNames(.,nm = c("SS","J","S","O","SB"))


# S_Metrics: Observed values
S_SokalSneath_Obs <- as.data.frame(as.matrix(dsimcom(Sampled_Abund_Data_Inter[,-c(1:2)], Dissim_Obs$SS, method = 1, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = "S_SokalSneath", Rep = Rep, Diag = F)   

S_Jaccard_Obs <- as.data.frame(as.matrix(dsimcom(Sampled_Abund_Data_Inter[,-c(1:2)], Dissim_Obs$J, method = 2, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = "S_Jaccard", Rep = Rep, Diag = F)

S_Sorensen_Obs <- as.data.frame(as.matrix(dsimcom(Sampled_Abund_Data_Inter[,-c(1:2)], Dissim_Obs$S, method = 3, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = "S_Sorensen", Rep = Rep, Diag = F)

S_Ochiai_Obs <-  as.data.frame(as.matrix(dsimcom(Sampled_Abund_Data_Inter[,-c(1:2)], Dissim_Obs$O, method = 4, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = "S_Ochiai", Rep = Rep, Diag = F)

S_Beta_Obs <- as.data.frame(as.matrix(dsimcom(Sampled_Abund_Data_Inter[,-c(1:2)], Dissim_Obs$SB, method = 5, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = "S_Beta", Rep = Rep, Diag = F)

# Message
cat(col_yellow("DONE:"),"S_Metrics observed", "\n") ; toc() 

# --- NULL-MODEL ------------------------------------------------------------- #
tic()

# Compute the dissimilarity values. 
Dissim_Null <- 
  foreach(j = 1:NNM, .verbose = F) %do% { 

    # We need, as the function require that the dissim species as rows and columns in the same order as in comm,
    # Therefore, we will need to reorder all the columns and rows. 
    sp_order <- colnames(Sampled_Pre_Abs_Data_Inter)
    sp_order <- sp_order [! sp_order  %in% c("rep","sample")]
      
    # Sokal-sneath
    SS <- dsimTree(phyl = List_Phylo_NULL[[j]], rootedge = NULL, method = 1, type = "similarity")
    SS <- SS[sp_order,sp_order]
    # Jaccard
    J <- dsimTree(phyl = List_Phylo_NULL[[j]], rootedge = NULL, method = 2, type = "similarity")
    J <- J[sp_order,sp_order]
    # Sorensen
    S <- dsimTree(phyl = List_Phylo_NULL[[j]], rootedge = NULL, method = 3, type = "similarity")
    S <- S[sp_order,sp_order]
    # Ochiai
    O <- dsimTree(phyl = List_Phylo_NULL[[j]], rootedge = NULL, method = 4, type = "similarity")
    O <- O[sp_order,sp_order]
    # S_Beta
    SB <- dsimTree(phyl = List_Phylo_NULL[[j]], rootedge = NULL, method = 5, type = "similarity")
    SB <- SB[sp_order,sp_order]
    # Create a list 
    Dissim_List <- list(SS,J,S,O,SB) %>% setNames(.,nm = c("SS","J","S","O","SB"))
} 
  
# S_Metrics: Observed values
S_SokalSneath_Null <-
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(dsimcom(comm = Sampled_Abund_Data_Inter[,-c(1:2)], Sigma = Dissim_Null[[j]]$SS, method = 1, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = paste0("S_SokalSneath_NM",j), Rep = Rep, Diag = F)
  }

S_Jaccard_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(dsimcom(comm = Sampled_Abund_Data_Inter[,-c(1:2)], Sigma = Dissim_Null[[j]]$J, method = 2, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = paste0("S_Jaccard_NM",j), Rep = Rep, Diag = F)
  } 
    
S_Sorensen_Null <-  
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(dsimcom(comm = Sampled_Abund_Data_Inter[,-c(1:2)], Sigma = Dissim_Null[[j]]$S, method = 3, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = paste0("S_Sorensen_NM",j), Rep = Rep, Diag = F)
  }

S_Ochiai_Null <-
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(dsimcom(comm = Sampled_Abund_Data_Inter[,-c(1:2)], Sigma = Dissim_Null[[j]]$O, method = 4, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = paste0("S_Ochiai_NM",j), Rep = Rep, Diag = F)
  }

S_Beta_Null <-
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.data.frame(as.matrix(dsimcom(comm = Sampled_Abund_Data_Inter[,-c(1:2)], Sigma = Dissim_Null[[j]]$SB, method = 5, option = "absolute", type = "similarity"))) %>%
      MetaData(Type = "Beta", Name = paste0("S_Beta_NM",j), Rep = Rep, Diag = F)
  }

# Message
cat(col_yellow("DONE:"),"Null_Model S_Metrics","\n")
toc()  

# --- COMBINE AND SAVE ------------------------------------------------------- #

S_SokalSneath_Total <-  full_join(S_SokalSneath_Obs,S_SokalSneath_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(S_SokalSneath_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S/Scenario",Scenario,"_S_SokalSneath_Pavoine_",paste(Parameters$Rep,collapse = "-"),".csv"))

S_Jaccard_Total <-  full_join(S_Jaccard_Obs,S_Jaccard_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(S_Jaccard_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S/Scenario",Scenario,"_S_Jaccard_Pavoine_",paste(Parameters$Rep,collapse = "-"),".csv"))

S_Sorensen_Total <-  full_join(S_Sorensen_Obs,S_Sorensen_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(S_Sorensen_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S/Scenario",Scenario,"_S_Sorensen_Pavoine_",paste(Parameters$Rep,collapse = "-"),".csv"))

S_Ochiai_Total <-  full_join(S_Ochiai_Obs,S_Ochiai_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(S_Ochiai_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S/Scenario",Scenario,"_S_Ochiai_Pavoine_",paste(Parameters$Rep,collapse = "-"),".csv"))

S_Beta_Total <-  full_join(S_Beta_Obs,S_Beta_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(S_Beta_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S/Scenario",Scenario,"_S_Beta_Pavoine_",paste(Parameters$Rep,collapse = "-"),".csv"))

# --- #
cat(rule(left = "DONE: S_Metrics", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of S_Metrics

# --------------------------------------------------------------------------------------------------------------------------------------------- #  
    
# -- Pst, Bst, PIst (Hardy & Senterre, 2007) ----
  
# Do we want to compute the PIst ? 
if (PIst == FALSE) {
  cli_alert_warning("PIst not computed")
} else {
  
  #  "Pst"  = MPDab BASED BETA DIVERSITY       = Additive decomposition of MPDab
  #  "Bst"  = InterMPDab BASED BETA DIVERSITY  = Additive decomposition of InterMPDab
  #  "PIst" = MPD BASED PROPORTIONAL DIVERSITY = Additive decomposition of MPD / Proportion of overall phylogenetic distinctness expressed among plots
  
cat(rule(left = "PST, PIst and BST", line_col = "yellow", line = "-", col = "br_yellow"))
  
  # Load the needed function
  source(paste0(getwd(),"/Utilities/Diversity_Metrics/SpacodirR.R"))  # Function to compute the F quadratic entropy
  
# --- OBSERVED --------------------------------------------------------------- #
tic()    

Hardy_Metrics <- spacodi.calc(
      sp.plot = t(Sampled_Abund_Data_Inter[,-c(1:2)]),   # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
      phy = List_Phylo,            # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
      all.together = TRUE,             # whether to treat all traits together or separately
      prune = TRUE,
      pairwise = TRUE)

Pst_Obs <- Hardy_Metrics$pairwise.Pst %>%
    MetaData(Type = "Beta", Name = "Pst", Rep = Rep, Diag = F)

Bst_Obs <- Hardy_Metrics$pairwise.Bst %>%
    MetaData(Type = "Beta", Name = "Bst", Rep = Rep, Diag = F)
    
PIst_Obs <- Hardy_Metrics$pairwise.PIst %>%
    MetaData(Type = "Beta", Name = "PIst", Rep = Rep, Diag = F)
      
# Message
cat(col_yellow("DONE:"),"Pst, Bst and PIst observed", "\n") ; toc() 

# --- NULL-MODEL ------------------------------------------------------------- #
tic()
    
Hardy_Metrics <-
  foreach(j = 1:NNM, .verbose = F) %do% { 
      
  # Compute all the metrics
    Metric <- spacodi.calc(
      sp.plot = t(Sampled_Abund_Data_Inter[,-c(1:2)]),   # sp.plot =  a community dataset in spacodiR format (see as.spacodi) i.e species in rows and plots in columns
      phy = List_Phylo_NULL[[j]],            # phy a phylogenetic tree of class phylo or evolutionary distance matrix between species (see cophenetic.phylo)                # sp.traits a species-by-trait(s) dataframe or a species traits distance matrix (see dist)
      all.together = TRUE,             # whether to treat all traits together or separately
      prune = TRUE,
      pairwise = TRUE)
    
}
    
Pst_Null <-   
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- Hardy_Metrics[[j]]$pairwise.Pst %>%
      MetaData(Type = "Beta", Name = paste0("Pst_NM",j), Rep = Rep, Diag = F)
  }

Bst_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- Hardy_Metrics[[j]]$pairwise.Bst  %>%
      MetaData(Type = "Beta", Name = paste0("Bst_NM",j), Rep = Rep, Diag = F)
  } 
      
PIst_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- Hardy_Metrics[[j]]$pairwise.PIst %>%
      MetaData(Type = "Beta", Name = paste0("PIst_NM",j), Rep = Rep, Diag = F)
  }
# Message
cat(col_yellow("DONE:"),"Null_Model Pst, Bst and PIst","\n")
toc()  

# --- COMBINE AND SAVE ------------------------------------------------------- #

Pst_Total <-  full_join(Pst_Obs,Pst_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(Pst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST/Scenario",Scenario,"_Pst_Hardy_",Rep,".csv"))

Bst_Total <-  full_join(Bst_Obs,Bst_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B)) 
write.csv(Bst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST/Scenario",Scenario,"_Bst_Hardy_",Rep,".csv"))

PIst_Total <-  full_join(PIst_Obs,PIst_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B)) 
write.csv(PIst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST/Scenario",Scenario,"_PIst_Hardy_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: Pst, Bst and PIst", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Pst, Bst and PIst

# --------------------------------------------------------------------------------------------------------------------------------------------- #  

# -- PCD (Helmus, 2007) ----

# Do we want to compute the PCD ? 
if (PCD == FALSE) {
  cli_alert_warning("PCD not computed")
} else {
  
  # "PCD" = PHYLOGENETIC COMMUNITY DISSIMILARITY = The pairwise differences between communities derived by asking how much of the variance among species in the values of a hypothetical non selected trait in one community can be predicted by the known trait values of species in another community.
  
  cat(rule(left = "PCD", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic()    
  
PCD_Obs <- as.matrix(pcd(comm = Sampled_Pre_Abs_Data_Inter[,-c(1:2)] ,tree = List_Phylo)$PCD) %>%
    MetaData(Type = "Beta", Name = "PCD", Rep = Rep, Diag = F)

# Message
cat(col_yellow("DONE:"),"PCD observed", "\n") ; toc() 

# --- NULL-MODEL ------------------------------------------------------------- #
tic()

PCD_Null <- 
  foreach(j = 1:NNM,.combine = merge, .verbose = F) %do% { 
    Metric <- as.matrix(pcd(comm = Sampled_Pre_Abs_Data_Inter[,-c(1:2)] ,tree = List_Phylo_NULL[[j]])$PCD) %>%
      MetaData(Type = "Beta", Name = paste0("PCD_NM",j), Rep = Rep, Diag = F)
  }
  
# Message
cat(col_yellow("DONE:"),"Null_Model PCD","\n")
toc()  

# --- COMBINE AND SAVE ------------------------------------------------------- #

PCD_Total <-  full_join(PCD_Obs,PCD_Null, by = join_by(Rep, Sample, Group, Sample_A, Habitat_A, Spatial_A, Sample_B, Habitat_B, Spatial_B))
write.csv(PCD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PCD/Scenario",Scenario,"_PCD_Helmus_",Rep,".csv"))

# --- #
cat(rule(left = "DONE: PCD", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of PCD

# --------------------------------------------------------------------------------------------------------------------------------------------- #  

# -- Simpson Diversity ----

#  "SIMPSON INDEX" = 

# Do we want to compute the SIMPSON INDEX? 
if (Simpson == FALSE) {
  cli_alert_warning("SIMPSON not computed")
} else {
  
  cat(rule(left = "SIMPSON", line_col = "yellow", line = "-", col = "br_yellow"))
  
# --- OBSERVED --------------------------------------------------------------- #
tic() 

  # Create a lite version of Sites_Data 
  Sites <- Sites_Data %>% 
    filter(rep == Rep) %>% # Select the replicate
    select(!c(X,Y)) # Remove the coordinates of the plot  
  
  Simpson_Obs <- as.data.frame(vegan::diversity(Sampled_Abund_Data_Inter[,-(1:2)],index = "simpson")) %>%
    setNames("Simpson") %>%
    # Add the sample number
    tibble::rownames_to_column("Sample") %>%
    # Get the MetaData
    left_join(y = Sites, join_by("Sample" == "sample")) %>% dplyr::rename(Rep = rep, Habitat = hab, Spatial = SampleType) %>% # Get the MetaData for the Sample
    # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
    mutate(Spatial = case_when(Spatial == 1 ~ "C", Spatial == 2 ~ "B")) %>%
    # Create a last column "Group" that is a combination of Habitat and Spatial and "Type" for the intra comparisons (Same Habitat and Spatial)
    mutate("Group" = paste0(Spatial,Habitat)) %>%
    # Transform what needs to be stransformed into factors
    mutate(across(c(Rep,Sample,Group,Habitat,Spatial),as.factor)) %>%
    # Reorder the columns in a more friendly way 
    select(Rep,Sample,Group,Habitat,Spatial,all_of("Simpson"))

# --- COMBINE AND SAVE ------------------------------------------------------- #
  
# Save the dataframe
write.csv(Simpson_Obs,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/Simpson/Scenario",Scenario,"_Simpson.csv"))
  
# --- #
  
toc() ; cat(rule(left = "DONE: SIMPSON INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     
  
} # End of Simpson

# END OF THE SCRIPT
cat(rule(left = "SCRIPT COMET_METRICS.R ENDING", line_col = "red", line = "-", col = "br_red")) 

