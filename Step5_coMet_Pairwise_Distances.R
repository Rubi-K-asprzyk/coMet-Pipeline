#!/usr/bin/env Rscript

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");require(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_PAIRWISE_DISTANCES.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

#----------------------------------------------------------------------------------#
##### coMet_Pairwise_Distances: Computation of phylogenetic pairwise distances #####
#----------------------------------------------------------------------------------#

# This script computes pairwise phylogenetic patristic distances from phylogenetic trees. It takes as input a simple tree or multiple trees (all treated as multiphylo objects as a multiphylo object to have a generalized script). 
# It outputs two matrices containing or not the intra specific distance (equal to 0) but important for the metrics using the mean phylogenetic distance. 

# Null-PW Distances are computed inside the coMet_Metric.R script by using the Obs_PW_distance and shuffling the tips_label from the list_Phylo_NULL trees. This allow to save of a huge amount of memory. 

# ------------------------------ #
##### STEP 0: INITIALISATION #####
# ------------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))

    # ----- #

  #### . Packages . ####

# Print a message
cli_alert_info("Loading the needed packages ... ")

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})
p_load("foreach",
       "phytools",
       "doParallel",
       "adephylo",
       "dplyr",
       "argparser")

    # ----- #

  #### . Local Mode . ####

Parameters <- list()
Parameters$Config <- "Foo.R"

    # ----- #

  #### . Argument Parser . ####

cli_alert_info("Reading the coMet_ConfigFile.R ...")

# Create the parser
arg_parser <- arg_parser("Check the communities computed by coMet_ComSim and create an output pdf.", name = "coMet_ComSim_Check.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = 1, help = "Configuration file for coMet_ComSim.R")

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Print a message to inform that the right configuration file was correctly loaded. 
cli::cli_alert_info(paste0("Working on configuration file: ",Parameters$Config))

# Load them into the environnement
source(Parameters$Config)

# Print a message
cli_alert_success("coMet_ConfigFile.R correctly read !")
  
# ------------------------------------------------------ #
##### STEP 1: COMPUTATION OF THE PAIRWISE DISTANCES  #####
# ------------------------------------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

# Display a message
cat(rule(left = "COMPUTATION OF THE PAIRWISE DISTANCES", line_col = "green", line = "~", col = "br_green"))  

# Keep in memory the time already used since the opening of the script to later only have the time of the computation. 
ptm <- proc.time() 

    # ----- #

  #### . Loading files . ####

# Load the tree
Phylo_tree <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_PhyloTrees.tree"), keep.multi = TRUE) 
NbTree <- length(Phylo_tree)
  
# Load one of the MultiPhylo file for one of the replicates.
# Phylo_Null <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step3_MultiPhylo_NullModel1.tree"), keep.multi = TRUE) 


# Print 
cli_alert_info(paste0("Scenario" ,Scenario,": Observed Multiphylo file correctly loaded (", length(Phylo_tree)," Replicate(s))."))

    # ----- #

##### ------------ 2.A: TOTAL PAIRWISE DISTANCE -------------------- #####
  
# Create a loop to loop on each trees (to take into account multiphylo objects) and compute the pairwise dist

# Total Pairwise Distance (Matrix Object)
Patristic_PW_M1 <- foreach(j = 1:NbTree) %dopar% {
  
  cli::cli_alert_info(paste0(Scenario,": Observed tree (",j,"/",NbTree," tree.)"))
  # Total Pairwise Distance (Distance Object)
  Patristic_PW <- distTips(Phylo_tree[[j]], tips = Phylo_tree[[j]]$tip.label, method = "patristic")
  # Total Pairwise Distance (Matrix Object)
  Patristic_PW <- as.matrix(Patristic_PW)

} %>% 
  setNames(nm = names(Phylo_tree)) # Set the names of the tree_NM

# Display a message
cat(rule(left = paste0("DONE: Pairwise distances of Observed tree."), line_col = "green", line = "~", col = "br_green"))  

    # ----- #
  
##### ------------ 2.B: INTRA PAIRWISE DISTANCE -------------------- #####

# Lower matrix with Diagonal - With Intra-Specific Distance
Patristic_PW_M2 <- foreach(j = 1:NbTree) %dopar% {
  
  # Extract the last number of the string file that is the replicate number to correctly identify them
  # rep <- stringi::stri_extract_last_regex(names(Patristic_PW_M1)[[j]], '\\d+')
  rep <- j
  
  # Create a copy of Patristic_PW_M1
  PW <- Patristic_PW_M1[[j]]
  # Make the matrix triangular
  PW[upper.tri(PW,diag = F)] <- NA  
  # Adding column identifying a tree ID, especially usefull for null models
  PW <- as.data.frame(PW) %>%
    `row.names<-`(colnames(PW)) %>% # Change the rownames with the colnames
    mutate(tree_NM = rep,.before = 1) %>% # Add the Tree_NM column and number
    mutate(Sp = rownames(.),.before = 1) # Add a column "sp" to easily identify the species
  
} %>% 
  setNames(nm = names(Phylo_tree)) # Set the names of the tree_NM

# Display a message
cat(rule(left = paste0("DONE: Intra pairwise distances of Observed tree."), line_col = "green", line = "~", col = "br_green"))

    # ----- #
  
##### ------------ 2.C: INTER PAIRWISE DISTANCE -------------------- #####

# Lower matrix without Diagonal - Without Intra-Specific Distance
Patristic_PW_M3 <- foreach(j = 1:NbTree) %dopar% {
  # Extract the last number of the string file that is the replicate number to correctly identify them
  rep <- stringi::stri_extract_last_regex(names(Patristic_PW_M1)[[j]], '\\d+')
  rep <- j 
  # Create a copy of Patristic_PW_M1
  PW <- Patristic_PW_M1[[j]]
  # Make the matrix triangular
  PW[upper.tri(PW,diag = T)] <- NA  
  # Adding column identifying a tree ID, especially usefull for null models
  PW <- as.data.frame(PW) %>%
    `row.names<-`(colnames(PW)) %>% # Change the rownames with the colnames
    mutate(tree_NM = rep,.before = 1) %>% # Add the Tree_NM column and number
    mutate(Sp = rownames(.),.before = 1) # Add a column "sp" to easily identify the species
  
} %>% 
  setNames(nm = names(Phylo_tree)) # Set the names of the tree_NM
  
# ----- #

# ---------------------- #
##### STEP 2: SAVING #####
# ---------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Create a list to return
PW_Obs <- list("Total_PW" = Patristic_PW_M1, "Intra_PW" = Patristic_PW_M2, "Inter_PW" = Patristic_PW_M3)

# Save the data
for (i in 1:length(PW_Obs$Intra_PW)) {
  
  # Identify the replicate
  R <- stringi::stri_extract_last_regex(names(PW_Obs$Intra_PW[i]) , '\\d+')
  R <- i 

  write.csv(as.data.frame(PW_Obs$Intra_PW[i], 
                          col.names = colnames(PW_Obs$Intra_PW[i])),
                          paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/TreeDistances/Scenario",Scenario,"_Step4_Intra_PW_Distances_Observed",R,".csv")) 
  
}

# Show the elapsed time of the simulation
print(proc.time() - ptm)
  
# Print a message 
cli_alert_success("Pairwise distances correctly saved !")   

# END OF THE SCRIPT
cat(rule(left = "SCRIPT COMET_PAIRWISE_DISTANCES.R ENDING", line_col = "red", line = "-", col = "br_red")) 

