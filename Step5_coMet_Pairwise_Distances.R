#!/usr/bin/env Rscript

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");require(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_PAIRWISE_DISTANCES.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

#----------------------------------------------------------------------------------#
##### coMet_Pairwise_Distances: Computation of phylogenetic pairwise distances #####
#----------------------------------------------------------------------------------#

# This script computes pairwise phylogenetic patristic distances from phylogenetic trees. It takes as input a simple tree or multiple trees (all treated as multiphylo objects as a multiphylo object to have a generalized script). 
# It outputs two matrices containing or not the intra specific distance (equal to 0) but important for the metrics using the mean phylogenetic distance. 

# Null-PW Distances are computed inside the coMet_Metric.R script by using the Obs_PW_distance and shuffling the tips_label from the list_Phylo_NULL trees. This avoid the save of a huge amount of memory. 

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
Parameters$Config <- "coMet_ConfigFiles/Equalizing100_Sig5.R"

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
  rep <- stringi::stri_extract_last_regex(names(Patristic_PW_M1)[[j]], '\\d+')
  
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


# THIS WAS THE OLD WAY, BEFORE THE TIP-SHUFFLING IN THE OBS-PW MATRIX. 

# ---------------------------- #
##### ---- NULL-MODEL ---- #####
# ---------------------------- # 

#   # Keep in memory the time already used since the opening of the script to later only have the time of the computation. 
#   ptm <- proc.time() 
# 
#   # Display a message
#   cat(rule(left = "Computing the NULL_MODEL pairwise distances", line_col = "green", line = "~", col = "br_green"))  
# 
#   # Create a loop to load all multiphylo trees for each simulations
#   PW_NM <-  foreach(i = 1:Nrep) %do% {
#     
#   # Load the tree(s)
#   Phylo_tree_NM <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"Step3_MultiPhylo_NullModel",i,".tree"), keep.multi = TRUE) 
#   NbTree <- length(Phylo_tree_NM)
# 
#   # Print 
#   cli_alert_info(paste0("Scenario" ,Scenario," Multiphylo file ",i,"/",Nrep," Rep correctly loaded (", NbTree," Tree(s))."))
#   
#   # Create a loop to loop on each trees (to take into account multiphylo objects) and compute the pairwise dist
#   # Total Pairwise Distance (Matrix Object)
#   Patristic_PW_M1 <- foreach(j = 1:NbTree) %dopar% {
#       cli::cli_alert_info(paste0("Pairwise distances of file Scenario",Scenario,"Step3_MultiPhylo_NullModel",i,".tree : ",j,"/",NbTree," tree."))
#       # Total Pairwise Distance (Distance Object)
#       Patristic_PW <- distTips(Phylo_tree_NM[[j]], tips = Phylo_tree_NM[[j]]$tip.label, method = "patristic")
#       # Total Pairwise Distance (Matrix Object)
#       Patristic_PW <- as.matrix(Patristic_PW)
#     } %>% setNames(nm = paste0("Tree_NM",1:length(.))) # Set the names of the tree_NM
#   
#   # Display a message
#   cat(rule(left = paste0("DONE: Pairwise distances of replicate ",i,"/",Nrep), line_col = "green", line = "~", col = "br_green"))  
# 
#   # Lower matrix with Diagonal - With Intra-Specific Distance
#   Patristic_PW_M2 <- foreach(j = 1:length(Patristic_PW_M1)) %dopar% {
#       # Create a copy of Patristic_PW_M1
#       PW <- Patristic_PW_M1[[j]]
#       # Make the matrix triangular
#       PW[upper.tri(PW,diag = F)] <- NA  
#       # Adding column identifying a tree ID, especially usefull for null models
#       PW <- as.data.frame(PW) %>%
#         `row.names<-`(colnames(PW)) %>% # Change the rownames with the colnames
#         mutate(tree_NM = j,.before = 1) %>% # Add the Tree_NM column and number
#         mutate(Sp = rownames(.),.before = 1) # Add a column "sp" to easily identify the species
#       
#     } %>% setNames(nm = paste0("Tree_NM",1:length(.))) # Set the names of the tree_NM
# 
#   # Lower matrix without Diagonal - Without Intra-Specific Distance
#   Patristic_PW_M3 <- foreach(j = 1:length(Patristic_PW_M1)) %dopar% {
#       # Create a copy of Patristic_PW_M1
#       PW <- Patristic_PW_M1[[j]]
#       # Make the matrix triangular
#       PW[upper.tri(PW,diag = T)] <- NA  
#       # Adding column identifying a tree ID, especially usefull for null models
#       PW <- as.data.frame(PW) %>%
#         `row.names<-`(colnames(PW)) %>% # Change the rownames with the colnames
#         mutate(tree_NM = j,.before = 1) %>% # Add the Tree_NM column and number
#         mutate(Sp = rownames(.),.before = 1) # Add a column "sp" to easily identify the species
#       
#     } %>% setNames(nm = paste0("Tree_NM",1:length(.))) # Set the names of the tree_NM
# 
#   # Create a list to return
#   PW <- list("Total_PW" = Patristic_PW_M1, "Intra_PW" = Patristic_PW_M2, "Inter_PW" = Patristic_PW_M3)
#   # Return it
#   return(PW)
#   
# } %>% setNames(nm = paste0("rep",1:length(.))) # Set the names of the replicates
# 
#   # Create global matrix with all replicates to save and work with if they were multiple trees.
#     # Split by replicates and then bind by NM_trees
#   Nihil <- foreach(i = 1:length(PW_NM)) %dopar%
#     if (length(PW_NM[[i]]$Total_PW) > 1){
#       # Get the intra value
#       Total_PW_M2  <- bind_rows(PW_NM[[i]]$Intra_PW) %>% 
#         write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"Step4_Intra_PW_Distances_Rep",i,".csv"))
#       # Total_PW_M3  <- bind_rows(Patristic_PW_M3) %>% 
#         # write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"Step4_Inter_PW_Distances.csv"))
#     } else {
#       # Get the intra value
#       write.csv(Patristic_PW_M2,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"Step4_Intra_PW_Distances.csv"))
#       # write.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"Step4_Inter_PW_Distances.csv"))
#   } 