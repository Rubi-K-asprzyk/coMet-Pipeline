#! /bin/bash

# ------------------------------------------#
# ----- WELCOME TO THE COMET PIPELINE ----- #
# ------------------------------------------#

# This bash script is designed to run the entire coMet_pipeline from one configuration file and the scenario name associated.
# By default, this configuration file is called Foo.R and the Scenario name is "Foo"

# This script will automatically create the directories needed for the complete pipeline to be run. 
# Be carefull that these directories will be created in the folder you are currently in. 

# It will after launch successively all R scripts of the pipeline.
# Each script will need the result from the previous steps. But, they can also be launched independantly from the console. 

# /!\ CAUTION /!\ : Step 1 and 6 are designed to be launched on an array on a cluster with slurm architecture to allow multiple replicates to be launched simultaneously. If you want to launch them on a local machine, you have to launch each replicate one after the other with the command line that is written. 

# -------------------------------------------------- #

# Save the name of the configuration file used and the option for the files to download.
Config=$1
Scenario=$2

# -----------------------------------#
# ----- SUB-DIRECTORY CREATION ----- #
# -----------------------------------#

# Create a dir for all outputs.
mkdir -p coMet_ComSim_Outputs

# Create a dir for the wanted scenario in coMet_ComSim_Outputs
mkdir -p coMet_ComSim_Outputs/$Scenario

# Create the sub-directories in coMet_ComSim_Outputs 
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics
mkdir -p coMet_ComSim_Outputs/$Scenario/Sampled_Communities
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities

# Create the sub-directories in coMet_ComSim_Outputs/Whole_Communities
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/DispRec
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/FitnessEsp
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/PhyloTree
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/summary_stats
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/TabFull
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/TraitsEsp
mkdir -p coMet_ComSim_Outputs/$Scenario/Whole_Communities/TreeDistances

# Create the sub-directories in coMet_ComSim_Outputs 
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PD
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/F
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PSV
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PDab
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/MPD
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/MNTD
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/ED
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PDb
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PhyloFroc
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/S
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PST
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/PCD
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/Simpson
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/InvSimpson
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/BC
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/SR
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/Sorensen
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/SorensenNest
mkdir -p coMet_ComSim_Outputs/$Scenario/Metrics/Temporary/SorensenTurn

# ------------------------------------ #
# ----- COMET-PIPELINE EXECUTION ----- #
# ------------------------------------ #

##### Step 1: coMet_ComSim #####

# Creation of the communities following the parameters and assembly rules wanted.
# This script is designed to be launched on a HPC to complete all replicates on an array. 
# However, if needed,  it is possible to launch each replicate "X" using this command line: Rscript ./Step1_coMet_ComSim.R Foo.R --Rep X
# Rscript ./Step1_coMet_ComSim.R $Config

##### Step 2: coMet_ComSim_Merge #####

# Combination of all the replicates from the first step into one file per "Scenario" for further processing. 
# Rscript ./Step2_coMet_ComSim_Merge.R $Config

##### Step 3: coMet_ComSim_Check #####

# Display useful informations from the outputs of the coMet_ComSim.R script
# Rscript ./Step3_coMet_ComSim_Check.R $Config

##### Step 4: coMet_Tree_NullModel #####

# Display useful informations from the outputs of the coMet_ComSim.R script
# Rscript ./Step4_coMet_Tree_NullModel.R $Config

##### Step 5: coMet_Pairwise_Distances #####

# Display useful informations from the outputs of the coMet_ComSim.R script
# Rscript ./Step5_coMet_Pairwise_Distances.R $Config

##### Step 6: coMet_Metrics #####

# This script is designed to be launched on a HPC to complete all replicates on an array. 
# However, if needed,  it is possible to launch each replicate "X" using this command line: Rscript ./Step6_coMet_Metrics.R Foo.R --Rep X
# Rscript ./Step6_coMet_Metrics.R $Config

##### Step 7: coMet_Metrics_Merge #####

# Merge all the replicates for each Metrics
# Rscript ./Step7_coMet_Metrics_Merge.R $Config

# - Remove temporary files - #


##### Step 8: coMet_Metrics_Check #####

# This script is used in the coMet pipeline to draw the plots of the results for each metric and their significance tests for one scenario. 
# It contains density plots of the values of the plots as well as plots of the metric significance.  
# Rscript ./Step8_coMet_Metric_Check.R $Config

# --------------------------------------------- #

# /!\ Both Step 9 and 10 requires outputs from multiple scenarios to compare the results /!\ 
# /!\ Therefore, they have to be launched after the computation of multiple scenarios. /!\

##### Step 9: coMet_ComSim_Comparison #####

# This script is used in the coMet pipeline to compare and draw plots of the simulated communities.  
# Rscript ./Step9_coMet_ComSim_Comparison.R Foo.R Foo2.R ...

##### Step 10: coMet_Metric_Comparison #####

# This script is used in the coMet pipeline to compare and draw plots of the metrics results.  
# Rscript ./Step10_coMet_Metric_Comparison.R Foo.R Foo2.R ...
