#!/usr/bin/env Rscript

##### ______________________________________________________________________________________________________________________________________ #####

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_METRICS_MERGE.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

# ---------------------------------------------------------------- #
# coMet_Metrics_Merge : Merge the pĥylogenetic diversity metrics.  #
# ---------------------------------------------------------------- #

# Merging the results of each metric for each scenario and computing the significativity of each metrics VS the Null Model. 

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
       "argparser",
       "ggplotify",
       "tidyr")

cli_alert_success("Packages correctly loaded !")

# ----- #

  #### . Local Mode . ####

Parameters <- list()
Parameters$Config <- "Foo.R"
registerDoParallel(10)

# ----- #

  #### . Argument Parser . ####

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
      filter(rep == Rep) %>% # Select the replicate
      select(!c(X,Y)) # Remove the coordinates of the plot  
      
    Metric <- Metric %>%
      # Add the sample number
      tibble::rownames_to_column("Sample") %>%
      # Get the MetaData
      left_join(y = Sites, join_by("Sample" == "sample")) %>% rename(Rep = rep, Habitat = hab, Spatial = SampleType) %>% # Get the MetaData for the Sample
      # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
      mutate(Spatial = case_when(Spatial == 1 ~ "C", Spatial == 2 ~ "B")) %>%
      mutate("Group" = paste0(Spatial,Habitat)) %>%
      # Transform what needs to be transformed into factors
      mutate(across(c(Rep,Sample,Group,Habitat,Spatial),as.factor)) %>% 
      # Reorder the columns in a more friendly way
      select(Rep,Sample,Group,Habitat,Spatial,all_of(Name))
  
    } # End of Type : Alpha 
    
  if (Type == "Beta"){
      
    # -- Transformation / MetaData -- # 
      
    # Create a lite version of Sites_Data 
    Sites <- Sites_Data %>% 
      filter(rep == Rep) %>% # Select the replicate
      select(!c(X,Y)) # Remove the coordinates of the plot
      
    Metric <- PW_to_Vector(Metric, Colname = Name, Type = "Intra", Diag = Diag) %>% # Transform the pairwise matrix into a vector
        # Get the MetaData
        left_join(y = Sites, join_by("Sample_A" == "sample")) %>% rename(Rep = rep, Habitat_A = hab, Spatial_A = SampleType) %>% # Get the MetaData for Sample_A
        left_join(y = Sites, join_by("Rep" == "rep","Sample_B" == "sample")) %>% rename(Habitat_B = hab, Spatial_B = SampleType) %>% # Get the MetaData for Sample_B
        # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
        mutate(Spatial_A = case_when(Spatial_A == 1 ~ "C", Spatial_A == 2 ~ "B")) %>%
        mutate(Spatial_B = case_when(Spatial_B == 1 ~ "C", Spatial_B == 2 ~ "B")) %>% 
        # Create a last column "Group" that is a combination of Habitat and Spatial and "Type" for the intra comparisons (Same Habitat and Spatial)
        mutate("Group" = paste0(Spatial_A,Habitat_A,"_",Spatial_B,Habitat_B)) %>%
        # Transform what needs to be stransformed into factors
        mutate(across(c(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B),as.factor)) %>% 
        # Reorder the columns in a more friendly way
        select(Rep,Sample,Group,Sample_A,Habitat_A,Spatial_A,Sample_B,Habitat_B,Spatial_B,all_of(Name))
      
    } # End of Type : Alpha 
 
  # Return the Results   
  return(Metric)
} 
  
  cli_li("DONE: MetaData ")

# -- Create of a function to compute the mean and the standard-deviation of each metric grouped by categories of SampleType-Habitat - #

  # Data is the Metric to Parse
  # Name is the name of the Metric

Stat_signif_mean <- function(Data,Name){
  
  # Compute the mean of the metric for each categories
  Mean_Metric <- Data %>%
    group_by(Rep,Group) %>% # Group the data by rep, Habitat and Sampletype
    dplyr::summarise(across(where(is.numeric),function(x) Mean = mean(x, na.rm = TRUE)), .groups = "drop")

  # Save the data
  write.csv(Mean_Metric,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_",Name,"_Mean.csv"))

  # Return the wanted values
  return(Mean_Metric)
  
}

  cli_li("DONE: Stat_signif_mean")
  
  # --- CREATION OF THE FUNCTION OBS_VS_RANDOM  --- #
  
  # obs_vs_random is a function created to determine if a unique observed value is different, inferior or superior to a distribution of other values (Null model for example)
  # Give back the result of the test, as a true or false as well as a p-value
  # The function return NA when all values are equal, like when Beta-diversity is computed intra-plot.
  
  # - Vec is a line (a sample) of a metric result file. It contains in order the columns.
  # rep / sample / type / habitat / sampletype / Metric / Metric_NMX ... 
  # - Metric is the Metric name
  # - Type is used to know if we want to test for a greater "g", less "l" , both "b" or different "t" for two sided, the default.
  
  # WARNING ! For now, only the Type ="b" and "All" was verified. 
  
  Obs_vs_Random <- function(Vec, Metric = Metric, Type ="All", AB) {     
    
    # Split "Vec" between the metadata and the actual values of the metric.
    
    # Create a vector of metadata names
    MD <- c("Rep","Group","Sample","Sample_A","Sample_B","Type","Habitat","Habitat_A","Habitat_B","Spatial","Spatial_A","Spatial_B")
    # Get the MetaData  
    MetaData <- Vec[(names(Vec) %in% MD)]
    # Get the metric Values
    MetricValues <- Vec[!(names(Vec) %in% MD)]
    
    # Transform it into a vector (with unlist to keep the names)
    MetricValues <- unlist(MetricValues)
    
    # Check if all the values are identical (Like intra-plot, All 0 or All 1)
    if (length(unique(MetricValues)) == 1) {
      if (Type == "All" ) {
        Result <- rep(NA,6)
        names(Result) <- c("Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower")
        # Bind back the metadata
        Result <- c(MetaData,Result) %>%
          unlist()
        
      }
      if (Type == "t" ) {Result <- as.double(data.frame(T_Different = NA, P_value_Different = NA))}
      # Break from this loop because the results are already made
      return(Result)
    } 
    
    # Get the observed value
    OBS <- MetricValues[Metric]
    # Sort all the data
    All <- sort(MetricValues)
    # Remove the observed value from the Vec to only keep the values from the Null Model
    MetricValues <- All[!(names(All) %in% Metric)]
    # Calculate the quantile of the Vector
    Q <- quantile(MetricValues,probs = c(0.025,0.05,0.95,0.975))
    
    # Check if the observed value is NA
    if (is.na(OBS) == T) {Result <- NA
    
    } else { 
      
      #### . Two-Sided (Type "t") . ####
      if (Type == "t") {
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(All == OBS)
        # Found the "two" possible p_value (One at the beginning, One at the end)
        P1 <- (P[length(P)]-1)/length(MetricValues)      # The first p-value, if the OBS is lower than the Null model
        
        # Sort again Vec but decreasing
        All <- sort(All, decreasing = T)
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(All == OBS)
        P2 <- (P[length(P)]-1)/length(MetricValues)    # The second p-value, if the OBS is greater than the Null model
        
        # Look if the test is passed or not for each side of the distribution
        if (OBS < Q[1]) {T1 <- TRUE} else {T1 <- FALSE}
        if (OBS > Q[4]) {T2 <- TRUE} else {T2 <- FALSE}
        
        # Resume the tests (If one of the test is TRUE, the OBS value is significantly different from the NULL model)
        if (T1 | T2 == TRUE) {Result <- TRUE} else {Result <- FALSE}
        
        # Keep the lowest p_value
        P <- min(c(P1,P2))
      }
      
      #### . Lower (Type "t") . ####
      if (Type == "l") {
        
        # Look for the position of "OBS" in this vector
        P <- which(All == OBS)
        # Fount the p_value
        P <- (P[length(P)]-1)/length(MetricValues) 
        # Make the test
        if (OBS < Q[2]) {Result <- TRUE} else {Result <- FALSE}}
      
      #### . Greater (Type "g") . ####
      if (Type == "g") {
        
        # Sort all the data
        All <- sort(c(OBS,MetricValues), decreasing = T)
        # Look for the position of "OBS" in this vector
        P <- which(All == OBS)
        # Fount the p_value
        P <- (P[length(P)])/length(MetricValues) 
        # Make the test 
        if (OBS > Q[3]) {Result <- TRUE} else {Result <- FALSE}}
      
      #### . Both lower and greater (Type "b") . ####
      if (Type == "b") {
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(All == OBS)
        # Found the "two" possible p_value (One at the beginning, One at the end)
        P_lower <- (P[length(P)]-1)/length(MetricValues)      # The first p-value, if the OBS is lower than the Null model
        
        # Sort again Vec but decreasing
        All <- sort(All, decreasing = T)
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(All == OBS)
        P_greater <- (P[length(P)]-1)/length(MetricValues)    # The second p-value, if the OBS is greater than the Null model
        
        # Look if the test is passed or not for each side of the distribution
        if (OBS < Q[2]) {T_lower <- 1} else {T_lower <- 0}
        if (OBS > Q[3]) {T_greater <- 1} else {T_greater <- 0}
        
        # --- Look at the results 
        Result <- c(T_greater, P_greater, T_lower, P_lower)
        # Name it
        names(Result) <- c("Test_Greater","P_value_Greater","Test_Lower","P_value_Lower")
        # Bind back the metadata
        Result <- c(MetaData,Result)
        
        # Break from this loop because the results are already made
        return(Result)
      }
      
      #### . All of them (Type "All") . ####
      if (Type == "All") {
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(names(All) == Metric)
        
        # Found the "two" possible p_value (One at the beginning, One at the end)
        P_lower <- (P[length(P)]-1)/length(MetricValues)      # The first p-value, if the OBS is lower than the Null model
        
        # Sort again Vec but decreasing
        All <- sort(All, decreasing = T)
        
        # Look for the position of "OBS" in the "All" vector
        P <- which(names(All) == Metric)
        P_greater <- (P[length(P)]-1)/length(MetricValues)    # The second p-value, if the OBS is greater than the Null model
        
        # Look if the test is passed or not for each side of the distribution
        if (OBS < Q[1]) {T_different1 <- 1} else {T_different1 <- 0}
        if (OBS < Q[2]) {T_lower <- 1} else {T_lower <- 0}
        if (OBS > Q[3]) {T_greater <- 1} else {T_greater <- 0}
        if (OBS > Q[4]) {T_different2 <- 1} else {T_different2 <- 0}
        
        # Resume the tests (If one of the test is TRUE, the OBS value is significantly different from the NULL model)
        if (T_different1 | T_different2 == TRUE) {T_different <- TRUE} else {T_different <- FALSE}
        
        # Keep the lowest p_value
        P_different <- min(c(P_greater,P_lower))
        
        ####  .. Standardized Effect Size .. ####
        
        # Compute the standardized effect size
        # SES <- (OBS - mean(MetricValues)) / sd(MetricValues)
        
        # --- Look at the results 
        Result <- c(T_different, P_different ,T_greater, P_greater, T_lower, P_lower)
        # Name it
        names(Result) <- c("Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower")
        # Bind back the metadata
        Result <- c(MetaData,Result) %>%
          unlist()
        # Change the wanted columns to factor
        
        # Break from this loop because the results are already made
        return(Result)
      }  
      
      # Get the results
      # Result <- t(c(Test = Result, P_Value = P))
      
      # Return the result of the test
      return(Result) 
    }
  }      
  
  cli_li("DONE: Obs_vs_Random")
  
  # --- CREATION OF THE FUNCTION STAT-SIGNIF --- #
  
  # Creation of a function to compute the significance of the metrics by rep, type, habitat and sampletype.  
  # It needs the function obs_vs_random to be run. 
  
  # It takes as input the results of the coMet_Metrics script. A data frame with metadatas, one observed value and N replicates. 
  # It takes also the name of the metric that will be computed. This name will be used to create the name of the output: Example: PD_Faith
  
  Stat_signif <- function(data, Name, Metric, Grouped = T, AB = "Alpha", Type = "All"){
    
    # Apply obs_vs_random on the dataset 
    OBSvsRND <- foreach(i = 1:nrow(data), .combine = rbind) %do% {
      Signif_Test <- unlist(Obs_vs_Random(Vec = data[i,], Type = "All", Metric = Metric))} %>%
      # Make it a dataframe
      as.data.frame() %>%
      # Transform the wanted columns into the good object
      mutate(across(c("Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"),as.numeric)) %>%
      mutate_if(is.numeric, round, digits = 2)
    
    # Group the data for further analysis depending on the type "Alpha" or "Beta" of the metric. 
    
    # -- ALPHA METRICS -- # 
    
    if (AB == "Alpha"){
      OBSvsRND <- OBSvsRND %>%
        group_by(Rep,Habitat,Spatial) 
      
      # Compute the mean percentage of significant tests, splitted by group
      if(Grouped == T) {
        Grouped_Signif <- stats::aggregate(
          x = subset(OBSvsRND,select = -c(Rep,Sample,Group,Habitat,Spatial)),       # Only select the rows we effectively want to compute the mean
          by = list(OBSvsRND$Rep,OBSvsRND$Group), # Select the grouping columns (The order of the grouping columns is important !)
          function(x) {mean(x, na.rm = TRUE) %>% round(digits = 2)},                          # Take the mean and round it
          simplify = TRUE, drop = TRUE) %>%
          setNames(nm = c("Rep","Group","Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"))
        
        # If we want to compute the means per replicates, not splitted by group  
      } else if (Grouped == F) {
        Grouped_Signif <- stats::aggregate(
          x = subset(OBSvsRND,select = -c(Rep,Sample,Group,Habitat,Spatial)),       # Only select the rows we effectively want to compute the mean
          by = list(OBSvsRND$Rep), # Select the grouping columns (The order of the grouping columns is important !)
          function(x) {mean(x, na.rm = TRUE) %>% round(digits = 2)},                          # Take the mean and round it
          simplify = TRUE, drop = TRUE) %>%
          setNames(nm = c("Rep","Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"))
        
      } else {return(OBSvsRND)}
      
    } # End of Alpha
    
    # -- BETA METRICS -- # 
    
    if (AB == "Beta"){
      OBSvsRND <- OBSvsRND %>%
        group_by(Rep,Group) 
      
      # Compute the mean percentage of significant tests
      if(Grouped == T) {
        
        Grouped_Signif <- stats::aggregate(
          x = subset(OBSvsRND,select = c(Test_Different,P_value_Different,Test_Greater,P_value_Greater,Test_Lower,P_value_Lower)),       # Only select the rows we effectively want to compute the mean
          by = list(OBSvsRND$Rep,OBSvsRND$Group), # Select the grouping columns (The order of the grouping columns is important !)
          function(x) {mean(x, na.rm = TRUE) %>% round(digits = 2)},                          # Take the mean and round it
          simplify = TRUE, drop = TRUE) %>%
          setNames(nm = c("Rep","Group","Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"))
        
      } else if (Grouped == F) {
        
        Grouped_Signif <- stats::aggregate(
          x = subset(OBSvsRND,select = c(Test_Different,P_value_Different,Test_Greater,P_value_Greater,Test_Lower,P_value_Lower)),
          by = list(OBSvsRND$Rep,OBSvsRND$Group),  # Select the grouping columns (The order of the grouping columns is important !)
          function(x) {mean(x, na.rm = TRUE) %>% round(digits = 2)},                          # Take the mean and round it
          simplify = TRUE, drop = TRUE) %>%
          setNames(nm = c("Rep","Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"))
        
      } else {return(OBSvsRND)}
      
    } # End of Beta
    
    return(Grouped_Signif)
  }
  
  # Message
  cli_li("DONE: Stat_signif")
  
  # --- CREATION OF THE FUNCTION STAT-SIGNIF-MEAN_2 --- #
  
  # Function to compute the significativity of the Mean observed values per group against the Ntree mean values of the Null Model
  
  Stat_signif_mean_2 <- function(data, Metric){
    
    # Apply Obs_VS_Random on the data
    OBSvsRND <- foreach(i = 1:nrow(data), .combine = rbind) %dopar% {
      Signif_Test <- Obs_vs_Random(Vec = data[i,], Type = "All", Metric = Metric)} %>%
      as.data.frame() %>%
      # Transform the wanted columns into the good object
      mutate(across(c("Test_Different","P_value_Different","Test_Greater","P_value_Greater","Test_Lower","P_value_Lower"),as.numeric)) %>%
      mutate_if(is.numeric, round, digits = 2)
    
    return(OBSvsRND)
    
  }
  
  # Message
  cli_li("DONE: Stat_signif_mean_2") ; cli_end()

# Print a message
cli_alert_success("Functions correctly created !")

# -------------------------------------------- #
##### STEP 1: MERGING THE METRICS RESULTS  #####
# -------------------------------------------- #


##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "MERGING OF THE METRICS RESULTS", line_col = "green", line = "-", col = "br_green"))

# ----- #

# -- PD (Faith 1992) and avPD (Clarke 2001) ----

# Do we want to compute PD and avPD ? 
if (PD == FALSE) {
  cli_alert_warning("PD and avPD not computed")
} else {

# ----- # 
cat(rule(left = "Merging of PD and avPD", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PD"),pattern = "_PD_",recursive = T, full.names = T)
suppressMessages(PD_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PD"),pattern = "_avPD_",recursive = T, full.names = T)
suppressMessages(avPD_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE -------------------------------------------------------------- #

write.csv(PD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PD_Faith.csv"))
write.csv(avPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPD_Faith.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_PD <- Stat_signif_mean(PD_Total,"PD_Faith")
Mean_avPD <- Stat_signif_mean(avPD_Total,"avPD_Faith")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = PD_Total, Name = "PD", Metric = "PD")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PD_Faith_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PD, Metric = "PD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PD_Faith_Mean_Signif.csv"))

Signif <- Stat_signif(data = avPD_Total, Name = "avPD", Metric = "avPD")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPD_Faith_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_avPD, Metric = "avPD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPD_Faith_Mean_Signif.csv"))

# --- #
cat(rule(left = "DONE: PD and avPD", line_col = "yellow", line = "-", col = "br_yellow")) 
# --- #

}

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# -- F (Izsak 2000) and J (Izsak 2000) ----

# Do we want to compute F and J ? 
if (FJ == FALSE) {
  cli_alert_warning("F and J not computed")
} else {

# ----- # 
cat(rule(left = "Merging of F and J", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/F"),pattern = "_F_",recursive = T, full.names = T)
suppressMessages(F_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/F"),pattern = "_J_",recursive = T, full.names = T)
suppressMessages(J_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE -------------------------------------------------------------- #

write.csv(F_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_F_Izsak.csv"))
write.csv(J_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_J_Izsak.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_F <- Stat_signif_mean(F_Total,"F_Izsak")
Mean_J <- Stat_signif_mean(J_Total,"J_Izsak")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = F_Total, Name = "F_", Metric = "F_")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_F_Izsak_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_F, Metric = "F_")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_F_Izsak_Mean_Signif.csv"))

Signif <- Stat_signif(data = J_Total, Name = "J_", Metric = "J_")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_J_Izsak_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_J, Metric = "J_")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_J_Izsak_Mean_Signif.csv"))

# --- #
cat(rule(left = "DONE: F and J", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of F and J

# ----- #

# -- PSV, PSR, and PSE (Helmus, 2007) ----

# Do we want to compute PSV, PSR and PSE ? 
if (PSV == FALSE) {
  cli_alert_warning("PSV, PSR and PSE not computed")
} else {
  
# ----- # 
cat(rule(left = "Merging of PSV, PSR and PSE", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV"),pattern = "_PSV_",recursive = T, full.names = T)
suppressMessages(PSV_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV"),pattern = "_PSR_",recursive = T, full.names = T)
suppressMessages(PSR_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PSV"),pattern = "_PSE_",recursive = T, full.names = T)
suppressMessages(PSE_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(PSV_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSV_Helmus.csv"))
write.csv(PSR_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSR_Helmus.csv"))
write.csv(PSE_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSE_Helmus.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_PSV <- Stat_signif_mean(PSV_Total,"PSV_Helmus")
Mean_PSR <- Stat_signif_mean(PSR_Total,"PSR_Helmus")
Mean_PSE <- Stat_signif_mean(PSE_Total,"PSE_Helmus")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = PSV_Total, Name = "PSV", Metric = "PSVs")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSV_Helmus_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PSV, Metric = "PSVs")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSV_Helmus_Mean_Signif.csv"))

Signif <- Stat_signif(data = PSR_Total, Name = "PSR", Metric = "PSR")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSR_Helmus_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PSR, Metric = "PSR")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSR_Helmus_Mean_Signif.csv"))

Signif <- Stat_signif(data = PSE_Total, Name = "PSE", Metric = "PSEs")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSE_Helmus_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PSE, Metric = "PSEs")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PSE_Helmus_Mean_Signif.csv"))

# --- #
cat(rule(left = "DONE: PSV, PSR and PSE", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of PSV, PSR and PSE

# ----- #

# -- PDab (Vellend 2010), ΔnPD (Barker 2002) and avPDab (Tucker 2017) ----

if (PDab == FALSE) {
  cli_alert_warning("PDab, DeltaNPD and avPDab not computed")
} else {
  
# ----- # 
cat(rule(left = "Merging of PDab, DnPD and avPDab", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab"),pattern = "_PDab_",recursive = T, full.names = T)
suppressMessages(PDab_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab"),pattern = "_DnPD_",recursive = T, full.names = T)
suppressMessages(DnPD_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDab"),pattern = "_avPDab_",recursive = T, full.names = T)
suppressMessages(avPDab_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(PDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDab_Vellend.csv"))
write.csv(DnPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_DnPD_Barker.csv"))
write.csv(avPDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPDab_Tucker.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_PDab <- Stat_signif_mean(PDab_Total,"PDab_Vellend")
Mean_DnPD <- Stat_signif_mean(DnPD_Total,"DnPD_Barker")
Mean_avPDab <- Stat_signif_mean(avPDab_Total,"avPDab_Tucker")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = PDab_Total, Name = "PDab", Metric = "PDab")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDab_Vellend_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PDab, Metric = "PDab")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDab_Vellend_Mean_Signif.csv"))

Signif <- Stat_signif(data = DnPD_Total, Name = "DnPD", Metric = "DnPD")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_DnPD_Barker_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_DnPD, Metric = "DnPD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_DnPD_Barker_Mean_Signif.csv"))

Signif <- Stat_signif(data = avPDab_Total, Name = "avPDab", Metric = "avPDab")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPDab_Tucker_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_avPDab, Metric = "avPDab")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_avPDab_Tucker_Mean_Signif.csv"))

# --- #
cat(rule(left = "DONE: PDab, DnPD and avPDab", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of PDab, DnPD and avPDab 


# -- MPD, MPDab and InterMPDab (Clarke & Warwick 1998) ----

if (MPD == FALSE) {
  cli_alert_warning("MPD, MPDab and InterMPDab not computed")
} else {
  
# ----- # 
cat(rule(left = "Merging of MPD and Dp", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MPD"),pattern = "_MPD_",recursive = T, full.names = T)
suppressMessages(MPD_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MPD"),pattern = "_Dp_",recursive = T, full.names = T)
suppressMessages(Dp_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(MPD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MPD_Clarke.csv"))
write.csv(Dp_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Dp_Hardy.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_MPD <- Stat_signif_mean(MPD_Total,"MPD_Clarke")
Mean_Dp <- Stat_signif_mean(Dp_Total,"Dp_Hardy")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = MPD_Total, Name = "MPD", Metric = "MPD")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MPD_Clarke_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_MPD, Metric = "MPD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MPD_Clarke_Mean_Signif.csv"))

Signif <- Stat_signif(data = Dp_Total, Name = "Dp", Metric = "Dp")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Dp_Hardy_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_Dp, Metric = "Dp")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Dp_Hardy_Mean_Signif.csv"))

# ---------------------------- #

# --- #
cat(rule(left = "DONE: MPD and Dp", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of MPD

# -- MNTD and MNTDab (Webb et al, 2002,2008) ----

if (MNTD == FALSE) {
  cli_alert_warning("MNTD and MNTDab not computed")
} else {
  
# ----- # 
cat(rule(left = "Merging of MNTD and MNTDab", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MNTD"),pattern = "_MNTD_",recursive = T, full.names = T)
suppressMessages(MNTD_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/MNTD"),pattern = "_MNTDab_",recursive = T, full.names = T)
suppressMessages(MNTDab_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(MNTD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTD_Webb.csv"))
write.csv(MNTDab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTDab_Webb.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_MNTD <- Stat_signif_mean(MNTD_Total,"MNTD_Webb")
Mean_MNTDab <- Stat_signif_mean(MNTDab_Total,"MNTDab_Webb")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = MNTD_Total, Name = "MNTD", Metric = "MNTD")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTD_Webb_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_MNTD, Metric = "MNTD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTD_Webb_Mean_Signif.csv"))

Signif <- Stat_signif(data = MNTDab_Total, Name = "MNTDab", Metric = "MNTDab")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTDab_Webb_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_MNTDab, Metric = "MNTDab")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MNTDab_Webb_Mean_Signif.csv"))

# ---------------------------- #

# --- #
cat(rule(left = "DONE: MNTD and MNTDab", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of MNTD and MNTDab

# --------------------------------------------------------------------------------------------------------------------------------------------- #    

# -- ED, AED and MED (Tucker 2017) ----

if (ED == FALSE) {
  cli_alert_warning("ED, AED and MED not computed")
} else {
  
# ----- # 
cat(rule(left = "Merging of ED, AED and MED", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED"),pattern = "_ED_",recursive = T, full.names = T)
suppressMessages(ED_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED"),pattern = "_AED_",recursive = T, full.names = T)
suppressMessages(AED_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/ED"),pattern = "_MED_",recursive = T, full.names = T)
suppressMessages(MED_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(ED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_ED_Tucker.csv"))
write.csv(MED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MED_Tucker.csv"))
write.csv(AED_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_AED_Tucker.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_ED <- Stat_signif_mean(ED_Total,"ED_Tucker")
Mean_AED <- Stat_signif_mean(AED_Total,"AED_Tucker")
Mean_MED <- Stat_signif_mean(MED_Total,"MED_Tucker")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = ED_Total, Name = "ED", Metric = "ED")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_ED_Tucker_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_ED, Metric = "ED")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_ED_Tucker_Mean_Signif.csv"))

Signif <- Stat_signif(data = AED_Total, Name = "AED", Metric = "AED")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_AED_Tucker_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_AED, Metric = "AED")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_AED_Tucker_Mean_Signif.csv"))

Signif <- Stat_signif(data = MED_Total, Name = "MED", Metric = "MED")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MED_Tucker_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_MED, Metric = "MED")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_MED_Tucker_Mean_Signif.csv"))

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

# ----- # 
cat(rule(left = "Merging of PDb", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PDb"),pattern = "_PDbeta_",recursive = T, full.names = T)
suppressMessages(PDBeta_Total <- foreach(j = 1:length(File), .combine = full_join) %do% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(PDBeta_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDbeta_Nipperess.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_PDBeta <- Stat_signif_mean(PDBeta_Total ,"PDbeta_Nipperess")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = PDBeta_Total, Name = "PDBeta", Metric = "PDb", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDbeta_Nipperess_Signif.csv"))

Mean_Signif <- Stat_signif_mean_2(data = Mean_PDBeta, Metric = "PDb")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PDbeta_Nipperess_Mean_Signif.csv"))

# ---------------------------- #

# --- #
cat(rule(left = "DONE: PD_Beta", line_col = "yellow", line = "-", col = "br_yellow")) 

} # End of PDbeta

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# -- PHYLOSOR and UNIFRAC (Leprieur 2012) ----

# Do we want to compute the PD (as well as avPD) ? 
if (UniPhylo == FALSE) {
  cli_alert_warning("PhyloSor and UniFrac not computed")
} else {

# ----- # 
cat(rule(left = "Merging of PhyloSor and UniFrac", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_PhyloSor_Bryant",recursive = T, full.names = T)
suppressMessages(PhyloSor_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(PhyloSor_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Bryant.csv"))
Mean_PhyloSor <- Stat_signif_mean(PhyloSor_Total ,"PhyloSor_Bryant")
Signif <- Stat_signif(data = PhyloSor_Total, Name = "PhyloSor", Metric = "PhyloSor", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Bryant_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_PhyloSor, Metric = "PhyloSor")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Bryant_Mean_Signif.csv"))
rm(PhyloSor_Total) ; gc()
cli_alert_success("PhyloSor")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_PhyloSor_Turn_",recursive = T, full.names = T)
suppressMessages(PhyloSor_Turn_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(PhyloSor_Turn_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Turn_Leprieur.csv"))
Mean_PhyloSor_Turn <- Stat_signif_mean(PhyloSor_Turn_Total ,"PhyloSor_Turn_Leprieur")
Signif <- Stat_signif(data = PhyloSor_Turn_Total, Name = "PhyloSor_Turn", Metric = "PhyloSor_Turn", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Turn_Leprieur_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_PhyloSor_Turn, Metric = "PhyloSor_Turn")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Turn_Leprieur_Mean_Signif.csv"))
rm(PhyloSor_Turn_Total) ; gc()
cli_alert_success("PhyloSor_Turn")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_PhyloSor_Ab_",recursive = T, full.names = T)
suppressMessages(PhyloSor_Ab_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(PhyloSor_Ab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Ab_Nipperess.csv"))
Mean_PhyloSor_Ab <- Stat_signif_mean(PhyloSor_Ab_Total ,"PhyloSor_Ab_Nipperess")
Signif <- Stat_signif(data = PhyloSor_Ab_Total, Name = "PhyloSor_Ab", Metric = "PhyloSor_Ab", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Ab_Nipperess_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_PhyloSor_Ab, Metric = "PhyloSor_Ab")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PhyloSor_Ab_Nipperess_Mean_Signif.csv"))
rm(PhyloSor_Ab_Total) ; gc()
cli_alert_success("PhyloSor_Ab")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_UniFrac_Lozupone",recursive = T, full.names = T)
suppressMessages(UniFrac_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(UniFrac_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Lozupone.csv"))
Mean_UniFrac <- Stat_signif_mean(UniFrac_Total ,"UniFrac_Lozupone")
Signif <- Stat_signif(data = UniFrac_Total, Name = "UniFrac", Metric = "UniFrac", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Lozupone_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_UniFrac, Metric = "UniFrac")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Lozupone_Mean_Signif.csv"))
rm(UniFrac_Total) ; gc()
cli_alert_success("UniFrac")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_UniFrac_Turn_",recursive = T, full.names = T)
suppressMessages(UniFrac_Turn_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(UniFrac_Turn_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Turn_Leprieur.csv"))
Mean_UniFrac_Turn <- Stat_signif_mean(UniFrac_Turn_Total ,"UniFrac_Turn_Leprieur")
Signif <- Stat_signif(data = UniFrac_Turn_Total, Name = "UniFrac_Turn", Metric = "UniFrac_Turn", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Turn_Leprieur_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_UniFrac_Turn, Metric = "UniFrac_Turn")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Turn_Leprieur_Mean_Signif.csv"))
rm(UniFrac_Turn_Total) ; gc()
cli_alert_success("UniFrac_Turn")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PhyloFroc"),pattern = "_UniFrac_Ab_",recursive = T, full.names = T)
suppressMessages(UniFrac_Ab_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(UniFrac_Ab_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Ab_Leprieur.csv"))
Mean_UniFrac_Ab <- Stat_signif_mean(UniFrac_Ab_Total ,"UniFrac_Ab_Leprieur")
Signif <- Stat_signif(data = UniFrac_Ab_Total, Name = "UniFrac_Ab", Metric = "UniFrac_Ab", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Ab_Leprieur_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_UniFrac_Ab, Metric = "UniFrac_Ab")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_UniFrac_Ab_Leprieur_Mean_Signif.csv"))
rm(UniFrac_Ab_Total) ; gc()
cli_alert_success("UniFrac_Ab")

# ---------------------------- #

# --- #
cat(rule(left = "DONE: PhyloSor and UniFrac", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Phylosor and Unifrac

# --------------------------------------------------------------------------------------------------------------------------------------------- #        

# -- S_METRICS (Pavoine & Ricotta, 2014) ----

# Do we want to compute the PD (as well as avPD) ? 
if (S == FALSE) {
  cli_alert_warning("S_Metrics not computed")
} else {

# ----- # 
cat(rule(left = "Merging of S_Metrics", line_col = "yellow", line = "-", col = "br_yellow"))

 File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S"),pattern = "_S_SokalSneath_",recursive = T, full.names = T)
 suppressMessages(S_SokalSneath_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
 write.csv(S_SokalSneath_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_SokalSneath_Pavoine.csv"))
 Mean_S_SokalSneath <- Stat_signif_mean(S_SokalSneath_Total ,"S_SokalSneath_Pavoine")
 Signif <- Stat_signif(data = S_SokalSneath_Total, Name = "S_SokalSneath", Metric = "S_SokalSneath", AB = "Beta")
 write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_SokalSneath_Pavoine_Signif.csv"))
 Mean_Signif <- Stat_signif_mean_2(data = Mean_S_SokalSneath, Metric = "S_SokalSneath")
 write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_SokalSneath_Pavoine_Mean_Signif.csv"))
 rm(S_SokalSneath_Total) ; gc()
 cli_alert_success("S_SokalSneath_Total")
File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S"),pattern = "_S_Jaccard_",recursive = T, full.names = T)
suppressMessages(S_Jaccard_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(S_Jaccard_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Jaccard_Pavoine.csv"))
Mean_S_Jaccard <- Stat_signif_mean(S_Jaccard_Total ,"S_Jaccard_Pavoine")
Signif <- Stat_signif(data = S_Jaccard_Total, Name = "S_Jaccard", Metric = "S_Jaccard", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Jaccard_Pavoine_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_S_Jaccard, Metric = "S_Jaccard")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Jaccard_Pavoine_Mean_Signif.csv"))
rm(S_Jaccard_Total) ; gc()
cli_alert_success("S_Jaccard_Total")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S"),pattern = "_S_Sorensen_",recursive = T, full.names = T)
suppressMessages(S_Sorensen_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(S_Sorensen_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Sorensen_Pavoine.csv"))
Mean_S_Sorensen <- Stat_signif_mean(S_Sorensen_Total ,"S_Sorensen_Pavoine")
Signif <- Stat_signif(data = S_Sorensen_Total, Name = "S_Sorensen", Metric = "S_Sorensen", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Sorensen_Pavoine_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_S_Sorensen, Metric = "S_Sorensen")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Sorensen_Pavoine_Mean_Signif.csv"))
rm(S_Sorensen_Total) ; gc()
cli_alert_success("S_Sorensen_Total")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S"),pattern = "_S_Ochiai_",recursive = T, full.names = T)
suppressMessages(S_Ochiai_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(S_Ochiai_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Ochiai_Pavoine.csv"))
Mean_S_Ochiai <- Stat_signif_mean(S_Ochiai_Total ,"S_Ochiai_Pavoine")
Signif <- Stat_signif(data = S_Ochiai_Total, Name = "S_Ochiai", Metric = "S_Ochiai", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Ochiai_Pavoine_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_S_Ochiai, Metric = "S_Ochiai")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Ochiai_Pavoine_Mean_Signif.csv"))
rm(S_Ochiai_Total) ; gc()
cli_alert_success("S_Ochiai_Total")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/S"),pattern = "_S_Beta_",recursive = T, full.names = T)
suppressMessages(S_Beta_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(S_Beta_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Beta_Pavoine.csv"))
Mean_S_Beta <- Stat_signif_mean(S_Beta_Total ,"S_Beta_Pavoine")
Signif <- Stat_signif(data = S_Beta_Total, Name = "S_Beta", Metric = "S_Beta", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Beta_Pavoine_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_S_Beta, Metric = "S_Beta")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_S_Beta_Pavoine_Mean_Signif.csv"))
rm(S_Beta_Total) ; gc()
cli_alert_success("S_Beta_Total")
# ---------------------------- #

# --- #
cat(rule(left = "DONE: S_Metrics", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of S_Metrics

# --------------------------------------------------------------------------------------------------------------------------------------------- #  

# -- Pst, Bst, PIst (Hardy & Senterre, 2007) ----

# Do we want to compute the PIst ? 
if (PIst == FALSE) {
  cli_alert_warning("PIst not computed")
} else {

# ----- # 
cat(rule(left = "Merging of Pst, Bst and PIst", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST"),pattern = "_Pst_",recursive = T, full.names = T)
suppressMessages(Pst_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(Pst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Pst_Hardy.csv"))
Mean_Pst <- Stat_signif_mean(Pst_Total ,"Pst_Hardy")
Signif <- Stat_signif(data = Pst_Total, Name = "Pst", Metric = "Pst", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Pst_Hardy_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_Pst, Metric = "Pst")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Pst_Hardy_Mean_Signif.csv"))
rm(Pst_Total) ; gc()
cli_alert_success("Pst_Total")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST"),pattern = "_Bst_",recursive = T, full.names = T)
suppressMessages(Bst_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(Bst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Bst_Hardy.csv"))
Mean_Bst <- Stat_signif_mean(Bst_Total ,"Bst_Hardy")
Signif <- Stat_signif(data = Bst_Total, Name = "Bst", Metric = "Bst", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Bst_Hardy_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_Bst, Metric = "Bst")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Bst_Hardy_Mean_Signif.csv"))
rm(Bst_Total) ; gc()
cli_alert_success("Bst_Total")

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PST"),pattern = "_PIst_",recursive = T, full.names = T)
suppressMessages(PIst_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})
write.csv(PIst_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PIst_Hardy.csv"))
Mean_PIst <- Stat_signif_mean(PIst_Total ,"PIst_Hardy")
Signif <- Stat_signif(data = PIst_Total, Name = "PIst", Metric = "PIst", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PIst_Hardy_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_PIst, Metric = "PIst")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PIst_Hardy_Mean_Signif.csv"))
rm(PIst_Total) ; gc()
cli_alert_success("PIst_Total")

# ---------------------------- #

# --- #
cat(rule(left = "DONE: Pst, Bst and PIst", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Pst, Bst and PIst

# --------------------------------------------------------------------------------------------------------------------------------------------- #  

# -- PCD (Helmus, 2007) ----

# Do we want to compute the PCD ? 
if (PCD == FALSE) {
  cli_alert_warning("PCD not computed")
} else {

# ----- # 
cat(rule(left = "Merging of PCD", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/PCD"),pattern = "_PCD_",recursive = T, full.names = T)
suppressMessages(PCD_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

write.csv(PCD_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PCD_Helmus.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_PCD <- Stat_signif_mean(PCD_Total ,"PCD_Helmus")

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

Signif <- Stat_signif(data = PCD_Total, Name = "PCD", Metric = "PCD", AB = "Beta")
write.csv(Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PCD_Helmus_Signif.csv"))
Mean_Signif <- Stat_signif_mean_2(data = Mean_PCD, Metric = "PCD")
write.csv(Mean_Signif,paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PCD_Helmus_Mean_Signif.csv"))

# ---------------------------- #

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

# ----- # 
cat(rule(left = "Merging of Simpson", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/Simpson"),pattern = "_Simpson",recursive = T, full.names = T)
suppressMessages(Simpson_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Save the dataframe
write.csv(Simpson_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Simpson.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_Simpson <- Stat_signif_mean(Simpson_Total ,"Simpson")

# ! Remove lines that contains NA for now 
Simpson_Total <- Simpson_Total %>% drop_na()
Mean_Simpson <- Mean_Simpson %>% drop_na()

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

# TEMA LA TAILLE DU GLAND, ON PEUT PAS CALCULER LA SIGNIFICATIVTE SANS MODÈLE NUL BAAAHAHA

cat(rule(left = "DONE: SIMPSON INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Simpson

# Do we want to compute the SIMPSON INDEX? 
if (InvSimpson == FALSE) {
  cli_alert_warning("InvSimpson not computed")
} else {

# ----- # 
cat(rule(left = "Merging of InvSimpson", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/InvSimpson"),pattern = "_InvSimpson",recursive = T, full.names = T)
suppressMessages(InvSimpson_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Save the dataframe
write.csv(InvSimpson_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_InvSimpson.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_InvSimpson <- Stat_signif_mean(InvSimpson_Total ,"InvSimpson")

# ! Remove lines that contains NA for now 
InvSimpson_Total <- InvSimpson_Total %>% drop_na()
Mean_InvSimpson <- Mean_InvSimpson %>% drop_na()

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

# TEMA LA TAILLE DU GLAND, ON PEUT PAS CALCULER LA SIGNIFICATIVTE SANS MODÈLE NUL BAAAHAHA

cat(rule(left = "DONE: INVSIMPSON INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of InvSimpson

# ------------------------------------------------------ #

# Do we want to compute the SIMPSON INDEX? 
if (BC == FALSE) {
  cli_alert_warning("Bray-Curtis distance not computed")
} else {

# ----- # 
cat(rule(left = "Merging of BrayCurtis", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/BC"),pattern = "_BC",recursive = T, full.names = T)
suppressMessages(BC_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Save the dataframe
write.csv(BC_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_BC.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_BC <- Stat_signif_mean(BC_Total ,"BC")

# ! Remove lines that contains NA for now 
BC_Total <- BC_Total %>% drop_na()
Mean_BC <- Mean_BC %>% drop_na()

# --- SIGNIFICATIVITY --------------------------------------------------------------- #

# TEMA LA TAILLE DU GLAND, ON PEUT PAS CALCULER LA SIGNIFICATIVTE SANS MODÈLE NUL BAAAHAHA

cat(rule(left = "DONE: BRAY-CURTIS INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Bray-Curtis

# ------------------------------------------------------ #

# Do we want to compute the SIMPSON INDEX? 
if (SR == FALSE) {
  cli_alert_warning("Species Richness not computed")
} else {

# ----- # 
cat(rule(left = "Merging of Species Richness", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/SR"),pattern = "_SR",recursive = T, full.names = T)
suppressMessages(SR_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Save the dataframe
write.csv(SR_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_SR.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_SR <- Stat_signif_mean(SR_Total ,"SR")

# ! Remove lines that contains NA for now 
SR_Total <- SR_Total %>% drop_na()
Mean_SR <- Mean_SR %>% drop_na()


cat(rule(left = "DONE: BRAY-CURTIS INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Bray-Curtis

# Do we want to compute the SORENSEN INDEX? 
if (Sorensen == FALSE) {
  cli_alert_warning("Sorensen distance not computed")
} else {

# ----- # 
cat(rule(left = "Merging of Sorensen", line_col = "yellow", line = "-", col = "br_yellow"))

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/Sorensen"),pattern = "_Sorensen",recursive = T, full.names = T)
suppressMessages(Sorensen_Total <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/SorensenNest"),pattern = "_SorensenNest",recursive = T, full.names = T)
suppressMessages(Sorensen_Nest <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

File <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Temporary/SorensenTurn"),pattern = "_SorensenTurn",recursive = T, full.names = T)
suppressMessages(Sorensen_Turn <- foreach(j = 1:length(File), .combine = full_join) %dopar% {Raw <- read.csv(File[j],row.names = 1)})

# --- COMBINE AND SAVE ------------------------------------------------------- #

# Save the dataframe
write.csv(Sorensen_Total,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_Sorensen.csv"))
write.csv(Sorensen_Nest,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_SorensenNest.csv"))
write.csv(Sorensen_Turn,paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_SorensenTurn.csv"))

# --- MEAN AND SD VALUES------------------------------------------------------------- #

Mean_Sorensen <- Stat_signif_mean(Sorensen_Total ,"Sorensen")
Mean_SorensenNest <- Stat_signif_mean(Sorensen_Nest ,"SorensenNest")
Mean_SorensenTurn <- Stat_signif_mean(Sorensen_Turn ,"SorensenTurn")

cat(rule(left = "DONE: SORENSEN INDEX", line_col = "yellow", line = "-", col = "br_yellow"))     

} # End of Sorensen

# ------------------------------------------------------ #


# END OF THE SCRIPT
cat(rule(left = "SCRIPT COMET_METRICS_JOIN.R ENDING", line_col = "red", line = "-", col = "br_red")) 


