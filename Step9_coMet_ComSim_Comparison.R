#!/usr/bin/env Rscript 

##### ______________________________________________________________________________________________________________________________________ #####

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");require(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_COMSIM_COMPARISON.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

#-----------------------------------------------------------------------#
##### coMet_Metric_Plots: Drawing of significancy plots for analysis  ###
#-----------------------------------------------------------------------#

# This script is used in the coMet pipeline to compare and draw plots of the simulated communities.  
# It is meant to be used in the command_line using R environnement.

# This script will all the results for all the input scenarios to allow the comparison between all the scenarios. 

# Because this script is separated between types of data, data frames will be loaded and unloaded each time to avoid memory problems. 

# ------------------------------ #
##### STEP 0: INITIALISATION #####
# ------------------------------ #

cat(rule(left = "INITIALISATION", line_col = "red", line = "-", col = "br_red"))

##### ______________________________________________________________________________________________________________________________________ #####

##### . Packages . #####

# Install/load pacman. 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");library(pacman)})
p_load("foreach",
       "dplyr",
       "doParallel",
       "ggplot2",
       "ggpubr",
       "tidyr",
       "stringr",
       "gridExtra",
       "ggplotify",
       "patchwork",
       "ggforce",
       "argparser",
       "purrr",
       "ggridges",
       "ggsignif",
       "ggrepel",
       "rstatix",
       "tibble",
       "vegan",
       "forcats",
       "ape",
       "stringi",
       "viridis",
       "stringr",
       "forcats",
       "argparser",
       "rlist",
       "ggpubr",
       "gridExtra",
       "ggplotify",
       "viridis",
       "paletteer",
       "ggtext",
       "patchwork",
       "scales",
       "pheatmap",
       "pracma",        # Create a logarithmic spaced vector
       "RColorBrewer",
       "vegan",
       "tibble",
       "GGally",
       "ggtree",
       "DescTools"
       )

# Avoid summarise() messages 
options("dplyr.summarise.inform" = FALSE)
# Set lifecycle_verbosity to 'warning', ensuring deprecated messages are shown.
options(lifecycle_verbosity = "warning")

#### . Local Mode . ####

Parameters <- list()
Parameters$Config <- c("coMet_ConfigFiles/Equalizing100_Sig2.R","coMet_ConfigFiles/Equalizing100_Sig5.R","coMet_ConfigFiles/Stabilizing100_Sig5.R","coMet_ConfigFiles/Stochastic100_Sig5.R",
"coMet_ConfigFiles/Equalizing100_Sig2_SF1.R","coMet_ConfigFiles/Equalizing100_Sig2_SF2.R","coMet_ConfigFiles/Equa50Stab50_Sig2.R","coMet_ConfigFiles/Equa50Stab50_Sig5.R")

#### . Argument Parser . ####

# Create the parser
arg_parser <- arg_parser("Check the communities computed by coMet_ComSim and create an output pdf.", name = "coMet_ComSim_Check.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = Inf, help = "Configuration files for coMet_ComSim.R")

# Add the name of the output file
arg_parser <- add_argument(arg_parser, arg = "--output", short = "-o", nargs = 1 , help = "Name of the output file.", default = "Metric_Total.pdf")

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Save the configuration files wanted
Config <- Parameters$Config

# Print a message to inform that the right configuration file was correctly loaded.
cli::cli_alert_info(paste0("Working on configuration files: ",paste(Config,collapse=' ')))

# Print a message
cli_alert_success("coMet_ConfigFile(s).R correctly read !")

    # ----- #

  #### . Plot Theme . ####

# Color Palettes
Sigma_Col <- c("#D72638","#3F88C5","#F49D37")
SampleType_Col <- c("#2A849D","#932F6D")
Habitat_Col <- c("#35C9AB","#6564DB")

# Print a message
cli_alert_info("Setting the theme ... ")

# This theme extends the 'theme_light' that comes with ggplot2.
# The "Lato" font is used as the base font. This is similar
# to the original font in Cedric's work, Avenir Next Condensed.
# theme_set(theme_light(base_family = "Lato"))
theme_set(theme_light())

theme_update(
  # Remove title for both x and y axes
  # axis.title = element_blank(),
  # Axes labels are grey
  axis.text = element_text(color = "grey40"),
  # The size of the axes labels are different for x and y.
  axis.text.x = element_text(size = 10, margin = margin(t = 5)),
  axis.text.y = element_text(size = 10, margin = margin(r = 5)),
  # Also, the ticks have a very light grey color
  axis.ticks = element_line(color = "grey91", linewidth = .5),
  # The length of the axis ticks is increased.
  axis.ticks.length.x = unit(.5, "lines"),
  axis.ticks.length.y = unit(.5, "lines"),
  # Remove the grid lines that come with ggplot2 plots by default
  # panel.grid = element_blank(),
  # Customize margin values (top, right, bottom, left)
  plot.margin = margin(10, 10, 10, 10),
  # Make the background transparent
  panel.background = element_rect(fill = "grey98", color = "grey98"),
  plot.background = element_rect(fill='transparent', color=NA),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  legend.background = element_rect(fill='transparent'),
  legend.box.background = element_rect(fill='transparent'),
  # Customize title appearence
  plot.title = element_text(
    color = "grey10",
    size = 10,
    face = "bold",
    margin = margin(t = 10)
  ),
  # Customize subtitle appearence
  plot.subtitle = element_text(
    color = "grey30",
    size = 8,
    lineheight = 1.35,
    margin = margin(t = 10, b = 20)
  ),
  # Title and caption are going to be aligned
  plot.title.position = "plot",
  plot.caption.position = "plot",
  plot.caption = element_text(
    color = "grey30",
    size = 8,
    lineheight = 1.2,
    hjust = 0,
    margin = margin(t = 20) # Large margin on the top of the caption.
  ),
  # Facet theme
  strip.background = element_rect(
    color="grey40", fill="grey98", linewidth=1, linetype="solid"),
  strip.text = element_text(color = "grey10")
  # Remove legend
  # legend.position = "none"
  # Change the background of the legend
  # legend.background = element_rect(fill = "grey98", color = "grey98"),
  # Add a square around the legend
  # legend.box.background = element_rect(colour = "black")
)

# Create a second theme for the ggarrange plots
arrange_theme <- function() {
  theme(panel.background = element_rect(fill = "grey98", color = "grey98"),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.margin = margin(5, 5, 5, 5),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10",
          size = 12,
          face = "bold",
          margin = margin(t = 10)
          # hjust = 0.5
        ),
        # Customize subtitle appearence
        plot.subtitle = element_text(
          color = "grey30",
          size = 8,
          lineheight = 1.35,
          margin = margin(l = 10 ,t = 10, b = 20)
          # hjust = 0.5
        ),
        # Title and caption are going to be aligned
        plot.caption = element_text(
          color = "grey30",
          size = 8,
          lineheight = 1.2,
          hjust = 0.1,
          margin = margin(t = 10, b = 10) # Large margin on the top of the caption.
        ))

} # End of arrange_theme

    # ----- # 

  ##### . Function Creation. #####

# Create a function "Add_Data" to add the wanted parameters of each scenario into the data frames.
  # - Data is the dataframe on which we want to append the parameters
  # - Parameters is an optional vector of parameters we want to append. 
Add_Data <- function(Data,Parameters = NULL) {
  
  Data <- Data %>% 
    group_by(Scenario) %>%
    # Create a column "Processes" by splitting the scenario name
    mutate(Processes = unlist(stri_split_fixed(str = Scenario, pattern = "_", n =2 ))[2], .after = Scenario) %>% # Split the Scenario name on the "_" and take the second field.
    # Add the wanted parameters at the end of the dataframe
    left_join(Parameters_Total %>% select(all_of(c("Scenario","SigmaValues",Parameters))), by = c("Scenario"))
  
  # Add ! ALL ! the other parameters at the end of the dataframe 
  # left_join(Parameters_Total, by = c("Scenario"))
  
}

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

# ----- #  

# ------------------------------------------- #
##### STEP 1: LOADING THE COMMUNITY FILES #####
# ------------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "LOADING OF THE COMMUNITY FILES", line_col = "green", line = "-", col = "br_green"))

# Read and combine all the simulation results into global summary files between the wanted scenarios. 

    # ----- #

  #### . Loading global files . ####

# Loop on all scenarios to retrieve the wanted datasets and combine them into one

  # - Parameters - # 
suppressMessages(Parameters_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- as.data.frame(t(read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Scenario",Scenario,"_Parameters.csv"),row.names = 1)))
  
}) ; cli_li(cat(col_yellow("LOADED"),": Parameters data.","\n"))

  # - Sampled Sites Data - # 
suppressMessages(Sites_Sampled <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Sites.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1) %>% dplyr::rename(Sample = sample, Rep = rep, Habitat = hab, Spatial = SampleType) %>% # Get the MetaData for the Sample
    # Transform the values in "Spatial_X" column into "C" and "B" for Core (1) and Border (2)
    mutate(Spatial = case_when(Spatial == 1 ~ "C", Spatial == 2 ~ "B")) %>%
    mutate("Group" = paste0(Spatial,Habitat)) %>%
    # Transform what needs to be transformed into factors
    mutate(across(c(Rep,Sample,Group,Habitat,Spatial),as.factor))
  
}) ; cli_li(cat(col_yellow("DONE"),": Sites data.","\n"))

  # - Modification of Data - # 
Parameters_Total <- Parameters_Total %>%
  # Transform some columns into numeric
  mutate(across(c(HabOpt1,HabOpt2), as.numeric))

# Print a message
cli_alert_success("Community files correctly loaded !")

# ------------------------ #
##### STEP 2: PLOTTING #####
# ------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "PLOTTING", line_col = "green", line = "-", col = "br_green"))

##### ------------ 2.A: RECRUITMENT -------------------- #####

# Print a message
cat(rule(left = "// Recruitment Plots //", line_col = "white", line = "-", col = "br_green"))

  #### . Loading files . ####

# - Total Recruitment results of the simulations - # 
suppressMessages(Result_Count_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Result_Count.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1)
  
} %>% Add_Data()) ; cli_li(cat(col_yellow("LOADED"),": Recruitment data.","\n"))

# - Total Recruitment by distance - # 
suppressMessages(Rec_by_Dist_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Rec_by_Dist.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1)
  
} %>% Add_Data()) ; cli_li(cat(col_yellow("LOADED"),": Recruitment by distance data.","\n"))

    # ----- #  

##### ------------ 2.A.a : % Success Total  ------------ #####

  #### . Statistics . ####

# -- Transform wanted columns into factors
Result_Count_Total <- Result_Count_Total %>% mutate(across(c(rep,gen),as.factor))
Rec_by_Dist_Total <- Rec_by_Dist_Total %>% mutate(across(rep,as.factor))

    # ----- #

  ##### .. Ratio Success .. #####

# -- Compute the number of successful recruitment for each category.
sum_success <- Result_Count_Total %>% 
  dplyr::select(contains("nbe_success"),Scenario, gen, rep) %>% 
  pivot_longer(cols = contains("nbe_success"), values_to = "value_success") %>% 
  mutate(name = str_remove(name, "_success"))

# -- Compute the number of failed recruitment for each category.
sum_fail <- Result_Count_Total %>% 
  dplyr::select(contains("nbe_failed"), Scenario, gen, rep) %>% 
  pivot_longer(cols = contains("nbe_failed"), values_to = "value_failed") %>% 
  mutate(name = str_remove(name, "_failed"))

# -- Compute the mean and SD ratio of failed recruitment for each category across replicate.
ratio_success <- 
  left_join(sum_success, sum_fail, by = c("Scenario","gen", "rep", "name")) %>% 
  mutate(prop = value_failed/(value_success + value_failed)*100) %>% 
  group_by(Scenario, gen, name) %>% 
  dplyr::summarise(across(prop, list(mean = mean, sd = sd))) %>%
  # Join again the parameters values
  Add_Data()

# -- Compute the Wilcoxon tests between the values from the same or a different habitat.

# WT_ratio_success <- ratio_success %>%
#   # Filter for habitat difference
#   filter(grepl(paste(c("_diff_hab", "_same_hab"), collapse = "|"), name)) %>%
#   # Group the data
#   group_by(Scenario) %>%
#   # Remove the groups with 0 variance (When there is 100% stochastic processes)
#   mutate(Var_Value = var(prop_mean)) %>%
#   filter(Var_Value != 0) %>%
#   # Compute the kruskall-test
#   wilcox_test(formula = prop_mean ~ name, p.adjust.method = "bonferroni") %>%
#   # Add the y.position for the significance brackets
#   mutate(y.position = 105)

    # ----- #

  ##### .. Stats Mean .. #####

# -- Compute the mean and SD for each numbers across replicate.
summary_stats_sumarized <- 
  Result_Count_Total %>% 
  group_by(Scenario, gen) %>% 
  dplyr::summarise(across(.cols = starts_with("nbe_"), list(mean = mean, sd = sd)))

# -- Compute the mean and SD across replicate of the numbers of recruitment for each category
stats_mean <- 
  summary_stats_sumarized %>% 
  pivot_longer(cols = contains("_mean"), names_to = "stats", values_to = "mean_value") %>% 
  dplyr::select(Scenario, gen, stats, mean_value) %>% 
  mutate(stats = str_remove(stats, "_mean"))

stats_sd <-
  summary_stats_sumarized %>% 
  pivot_longer(cols = contains("_sd"), names_to = "stats", values_to =  "sd_value") %>%
  dplyr::select(Scenario, gen, stats, sd_value) %>% 
  mutate(stats = str_remove(stats, "_sd"))

# -- Joining mean and SD by gen.
Result_Count <- left_join(stats_mean, stats_sd,  by = c("Scenario","gen", "stats")) %>%
  # Join again the parameters values
  Add_Data()

    # ----- #

  #### . Plots . ####

# # New facet label names for processes variable
# Processes <- c("100% Egalisateur", "100% Stabilisateur", "100% Stochastique")
# names(Processes) <- c("100eq", "100stab", "100stoc")
# 
# # New facet label names for sigma variable
# Sigma <- c("Sigma 1", "Sigma 2", "Sigma 5")
# names(Sigma) <- c("1", "2", "5")

# LINE PLOT: Proportion of failed recruitment / habitats ---
Rec_1a <- ggplot(
    # Keep only the columns concerning the recruitment from habitats
    data = ratio_success %>% filter(grepl(paste(c("_diff_hab", "_same_hab"), collapse = "|"), name)), # "name" being a character class for knowing if the recruitment is from the same or a different habitat  
    aes(x = gen, group = name)) +
  # Draw the line plot
  geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
  # Draw the confidence interval
  geom_ribbon(aes(
    y = prop_mean,
    ymin = prop_mean - prop_sd,
    ymax = prop_mean + prop_sd,
    fill = name),
    alpha = .2) +
  # Change the x-axis ticks 
  scale_x_discrete(breaks = pretty(1:as.numeric(max(Parameters_Total$Ngen)),n = 5)) +
  # scale_y_discrete(breaks = pretty(1:100,n = 10)) +
  # Make a facetted plot
  facet_grid(. ~ Scenario ) +
             # labeller = labeller(Processes = Processes, SigmaValues = Sigma)) +
             # labeller = label_both

  # Legends
  scale_color_discrete(name = " ", labels = c("Different habitat","Same habitat")) +
  guides(fill = "none") +
  theme(legend.position="bottom") +
  # Labels
  ylab(label = "Proportion of failed recruitment") +
  xlab(label = "Generation")
  # labs(
  #  title = "Proportion of failed recruitment between habitats."
  #  subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited that could be from a plot of the same or a different habitat.\nEach one of the squares represent one of our scenarios by the combination of both the processes and the sigma values."
  # )

ggsave(filename = "Rec_Habitat_G1.png",
       plot = Rec_1a,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 50,
       height = 20,
       units = "cm",
       dpi = "retina",
       limitsize = F)


# LINE PLOT: Proportion of failed recruitment / plots ---
Rec_2a <- ggplot(
  # Keep only the columns concerning the recruitment from plots
    data = ratio_success %>% filter(grepl(paste(c("_within", "_between"), collapse = "|"), name)),
    aes(x = gen, group = interaction(Scenario,name))) + 
  # Draw the line plot
  geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
  # Draw the confidence interval
  geom_ribbon(aes(
    y = prop_mean,
    ymin = prop_mean - prop_sd,
    ymax = prop_mean + prop_sd,
    fill = name),
    alpha = .2) +
  # Change the x-axis ticks 
  scale_x_discrete(breaks = pretty(1:max(Parameters_Total$Ngen),n = 5)) +
  # Make a faceted plot
  # Make a facetted plot
  facet_grid(. ~ Scenario ) +
            # labeller = labeller(Processes = Processes, SigmaValues = Sigma)) +
            # labeller = label_both
  # Legends
  scale_color_discrete(name = " ", labels = c("Different plot","Same plot")) +
  guides(fill = "none") +
  theme(legend.position="bottom")+
  # Labels
  ylab(label = "Proportion of failed recruitment") +
  xlab(label = "Generation") 
# labs(
#    title = "Proportion of failed recruitment between plots || SigmaValues ~ Processes",
#    subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited that could be from the same or a different plot\nEach one of the squares represent one of our scenarios by the combination of both the processes and the sigma values."
#  )
 
ggsave(filename = "Rec_Plot_G1.png",
       plot = Rec_2a,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 50,
       height = 20,
       units = "cm",
       dpi = "retina",
       limitsize = F)

    # ----- #

  ##### . Patchworking . #####

# # Recruitment by Habitat
# Rec_Habitat <- ggarrange(
#   Rec_1a, Rec_1b, ncol = 1, common.legend = T, legend = "right", heights = c(3,1)) +
#   arrange_theme()
# 
# # Recruitment by Plot
# Rec_Plot <- ggarrange(
#   Rec_2a, Rec_2b, ncol = 1, common.legend = T, legend = "right", heights = c(3,1)) +
#   arrange_theme()

##### ------------ 2.A.b : % Success Distance ------------ #####

  # /!\ For the geom-text, the text is overlapping because one annotation is created per X value. Something has to be found to coreect this. 

options(scipen = 1)
options(digits = 2)

#### . Plots . ####

# First, we summarize the values across the replicates to get the mean and the sd for each distance. For 1 dist, the mean across all the replicates
Rec_Summarize <- Rec_by_Dist_Total %>% group_by(Dist, Scenario) %>% 
  summarise_at(c("Prop_Rec","Prop_success_Rec"),list(Mean = mean, Sd = sd), na.rm = TRUE) %>%
  # Add the MetaData
  Add_Data(Parameters = c("Mno","betadisp","deltadisp"))

# Modify the "Dist" argument to be -1. So intra-plot dispersal corresponds to Dist=0
Rec_Summarize$Dist <- Rec_Summarize$Dist -1

# Version bis with all plots combined. 
Rec_Summarize_Bis <- Rec_by_Dist_Total %>% group_by(Dist) %>% 
  summarise_at(c("Prop_Rec","Prop_success_Rec"),list(Mean = mean, Sd = sd), na.rm = TRUE)

Rec_Summarize_Bis$Prop_Rec_Mean <- round(Rec_Summarize_Bis$Prop_Rec_Mean, 3)

# Modify the "Dist" argument to be -1. So intra-plot dispersal corresponds to Dist=0
Rec_Summarize_Bis$Dist <- Rec_Summarize_Bis$Dist -1

##### .. Recruitment ~ Distance .. #####

# Create a dataset of annotations.

#   # Select the variables used to further create the grid and the graph
# Ann <- select(Parameters_Total,c(Scenario)) %>%
#   # Add a column that is the annotation
#   mutate(
#            # Create the annotation by a combination of the wanted parameters
#            imap(select(Parameters_Total,c(deltadisp, betadisp,Mno)),
#                 ~ paste(.y, .x, sep = " = ")) %>%
#                 as.data.frame() %>%
#                 unite(Annotation, sep=" / ")
#          ) %>%
#   # Add a Dist and Prop_Rec_Mean column to have only one annotation instead of a bunch of overlapping ones, as x and y values for the annotations
#     # max(Dist) is the max distance between to plots on the Tore Grid, we will use half this distance.
#   mutate(x = max(Rec_Summarize$Dist)/2) %>%
#   mutate(y = max(Rec_Summarize$Prop_Rec_Mean))
  
# Plot the proportion of recruitment in function of the distance
# Rec_3 <- ggplot(
#     data = Rec_Summarize, 
#     aes(x = Dist, y = Prop_Rec_Mean, group = 1)) +
#   # Make a faceted plot
#   facet_grid(. ~ Scenario, 
#              labeller = label_both) +
#   # Draw the line
#   geom_line() +
#   # Draw points
#   geom_point(aes(
#     colour = Prop_Rec_Mean, 
#     size = Prop_Rec_Mean))
#   # Repel the labels
#   geom_text_repel(
#     # Only keep the positive mean of successful recruitment
#     data = . %>% filter(Prop_Rec_Mean > 0),
#     aes(label = Prop_Rec_Mean),
#     force_pull   = 0, # do not pull toward data points
#     nudge_x      = 2,
#     direction    = "y",
#     xlim  = c(0.2,NA),         # Define a zone where the labels should be present
#     ylim  = c(0.1,NA),
#     min.segment.length = 0) +  # Always draw line segments
#   # Legend
#   theme(legend.position = "none") +
#   # Labels
#   ylab(label = "Proportion of Recruitment") +
#   xlab(label = "Distance") +
#   labs(
#     title = "Proportion of recruitment by distance. ",
#     subtitle = "The values represent the proportion of recruitments between plots in function of the distance.\nThe distance is in plot unit between the plot receiving an individual and the plot giving this individual accross all the meta-community. Values of the proportion of recruitment were averaged across the replicates."
#   )

# Version with all scenarios combined. 
Rec_3bis <- ggplot(
  data = Rec_Summarize_Bis, 
  aes(x = Dist, y = Prop_Rec_Mean, group = 1)) +
  # Make a faceted plot
  # facet_grid(SigmaValues ~ Processes, 
  #            labeller = label_both) +
  # Draw the line
  geom_line() +
  # Draw points
  geom_point(aes(
    colour = Prop_Rec_Mean, 
    size = Prop_Rec_Mean)) +
  # Repel the labels
  geom_text_repel(
    # Only keep the positive mean of successful recruitment
    data = . %>% filter(Prop_Rec_Mean > 0),
    aes(label = round(Prop_Rec_Mean,4)),
    force_pull   = 0, # do not pull toward data points
    nudge_x      = 2,
    direction    = "y",
    xlim  = c(0.2,NA),         # Define a zone where the labels should be present
    ylim  = c(0.1,NA),
    min.segment.length = 0) +  # Always draw line segments
  # Add the values of Mno, deltadisp and betadisp on the graphs
  #geom_text(data = Ann, 
  #           aes(x = x, y = y, label = Annotation), 
  #          vjust = "inward", hjust = "inward",
  #          inherit.aes=FALSE) +
  # Legend
  theme(legend.position = "none") +
  # Labels
  ylab(label = "Proportion of Recruitment") +
  xlab(label = "Distance")
  # labs(
  #  title = "Proportion of recruitment by distance. ",
  #  subtitle = "The values represent the proportion of recruitments between plots in function of the distance.\nThe distance is in plot unit between the plot receiving an individual and the plot giving this individual accross all the meta-community. Values of the proportion of recruitment were averaged across the replicates."
  # )

ggsave(filename = "Recruitment_Distance.png",
       plot = Rec_3bis,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 40,
       height = 20,
       units = "cm")

  ##### .. Successful recruitment ~ Distance .. #####

# Create a dataset of annotations.

# Select the variables used to further create the grid and the graph
# Ann <- select(Parameters_Total,c(Scenario)) %>%
#   # Add a column that is the annotation
#   mutate(
#     # Create the annotation by a combination of the wanted parameters
#     imap(select(Parameters_Total,c(deltadisp, betadisp,Mno)),
#          ~ paste(.y, .x, sep = " = ")) %>%
#       as.data.frame() %>%
#       unite(Annotation, sep=" / ")
#   ) %>%
#   # Add a Dist and Prop_Rec_Mean column to have only one annotation instead of a bunch of overlapping ones, as x and y values for the annotations
#   # max(Dist) is the max distance between to plots on the Tore Grid, we will use half this distance.
#   mutate(x = max(Rec_Summarize$Dist)/2) %>%
#   mutate(y = max(Rec_Summarize$Prop_success_Rec_Mean))

# Plot the proportion of successful recruitment in function of the distance
# Rec_4 <-ggplot(
#   data = Rec_Summarize, 
#   aes(x = Dist, y = Prop_success_Rec_Mean, group = 1)) +
#   # Make a faceted plot
#   facet_grid(. ~ Scenario, 
#              labeller = label_both) +
#   # Draw the line
#   geom_line() +
#   # Draw points
#   geom_point(aes(
#     colour = Prop_success_Rec_Mean, 
#     size = Prop_success_Rec_Mean)) +
#   # Repel the labels
#   geom_text_repel(
#     # Only keep the positive mean of successful recruitment
#     data = . %>% filter(Prop_success_Rec_Mean > 0),
#     aes(label = Prop_success_Rec_Mean),
#     force_pull   = 0, # do not pull toward data points
#     nudge_x      = 2,
#     direction    = "y",
#     xlim  = c(0.2,NA),         # Define a zone where the labels should be present
#     ylim  = c(0.1,NA),
#     min.segment.length = 0) +  # Always draw line segments
#   # Add the values of Mno, deltadisp and betadisp on the graphs
#   geom_text(data = Ann, 
#             aes(x = x, y = y, label = Annotation), 
#             vjust = "inward", hjust = "inward",
#             inherit.aes=FALSE) +
#   # Legend
#   theme(legend.position = "none") +
#   # Labels
#   ylab(label = "Proportion of successful Recruitment") +
#   xlab(label = "Distance")
#   #labs(
#   #  title = "Proportion of successful recruitment by distance. ",
#   #  subtitle = "The values represent the proportion of successful recruitments between plots in function of the distance.\nThe distance is in plot unit between the plot receiving an individual and the plot giving this individual accross all the meta-community. Values of the proportion of recruitment were averaged across the replicates.",
#   #)

# Version with all scenarios combined. 
Rec_4bis <- ggplot(
  data = Rec_Summarize_Bis, 
  aes(x = Dist, y = Prop_success_Rec_Mean, group = 1)) +
  # Make a faceted plot
  # facet_grid(SigmaValues ~ Processes, 
  #            labeller = label_both) +
  # Draw the line
  geom_line() +
  # Draw points
  geom_point(aes(
    colour = Prop_success_Rec_Mean, 
    size = Prop_success_Rec_Mean)) +
  # Repel the labels
  geom_text_repel(
    # Only keep the positive mean of successful recruitment
    data = . %>% filter(Prop_success_Rec_Mean > 0),
    aes(label = round(Prop_success_Rec_Mean,4)),
    force_pull   = 0, # do not pull toward data points
    nudge_x      = 2,
    direction    = "y",
    xlim  = c(0.2,NA),         # Define a zone where the labels should be present
    ylim  = c(0.1,NA),
    min.segment.length = 0) +  # Always draw line segments
  # Add the values of Mno, deltadisp and betadisp on the graphs
  #geom_text(data = Ann, 
  #           aes(x = x, y = y, label = Annotation), 
  #          vjust = "inward", hjust = "inward",
  #          inherit.aes=FALSE) +
  # Legend
  theme(legend.position = "none") +
  # Labels
  ylab(label = "Proportion of successful recruitment") +
  xlab(label = "Distance") 

ggsave(filename = "RecruitmentSuccess_Distance.png",
       plot = Rec_4bis,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

# Print a message
cli_alert_success("DONE: Recruitment plots.")

# Clean the work space

    # ----- #

  ##### . Saving . #####

# Recruitment by Distance
Rec_Dist <- ggarrange(
  Rec_3bis, Rec_4bis, ncol = 1, legend = NULL) +
  arrange_theme()


# Save the plots
# ggsave(filename = "Recruitment_Distance.png",
#        plot = Rec_3,
#        device = "png",
#        path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
#        width = 30,
#        height = 20,
#        units = "cm")

ggsave(filename = "Recruitment_Distance_Total.png",
       plot = Rec_Dist,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")


##### ------------ 2.B: ABUNDANCE / OCCURENCE -------------------- #####

# Print a message
cat(rule(left = "// Abundance and Occurence Plots //", line_col = "white", line = "-", col = "br_green"))

#### . Loading files . ####

# - Total Abundance - # 
suppressMessages(Abundance_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Abundance.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1) %>%
    # Rename
    dplyr::rename(Rep = rep, Sample = sample) %>%
    # Transform what needs to be transformed into factors
    mutate(across(c(Rep,Sample),as.factor)) %>%
    # Get the MetaData
    left_join(y = Sites_Sampled, join_by("Scenario","Sample","Rep")) %>%
    # Change the columns order
    relocate(starts_with("t"), .after = last_col())
  
  
}) ; cli_li(cat(col_yellow("DONE"),": Abundance data.","\n"))

# - Total Occurrence - # 
suppressMessages(Occurence_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Occurence.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1) %>%
    # Rename
    dplyr::rename(Rep = rep, Sample = sample) %>%
    # Transform what needs to be transformed into factors
    mutate(across(c(Rep,Sample),as.factor)) %>%
    # Get the MetaData
    left_join(y = Sites_Sampled, join_by("Scenario","Sample","Rep")) %>%
    # Change the columns order
    relocate(starts_with("t"), .after = last_col())
  
  
  }) ; cli_li(cat(col_yellow("DONE"),": Occurence data.","\n"))

    # ----- #  

  #### . Statistics . ####

# Compute total species abundance per replicates.
Abundance_Total_Plot <- Abundance_Total %>%
  # Group 
  group_by(Scenario,Rep,Group) %>%
  # Compute the sum of each species across all samples for each replicates
  dplyr::summarise_if(is.numeric,sum, na.rm = TRUE) %>%
  # Compute the total number of individuals for each replicates. It is meant to be the same for each replicate of a Scenario, because we have Nip x L x W individuals. 
  mutate(Total_Sum = rowSums(pick(starts_with("t"))),.after = Rep) %>% 
  # Compute the Simpson diversity of each of the replicates.
  mutate(Simpson = diversity(pick(starts_with("t")),index = "simpson"),.after = Total_Sum)

# Compute the sum of each species abundance across each plots for each replicates. 
 # Sum by columns. This is not plot abundance but species abundance across all the plots.
Abundance_Total_Species <- Abundance_Total %>%
  # Remove sample column
  dplyr::select(-c(Sample,X,Y,Habitat,Spatial)) %>% 
  # Group
  group_by(Rep,Scenario,Group) %>%
  # Compute the sum of each species across all samples for each replicates
  dplyr::summarise_if(is.numeric,sum, na.rm = TRUE) %>%
  # Get the species as a column and the values as another
  pivot_longer(cols = starts_with("t"),
               names_to = "Species",
               values_to = "Abundance")


  ##### .. BP: Diversity .. #####

# Plot the diversity of the replicates (With ggpubr)
BP_Diversity <- ggboxplot(
               Abundance_Total_Plot, 
               x = "Group", y = "Simpson",
               color = "Group", 
               # palette = Sigma_Col,
               add = "jitter",
               facet.by = "Scenario",
               ggtheme = arrange_theme()) + 
  # Add the results of the kruskall Wallis tests.
  # geom_pwc(label = "P-value = {p.format} / {p.signif}",
  #         hide.ns = F) +
  # Legend
  theme(legend.position = "none") +
  # Labels
  ylab(label = "Simpson Diversity") +
  xlab(label = "SigmaValues")
  # labs(
  #   title = "Boxplot of the Simpson diversity || . ~ Processes",
  #   subtitle = "Kruskall-Wallis tests are realized between the Simpson Diversity computed for different SigmaValues for each Processes.\nEach point represent one replicate of the Scenario."
  # )

    # ----- #

  ##### . Saving . #####

# Save the plots
ggsave(filename = "BP_Diversity_G1.png",
       plot = BP_Diversity,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

##### ------------ 2.C: RICHNESS -------------------- #####

# Print a message
cat(rule(left = "// Richness Plots //", line_col = "white", line = "-", col = "br_green"))

#### . Loading files . ####

# - Total Species Richness - # 
suppressMessages(SR_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_SR.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1) %>%
    # Rename
    dplyr::rename(Rep = rep, Sample = sample) %>%
    # Transform what needs to be transformed into factors
    mutate(across(c(Rep,Sample),as.factor)) %>%
    # Get the MetaData
    left_join(y = Sites_Sampled, join_by("Scenario","Sample","Rep")) %>%
    # Change the columns order
    relocate(SR_Data, .after = last_col())
  
}) ; cli_li(cat(col_yellow("DONE"),": Species richness data.","\n"))

    # ----- #  

  #### . Statistics . ####

# Compute the mean SR by Processes ~ SigmaValues
Summ_SR_Total <- SR_Total %>%
  # Group
  group_by(Rep,Scenario,Spatial) %>%
  # Compute the mean Abundance ~ Replicate ~ Scenario
  dplyr::summarise(across(SR_Data, mean, na.rm = TRUE))

# Save the value for each scenario
Summ_SR_Total2 <- SR_Total %>%
  # Group
  group_by(Scenario) %>%
  # Compute the mean Abundance ~ Replicate ~ Scenario
  dplyr::summarise(across(SR_Data, list(mean = mean, sd = sd), na.rm = TRUE)) %>%
  # round the results
  round(digits = 2)

  #### .. BP: Richness .. ####

# Plot the abundance of the replicates (With ggpubr)
BP_Richness <- ggboxplot(
  Summ_SR_Total, 
  x = "Group", y = "SR_Data",
  color = "Group", 
  # palette = Sigma_Col,
  add = "jitter",
  facet.by = "Scenario",
  ggtheme = arrange_theme()) + 
  # Add the results of the kruskall Wallis tests.
  geom_pwc(label = "P-value = {p.format} / {p.signif}",
            hide.ns = F) +
  # Legend
  theme(legend.position = "none") +
  # Labels
  ylab(label = "Species Richness") +
  xlab(label = "Groups") 
  # labs(
  #  title = "Boxplot of the species richness || . ~ Processes",
  #  subtitle = "Kruskall-Wallis tests are realized between the species richness computed for different SigmaValues for each Processes.\nEach point represent one replicate of the Scenario."
  #)

BP_Richness2 <- ggviolin(
  Summ_SR_Total, 
  x = "Scenario", y = "SR_Data",
  color = "Scenario", 
  # palette = Sigma_Col,
  # add = "jitter",
  # facet.by = "Scenario",
  ggtheme = arrange_theme()) + 
  # Add the results of the kruskall Wallis tests.
  # geom_pwc(label = "P-value = {p.format} / {p.signif}",
  #          hide.ns = F) +
  # Legend
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  # Labels
  ylab(label = "Species Richness") +
  xlab(label = "Scenario") +
  stat_summary(fun = mean, geom = "text", 
               col = "black",     # Add text to plot
               # vjust = 1.5, 
               aes(label = paste("Mean:", round(after_stat(y), digits = 1))))
# labs(
#  title = "Boxplot of the species richness || . ~ Processes",
#  subtitle = "Kruskall-Wallis tests are realized between the species richness computed for different SigmaValues for each Processes.\nEach point represent one replicate of the Scenario."
#)



    # ----- #

  #### .. DP: Richness ~ Group .. ####

# Do the wilcoxon tests  
KT <- compare_means(data = SR_Total, SR_Data ~ Group, group.by = c("Scenario")) %>%
  # Add dummy columns for the next plot to understand where to plot. 
  mutate("SR_Data" = 1, "Group" = 1)

# The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 
DP_Richness <- SR_Total %>%
  ggplot(aes(x=SR_Data, group = Group)) +
  # Draw the plot
  geom_density(aes(color = Group, fill = Group, alpha = 0.2),na.rm = T) + 
  # Draw a facet
  facet_grid(. ~ Scenario) +
  # Significant
  # geom_label(data = KT, 
  #           label = paste0("p-value = ",KT$p.adj," / ",KT$p.signif), 
  #           x = Inf,  y = Inf, hjust = "inward", vjust = "inward") + 
  # Aesthetics
  # scale_color_manual(values = SampleType_Col, labels = c("Core","Ecotone")) +
  # scale_fill_manual(values = SampleType_Col) +
  # Legends
  guides(alpha = "none", fill = "none") +
  theme(legend.position = "bottom")+
  # Labels
  ylab(label = "Density") +
  xlab(label = "Species Richness")
  # labs(
  #  title = "Density of the plot species richness || SampleType ~ SigmaValues + Processes",
  #  subtitle = " "
  #) ; DP_Richness

  # ----- #

    ##### . Saving . #####

# Save the plots
ggsave(filename = "BP_Richness.png",
       plot = BP_Richness,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

ggsave(filename = "BP_Richness2.png",
       plot = BP_Richness2,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

ggsave(filename = "DP_Richness.png",
       plot = DP_Richness,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")



# ----- #

##### ------------ 2.D: TRAITS -------------------- #####

# Print a message
cat(rule(left = "// Traits Plots //", line_col = "white", line = "-", col = "br_green"))

#### . Loading files . ####

# - Total Species Traits - # 
suppressMessages(Traits_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Traits.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1)
  
} %>% Add_Data()) ; cli_li(cat(col_yellow("DONE"),": Traits data.","\n"))

    # ----- #  

  #### . Statistics . ####

# Prepare the data
Traits_Total <- Traits_Total %>%
  # Pivot the data
  pivot_longer( cols = starts_with("Trait"),
                names_to = "Traits",
                values_to = "TraitValues") %>%
  # Add the values of habitat optimums for each scenario
  left_join(select(Parameters_Total, Scenario, HabOpt1, HabOpt2), by = "Scenario") %>%
  # Create a new column to have the trait number
  mutate("Trait_Number" = str_split(Traits,pattern = "_", simplify = T)[,2]) %>% 
  # Change the processes concerned by the traits
  mutate(Processes = case_when(Trait_Number == 1:3 ~ "Stabilizing", 
                               Trait_Number == 4:6 ~ "Equalizing",
                               Trait_Number == 7:9 ~ "Stochastic")) %>%
  # Transform what needs to be transformed into factors
  mutate(across(c(rep,sp),as.factor))
  
  

# Combine the data of traits with only the species that are present (Occurrence Data)
Occurence_Traits <- Occurence_Total %>%
  # Get the species as a column and the values as another
  pivot_longer(cols = starts_with("t"),
               names_to = "Species",
               values_to = "Occurence") %>%
  # Remove all the species that are absent
  filter(!Occurence == 0) %>%
  # Combine with the trait values
  left_join(Traits_Total,join_by("Rep" == "rep","Scenario","Species" == "sp"), relationship = "many-to-many")

  #### . Plots . #### 

# Density plot total ---
DP_Traits <- Traits_Total %>%
  ggplot(aes(x=TraitValues, color = Scenario,  group = Scenario)) +
  # Draw the plot
  geom_density() +
  geom_vline(data = Parameters_Total, aes(xintercept = c(HabOpt1))) +
  geom_vline(data = Parameters_Total, aes(xintercept = c(HabOpt2))) +
  # Legends
  guides(alpha = "none") +
  # Labels
  ylab(label = "Density") +
  xlab(label = "Trait Values") +
  labs(
    title = "Distribution of the traits values for each Scenario for all species.",
    subtitle = " "
  )
  
# Boxplots Trait ~ Group + Processes ---
BP_Traits <- Occurence_Traits %>%
  ggplot(aes(x=Group, y= TraitValues, color = Processes, fill = Processes, alpha = 0.2)) +
  # Draw the plot
  geom_boxplot(outlier.shape = NA) +
  facet_grid(. ~Scenario) +
  # Legends
  guides(alpha = "none") +
  # Labels
  ylab(label = "Trait Values") +
  xlab(label = "Groups") +
  labs(
    title = "Distribution of the traits values for each Scenario.",
    subtitle = "Distribution of the 3x3 traits of the species present in the samples, splitted by group and type of traits (processus involved). Outliers are masked."
  )

    # ----- #

  ##### . Saving . #####

# Save the plots
ggsave(filename = "DP_Traits.png",
       plot = DP_Traits,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

# Save the plots
ggsave(filename = "BP_Traits.png",
       plot = BP_Traits,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 40,
       height = 20,
       units = "cm")


##### ------------ 2.F: FITNESS -------------------- #####

# Print a message
cat(rule(left = "// Fitness Plots //", line_col = "white", line = "-", col = "br_green"))

  #### . Loading files . ####

# - Total Species Fitness - # 
suppressMessages(Fitness_Total <- foreach(Conf = 1:length(Config),.combine = full_join) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Retrieve the number of metrics already computed for each Scenario
  Data <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Fitness.csv"),row.names = 1) %>%
    # Add the scenario name 
    mutate("Scenario" = Scenario, .before = 1)
  
}) ; cli_li(cat(col_yellow("DONE"),": Fitness data.","\n"))

    # ----- #  

  #### . Statistics . ####

# Add the difference of fitness between the two habitats as a new column
Fitness_Total <- Fitness_Total %>%
      ungroup() %>%
      mutate(Abs_Difference = abs(Fitness_Total$habitat_1 - Fitness_Total$habitat_2)) %>%
      # Transform the wanted columns into factor
      mutate_at(c("Scenario","rep"),as.factor)

# Rank the 5 most abundant species per scenario and Group
Abundance_Species_20R <- Abundance_Total_Species %>%
  group_by(Scenario, Group) %>%
  # Keep the 20 most abundant species
  slice_max(order_by = Abundance, n = 5)

# Create a dataset with only these species
Fitness_20R <- right_join(Fitness_Total, Abundance_Species_20R, by = c("rep" = "Rep", "sp" = "Species","Scenario")) %>%
  # Group
  group_by(Scenario, Group) %>%
  # Sort the abundance decreasing
  dplyr::arrange(desc(Abundance), .by_group = TRUE)

    # ----- #

  #### . DP : Absolute difference . ####

DP_Fitness <- 
  ggplot(Fitness_Total, aes(x = Abs_Difference)) +
  geom_density(alpha = 0.2, fill = "grey56") +
  # Draw a facet
  facet_grid(. ~ Scenario, labeller = "label_both") +
  # Labs
  labs(title = "Distribution of the fitness difference. || . ~ SigmaValues + Processes",
       subtitle = "Density plot of the absolute difference of fitness for all species.",
       x ="Absolute difference of fitness",
       y = "Density") 

    # ----- #

  #### . BP : Absolute difference . ####

BP_Fitness <- ggboxplot(
  Fitness_Total, 
  x = "Scenario", y = "Abs_Difference",
  color = "Scenario",
  ggtheme = arrange_theme()) + 
  # Legend
  theme(legend.position = "none") +
  # Labels
  ylab(label = "Absolute fitness difference") +
  xlab(label = "SigmaValues") +
# Aesthetic
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 0.5)) +
  guides(size = "none", fill = "none") +
  arrange_theme()
  # labs(
  #  title = "Boxplot of the absolute fitness difference || . ~ Processes",
  #  subtitle = "Kruskall-Wallis tests are realized between the absolute fitness difference for the two habitats computed for different SigmaValues for each Processes."
  #)


  #### . LP: Absolute difference . ####

# Plot the fitness of the 20 most abundant species for both habitats.
LP_Fitness <- Fitness_20R %>%
  # Draw the plot
  ggplot(aes(x = fct_inorder(sp))) +
  geom_segment(aes(x = fct_inorder(sp), xend = fct_inorder(sp), y=habitat_1, yend=habitat_2), color="grey") +
  geom_point(aes(y = habitat_1, color= "Habitat 1")) +
  geom_point(aes(y = habitat_2, color= "Habitat 2")) +
  # Draw a facet
  facet_grid(Group ~ Scenario, labeller = "label_both") +
  # Scales
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_color_manual(values = c("Habitat 1" = rgb(0.2,0.7,0.1,0.5),"Habitat 2" = rgb(0.7,0.2,0.1,0.5))) +
  labs(title = paste0("Lollipop plot of Species fitness."),
       subtitle = "Relative fitness of the 5 ranked most abundant species in the metacommunity for each scenarios and groups.",
       x = "Species",
       y = "Fitness",
       color = "Habitat")

    # ----- # 

  ##### . Saving . #####

# Save the plots
ggsave(filename = "DP_Fitness.png",
       plot = DP_Fitness,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

ggsave(filename = "BP_Fitness.png",
       plot = BP_Fitness,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 30,
       height = 20,
       units = "cm")

# ----- #

##### ------------ 2.G: HEATMAPS -------------------- #####

# Print a message
cat(rule(left = "// Heatmap Plots //", line_col = "white", line = "-", col = "br_green"))

### .. HM Data .. ###

# Only keep the species that are the most abundant in each habitat for each scenario
Total_Grid <- foreach(Conf = 1:length(Config)) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])
  
  # Load the Whole communities
  Whole_Com <- read.csv(paste0("coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Community_Abundance.csv"), row.names = 1) %>%
    # Select only the first replicate
    filter(rep == unique(.$rep)[1]) %>%
    # In addition, remove the columns of the species not at all presen in that replicate
    select(where(~sum(!is.na(.)) > 0))
  
  # Load the Sites
  Sites <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step2_Sites_Data.csv"),row.names = 1) %>%
    mutate(across(c("hab","SampleType","sample"),as.factor))
  
  cli_alert_info(paste0(Config[Conf]," : Launched"))
  
  # Compute the abundance of each species per replicate. 
  HM_Data <- Whole_Com %>%
    # Add the Sites Data
    cbind(Sites) %>%
    # Transform wanted columns into factors
    mutate_at(vars(hab), factor) %>% # Transform all the stuff into factor
    # Group_by Scenario and Hab
    group_by(hab) %>%
    # Sum the data
    dplyr::summarise(across(starts_with("t"), ~ sum(.x))) %>%
    # Pivot longer the species
    pivot_longer(
      cols = starts_with("t"), names_to = "Species",values_to = "Abundance") %>%
    # Sort them decreasingly
    arrange(desc(Abundance))
  
  # Create a data set with the most abundant species ~ Scenario + Hab
  HM_Data_Best <- HM_Data %>%
    # Group_by Scenario and Hab
    group_by(hab) %>%
    # Keep the most abundant species per group
    slice_head(n = 1)
  

  HM_Data_Lite <- Whole_Com %>%
    # select by the most abundant species
    select(HM_Data_Best$Species) %>%
    # Add the Sites Data
    cbind(.,Sites) %>%
    # Pivot longer the species
    pivot_longer(
      cols = starts_with("t"), names_to = "Species",values_to = "Abundance") %>%
    mutate_at(vars(!Abundance),factor) 
  
  # Select the data on which we are gonna work
  HM_Filter <- HM_Data_Lite %>%
    subset(Species %in% HM_Data_Best$Species)
  
  # Create the labels for the plot
  Facet_Names <- c(paste0("Hab_1 / ",HM_Data_Best$Species[1]," / ",HM_Data_Best$Abundance[1]),paste0("Hab_2 / ",HM_Data_Best$Species[2]," / ",HM_Data_Best$Abundance[2]))
  names(Facet_Names) <- c(HM_Data_Best$Species[1:2])
  
  # Change the color of the plotting window
  par(bg = "grey10")
  
  # Sometimes, the same species is the most abundant for both habitats, therefore only one plot is produced
  HM <- ggplot(HM_Filter, aes(X,Y)) +
    # Grid Color
    geom_point(aes(color = Abundance, fill = Abundance, size = 1, shape = hab)) + # Try t59 for Habitat 2, T419 for hab 1
    scale_shape_manual(values = c(15,22)) +
    scale_color_viridis_c(na.value = "grey80", 
                          option = "plasma",
                          limits = c(0, 50)
    ) +
    scale_fill_viridis_c(na.value = "grey80", 
                         # option = "plasma",
                         limits = c(0, 50)
    ) +
    
    coord_fixed() +
    guides(size = "none", shape = "none") +
    # Facets
    facet_grid(. ~ Species, labeller = labeller(Species=Facet_Names)) +
    # Label
    labs(title = paste0(Scenario)) +
    # Theme
    theme_void() +
    theme(
      # Title
      plot.title = element_text(color="grey80", size=14, face="bold",hjust = 0.5),
      # Background
      plot.background = element_rect(fill = "grey10", color = "grey10"),
      panel.background = element_rect(fill = "grey10", color = "grey10"),
      # Legend
      legend.title = element_text(color = "grey80", size = 10),
      legend.text = element_text(color = "grey80"),
      # aspect.ratio = 1/1,
      # Facet
      strip.text.x =  element_text(size = 10, color = "grey80"),
      strip.text.y =  element_text(size = 10, color = "grey80")
    )
  
  ggsave(filename = paste0("Abundance_HeatMaps_",Scenario,".png"),
         plot = HM,
         device = "png",
         path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots/Heatmaps"),
         width = 40,
         height = 25,
         units = "cm")
    
  # Return the data
  return(HM)

}

# ---- #

# Append them together
pdf(file = NULL)
pl <- lapply(Total_Grid, as.grob)
dev.off()
# Create the pdf
ml <- marrangeGrob(pl,nrow=4,ncol =2 , top = NULL)
dev.size()
# Save the plots
ggsave(filename = "Abundance_HeatMaps.png",
       plot = ml,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/Plots"),
       width = 40,
       height = 120,
       units = "cm")
