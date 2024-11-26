#!/usr/bin/env Rscript 

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");require(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_METRIC_CHECK BEGINNING", line_col = "red", line = "-", col = "br_red")) 

#----------------------------------------------------------------------------------#
##### coMet_Metric_Check: Drawing of metric and significance plots for analysis  ###
#----------------------------------------------------------------------------------#

# This script is used in the coMet pipeline to draw the plots of the results for each metric and their significance tests for one scenario. 
# It contains density plots of the values of the plots as well as plots of the metric significance.  
# It is meant to be used in the command_line using R environnement.

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
       "tidyr",
       "doParallel",
       "ggplot2",
       "ggpubr",
       "stringr",
       "gridExtra",
       "ggplotify",
       "patchwork",
       "ggforce",
       "argparser",
       "purrr")

# Avoid summarise() messages 
options("dplyr.summarise.inform" = FALSE)
# Set lifecycle_verbosity to 'warning', ensuring deprecated messages are shown.
options(lifecycle_verbosity = "warning")

# ----- #

#### . Local Mode . ####

Parameters <- list()
Parameters$Config <- "coMet_ConfigFiles/Stochastic100_Sig5.R"

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

# Print a message
cli_alert_success("coMet_ConfigFile.R correctly read !")

    # ----- #

  #### . Plot Theme . ####

# Print a message
cli_alert_info("Setting the theme ... ")

# This theme extends the 'theme_light' that comes with ggplot2.
# The "Lato" font is used as the base font. This is similar
# to the original font in Cedric's work, Avenir Next Condensed.
# theme_set(theme_light(base_family = "Lato"))
  
# Set one known theme
theme_set(theme_light())
  # Modify it
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
  plot.margin = margin(5, 5, 5, 5),
  # Use a light grey color for the background of both the plot and the panel
  plot.background = element_rect(fill = "grey98", color = "grey98"),
  panel.background = element_rect(fill = "grey98", color = "grey98"),
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
  plot.title.position = "panel",
  plot.caption.position = "panel",
  plot.caption = element_text(
    color = "grey30", 
    size = 8,
    lineheight = 1.2, 
    hjust = 0,
    margin = margin(t = 20) # Large margin on the top of the caption.
  ),
  # Remove legend
  # legend.position = "none"
  # Change the background of the legend
  legend.background = element_rect(fill = "grey98", color = "grey98"),
  # Add a square around the legend
  # legend.box.background = element_rect(colour = "black")
)

# Create a second theme for the ggarrange plots
arrange_theme <- function() {
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.margin = margin(5, 5, 5, 5),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ),
        # Customize subtitle appearence
        plot.subtitle = element_text(
          color = "grey30", 
          size = 8,
          lineheight = 1.35,
          margin = margin(l = 10 ,t = 10, b = 20),
          hjust = 0.5
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

# Print a message
cli_alert_success("Plot theme correctly set !")

#### . Loading files . ####

cat(rule(left = "Loading the community files ...", line_col = "green", line = "-", col = "br_green"))

# Read and combine all the metrics result into two summary files. One for the metric result, the other for the significant result

# Find all the result file for the metric
  # Find the Mean_Signif files
File_Mean_Signif <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics"), pattern = "Mean_Signif\\.csv$" ,recursive = T, full.names = T)
  # Find the Signif files by removing the "_Mean" string
# File_Signif <- gsub(pattern = "_Mean", replacement = "", File_Mean_Signif)
  # Find the Raw_Mean file by removing the "Signif" string
File_Mean_Raw <- gsub(pattern = "_Signif", replacement = "", File_Mean_Signif)
  # Find the Raw file by removing the "Mean_Signif" string
# File_Raw <- gsub(pattern = "_Mean_Signif", replacement = "", File_Mean_Signif)

  # Extract the Metric Name
Metrics <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/"), pattern = "Mean_Signif\\.csv$")
# Remove the name of the scenario
Metrics <- gsub(pattern = paste0("Scenario",Scenario,"_"), replacement = "", Metrics)
# Remove the last parts of the file
Metrics <- gsub(pattern = "_Mean_Signif.csv", replacement = "", Metrics)


# Loop across the found files 
Mean_Signif_Data <- foreach(j = 1:length(File_Mean_Signif)) %do% { Raw <- read.csv(File_Mean_Signif[j],row.names = 1) %>% mutate("Metric" = Metrics[j], .before = 1)} %>% setNames(Metrics)
# Loop across the found files 
# Signif_Data <- foreach(j = 1:length(File_Signif)) %do% { Raw <- read.csv(File_Signif[j],row.names = 1) %>% mutate("Metric" = Metrics[j], .before = 1)} %>% setNames(Metrics)
# Loop across the found files 
Mean_Raw_Data <- foreach(j = 1:length(File_Mean_Raw)) %do% { Raw <- read.csv(File_Mean_Raw[j], row.names = 1) %>% mutate("Metric" = Metrics[j], .before = 1)} %>% setNames(Metrics)
# Loop across the found files 
# Raw_Data <- foreach(j = 1:length(File_Raw)) %do% { Raw <- read.csv(File_Raw[j],row.names = 1) %>% mutate("Metric" = Metrics[j], .before = 1)} %>% setNames(Metrics)

# Print a message
cli_alert_success("Files correctly loaded !")

# ------------------------ #
##### STEP 1: PLOTTING #####
# ------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

# --- #
cat(rule(left = "PLOTTING", line_col = "green", line = "-", col = "br_green"))
# --- #

##### ------------ 1.A: RAW DATA -------------------- #####

cat(rule(left = "Plotting of the raw data ...", line_col = "yellow", line = "-", col = "br_yellow"))

# #### . Function creation . ####
# 
# # Create a function to create the plots for each individual metric.
#   # Raw_Data is the list of results of raw data for the metrics
# Raw_Fct  <- function(Raw_Data) {
#   
#   # Remove any unwanted columns
#   Raw_Data <- Raw_Data %>%
#     select(-any_of(c("Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial")))
#   
#   # Extract the Metric Name
#   Name <- colnames(Raw_Data)[5]
#   
#   # Transform the data to add a column for the null vs obs data
#   Raw_Trans <- Raw_Data %>%
#     pivot_longer(
#       cols = contains(Name, ignore.case = TRUE, vars = NULL),
#       names_to = "OBSvsNM", values_to = "Value")
#   
#   # Change the values contained in the new column
#   Raw_Trans$OBSvsNM[grep("NM",Raw_Trans$OBSvsNM)] <- paste0(Name,"_NM")
#   
#   # Transform all the stuff into factor
#   Raw_Trans <-  Raw_Trans %>% 
#     mutate_at(vars(Metric,Rep,Sample,Group,), factor) %>%
#     as.data.frame()
#   
#   Raw_Data <- Raw_Data %>%
#     mutate_at(vars(Metric,Rep,Sample,Group,), factor) %>%
#     as.data.frame()
# 
#   # ----- #
#   
#   #### . Observed VS Null Model . ####
#   
# 
#   # Compute summary statistics
#   Sum_Raw <- Raw_Trans %>%
#     group_by(OBSvsNM) %>%  # Group_by to have two groups in a tibble for faceting
#     get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups 
#     # Draw the summary statistics
#     ggsummarytable(x = "OBSvsNM",           # Sub scenario 
#                    y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
#                    digits = 4,                        # Number of digits 
#                    size = 4,                          # Size of the text     
#                    ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
#                      
#   # Create the density plot
#   Dens_Raw <- Raw_Trans %>%
#     # Plotting
#     ggplot(aes(x = Value, group = OBSvsNM, color = OBSvsNM, fill = OBSvsNM)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
#     geom_density(alpha = 0.3) +
#     theme(legend.position= "right",
#           axis.title.x=element_blank(),
#           legend.title =element_blank())
#   
#   # Put the two graphs together
#   Plot_Raw_Data_1 <- as.ggplot(ggarrange(
#     Dens_Raw,Sum_Raw, nrow = 2, heights = c(2,1))) +
#     labs(
#       title = paste0(Name," : Distribution of observed and NM values."),
#       subtitle = paste0("Density plot and statistics of the observed and NM values of ",Name," for each scenario.\nAll replicates are mixed.")
#     ) + arrange_theme()
#   
#   # ----- #
#   
#   # #### . Observed ~ Replicates . ####
#   # 
#   # # Compute summary statistics
#   # Sum_Raw <- Raw_Data %>%
#   #   as.data.frame() %>%
#   #   # group_by(rep) %>%  # Group_by to have two groups in a tibble for faceting
#   #   # mutate_at(vars(rep,type,Habitat,SampleType), factor) %>%  # Transform all the stuff into factor
#   #   get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups 
#   #   # Draw the summary statistics
#   #   ggsummarytable(x = "rep",           # Sub scenario 
#   #                  y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
#   #                  labeller = "label_both",
#   #                  digits = 4,                        # Number of digits 
#   #                  size = 4,                          # Size of the text     
#   #                  ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
#   # 
#   # # Create the density plot
#   # Dens_Raw <- Raw_Data %>%
#   #   # Preparation
#   #   as.data.frame() %>%
#   #   group_by(Rep) %>%
#   #   mutate_at(vars(Rep,Group), factor) %>%  # Transform all the stuff into factor
#   #   # Plotting
#   #   ggplot(aes(x = .data[[Name]], group = Rep, color = Rep, fill = Rep)) + 
#   #   geom_density(alpha = 0.3) +
#   #   theme(legend.position= "right",
#   #         legend.title =element_blank()) +
#   #   xlab(" ")
#   # 
#   # # Create a violin plot
#   # Violin_Raw <- Raw_Data %>%
#   #   # Preparation
#   #   as.data.frame() %>%
#   #   group_by(Rep) %>%
#   #   mutate_at(vars(Rep,Group), factor) %>%  # Transform all the stuff into factor
#   #   # Plotting
#   #   ggplot(aes(x = Rep, y = .data[[Name]], group = Rep, color = Rep, fill = Rep)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
#   #   geom_violin(alpha = 0.3) +
#   #   theme(legend.position= "right",
#   #         legend.title =element_blank()) +
#   #   xlab(" ")
#   # 
#   # # Put the two graphs together
#   # Plot_Raw_Data_2 <- as.ggplot(ggarrange(
#   #   ggarrange(Dens_Raw,Violin_Raw, common.legend = T, legend = "right"),
#   #   Sum_Raw,nrow = 2, heights = c(2,1))) +
#   #   labs(
#   #     title = paste0(Name," : Distribution of observed values."),
#   #     subtitle = paste0("Density plot and statistics of the observed values of ",Name," for each scenario.\nObserved values are splitted between replicates.")
#   #   ) + arrange_theme()
#   # 
#   # # ----- #
#   # 
#   #### . Observed ~ Group . ####
#   
#   # Compute summary statistics
#   Sum_Raw <- Raw_Data %>%
#     group_by(Group) %>%  # Group_by to have two groups in a tibble for faceting
#     get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups 
#     # Draw the summaru statistics
#     ggsummarytable(x = "Group",           # Sub scenario 
#                    y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
#                    labeller = "label_both",
#                    digits = 4,                        # Number of digits 
#                    size = 4,                          # Size of the text     
#                    ggtheme = arrange_theme()) %>% as.ggplot()
# 
#   # Create a dataset of mean values
#   DataMean <- Raw_Data %>%
#     group_by(Group) %>%
#     reframe(mean_x = mean(get(Name)))
#   
#   # Create the density plot
#   Dens_Raw <- Raw_Data %>%
#     # Plotting
#     ggplot(aes(x = .data[[Name]], group = Group, fill = Group)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
#     geom_density(alpha = 0.3) +
#     geom_vline(data  = DataMean, aes(xintercept = mean_x, colour = Group)) +
#     theme(legend.position="right") +
#     xlab(" ")
# 
#   # Put the two graphs together
#   Plot_Raw_Data_3 <- as.ggplot(ggarrange(
#   Dens_Raw,Sum_Raw, nrow = 2, heights = c(2,1))) +
#   labs(
#     title = paste0(Name," : Distribution of observed values."),
#     subtitle = paste0("Density plot and statistics of the observed values of ",Name," for each scenario.\nObserved values are splitted between groups. Drawn lines represent the mean values of each ditribution."),
#     caption = "Habitat 1: Continuous / Habitat 2: Fragmented.\n SampleType 1: Core / SampleType 2: Ecotone "
#   ) + arrange_theme()
#   
#     # ----- #
#   
#   #### . Patchworking . ####
#   
# # Finally, bind the three plots alltogether on different pages using marrange grob
# pl <- list(Plot_Raw_Data_1,
#            # Plot_Raw_Data_2,
#            Plot_Raw_Data_3)
# # Transform them into grobs
# pl <- lapply(pl, as.grob)
# # Create the pdf
# Plot_Raw_Data_Total <- marrangeGrob(pl,nrow=1,ncol = 1, top = NULL)
# 
# # Print a message
# cli_alert_success(paste0(Name))
# 
# # Return the plot
# return(Plot_Raw_Data_Total)
# 
# } # END OF RAW FUNCTION
# 
# #### . Apply Raw_Function . ####
# 
# Raw_Total <- lapply(Raw_Data,Raw_Fct)
# 
# cli_alert_success("Raw data correctly plotted !")

# ##### ------------ 1.B: SIGNIFICANT RAW DATA -------------------- #####
# 
# cat(rule(left = "Plotting of the significant data", line_col = "yellow", line = "-", col = "br_yellow"))
# 
# #### . Function creation . ####
# 
# # Create a function to create the plots for each individual metric.
# # Signif Raw_Data is the list of results of raw data for the metrics
# Signif_Raw_Fct  <- function(Signif_Data) {
#   
#   # Get the name
#   Name <- unique(Signif_Data$Metric)
#   
#   # Transform all the stuff into factor
#   Signif_Data <- Signif_Data %>%
#     mutate_at(vars(Metric,Rep,Group), factor) %>%
#     as.data.frame() %>%
#     # Remove the "P-value" columns that we don't need
#     select(-starts_with("P_")) %>%
#     # Pivot the dataframe to have the different tests as a variable
#     pivot_longer(
#       cols = starts_with("Test"),
#       names_to = "Test", values_to = "Result")
# 
#     # ----- # 
#   
#   #### . Total . ####
#   
#   # --- Compute summary statistics --- #
#   Summ_Signif <- Signif_Data %>%
#     group_by(Metric, Test) %>%
#     get_summary_stats(type = "common")
#   
#   # Are the tests significantly different ? 
#   S1 <- Summ_Signif %>%
#     ggsummarytable( x = "Test",
#                     y = c("n","min", "max", "mean","sd"),
#                     labeller = "label_both",
#                     digits = 4,
#                     size = 4, 
#                     ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
#       
#   # --- Create violin plots --- # 
#   
#   # Are the tests significantly different ? 
#    P1 <- Signif_Data %>%
#     group_by(Metric,Test) %>%  # Group_by to have two groups in a tibble for faceting
#     ggplot(aes(x=Test, y=Result, fill = Test)) +
#     scale_fill_manual(values = c("grey","forestgreen","red3")) +
#     geom_violin() +
#     xlab(" ") +
#     scale_y_continuous(name="% Significant tests", limits=c(0, 1))
#   
# 
#   # Put the two graphs together
#   Plot_Signif_Data_1 <- as.ggplot(ggarrange(
#     P1,S1,nrow = 2, heights = c(4,1))) +
#     labs(
#       title = paste0(Name," : Significancy tests."),
#       subtitle = paste0("Percentage of significant tests (Observed values different from the Null Model) of ",Name,".\nAll replicates and groups are mixed.")
#     ) + arrange_theme()
#   
#   
#     # ----- # 
#   
#   # #### . Splitted by Replicates . ####
#   # 
#   # # --- Compute summary statistics --- #
#   # Summ_Signif <- Signif_Data %>%
#   #   group_by(Metric,rep) %>%
#   #   mutate_at(vars(Metric,rep), factor) %>%  # Transform all the stuff into factor
#   #   get_summary_stats(type = "common")
#   # 
#   # # Are the tests significantly different ? 
#   # S2 <- Summ_Signif[Summ_Signif$variable == "Test_Different",] %>%
#   #   ggsummarytable( x = "rep",
#   #                 y = c("n","min", "max", "mean","sd"),
#   #                 labeller = "label_both",
#   #                 digits = 4,
#   #                 size = 4, 
#   #                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
#   # 
#   # # --- Create violin plots --- # 
#   # 
#   # # Are the tests significantly different ? 
#   # P2 <- Signif_Data %>%
#   #   group_by(Metric,rep) %>%  # Group_by to have two groups in a tibble for faceting
#   #   mutate_at(vars(rep,type,Habitat,SampleType), factor) %>%  # Transform all the stuff into factor
#   #   ggplot(aes(x=rep, y=Test_Different ,color = rep, fill = rep)) +
#   #   #scale_fill_manual(values = c("#CA5310","#FBBA72")) +
#   #   geom_violin() +
#   #   #facet_grid(Scenario ~ Habitat ,labeller = labeller(Scenario =label_value, Habitat = label_both))+
#   #   # stat_compare_means() + 
#   #   xlab(" ") +
#   #   theme(legend.position = "none") +
#   #   scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
#   #   theme(
#   #     plot.margin = margin(b = 0),
#   #     # The size of the axes labels are different for x and y.
#   #     axis.text.x = element_blank(),
#   #     # Also, the ticks have a very light grey color
#   #     axis.ticks = element_line(color = "grey91", linewidth = .5),
#   #     # The length of the axis ticks is increased.
#   #     axis.ticks.length.x = unit(.5, "lines"))
#   # 
#   # # Put the two graphs together
#   # Plot_Signif_Data_2 <- as.ggplot(ggarrange(
#   #   P2,S2,nrow = 2, heights = c(2,1))) +
#   #   labs(
#   #     subtitle = paste0("Percentage of significant tests (Observed values different from the Null Model) of ",Name," splitted by replicates.")
#   #   ) + arrange_theme()
#   # 
#   #   # ----- #
#   # 
#   #### . Splitted by Groups . ####
#   
#   # --- Compute summary statistics --- #
#   Summ_Signif <- Signif_Data %>%
#     group_by(Metric, Group, Test) %>%
#     get_summary_stats(type = "common")
# 
#   # Are the tests significantly different ? 
#   S3 <- Summ_Signif %>%
#     ggsummarytable(x = "Group",
#                   y = c("n","min", "max", "mean","sd"),
#                   facet.by = "Test",
#                   labeller = "label_both",
#                   digits = 4,
#                   size = 3, 
#                   ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot()
# 
# 
#   # --- Create violin plots --- # 
# 
#   # Are the tests significantly different ? 
#   P3 <- Signif_Data %>%
#     group_by(Metric, Group, Test) %>%
#     ggplot(aes(x=Group, y=Result, fill = Test)) +
#     scale_fill_manual(values = c("grey","forestgreen","red3")) +
#     geom_violin() +
#     facet_grid(. ~ Test )+
#     # stat_compare_means() + 
#     xlab(" ") +
#     scale_y_continuous(name="% Significant tests", limits=c(0, 1)) 
# 
#   # Put the two graphs together
#   Plot_Signif_Data_3 <- as.ggplot(ggarrange(
#     P3,S3,nrow = 2, heights = c(1,1),common.legend = T,legend = "right")) +
#     labs(
#       title = paste0(Name," : Significancy tests."),
#       subtitle = paste0("Percentage of significant tests (Observed values different from the Null Model) of ",Name," splitted by Group.")
#     ) + arrange_theme()
#   
#     # ----- #
#   
#   #### . Patchworking . ####
# 
# 
#   # Finally, bind the three plots alltogether on different pages using marrange grob
#   pl <- list(Plot_Signif_Data_1,
#              # Plot_Signif_Data_2,
#              Plot_Signif_Data_3)
#   # Transform them into grobs
#   pl <- lapply(pl, as.grob)
#   # Create the pdf
#   Plot_Signif_Data_Total <- marrangeGrob(pl,nrow=1,ncol = 1, top = NULL)
#   
#   # Print a message
#   Name <- gsub(pattern = "_.*", replacement = "", Name, perl = T) # Just change the name 
#   cli_alert_success(paste0(Name))
#   
#   # Retrun the plot
#   return(Plot_Signif_Data_Total)
# 
# } # END OF SIGNIF FUNCTION
#   
# #### . Apply Signif_Function . ####
# 
# Signif_Total <- lapply(Signif_Data,Signif_Raw_Fct)
# 
# cli_alert_success("Significant Raw data correctly plotted !")

##### ------------ 1.C: MEAN DATA -------------------- #####

cat(rule(left = "Plotting of the mean raw data ...", line_col = "yellow", line = "-", col = "br_yellow"))

#### . Function creation . ####

# Create a function to create the plots for each individual metric.
# Raw_Data is the list of results of raw data for the metrics
Mean_Raw_Fct  <- function(Mean_Raw_Data) {

  # Extract the Metric Name in function of the type of metric (Alpha or Beta) using the number of different groups
  if(length(unique(Mean_Raw_Data$Group)) == 4){  Name <- colnames(Mean_Raw_Data)[6]} else {Name <- colnames(Mean_Raw_Data)[8]}

  # Add an "if" because of a wrong name for PSV and PSE
  if(Name == "PSVs"){Name <- "PSV"}
  if(Name == "PSEs"){Name <- "PSE"}
  
  # Transform the data to add a column for the null vs obs data
  Mean_Raw_Trans <- Mean_Raw_Data %>%
    pivot_longer(
      cols = contains(Name, ignore.case = TRUE, vars = NULL),
      names_to = "OBSvsNM", values_to = "Value")

  # Change the values contained in the new column
  Mean_Raw_Trans$OBSvsNM[grep("NM",Mean_Raw_Trans$OBSvsNM)] <- paste0(Name,"_NM")

  # Transform all the stuff into factor
  Mean_Raw_Trans <-  Mean_Raw_Trans %>%
    mutate_at(vars(Metric,Rep,Group), factor) %>%
    as.data.frame()

  Mean_Raw_Data <- Mean_Raw_Data %>%
    mutate_at(vars(Metric,Rep,Group), factor) %>%
    as.data.frame()

  # ----- #

  #### . Observed VS Null Model . ####

  # n = Number of replicates X Number of groups (4 or 10) X (Number of NullModel (99)).

  # Compute summary statistics
  Mean_Sum_Raw <- Mean_Raw_Trans %>%
    group_by(OBSvsNM) %>%  # Group_by to have two groups in a tibble for faceting
    get_summary_stats(Value,type = "common") %>% # Compute the summary statistics for each groups
    # Draw the summary statistics
    ggsummarytable(x = "OBSvsNM",           # Sub scenario
                   y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                   digits = 4,                        # Number of digits
                   size = 4,                          # Size of the text
                   ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot()

  # Create the density plot
  Mean_Dens_Raw <- Mean_Raw_Trans %>%
    # Plotting
    ggplot(aes(x = Value, group = OBSvsNM, color = OBSvsNM, fill = OBSvsNM)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
    geom_density(alpha = 0.3) +
    theme(legend.position= "right",
          axis.title.x=element_blank(),
          legend.title =element_blank())

  # Put the two graphs together
  Plot_Mean_Raw_Data_1 <- as.ggplot(ggarrange(
    Mean_Dens_Raw,Mean_Sum_Raw, nrow = 2, heights = c(2,1))) +
    labs(
      title = paste0(Name," : Distribution of mean observed and mean NM values."),
      subtitle = paste0("Density plot and statistics of the mean observed and mean NM values of ",Name," for each group.\nAll replicates are mixed."),
      caption = paste0("n = Number of replicates (",Nrep,") X Number of groups (",length(levels(Mean_Raw_Data$Group)),") X (Number of NullModel (",Ntree,") ).")
    ) + arrange_theme()

  # ----- #

  #### . Observed ~ Group . ####
  
  # Add an "if" because of a wrong name for PSV and PSE
  if(Name == "PSV"){Name <- "PSVs"}
  if(Name == "PSE"){Name <- "PSEs"}

  # Compute summary statistics
  Mean_Sum_Raw <- Mean_Raw_Data %>%
    group_by(Group) %>%  # Group_by to have two groups in a tibble for faceting
    get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups
    # Draw the summaru statistics
    ggsummarytable(x = "Group",           # Sub scenario
                   y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                   labeller = "label_both",
                   digits = 4,                        # Number of digits
                   size = 4,                          # Size of the text
                   ggtheme = arrange_theme()) %>% as.ggplot()

  # Create a dataset of mean values
  DataMean <- Mean_Raw_Data %>%
    group_by(Group) %>%
    reframe(mean_x = mean(get(Name)))

  # Create the density plot
  Mean_Dens_Raw <- Mean_Raw_Data %>%
    # Plotting
    ggplot(aes(x = .data[[Name]], group = Group, fill = Group)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
    # facet_grid(. ~Group) +
    geom_density(alpha = 0.3) +
    geom_vline(data  = DataMean, aes(xintercept = mean_x, colour = Group)) +
    theme(legend.position="right") +
    xlab(" ")

  # Put the two graphs together
  Plot_Mean_Raw_Data_3 <- as.ggplot(ggarrange(
    Mean_Dens_Raw,Mean_Sum_Raw, nrow = 2, heights = c(2,1))) +
    labs(
      title = paste0(Name," : Distribution of mean observed values."),
      subtitle = paste0("Density plot and statistics of mean the observed values of ",Name," for each group\nObserved values are splitted between groups. Drawn lines represent the mean values of each ditribution."),
      caption = "Habitat 1: Continuous / Habitat 2: Fragmented.\n SampleType 1: Core / SampleType 2: Border "
    ) + arrange_theme()

  # ----- #

  #### . Patchworking . ####

  # Finally, bind the three plots alltogether on different pages using marrange grob
  pl <- list(Plot_Mean_Raw_Data_1,
             # Plot_Raw_Data_2,
             Plot_Mean_Raw_Data_3)
  # Transform them into grobs
  pl <- lapply(pl, as.grob)
  # Create the pdf
  Plot_Mean_Raw_Data_Total <- marrangeGrob(pl,nrow=1,ncol = 1, top = NULL)

  # Print a message
  cli_alert_success(paste0(Name))

  # Return the plot
  return(Plot_Mean_Raw_Data_Total)

} # END OF RAW FUNCTION

#### . Apply Mean_Raw_Function . ####

Mean_Raw_Total <- lapply(Mean_Raw_Data,Mean_Raw_Fct)

cli_alert_success("Mean Raw data correctly plotted !")

##### ------------ 1.D: SIGNIFICANT MEAN DATA -------------------- #####

cat(rule(left = "Plotting of the significant mean data", line_col = "yellow", line = "-", col = "br_yellow"))

#### . Function creation . ####

# Create a function to create the plots for each individual metric.
# Signif Raw_Data is the list of results of raw data for the metrics
Signif_Mean_Fct  <- function(Mean_Signif_Data) {
  
  # Get the name
  Name <- unique(Mean_Signif_Data$Metric)
  
  # Transform all the stuff into factor
  Mean_Signif_Data <- Mean_Signif_Data %>%
    mutate_at(vars(Metric,Rep,Group), factor) %>%
    as.data.frame() %>%
    # Remove the "P-value" columns that we don't need and an eventual column samples
    select(-starts_with("P_")) %>%
    select(-any_of(c("Sample","Sample_A","Sample_B","Habitat","Habitat_A","Habitat_B"))) %>%
    # Pivot the dataframe to have the different tests as a variable
    pivot_longer(
      cols = starts_with("Test"),
      names_to = "Test", values_to = "Result")
  
  # ----- # 
  
  #### . Total . ####
  
  # --- Compute summary statistics --- #
  Mean_Summ_Signif <- Mean_Signif_Data %>%
    group_by(Metric, Test) %>%
    get_summary_stats(type = "common")
  
  # Are the tests significantly different ? 
  S1 <- Mean_Summ_Signif %>%
    ggsummarytable( x = "Test",
                    y = c("n","min", "max", "mean","sd"),
                    labeller = "label_both",
                    digits = 4,
                    size = 4, 
                    ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
  
  # --- Create violin plots --- # 
  
  # Are the tests significantly different ? 
  P1 <- Mean_Summ_Signif %>%
    group_by(Metric,Test) %>%  # Group_by to have two groups in a tibble for faceting
    ggplot(aes(x=Test, y=mean, color = Test)) +
    scale_color_manual(values = c("grey","forestgreen","red3")) +
    geom_point(size = 5) +
    geom_segment(aes(x=Test, xend=Test, y=0, yend=mean)) +
    xlab(" ") +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1))
  
  # Put the two graphs together
  Plot_Mean_Signif_Data_1 <- as.ggplot(ggarrange(
    P1,S1,nrow = 2, heights = c(4,1))) +
    labs(
      title = paste0(Name," : Significancy tests."),
      subtitle = paste0("Percentage of significant tests (Mean observed values different from the mean Null Models) of ",Name,".\nAll replicates and groups are mixed.")
    ) + arrange_theme()
  
  
  # ----- # 
  
  # #### . Splitted by Replicates . ####
  # 
  # # --- Compute summary statistics --- #
  # Summ_Signif <- Signif_Data %>%
  #   group_by(Metric,rep) %>%
  #   mutate_at(vars(Metric,rep), factor) %>%  # Transform all the stuff into factor
  #   get_summary_stats(type = "common")
  # 
  # # Are the tests significantly different ? 
  # S2 <- Summ_Signif[Summ_Signif$variable == "Test_Different",] %>%
  #   ggsummarytable( x = "rep",
  #                 y = c("n","min", "max", "mean","sd"),
  #                 labeller = "label_both",
  #                 digits = 4,
  #                 size = 4, 
  #                 ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot() 
  # 
  # # --- Create violin plots --- # 
  # 
  # # Are the tests significantly different ? 
  # P2 <- Signif_Data %>%
  #   group_by(Metric,rep) %>%  # Group_by to have two groups in a tibble for faceting
  #   mutate_at(vars(rep,type,Habitat,SampleType), factor) %>%  # Transform all the stuff into factor
  #   ggplot(aes(x=rep, y=Test_Different ,color = rep, fill = rep)) +
  #   #scale_fill_manual(values = c("#CA5310","#FBBA72")) +
  #   geom_violin() +
  #   #facet_grid(Scenario ~ Habitat ,labeller = labeller(Scenario =label_value, Habitat = label_both))+
  #   # stat_compare_means() + 
  #   xlab(" ") +
  #   theme(legend.position = "none") +
  #   scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
  #   theme(
  #     plot.margin = margin(b = 0),
  #     # The size of the axes labels are different for x and y.
  #     axis.text.x = element_blank(),
  #     # Also, the ticks have a very light grey color
  #     axis.ticks = element_line(color = "grey91", linewidth = .5),
  #     # The length of the axis ticks is increased.
  #     axis.ticks.length.x = unit(.5, "lines"))
  # 
  # # Put the two graphs together
  # Plot_Signif_Data_2 <- as.ggplot(ggarrange(
  #   P2,S2,nrow = 2, heights = c(2,1))) +
  #   labs(
  #     subtitle = paste0("Percentage of significant tests (Observed values different from the Null Model) of ",Name," splitted by replicates.")
  #   ) + arrange_theme()
  # 
  #   # ----- #
  # 
  #### . Splitted by Groups . ####
  
  # --- Compute summary statistics --- #
  Mean_Summ_Signif <- Mean_Signif_Data %>%
    group_by(Metric, Group, Test) %>%
    get_summary_stats(type = "common")
  
  # Are the tests significantly different ? 
  Mean_S3 <- Mean_Summ_Signif %>%
    ggsummarytable(x = "Group",
                   y = c("n","min", "max", "mean","sd"),
                   facet.by = "Test",
                   labeller = "label_both",
                   digits = 4,
                   size = 3, 
                   ggtheme = arrange_theme() + theme(axis.title.x=element_blank())) %>% as.ggplot()
  
  
  # --- Create lollipop plots --- # 
  
  # Are the tests significantly different ? 
  Mean_P3 <- Mean_Summ_Signif %>%
    group_by(Metric, Group, Test) %>%
    ggplot(aes(x=Group, y=mean, color= Test)) +
    scale_color_manual(values = c("grey","forestgreen","red3")) +
    geom_point() +
    geom_segment(aes(x=Group, xend=Group, y=0, yend=mean)) +
    facet_grid(. ~ Test )+
    # stat_compare_means() + 
    xlab(" ") +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) 
  
  # Put the two graphs together
  Plot_Mean_Signif_Data_3 <- as.ggplot(ggarrange(
    Mean_P3,Mean_S3,nrow = 2, heights = c(1,1),common.legend = T,legend = "right")) +
    labs(
      title = paste0(Name," : Significancy tests."),
      subtitle = paste0("Percentage of significant tests (Mean observed values different from the mean Null Model) of ",Name," splitted by Group.")
    ) + arrange_theme()
  
  # ----- #
  
  #### . Patchworking . ####
  
  
  # Finally, bind the three plots alltogether on different pages using marrange grob
  pl <- list(Plot_Mean_Signif_Data_1,
             # Plot_Signif_Data_2,
             Plot_Mean_Signif_Data_3)
  # Transform them into grobs
  pl <- lapply(pl, as.grob)
  # Create the pdf
  Plot_Signif_Mean_Data_Total <- marrangeGrob(pl,nrow=1,ncol = 1, top = NULL)
  
  # Print a message
  Name <- gsub(pattern = "_.*", replacement = "", Name, perl = T) # Just change the name 
  cli_alert_success(paste0(Name))
  
  # Retrun the plot
  return(Plot_Signif_Mean_Data_Total)
  
} # END OF SIGNIF FUNCTION

#### . Apply Signif_Function . ####

Signif_Mean_Total <- lapply(Mean_Signif_Data,Signif_Mean_Fct)

cli_alert_success("Significant Mean data correctly plotted !")
# ---------------------------- #
##### STEP 2: PATCHWORKING #####
# ---------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

cat(rule(left = "Saving the metric plot report ...", line_col = "yellow", line = "-", col = "br_yellow"))

Merged_List <- mapply(list, 
                      #Raw_Total, 
                      #Signif_Total, 
                      Mean_Raw_Total, 
                      Signif_Mean_Total, 
                      SIMPLIFY=FALSE) %>% list_flatten() %>% list_flatten()

# Transform them into grobs
pl <- lapply(Merged_List, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)

# Save the pdf
ggsave(filename = paste0("Scenario",Scenario,"_MetricReport.pdf"),
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario),
       width = 35,
       height = 25,
       units = "cm")

cli_alert_success("Report correctly saved !")

# ----------------------------------------------------------------------- #

cat(rule(left = "SCRIPT COMET_METRIC_CHECK.R ENDING", line_col = "red", line = "-", col = "br_red"))

# ----------------------------------------------------------------------- #
