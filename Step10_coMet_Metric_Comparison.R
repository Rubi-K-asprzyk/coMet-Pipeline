#!/usr/bin/env Rscript 

##### ______________________________________________________________________________________________________________________________________ #####

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");require(cli)})
# Display a beginning message. 
cat(rule(left = "SCRIPT COMET_FINAL_PLOTS.R BEGINNING", line_col = "red", line = "-", col = "br_red")) 

#-----------------------------------------------------------------------#
##### coMet_Metric_Plots: Drawing of significancy plots for analysis  ###
#-----------------------------------------------------------------------#

# This script is used in the coMet pipeline to draw the final plots of the results for each metric and to compare all the metrics all_together ! 
# It contains density plots of the values of the plots as well as plots of the metric significance.  

# This script will all the results for all the input scenarios to allow the comparison between all the scenarios

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
       "stringr",
       "gridExtra",
       "ggplotify",
       "patchwork",
       "ggforce",
       "argparser",
       "purrr",
       "ggridges",
       "ggsignif",
       "tidyr",
       "rstatix"
       )
# Avoid summarise() messages 
options("dplyr.summarise.inform" = FALSE)
# Set lifecycle_verbosity to 'warning' or 'quiet', ensuring deprecated messages are shown or not.
options(lifecycle_verbosity = "quiet")

#### . Local Mode . ####

Parameters <- list()
Parameters$Config <- c("coMet_ConfigFiles/Equalizing100_Sig2.R","coMet_ConfigFiles/Equalizing100_Sig5.R","coMet_ConfigFiles/Stabilizing100_Sig5.R","coMet_ConfigFiles/Stochastic100_Sig5.R",
                                            "coMet_ConfigFiles/Equalizing100_Sig2_SF1.R","coMet_ConfigFiles/Equalizing100_Sig2_SF2.R","coMet_ConfigFiles/Equa50Stab50_Sig2.R","coMet_ConfigFiles/Equa50Stab50_Sig5.R")
Parameters$output <- "Metric_Total.pdf"

#### . Argument Parser . ####

# Create the parser
arg_parser <- arg_parser("Check the communities computed by coMet_ComSim and create an output pdf.", name = "coMet_ComSim_Check.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = Inf, help = "Configuration file for coMet_ComSim.R")

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
  plot.title.position = "plot",
  plot.caption.position = "plot",
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

# ------------------------------------------- #
##### STEP 1: LOADING THE COMMUNITY FILES #####
# ------------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "LOADING OF THE COMMUNITY FILES", line_col = "green", line = "-", col = "br_green"))

# Read and combine all the metrics result into two summary files. One for the metric result, the other for the significant result

    # ----- #

  #### . Metrics Available . ####

# Loop on all scenarios to retrieve the metrics available
Metric.List <- foreach(Conf = 1:length(Config),.combine = c) %do% {
  
  # Get the Scenario Name
  source(Config[Conf])

  # Retrieve the number of metrics already computed for each Scenario
    # We use the files of "Signif" because they are the last step of the processing.
  Metrics <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Metrics/"))

  # Extract the name of the metric from the name of the file
    # Remove the name of the scenario
  Metrics <- gsub(pattern = paste0("Scenario",Scenario,"_"), replacement = "", Metrics)
    # Remove the "signif" parts of the file
  Metrics <- gsub(pattern = "_Signif*", replacement = "", Metrics)
  # Remove the "Mean" parts of the file
  Metrics <- gsub(pattern = "_Mean*", replacement = "", Metrics)
  # Remove the "csv" parts of the file
  Metrics <- gsub(pattern = ".csv", replacement = "", Metrics)

}

  #### . /!\ Remove metrics . #####

# Remove the DISC metric because it causes troubles
Metric.List <- Metric.List[!Metric.List %in% c("DISC_Rao","DISC_Rao_Mean")]

# Find which metrics are present in all scenarios.
  # Count the presence of the metrics
  Metric.Count <- plyr::count(Metric.List)

  # Search for the metrics present in all scenarios
  Metric.Count <- Metric.Count$x[Metric.Count$freq == length(Config)*4] # We multiply by 4 because we want Raw, Raw_Signif, Mean and Mean_Signif datas

  # Create a version without the author name
  Metric.Count.Lite <- gsub("\\_[^_]*$","",Metric.Count)
  
  # Replace F and J values
  Metric.Count.Lite <- Metric.Count.Lite %>%
    replace(Metric.Count.Lite == "F","F_") %>%
    replace(Metric.Count.Lite == "J","J_") %>%
    replace(Metric.Count.Lite == "PDbeta","PDb") %>%
    replace(Metric.Count.Lite == "PSE","PSEs") %>%
    replace(Metric.Count.Lite == "PSV","PSVs") 
  
  # Test if no common metrics were found
  if(is.null(Metric.Count) == T){stop("No metric results common to all the scenario were found, please verify that all the scenarios where correctly run.")
  } else {cat(col_green(paste0(length(Metric.Count)," Metrics were found in common in all the wanted Scenarios.")))}

# Metric.Count is therefore the list of metrics that were computed for all the wanted scenarios.

  # Retrieve the Scenario Names
Scenario.List <- foreach(Conf = 1:length(Config),.combine = c) %do% {
  # Get the Scenario Name
  source(Config[Conf])
  # Display the Scenario
  Scenario
}

    # ----- #

  #### . Prepare total dataset . ####

# /!\ WARNING WARNING /!\ 

# Scenario order in 'File_Mean_Raw" is not the same as in "Scenario_List" therefore leading to misleading results. "Scenario_List" will therefore 
# Be reordered alphabetically 

Scenario.List <- sort(Scenario.List)

# Equalizing scenarios are still not well reorganised, moving them. 
Scenario.List[c(3,5)] <- Scenario.List[c(5,3)]
Scenario.List[c(3,4)] <- Scenario.List[c(4,3)]

# Loop across the metrics available for each scenario.
Metric_Total <- foreach(i = 1:length(Metric.Count)) %do% {

  # Find all the result file for the metric
  # File_Raw <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario.List,"/Metrics"), pattern = paste0("_",Metric.Count[i],".csv"),recursive = T, full.names = T)
  # File_Signif <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario.List,"/Metrics"), pattern = paste0("_",Metric.Count[i],"_Signif\\.csv$"),recursive = T, full.names = T)
  File_Mean_Raw <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario.List,"/Metrics"), pattern = paste0("_",Metric.Count[i],"_Mean.csv"),recursive = T, full.names = T)
  File_Mean_Signif <- list.files(path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario.List,"/Metrics"), pattern = paste0("_",Metric.Count[i],"_Mean_Signif.csv"),recursive = T, full.names = T)
  
  # Loop across the found files of raw data
  # suppressMessages(Raw_Data <- foreach(j = 1:length(File_Raw), .combine = full_join) %do% {
  #   # Load the file
  #   Raw <- read.csv(File_Raw[j],row.names = 1) %>% # Load the file
  #     mutate("Scenario" = Scenario.List[j], .before = 1) # Append the scenario
  # })
  # Loop across the found files of Significant data
  # suppressMessages(Signif_Data <- foreach(j = 1:length(File_Signif), .combine = full_join) %do% {
  # 
  #   # Load the file
  #   Signif <- read.csv(File_Signif[j],row.names = 1) %>% # Load the file
  #     mutate("Scenario" = Scenario.List[j], .before = 1)
  # })
  # Loop across the found files of raw data
  
  suppressMessages(Mean_Raw_Data <- foreach(j = 1:length(File_Mean_Raw), .combine = full_join) %do% {
    # Load the file
    Raw <- read.csv(File_Mean_Raw[j],row.names = 1) %>% # Load the file
      mutate("Scenario" = Scenario.List[j], .before = 1) # Append the scenario
  })
  # Loop across the found files of raw data
  suppressMessages(Mean_Signif_Data <- foreach(j = 1:length(File_Mean_Signif), .combine = full_join) %do% {
    # Load the file
    Raw <- read.csv(File_Mean_Signif[j],row.names = 1) %>% # Load the file
      mutate("Scenario" = Scenario.List[j], .before = 1) # Append the scenario
  })
  
  # Combine the data
  Total <- list(Mean_Raw_Data,Mean_Signif_Data) %>%
    setNames(nm = c("Mean_Raw","Mean_Signif")) # Set the names of the replicates.

  # Return
  return(Total)
} %>% setNames(nm = Metric.Count)

# Print a message
cli_alert_success("coMet_ConfigFile(s).R correctly read !")

##### ______________________________________________________________________________________________________________________________________ #####

# ---------------------------------------- #
##### STEP 2: INDIVIDUAL METRIC REPORT #####
# ---------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "PLOTTING", line_col = "green", line = "-", col = "br_green"))

    # ----- #

  #### . Function Creation . ####

# Create a function to make all the plots needed for each metrics
# The argument "List' is a list of the raw and risgnif data for one metric. One object of the list Metric_Total.
# "Name" comes with the List object number and is the name of the metric to compute

ploturbo <- function(List) {

  # Extract the Metric Name
  Name <- colnames(List$Mean_Raw)[6]
  # If the metric is a beta metric, we need to change the name 
  if (Name == "Sample_B") {
    Name <- colnames(List$Mean_Raw)[8]
    }
  
  # ----- #
  
  #### _RAW_ ####
  
  #### . Observed ~ Scenario + Group . ####

  # --- Compute summary statistics --- #
  # Summ_Raw <- List$Raw %>%
  #       group_by(Scenario,Group) %>%  # Group_by to have two groups in a tibble for faceting
  #       get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups
  #       # Draw the summaru statistics
  #       ggsummarytable(x = "Group",           # Sub scenario
  #                  facet.by = "Scenario",
  #                  y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
  #                  digits = 4,                        # Number of digits
  #                  size = 4,                          # Size of the text
  #                  ggtheme = theme_gray()) %>%        # Theme of the plot
  #       ggpar(font.ytickslab = c(13, "plain", "black")) + # Aesthetics of the tick labs
  #       arrange_theme()
  # 
  # # --- Create density plot --- #
  # 
  # # Create a dataset of mean values
  # DataMean <- List$Raw %>%
  #   group_by(Scenario,Group) %>%
  #   mutate_at("Group", factor) %>%
  #   reframe(mean_x = mean(get(Name)))
  # 
  # Dens_Raw <- List$Raw %>%
  #   group_by(Rep,Group) %>%
  #   mutate_at(vars(Rep,Group), factor) %>%  # Transform all the stuff into factor
  #   ggplot(aes(x = .data[[Name]], group = Group, fill = Group)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
  #   geom_density(alpha = 0.3) +
  #   facet_grid(. ~ Scenario)+
  #   # scale_fill_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
  #   geom_vline(data  = DataMean, aes(xintercept = mean_x, colour = Group)) +
  #   # scale_color_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
  #   theme(legend.position="right") +
  #   labs(caption = paste0("n = Number of Samples (",length(unique(List$Raw$Sample)),") * Number of Replicates (",max(as.numeric(List$Raw$Rep)),").")) +
  #   xlab(" ")
  # 
  # 
  # # ----- #
  # 
  # #### . Observed ~ Scenario . ####
  # 
  # # --- Compute summary statistics --- #
  # Summ_Raw_2 <- List$Raw %>%
  #   group_by(Scenario) %>%  # Group_by to have two groups in a tibble for faceting
  #   get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups
  #   # Draw the summaru statistics
  #   ggsummarytable(x = "Scenario",           # Sub scenario
  #                  y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
  #                  labeller = "label_both",
  #                  digits = 4,                        # Number of digits
  #                  size = 4,                          # Size of the text
  #                  ggtheme = theme_gray()) %>%        # Theme of the plot
  #   ggpar(font.ytickslab = c(13, "plain", "black")) + # Aesthetics of the tick labs
  #   xlab(" ") + 
  #   labs(caption = paste0("n = Number of Samples (",length(unique(List$Raw$Sample)),") * Number of Replicates (",max(as.numeric(List$Raw$Rep)),").")) +
  #   arrange_theme()
  # 
  # # Create a dataset of mean values
  # DataMean_2 <- List$Raw %>%
  #   group_by(Scenario) %>%
  #   reframe(mean_x = mean(get(Name)))
  # 
  # Dens_Raw_2 <- List$Raw %>%
  #   group_by(Scenario) %>%
  #   ggplot(aes(x = .data[[Name]], y = Scenario, fill = Scenario)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
  #   geom_density_ridges(alpha = 0.5,quantile_lines = TRUE, quantiles = 2,
  #                       jittered_points = FALSE,
  #                       position = position_points_jitter(width = 0.05, height = 0),
  #                       point_shape = '/', point_size = 3, point_alpha = 1) +
  #   scale_fill_viridis_d() + 
  #   geom_signif(comparisons = combn(Scenario.List,2, simplify = F),
  #               map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
  #               margin_top = 0.2,
  #               step_increase = 0.1,
  #               tip_length = 0.01) +
  #   theme(legend.position="right") +
  #   labs(caption ="P_value: *** = 0.001 / ** = 0.01 / * = 0.05 ") +
  #   xlab(" ") + arrange_theme()
  # 
  # # __ Combine the plots __ # 
  # suppressMessages(P1 <- as.ggplot(ggarrange(
  #   Dens_Raw_2,Summ_Raw_2, nrow = 2, heights = c(2,1))) +
  #   labs(
  #     title = paste0(Name," : Distribution of observed values."),
  #     subtitle = paste0(" Density plot and statistics of the observed values of ",Name," for each scenario. Drawn lines represent the mean values of each ditribution.")
  #   ) + arrange_theme())
  # 
  # suppressMessages(P2 <- as.ggplot(ggarrange(
  #   Dens_Raw,Summ_Raw, nrow = 2, heights = c(2,1))) +
  #     labs(
  #       title = paste0(Name," : Distribution of observed values."),
  #       subtitle = paste0("Density plot and statistics of the observed values of ",Name," for each scenario splitted by group, all replicates being merged together.\nDrawn lines represent the mean values of each ditribution.")
  #     ) + arrange_theme())
  # 
  # # ----- #
  # 
  # #### . Significant Results ~ Scenario + Group. ####
  # 
  # # --- Compute summary statistics --- #
  # Summ_Signif <- List$Signif %>%
  #   group_by(Scenario,Group) %>%
  #   mutate_at(vars(Scenario,Group), factor) %>%  # Transform all the stuff into factor
  #   get_summary_stats(!Rep,type = "common") %>%
  #   # Remove the p-values
  #   filter(!str_detect(variable, 'P_value'))
  # 
  # # Are the tests significantly different ?
  # P4 <- Summ_Signif %>%
  #   ggsummarytable( x = "Group",
  #                   y = c("n","min", "max", "mean","sd"),
  #                   facet.by = c("Scenario","variable"),
  #                   digits = 4,
  #                   size = 4,
  #                   ggtheme = theme_gray()) %>%
  #   ggpar(font.ytickslab = c(13, "plain", "black")) +
  #   labs(
  #     title = paste0(Name," : Results of significancy tests."),
  #     subtitle = paste0("Boxplots and statistics of the significancy tests of ",Name," for each scenario splitted by group, all replicates being merged together.\nFor each replicates, the mean number of significant test is computed and plotted."),
  #     caption = paste0("For one replicate, the significance means the percentage of Observed values that are significantly different from the null model for each sample.\nn = Number of Replicates (",max(as.numeric(List$Raw$Rep)),")."),
  #   ) +
  #   arrange_theme()
  # 
  # # --- Create Box plots --- #
  # 
  # # Are the tests significantly different ?
  # P3 <- List$Signif %>%
  #       # Pivot the data of test results
  #       pivot_longer( cols = starts_with("Test"),
  #                 names_to = "variable",
  #                 values_to = "Results") %>%
  #       group_by(Scenario,Group) %>%  # Group_by to have two groups in a tibble for faceting
  #       ggplot(aes(x=Group, y=Results, group = Group, fill = Group)) +
  #       geom_boxplot() +
  #       facet_grid(Scenario ~ variable)+
  #       xlab(" ") +
  #       ylim(c(0,1)) +
  #   labs(
  #     title = paste0(Name," : Boxplots of significancy tests."),
  #     subtitle = paste0("Boxplots of the significancy tests of ",Name," for each scenario splitted by group, all replicates being merged together.\nFor each replicates, the mean number of significant test is computed and plotted."),
  #     caption = paste0("For one replicate, the significance means the percentage of Observed values that are significantly different from the null model for each sample.\nn = Number of Replicates (",max(as.numeric(List$Raw$Rep)),")."),
  #   ) + arrange_theme()
  #       
  # # ----- #
  
  #### _MEAN_ ####
  
  #### . Observed ~ Scenario + Group . ####
  
  # --- Compute summary statistics --- #
  Summ_Raw_Mean <- List$Mean_Raw %>%
    group_by(Scenario,Group) %>%  # Group_by to have two groups in a tibble for faceting
    get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups
    # Draw the summaru statistics
    ggsummarytable(x = "Group",           # Sub scenario
                   facet.by = "Scenario",
                   y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                   digits = 4,                        # Number of digits
                   size = 2,                          # Size of the text
                   ggtheme = theme_gray()) %>%        # Theme of the plot
    ggpar(font.ytickslab = c(8, "plain", "black")) + # Aesthetics of the tick labs
    arrange_theme()
  
  # --- Create density plot --- #
  
  # Create a dataset of mean values
  DataMean_Mean <- List$Mean_Raw %>%
    group_by(Scenario,Group) %>%
    mutate_at("Group", factor) %>%
    reframe(mean_x = mean(get(Name)))
  
  Dens_Raw_Mean <- List$Mean_Raw %>%
    group_by(Rep,Group) %>%
    mutate_at(vars(Rep,Group), factor) %>%  # Transform all the stuff into factor
    ggplot(aes(x = .data[[Name]], group = Group, fill = Group)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
    geom_density(alpha = 0.3) +
    facet_grid(. ~ Scenario)+
    # scale_fill_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
      # Draw a line at the mean value of each metric
    # geom_vline(data  = DataMean_Mean, aes(xintercept = mean_x, colour = Group)) +
    # scale_color_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
    theme(legend.position="top") +
    theme(strip.text.x.top = element_text(size = 8)) +
    xlab(" ") + arrange_theme()
  
  # ----- #
  
  #### . Observed ~ Scenario . ####
  
  # --- Compute summary statistics --- #
  Summ_Raw_Mean_2 <- List$Mean_Raw %>%
    group_by(Scenario) %>%  # Group_by to have two groups in a tibble for faceting
    get_summary_stats(all_of(Name),type = "common") %>% # Compute the summary statistics for each groups
    # Draw the summaru statistics
    ggsummarytable(x = "Scenario",           # Sub scenario
                   y = c("n","min", "max", "mean","sd"),  # Choose the metrics we want to display
                   labeller = "label_both",
                   digits = 4,                        # Number of digits
                   size = 4,                          # Size of the text
                   ggtheme = theme_gray()) %>%        # Theme of the plot
    ggpar(font.ytickslab = c(8, "plain", "black")) + # Aesthetics of the tick labs
    xlab(" ") + 
    theme(axis.text.x = element_text(size = 8)) +
    arrange_theme()
  
  # Create a dataset of mean values
  DataMean_Mean2 <- List$Mean_Raw %>%
    group_by(Scenario) %>%
    reframe(mean_x = mean(get(Name)))
  
  Dens_Raw_Mean2 <- List$Mean_Raw %>%
    group_by(Scenario) %>%
    ggplot(aes(x = .data[[Name]], y = Scenario, fill = Scenario)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
    geom_density_ridges(alpha = 0.5,quantile_lines = TRUE, quantiles = 2,
                        jittered_points = FALSE,
                        position = position_points_jitter(width = 0.05, height = 0),
                        point_shape = '/', point_size = 3, point_alpha = 1) +
    scale_fill_viridis_d() + 
    # geom_signif(comparisons = combn(Scenario.List,2, simplify = F),
    #             map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    #             margin_top = 0.2,
    #             step_increase = 0.1,
    #             tip_length = 0.01) +
    theme(legend.position="right") +
    labs(caption ="P_value: *** = 0.001 / ** = 0.01 / * = 0.05 ") +
    xlab(" ") + arrange_theme()

  # Combine the two plots
  
  suppressMessages(P5 <- as.ggplot(ggarrange(
    Dens_Raw_Mean2,Summ_Raw_Mean_2, nrow = 2, heights = c(2,1))) +
    labs(
      title = paste0(Name," : Mean observed values."),
      subtitle = paste0("Density plot and statistics of the mean observed values of ",Name," for each scenario, all replicates being merged together.\nDrawn lines represent the mean values of each ditribution."),
      caption =  paste0("Means are computed per group. n = Number of groups (",length(table(List$Mean_Signif$Group)),") * Number of Replicates (100).")) +
    arrange_theme())
  
  suppressMessages(P6 <- as.ggplot(ggarrange(
    Dens_Raw_Mean,Summ_Raw_Mean, nrow = 2, heights = c(1,1))) +
    labs(
      title = paste0(Name," : Mean observed values ~ Groups."),
      subtitle = paste0("Density plot and statistics of the mean observed values per groups of ",Name," for each scenario, all replicates being merged together.\nDrawn lines represent the mean values of each ditribution."),
      caption =  paste0("Means are computed per group. n = Number of Replicates (100).")) +
    arrange_theme())
  
  # ----- #
  
  #### . Significant Results ~ Scenario + Group. ####
  
  # --- Compute summary statistics --- #
  Summ_Signif_Mean <- List$Mean_Signif %>%
    group_by(Scenario,Group) %>%
    mutate_at(vars(Scenario,Group), factor) %>%  # Transform all the stuff into factor
    get_summary_stats(!Rep,type = "common") %>%
    # Remove the p-values
    filter(str_detect(variable, 'Test'))
  
  # Are the tests significantly different ?
  P8 <- Summ_Signif_Mean %>%
    ggsummarytable( x = "Group",
                    y = c("n", "mean","sd"),
                    facet.by = c("Scenario","variable"),
                    digits = 3,
                    size = 3,
                    ggtheme = theme_gray()) %>%
    ggpar(font.ytickslab = c(13, "plain", "black")) +
    labs(
      title = paste0(Name," : Results of significancy tests."),
      subtitle = paste0("Statistics of the significancy tests of ",Name," for each scenario splitted by group, all replicates being merged together.\nFor each replicates, the mean number of significant test is computed and plotted."),
      caption = paste0("The mean represent the percentage of mean observed values that are significantly different from the null model for each sample.\nn = Number of Replicates (100)."),
    ) +
    theme(strip.text.y.right = element_text(size = 5)) +
    arrange_theme()
  
  # --- Create Box plots --- #
  
  # Are the tests significantly different ? 
  P7 <- Summ_Signif_Mean %>%
    group_by(Scenario, Group, variable) %>%
    ggplot(aes(x=Group, y=mean, color= variable)) +
    scale_color_manual(values = c("grey","forestgreen","red3")) +
    geom_point() +
    geom_segment(aes(x=Group, xend=Group, y=0, yend=mean)) +
    facet_grid(variable ~ Scenario )+
    # stat_compare_means() + 
    xlab(" ") +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
    theme(axis.text.x = element_text(size = 6, angle = 90)) +
    labs(
      title = paste0(Name," : Results of significancy tests."),
      subtitle = paste0("Statistics of the significancy tests of ",Name," for each scenario splitted by group, all replicates being merged together.\nFor each replicates, the mean number of significant test is computed and plotted."),
      caption = paste0("The mean represent the percentage of mean observed values that are significantly different from the null model for each sample.\nn = Number of Replicates (100)."),
    ) +
    arrange_theme()
  
  # Print a message
  cli_alert_success(paste0(Name))

  # Return the two plots
  return(list(P5,P6,P7,P8))
}

    # ----- #

  #### . Function Call . ####

# Apply ploturbo to all the metrics
All_Plot <- lapply(Metric_Total,ploturbo)

#### . Patchworking and savings . ####

cat(rule(left = "Saving the plot report ...", line_col = "yellow", line = "-", col = "br_yellow"))
# Save the plot
pl <- list_flatten(All_Plot)
# Transform them into grobs
pl <- lapply(pl, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)
# Save a global plot
ggsave(filename = paste0(Parameters$output),
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/coMet_ComSim_Outputs"),
       width = 35,
       height = 25,
       units = "cm")

# New version with a BoxPLot instead
boxploturbo <- function(List) {
  
  # Extract the Metric Name
  Name <- colnames(List$Mean_Raw)[6]
  # If the metric is a beta metric, we need to change the name 
  if (Name == "Sample_B") {
    Name <- colnames(List$Mean_Raw)[8]
  }

BoxPlot_Raw_Mean2 <- List$Mean_Raw %>%
  group_by(Scenario) %>%
  ggplot(aes(x = .data[[Name]], y = Scenario, fill = Scenario)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
  geom_boxplot(alpha = 0.5) +
  scale_fill_viridis_d() + 
  # geom_signif(comparisons = combn(Scenario.List,2, simplify = F),
  #             map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
  #             margin_top = 0.2,
  #             step_increase = 0.1,
  #             tip_length = 0.01) +
  theme(legend.position="right") +
  labs(caption ="P_value: *** = 0.001 / ** = 0.01 / * = 0.05 ") +
  xlab(" ") + arrange_theme() 
}

# Apply boxploturbo to all the metrics
All_BoxPlot <- lapply(Metric_Total,boxploturbo)

# Save the plot
pl <- list_flatten(All_BoxPlot)
# Transform them into grobs
pl <- lapply(pl, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=3,ncol = 11)
ml
##### ______________________________________________________________________________________________________________________________________ #####

# ---------------------------- #
##### STEP 3: GLOBAL PLOTS #####
# ---------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

##### . Total files . #####

# Split and combine together the raw data and the signif data.
# Raw_Total <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {
#   # Extract the raw data
#   Raw <- as.data.frame(flatten(Metric_Total[[l]][1])) %>%
#   # Remove the columns of Null-Model values
#   select(-contains("NM")) %>% 
#   # Pivot the metric value. 
#   pivot_longer(last_col(), names_to = "Metric", values_to = "Value") %>%
#   # Remove any unwanted columns to avoid differences between alpha and Beta metrics
#   select(-any_of(c("Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial"))) %>%
#   # Change every metadata column to factor
#   mutate(across(any_of(c("Scenario","Rep","Sample","Group","Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial")),as.factor))
#   
# }

# Signif_Total <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {
#   # Extract the raw data
#   Signif <- as.data.frame(flatten(Metric_Total[[l]][2])) %>%
#     # Append the metric
#     mutate("Metric"= Metric.Count.Lite[l], .before = 1) %>%
#     # Remove the p-values
#     select(!starts_with('P_value')) %>%
#     # Pivot the data of test results
#     pivot_longer( cols = starts_with("Test"),
#                   names_to = "variable",
#                   values_to = "Results") %>%
#     # Transform the wanted columns into character
#     mutate(across(c("Scenario","Metric","Rep","Group","variable"), as.factor))
# }

Mean_Raw_Total <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {
  # Extract the raw data
  Raw <- as.data.frame(Metric_Total[[l]]$Mean_Raw) %>%
    # Remove the columns of Null-Model values
    select(-contains("NM")) %>% 
    # Pivot the metric value. 
    pivot_longer(last_col(), names_to = "Metric", values_to = "Value") %>%
    # Remove any unwanted columns to avoid differences between alpha and Beta metrics
    select(-any_of(c("Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial"))) %>%
    # Change every metadata column to factor
    mutate(across(any_of(c("Scenario","Rep","Sample","Group","Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial")),as.factor)) %>%
    # Add in a ugly way the type of metric (alpha or beta)
    mutate(Type = ifelse(Metric %in% c("Bst","PCD","PDb","PhyloSor","PhyloSor_Ab","PhyloSor_Turn","PIst","Pst","S_Beta","S_Jaccard","S_Ochiai","S_SokalSneath","S_Sorensen","UniFrac","UniFrac_Turn","UniFrac_Ab"), "Beta", "Alpha" ))
}

Mean_Signif_Total <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {
  # Extract the raw data
  Signif <- as.data.frame(Metric_Total[[l]]$Mean_Signif) %>%
    # Append the metric
    mutate("Metric"= Metric.Count.Lite[l], .before = 1) %>%
    # Remove the p-values
    dplyr::select(!starts_with('P_value')) %>%
    # Pivot the data of test results
    pivot_longer( cols = starts_with("Test"),
                  names_to = "variable",
                  values_to = "Results") %>%
    # Remove any unwanted columns to avoid differences between alpha and Beta metrics
    dplyr::select(-any_of(c("Sample_A","Sample_B", "Habitat_A", "Habitat_B", "Spatial_A", "Spatial_B", "Habitat", "Spatial","Sample"))) %>%
    # Change every metadata column to factor
    mutate(across(-any_of(c("Results")),as.factor)) %>%
    # Add in a ugly way the type of metric (alpha or beta)
    mutate(Type = ifelse(Metric %in% c("Bst","PCD","PDb","PhyloSor","PhyloSor_Ab","PhyloSor_Turn","PIst","Pst","S_Beta","S_Jaccard","S_Ochiai","S_SokalSneath","S_Sorensen","UniFrac","UniFrac_Turn","UniFrac_Ab"), "Beta", "Alpha" )) %>%
    # Add in a ugly way the type of comparison (Intra or inter)
    mutate(Comp = ifelse(Group %in% c("B2_B1","C1_B2","C2_B1","C2_C1"), "Inter", "Intra" ))
}

  ### ------- ###

  #### _MEAN_TESTING_ ####

  # The goal of this analysis is to know for each metric, if the observed values are significantly different between scenarios. 
  # To do so, we have to split the data for each metric, compute the test between each scenario. 

Mean_Test <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {

Kruskall_Test <- Mean_Raw_Total %>%
  # EVENTUALLY, for the first axes of the article, we remove two scenarios SF1 and SF2, modifying the significativity of the metrics. 
  # filter(Scenario != "Equalizing100_Sig2_SF1" & Scenario != "Equalizing100_Sig2_SF2") %>%
  # Filter based on the metric
  filter(Metric == Metric.Count.Lite[l]) %>%
  # Compute the test
  dunn_test(Value ~ Scenario, p.adjust.method = "holm") %>%
  # Add the metric name
  mutate("Metric" = Metric.Count.Lite[l], .before = 1) %>%
  # Add a modal column ( Is the test significant or not) 
  mutate(Signif = ifelse(p.adj.signif %in% c("ns"), "ns", "s" ))

}

# ----- #

# Have a version with only non-significant comparisons
Mean_Test_NS <- Mean_Test %>%
  # Filter to only keep the non-significant values
  filter(p.adj.signif == "ns")

# ----- #

# Count the number of significant results for each metric to sort by the number of significant comparisons
# If the number of tests significatif is 28, it means that all tests were significatives. 
Mean_Test_Count <- Mean_Test %>% 
  # Group
  group_by(Metric) %>%
  # Count the number of significative tests
  count(Signif, sort = TRUE)
  
# Modify the data frame to have Metric, Number of significatives and number of not significant data
Mean_Test_Count <- full_join(filter(Mean_Test_Count,Signif == "s"),filter(Mean_Test_Count,Signif == "ns"), by = c("Metric")) %>%
  # Select and change the names of the columns
  select(Metric,"s" = n.x,"ns" = n.y) %>%
  # CHnage the NA's (if no non-significant tests) by 0.
  replace_na(list(s = 0, ns = 0))

# ----- #

# Extract only the comparisons with the stochastic processes

Mean_Test_Stochastic <- Mean_Test %>% 
  # Filter the comparison were the stochastic scenario is involved
  filter(group1 == "Stochastic100_Sig5" | group2 == "Stochastic100_Sig5") %>%
  # Count the number of significative tests
  count(Signif, sort = TRUE)


# Version without differentiation between the different scenario 
# (Therefore not a pairwise differentiation between scenarios but a global appreciation for each metric)
Mean_Test_Unsplitted <- foreach(l = 1:length(Metric_Total), .combine = full_join) %do% {
  
  Kruskall_Test <- Mean_Raw_Total %>%
    # Filter based on the metric
    filter(Metric == Metric.Count.Lite[l]) %>%
    # Compute the test
    kruskal_test(Value ~ Scenario) %>%
    # Add the metric name
    mutate("Metric" = Metric.Count.Lite[l], .before = 1)
  
}

# Compute the mean per comparison inter or intra plot.

Mean_InterIntra <- Mean_Signif_Total %>%
  # Compute the mean and sd by Scenario and by type of test
  group_by(Metric,Scenario,variable,Comp) %>%
  dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
  # Remove the results of "Test different" for better visibility.  
  filter(variable != "Test_Different") %>%
  # Mutate 
  mutate(across(c(Mean,Sd), round, 2))
  
# ----- #

# Save the data
write.csv(Mean_Test, paste0(getwd(),"/coMet_ComSim_Outputs/Obs_Mean_Comparison.csv"))
write.csv(Mean_Test_NS, paste0(getwd(),"/coMet_ComSim_Outputs/Obs_Mean_NS_Comparison.csv"))
write.csv(Mean_Test_Count, paste0(getwd(),"/coMet_ComSim_Outputs/Obs_Mean_Comparison_Count.csv"))
write.csv(Mean_InterIntra, paste0(getwd(),"/coMet_ComSim_Outputs/Obs_Mean_Comparison_InterIntra.csv"))

# ----- #

  ### ------- ###

  #### _RAW_ ####

#### . Significant Results ~ Scenario + Metric . ####

# GP_1 <- Signif_Total %>%
#   # Compute the mean and sd by Scenario and by type of test
#   group_by(Metric,Scenario,variable) %>%
#   dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
#   # Draw the plot 
#   ggplot(aes(x = Metric, y = Mean, fill = variable)) +
#   geom_bar(stat = "identity", position = "dodge")  +
#   geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = variable), width=0.4, alpha=0.9, linewidth=0.5,
#                  position=position_dodge(.9)) +
#   coord_flip() +
#   facet_grid( ~ variable  + Scenario , labeller = "label_both") +
#   # theme
#   theme(axis.title.x=element_blank(),
#         legend.position = "bottom") +
#   guides(colour = "none")

#### . Significant Results ~ Metric + Group . ####

  # We have to do one Scenario at the time or it will be .... unreadable. 

# GP_2 <- foreach(L = 1:length(Scenario.List)) %do% {
# 
#   # Filter the dataset to only keep one scenario
#   Data <- filter(Signif_Total,Scenario == Scenario.List[L])
#   
# Plot <- Data %>%
#   # Compute the mean and sd by Scenario and by type of test
#   group_by(Metric,Scenario,variable,Group) %>%
#   dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
#   # Draw the plot 
#   ggplot(aes(x = Metric, y = Mean, fill = Group)) +
#   geom_bar(stat = "identity", position = "dodge")  +
#   #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
#   #               position=position_dodge(.9)) +
#   coord_flip() +
#   facet_grid(. ~ Scenario + variable) +
#   # theme
#   theme(axis.title.x=element_blank(),
#         legend.position = "right") +
#   guides(colour = "none")
#  
# }

#### _MEAN_ ####

GP_2 <- Mean_Raw_Total %>%
  group_by(Scenario) %>%
  ggplot(aes(x = Value , y = Scenario, fill = Scenario)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
  geom_boxplot(alpha = 0.5) +
  scale_fill_viridis_d() + 
  # geom_signif(comparisons = combn(Scenario.List,2, simplify = F),
  #             map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
  #             margin_top = 0.2,
  #             step_increase = 0.1,
  #             tip_length = 0.01) +
  theme(legend.position="right") +
  facet_grid(Type ~ Metric, scales = "free") +
  labs(caption ="P_value: *** = 0.001 / ** = 0.01 / * = 0.05 ") +
  xlab(" ") + arrange_theme() 



#### . Significant Results ~ Scenario + Metric . ####

GP_3 <- Mean_Signif_Total %>%
  # Compute the mean and sd by Scenario and by type of test
  group_by(Metric,Scenario,variable,Type) %>%
  dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
  # Remove the results of "Test different" for better visibility.  
  filter(variable != "Test_Different") %>%
  # Draw the plot 
  ggplot(aes(x = Metric, y = Mean, fill = variable, group = variable, alpha = 0.5)) +
  scale_fill_manual(values = c("forestgreen","red3")) +
  geom_bar(stat = "identity", position = "dodge")  +
  geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = variable), width=0.4, alpha=0.9, linewidth=0.5,
                position=position_dodge(.9)) +
  scale_color_manual(values = c("forestgreen","red3")) +
  coord_flip() +
  facet_grid(Type ~ Scenario, scales = "free_y") +
  # theme
  theme(axis.title.x= element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8),
        legend.position = "bottom") +
  guides(colour = "none", alpha = "none") 

#### . Significant Results ~ Metric + Group . ####

# We have to do one Scenario at the time or it will be .... unreadable. 

GP_4 <- foreach(L = 1:length(Scenario.List)) %do% {
  
  # Filter the dataset to only keep one scenario
  Data <- filter(Mean_Signif_Total,Scenario == Scenario.List[L])
  
  # Get the Alpha and Beta Metrics splitted
  Alpha <- filter(Data, Type == "Alpha")
  Beta <- filter(Data, Type == "Beta")
  
  Plot_Alpha <- Alpha %>%
    # Compute the mean and sd by Scenario and by type of test
    group_by(Metric,Scenario,variable,Group,Type) %>%
    dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
    # Remove the results of "Test different" for better visibility.  
    filter(variable != "Test_Different") %>%
    # Draw the plot 
    ggplot(aes(x = Metric, y = Mean, fill = variable, alpha = 0.5)) +
    scale_fill_manual(values = c("forestgreen","red3")) +
    geom_bar(stat = "identity", position = "dodge")  +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
    #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
    #               position=position_dodge(.9)) +
    coord_flip() +
    facet_grid(Scenario ~ Group) +
    # theme
    theme(axis.title.x=element_blank(),
          legend.position = "top") +
    guides(colour = "none",alpha = "none")
  
  Plot_Beta <- Beta %>%
    # Compute the mean and sd by Scenario and by type of test
    group_by(Metric,Scenario,variable,Group) %>%
    dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
    # Remove the results of "Test different" for better visibility.  
    filter(variable != "Test_Different") %>%
    # Draw the plot 
    ggplot(aes(x = Metric, y = Mean, fill = variable, alpha = 0.5)) +
    scale_fill_manual(values = c("forestgreen","red3")) +
    geom_bar(stat = "identity", position = "dodge")  +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
    #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
    #               position=position_dodge(.9)) +
    coord_flip() +
    facet_grid(Scenario ~ Group) +
    # theme
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    guides(colour = "none",alpha = "none")
  
  # Return the two graphs
  Plot_Alpha + Plot_Beta + plot_layout(ncol = 1) # Plotted with patchwork to have the plots aligned (and not displaced because of the axis names length)
  
}

#### . Significant Results ~ Metric + Group (Inter or Intra) . ####

# We have to do one Scenario at the time or it will be .... unreadable. 

GP_5 <- foreach(L = 1:length(Scenario.List)) %do% {
  
  # Filter the dataset to only keep one scenario
  Data <- filter(Mean_Signif_Total,Scenario == Scenario.List[L])
  
  # Get the Alpha and Beta Metrics splitted
  Alpha <- filter(Data, Type == "Alpha")
  Beta <- filter(Data, Type == "Beta")
  
  Plot_Alpha <- Alpha %>%
    # Compute the mean and sd by Scenario and by type of test
    group_by(Metric,Scenario,variable,Comp,Type) %>%
    dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
    # Remove the results of "Test different" for better visibility.  
    filter(variable != "Test_Different") %>%
    # Draw the plot 
    ggplot(aes(x = Metric, y = Mean, fill = variable, alpha = 0.5)) +
    scale_fill_manual(values = c("forestgreen","red3")) +
    geom_bar(stat = "identity", position = "dodge")  +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
    #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
    #               position=position_dodge(.9)) +
    coord_flip() +
    facet_grid(Scenario ~ Comp) +
    # theme
    theme(axis.title.x=element_blank(),
          legend.position = "top") +
    guides(colour = "none",alpha = "none")
  
  Plot_Beta <- Beta %>%
    # Compute the mean and sd by Scenario and by type of test
    group_by(Metric,Scenario,variable,Comp) %>%
    dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
    # Remove the results of "Test different" for better visibility.  
    filter(variable != "Test_Different") %>%
    # Draw the plot 
    ggplot(aes(x = Metric, y = Mean, fill = variable, alpha = 0.5)) +
    scale_fill_manual(values = c("forestgreen","red3")) +
    geom_bar(stat = "identity", position = "dodge")  +
    scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
    #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
    #               position=position_dodge(.9)) +
    coord_flip() +
    facet_grid(Scenario ~ Comp) +
    # theme
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    guides(colour = "none",alpha = "none")
  
  # Return the two graphs
  Plot_Alpha + Plot_Beta + plot_layout(ncol = 1) # Plotted with patchwork to have the plots aligned (and not displaced because of the axis names length)
  
}

#### . Significant Results ~ Beta Metric + Group (Inter or Intra) + Equalizing Scenarios . ####

Scenario.List.Bis <- c("Equa50Stab50_Sig2","Equa50Stab50_Sig5","Equalizing100_Sig2_SF1","Equalizing100_Sig2_SF2","Equalizing100_Sig2","Equalizing100_Sig5")

# Filter the dataset to only keep one scenario
Data <- filter(Mean_Signif_Total,Scenario %in% Scenario.List.Bis) %>%
   filter(Type == "Beta")
  
GP_6 <- Data %>%
  # Compute the mean and sd by Scenario and by type of test
  group_by(Metric,Scenario,variable,Comp) %>%
  dplyr::summarize(Mean = mean(Results, na.rm=TRUE), Sd = sd(Results, na.rm = TRUE)) %>%
  # Remove the results of "Test different" for better visibility.  
  filter(variable != "Test_Different") %>%
  # Draw the plot 
  ggplot(aes(x = Metric, y = Mean, fill = variable, alpha = 0.5)) +
  scale_fill_manual(values = c("forestgreen","red3")) +
  geom_bar(stat = "identity", position = "dodge")  +
  scale_y_continuous(name="% Significant tests", limits=c(0, 1)) +
  #geom_errorbar(aes(ymin=Mean, ymax=Mean+Sd, color = Group), width=0.4, alpha=0.9, linewidth=0.5,
  #               position=position_dodge(.9)) +
  coord_flip() +
  # facet_grid(Scenario ~ Comp) +
  facet_grid(Comp ~ Scenario) +
  # theme
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  guides(colour = "none",alpha = "none")
  
ggsave(filename = "Metric_Signif_Beta.png",
       plot = GP_6,
       device = "png",
       path = paste0(getwd(),"/coMet_ComSim_Outputs"),
       width = 35,
       height = 25,
       units = "cm")

#### . Patchworking and savings . ####

cat(rule(left = "Saving the plot report ...", line_col = "yellow", line = "-", col = "br_yellow"))
# Save the plot
pl <- list_flatten(list(GP_3,GP_4,GP_5))
# Transform them into grobs
pl <- lapply(pl, as.grob)
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)
# Save a global plot
ggsave(filename = "Metric_Signif_Total.pdf",
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/coMet_ComSim_Outputs"),
       width = 35,
       height = 25,
       units = "cm")

# END OF THE SCRIPT
cat(rule(left = "SCRIPT COMET_METRIC_PLOTS.R ENDING", line_col = "red", line = "-", col = "br_red")) 


# ---------------------------- #
##### BONUS: VARIOUS PLOTS #####
# ---------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

  # ----- #

#### . Phylo-Trees and Traits . ####

#
# library(ape)
# library(untb)
# library(ggtree)
#
# # Create a false community
#   Comm <- MetaCom <- fisher.ecosystem(1000000, 1000, 800)
# # Load the phylo_tree associated
#   Phylo_tree <- read.tree(paste0(getwd(),"/coMet_ComSim_Outputs/Foo_a1/Whole_Communities/ScenarioFoo_a1_Step1_Total_PhyloTrees.tree"))$rep1
# # Load the trait values associated.
#   Trait <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/Foo_a1/Whole_Communities/ScenarioFoo_a1_Step1_Total_Species_Traits.csv"))
#   # Only keep the first replicate
#   Trait <- Trait[which(Trait$rep == 1),]
# # Make colored tree ---------------------------- #
#
# # Transform the tree into a tibble
#   Phylo_tbl <- tidytree::as_tibble(Phylo_tree)
# # It can be transformed back with as.phylo()
#
# # Make sure that the trait data contains a columns of the species (tips) names of the same name that in the phylo_tibble
#   colnames(Trait)[3] <- "label"
#
# # Join the data of the tree and the data of the traits
#   PhyData <- full_join(Phylo_tbl,Trait, by = 'label')
#
# # Transform it into a "treedata" object
#   TreeData <- tidytree::as.treedata(PhyData)
#
#   # Plot the tree with trait value
#   ggtree(TreeData)
#
#   # For display
#   ggsave(filename = paste0("Species_Traits1_Plot.tiff"), path = paste0(getwd(),"/Markdown/"), width = 20, height = 20, device='png', dpi=700)
#
# # Plot the tree with trait value
#   ggtree(TreeData, aes(color=Trait_2)) +
#     scale_color_continuous(low="blue", high="red") +
#     theme(
#       legend.position = c(.05, .95),
#       legend.justification = c("left", "top"),
#       legend.box.just = "left",
#       legend.margin = margin(6, 6, 6, 6),
#       legend.background = element_rect(fill="white", size=.5)
#     ) +
#     theme(plot.background = element_rect(colour = "black", size = 1))
#
#   ggsave(filename = paste0("Species_Traits2_Plot.tiff"), path = paste0(getwd(),"/Markdown/"), width = 20, height = 20, device='png', dpi=700)
#
#
#   # ---- Present ALL the graphs simultaneously ---- #
#
#   # Create a string with a list of color
#   Col <- rainbow(n = 20) # Here, we want n contiguous colors
#
#   # ---- #
#
#   # Initiate a list
#   TreeList <- vector(mode = "list", length = 9) # Number of variables/traits to test
#
#   # Initiate a loop to make the plots for all plots
#   for (i in 1:length(TreeList)){
#     # Create a tree
#     TreeList[[i]] <- ggtree(TreeData, aes_string(color=paste0("Trait_",i))) +   # Select the trait to map, it is important to use "aes_string" to use a string of character with paste0 during the loop
#       scale_color_continuous(low=sample(Col,1), high=sample(Col,1)) +           # Pick two colors randomly from the palette to represent the beginning and the end colors of the continuous scale
#       # Change legend things
#       theme(legend.position = c(.05, .95),legend.justification = c("left", "top"),legend.box.just = "left",legend.margin = margin(6, 6, 6, 6),legend.background = element_rect(fill="white", size=.5))
#   }
#
#   # Arrange them all.together and visualize them
#   # do.call(grid.arrange, c(TreeList, list(ncol=3)))
#
#   # Save the plot
#   # Save it with arrangeGrob
#   Plot_to_save <- do.call(arrangeGrob, c(TreeList, list(ncol=3)))
#   ggsave(filename = paste0("S0_Simul_Species_Traits_Plot.tiff"), plot = Plot_to_save, path = paste0(getwd(),"/Simulation_Outputs/Scenario",Scenario,"/"), width = 20, height = 20, device='tiff', dpi=700)
#   # For the markdown
#   ggsave(filename = paste0("S0_Simul_Species_Traits_Plot.png"), plot = Plot_to_save, path ="C:/Users/Administrateur/Desktop/Simulation_Doctorat/Markdown/figs", width = 20, height = 20, device='png', dpi=200)
#

  # Make a bar plot ------------------------------ #
  #
  #   # Create a phylo4d object
  #   Trait_Tree <- phylo4d(rtree,TraitsEsp)
  #
  #   # Create the bar plot
  #   BP_Trait <- barplot(Trait_Tree,
  #                       # trait = c("Trait_1","Trait_2"), # Choose the traits to show
  #                       show.tip = FALSE)    # Remove the tips
  #
  #

# ----- #

#### . Ungrouped Significativity  . ####

 # # Load the Data
 # Data <- read.csv(paste0("coMet_ComSim_Outputs/",Scenario,"/Metrics/Scenario",Scenario,"_PD_Faith_Signif_Ungrouped.csv"),row.names = 1)
 #
 # # Are the tests significantly different ?
 #  P1 <- Data %>%
 #    group_by(type,SampleType) %>%  # Group_by to have two groups in a tibble for faceting
 #    mutate_at(vars(rep,type,Habitat,SampleType), factor) %>%  # Transform all the stuff into factor
 #    ggplot(aes(x=SampleType, y=P_value_Different, group = SampleType, fill = SampleType)) +
 #    scale_fill_manual(values = c("#CA5310","#FBBA72")) +
 #    geom_boxplot() +
 #    facet_grid(. ~ Habitat ,labeller = labeller(Scenario = label_value, Habitat = label_both))+
 #    # stat_compare_means() +
 #    xlab(" ") +
 #    ylim(c(0,1.3)) ; P1

# ----- #

#### . Mean values of the Metrics . ####

 # The goal of these plots is to show the value of the Mean Observed metric compared to the density plot of the Null Metric

# --- 1 Metric at a time --- #

# # Apply ploturbo
# Data <- Metric_Total$PD_Faith_Mean$Raw
#
# # Make the plot for one line
# Data_Test <- Data[1,7:ncol(Data)]
#
# # Compute the density BEFORE
# dens <- density(as.numeric(Data_Test[2:101]))
#
# # Plot the density
# ggplot(data.frame(x = dens$x, y = dens$y)) +
#   aes(x = x, y = y) + geom_line() +
#   geom_vline(aes(xintercept = Data_Test[[1]])) # Add the mean
#
#
# Plot_Test <- ggplot(Data_Test,aes(x = PD_Mean)) +
#   geom_vline(aes(xintercept = PD_Mean)) +
#   geom_density(data = Data_Test); Plot_Test
#
# ggplot(data.frame(x = Data_Test)) +
#   aes(x = x) + geom_density()
#
# # Create the density plot for one line
# Dens_Mean <- Data[1,] %>%
#   group_by(Scenario,rep,sample,type,Habitat,SampleType) %>%
#   mutate_at(vars(Scenario,rep,sample,type,Habitat,SampleType), factor) %>%  # Transform all the stuff into factor
#   ggplot(aes(x = PD_Mean, group = rep, fill = rep)) +  # We need the strange ."as.name" use with "aes_" to resolve the problem of finding the column containing the metric values
#   geom_density(data = Data[1,7:ncol(Data)],alpha = 0.3) +
#   geom_vline(aes(xintercept = PD_Mean)) ; Dens_Mean
#
#
#
#
#   facet_grid(Scenario ~ . + SampleType, labeller = labeller(Scenario =label_value, SampleType = label_both))+
#   scale_fill_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
#   geom_vline(data  = DataMean, aes(xintercept = mean_x, colour = Habitat)) +
#   scale_color_manual(values = c("1" = "#5FAD56","2" ="#8B80F9","1-2" ="grey50")) +
#   theme(legend.position="right") +
#   xlab(" ") ; Dens_Raw



