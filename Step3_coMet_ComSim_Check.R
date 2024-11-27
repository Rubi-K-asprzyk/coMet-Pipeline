#!/usr/bin/env Rscript

##### ______________________________________________________________________________________________________________________________________ #####

suppressPackageStartupMessages(if(!require(cli)){install.packages("cli");library(cli)})
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");require(pacman)})
cat(rule(left = "SCRIPT COMET_COMSIM_CHECK.R BEGINNING", line_col = "red", line = "-", col = "br_red"))

# --------------------------------------------------------------------- #
# coMet_ComSim_Check: Exploratory analyses of the simulated communities #
# --------------------------------------------------------------------- #

# This script is meant to display useful informations from the outputs of the coMet_ComSim.R script
# to determine and verify if the simulated communities are suitable for the rest of the coMet pipeline.

# It takes as input the Config_File of the simulations we want to test to retrieve the corresponding outputs with the simulation name. 

# ------------------------------ #
##### STEP 0: INITIALISATION #####
# ------------------------------ #

##### ______________________________________________________________________________________________________________________________________ #####

cat(rule(left = "INITIALISATION", line_col = "green", line = "-", col = "br_green"))

    # ----- #

#### . Packages . ####

# Print a message
cli_alert_info("Loading the needed packages ... ")

# Packages
# Install/load pacman 
suppressPackageStartupMessages(if(!require(pacman)){install.packages("pacman");require(pacman)})
# Install/load tons of packages
p_load("dplyr",
       "plyr",
       "tidyr",
       "doParallel",
       "stringr",
       "ggplot2",
       "forcats",
       "foreach",
       "argparser",
       "rlist",
       "ggpubr",
       "gridExtra",
       "ggplotify",
       "viridis",
       "paletteer",
       "ggrepel",
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
       "ape",
       "DescTools"
) 

# Print a message
cli_alert_success("Packages correctly loaded !")

# Avoid summarise() messages 
options("dplyr.summarise.inform" = FALSE)

# Avoid undesired warning messages (set 0 to default)
options(warn = -1)

  #### . Local Mode . ####

# Parameters <- list()
# Parameters$Config <- "Foo.R"

    # ----- #

  #### . Argument Parser . ####

cli_alert_info("Reading the coMet_ConfigFile.R ...")

# Create the parser
arg_parser <- arg_parser("Check the communities computed by coMet_ComSim and create an output pdf.", name = "coMet_ComSim_Check.R", hide.opts = FALSE)

# Add the config_files.R as positional arguments
arg_parser <- add_argument(arg_parser, arg = "Configuration", nargs = 1, help = "Configuration file for coMet_ComSim.R")

# Parse the arguments
Parameters <- parse_args(parser = arg_parser, argv = commandArgs(trailingOnly = T))

# Load them into the environnement
source(Parameters$Config)

# Print a message to inform that the right configuration file was correctly loaded. 
cli::cli_alert_info(paste0("Working on configuration file: ",Parameters$Config))

# Print a message
cli_alert_success("coMet_ConfigFile.R correctly read !")

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
    axis.ticks = element_line(color = "grey91", size = .5),
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
      legend.box.background = element_rect(colour = "black")
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
  
# Print a message
cli_alert_success("Theme set !")

# ------------------------------------------- #
##### STEP 1: LOADING THE COMMUNITY FILES #####
# ------------------------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "LOADING OF THE COMMUNITY FILES", line_col = "green", line = "-", col = "br_green"))

    # ----- #

  #### . Loading files . ####

# Total Recruitment results of the simulations
Result_Count_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Result_Count.csv"),row.names = 1)
# Total Recruitment by distance
Rec_by_Dist_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Rec_by_Dist.csv"),row.names = 1)
# Total Abundance
Abundance_Total  <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Abundance.csv"),row.names = 1)
# Total Occurrence
Occurrence_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Occurence.csv"),row.names = 1)
# Total Species Richness
SR_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_SR.csv"),row.names = 1)
# Total Species Traits
Traits_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Traits.csv"),row.names = 1)
# Total species fitness
Fitness_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Fitness.csv"),row.names = 1)
# Sampled Sites data
Sites_Sampled <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Sites.csv"),row.names = 1)
# Load the global sites data
sites <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step2_Sites_Data.csv"),row.names = 1) %>%
  mutate(across(c("hab","SampleType","sample"),as.factor)) # Change into numeric and factor
    # ----- #

#   #### . Testing files . ####
# 
# # Abundance Data
# if(length(unique(Abundance_Total$rep)) != Nrep) {
#   stop("Number of Sampled_Abundance replicates different from the number of Replicates ; EXPECTED-NREP = ",Nrep," - SEEN-NREP = ",length(Abundance_Total))
# } else if (nrow(Abundance_Total) != nbe_sample*Nrep){
#   stop("Number of sampled plots in Abundance_Data different from expected ! EXPECTED: ",nbe_sample*Nrep," Rep - SEEN: ",nrow(Abundance_Total))
# } else {cli_li(cat(col_yellow("DONE"),": Abundance_Data","\n"))}
# 
# # Occurrence Data
# if(length(unique(Occurrence_Total$rep)) != Nrep){
#   stop("Number of Sampled_Occurence replicates different from the number of Replicates ; EXPECTED-NREP = ",Nrep," - SEEN-NREP = ",length(Occurrence_Total))
# } else if (nrow(Occurrence_Total) != nbe_sample*Nrep){
#   stop("Number of sampled plots in Pre_Abs_Data different from expected ! EXPECTED: ",nbe_sample*Nrep," Rep - SEEN: ",nrow(Occurrence_Total))
# } else {cli_li(cat(col_yellow("DONE"),": Pre_Abs_Data","\n"))}
# 
# # Sites Data
# if(max(Sites_Sampled$sample) != nbe_sample){
#   stop("The last sample number is not the same as the number of sample expected ; EXPECTED = ",nbe_sample," - SEEN = ", max(Sites_Sampled$sample))
# } else {cli_li(cat(col_yellow("DONE"),": Sites_Data","\n"))}
# 
# # Species Richness Data
# if(length(unique(SR_Total$rep)) != Nrep){
#   stop("Number of Sampled_SR Data replicates files different from the number of Replicates ; NREP = ",Nrep," - FILES = ",length(SR_Total))
# } else {cli_li(cat(col_yellow("DONE"),": Richness_Data","\n"))}
# 
# # Test each files for the number of sampled plots.
# Nihil <- foreach(i = 1:Nrep) %dopar% {
#   if (nrow(Abundance_Total[which(Abundance_Total$rep == i),]) != nbe_sample)
#   {stop("Number of sampled plots in Abundance_Data / Replicate ",i," is different from expected ! EXPECTED: ",nbe_sample," - SEEN: ",nrow(Abundance_Total[which(Abundance_Total$rep == i),]))}
#   if (nrow(Occurrence_Total[which(Occurrence_Total$rep == i),]) != nbe_sample)
#   {stop("Number of sampled plots in Pre_Abs_Data / Replicate ",i," is different from expected ! EXPECTED: ",nbe_sample," - SEEN: ",nrow(Occurrence_Total[which(Occurrence_Total$rep == i),]))}
#   if (nrow(SR_Total[which(SR_Total$rep == i),]) != nbe_sample)
#   {stop("Number of sampled plots in Richness_Data / Replicate ",i," is different from expected ! EXPECTED: ",nbe_sample," - SEEN: ",nrow(SR_Total[which(SR_Total$rep == i),]))}
# } ; rm(Nihil)

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
cli_alert_info("Recruitment plots ... ")

##### ------------ 2.A.a : % Success Total  ------------ #####

  #### . Statistics . ####

# -- Compute the number of successful recruitment for each category
sum_success <- Result_Count_Total %>%
  mutate(across(c(rep,gen),as.factor)) %>% # Transform what needs to be transformed into factors
  dplyr::select(contains("nbe_success"), gen, rep) %>% 
  pivot_longer(cols = contains("nbe_success"), values_to = "value_success") %>% 
  mutate(name = str_remove(name, "_success"))

# -- Compute the number of failed recruitment for each category
sum_fail <- Result_Count_Total %>% 
  mutate(across(c(rep,gen),as.factor)) %>% # Transform what needs to be transformed into factors
  dplyr::select(contains("nbe_failed"), gen, rep) %>% 
  pivot_longer(cols = contains("nbe_failed"), values_to = "value_failed") %>% 
  mutate(name = str_remove(name, "_failed")) 

# -- Compute the mean and SD ratio of failed recruitment for each category across replicate.
ratio_succes <- 
  left_join(sum_success, sum_fail, by = c("gen", "rep", "name")) %>% 
  mutate(prop = (value_failed/(value_success + value_failed))*100) %>% 
  group_by(gen, name) %>% 
  dplyr::summarise(across(prop, list(mean = mean, sd = sd)))

# -- Compute the mean and SD for each numbers across replicate.
summary_stats_sumarized <- 
  Result_Count_Total %>% 
  group_by(gen) %>% 
  dplyr::summarise(across(.cols = starts_with("nbe_"), list(mean = mean, sd = sd)))

# -- Compute the mean and SD across replicate of the numbers of recruitment for each category
stats_mean <- 
  summary_stats_sumarized %>% 
  pivot_longer(cols = contains("_mean"), names_to = "stats", values_to = "mean_value") %>% 
  dplyr::select(gen, stats, mean_value) %>% 
  mutate(stats = str_remove(stats, "_mean"))

stats_sd <-
  summary_stats_sumarized %>% 
  pivot_longer(cols = contains("_sd"), names_to = "stats", values_to =  "sd_value") %>%
  dplyr::select(gen, stats, sd_value) %>% 
  mutate(stats = str_remove(stats, "_sd"))

# -- Joining mean and SD by gen.
Result_Count <- left_join(stats_mean, stats_sd,  by = c("gen", "stats")) %>%
  mutate(across(c(stats,gen),as.factor))
# Save the data of recruitment
write.csv(Result_Count, paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step2_Results_Count.csv"))

    # ----- #

  #### . Plots . ####

# Number of failed recruitment / habitats ---
GG_Rec_1 <- ggplot(data = Result_Count %>%
                   filter(grepl(paste(c("_failed_diff", "_failed_same"), collapse = "|"), stats)),
                   aes(x = gen, group = stats)) +
 geom_line(aes(y = mean_value, color = stats), linewidth = 1) +
 geom_ribbon(aes(
   y = mean_value,
   ymin = mean_value - sd_value,
   ymax = mean_value + sd_value,
   fill = stats),
   alpha = .2) +
 scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
 ylab(label = "Number of failed recruitment") +
 xlab(label = "Generation") +
 scale_color_discrete(labels = c("Different","Same")) +
 guides(fill = "none") +
 labs(
    title = "Number of failed recruitment by habitat.",
    subtitle = "Number of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from a plot from the same or a different habitat.",
    color= NULL
  )

# Proportion of failed recruitment / habitats ---
GG_Rec_2 <- ggplot(data = ratio_succes %>%
                   filter(grepl(paste(c("_diff_hab", "_same_hab"), collapse = "|"), name)),
                   aes(x = gen, group = name)) +
 geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
 geom_ribbon(aes(
   y = prop_mean,
   ymin = prop_mean - prop_sd,
   ymax = prop_mean + prop_sd,
   fill = name),
   alpha = .2) +
 scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
 ylab(label = "Proportion of failed recruitment") +
 xlab(label = "Generation") +
 scale_color_discrete(name = " ", labels = c("Different","Same")) +
 guides(fill = "none") +
 labs(
    title = "Proportion of failed recruitment by habitat.",
    subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from a plot from the same or a different habitat.",
    color= NULL
  )

# Number of failed recruitment / plots ---
GG_Rec_3 <- ggplot(data = Result_Count %>%
                   filter(grepl(paste(c("_failed_within", "_failed_between"), collapse = "|"), stats)),
                   aes(x = gen, group = stats)) +
 geom_line(aes(y = mean_value, color = stats), linewidth  = 1) +
 geom_ribbon(aes(
 y = mean_value,
   ymin = mean_value - sd_value,
   ymax = mean_value + sd_value,
   fill = stats),
   alpha = .2) +
 scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
 ylab(label = "Number of failed recruitment") +
 xlab(label = "Generation") +
 scale_color_discrete(name = " ",  labels = c("Different","Same")) +
 guides(fill = "none") +
  labs(
    title = "Number of failed recruitment by plot.",
    subtitle = "Number of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from the same or a different plot, regardless of the habitat.",
    color= NULL
  )

# Proportion of failed recruitment / plots ---
GG_Rec_4 <- ggplot(data = ratio_succes %>%
                  filter(grepl(paste(c("_within", "_between"), collapse = "|"), name)),
                  aes(x = gen, group = name)) +
  geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
  geom_ribbon(aes(
    y = prop_mean,
    ymin = prop_mean - prop_sd,
    ymax = prop_mean + prop_sd,
    fill = name),
    alpha = .2) +
 scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
 ylab(label = "Proportion of failed recruitment") +
 xlab(label = "Generation") +
 scale_color_discrete(name = " ", labels = c("Different","Same")) +
 guides(fill = "none") +
 labs(
   title = "Proportion of failed recruitment by plot.",
   subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from the same or a different plot, regardless of the habitat.",
   color= NULL
  )

# Arrange plots ---
RecArr <- as.ggplot(ggarrange(
  GG_Rec_1,GG_Rec_2,GG_Rec_3,GG_Rec_4,common.legend = T, legend = "right")) +
  labs(title = paste0("Scenario ",Scenario,": Species recruitment data by generation.")) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
          # Customize title appearence
        plot.title = element_text(
            color = "grey10", 
            size = 12, 
            face = "bold",
            margin = margin(t = 10),
            hjust = 0.5
        ))

    # ----- #

# ##### ------------ 2.A.b : % Success Replicates  ------------ #####
# 
#   # We can only use that when the number of replicates is not to high ( Less than 10 )
# 
#   #### . Statistics . ####
# 
# # -- Compute the mean and SD ratio of failed recruitment for each category for each replicates. .
# ratio_succes_rep <- 
#   left_join(sum_success, sum_fail, by = c("gen", "rep", "name")) %>% 
#   mutate(prop = value_failed/(value_success + value_failed)*100) %>% 
#   group_by(gen, name, rep) %>% 
#   dplyr::summarise(across(prop, list(mean = mean)))
# 
# # -- Compute the mean for each numbers across replicate.
# summary_stats_sumarized_rep <- 
#   Result_Count_Total %>% 
#   group_by(rep,gen) %>% 
#   dplyr::summarise(across(.cols = starts_with("nbe_"), list(mean = mean)))
# 
# # -- Compute the mean replicate of the numbers of recruitment for each category
# stats_mean_rep <- 
#   summary_stats_sumarized_rep %>% 
#   pivot_longer(cols = contains("_mean"), names_to = "stats", values_to = "mean_value") %>%
#   dplyr::select(gen, stats,rep, mean_value) %>% 
#   mutate(stats = str_remove(stats, "_mean"))
# 
#     # ----- #
# 
#   #### . Plots . ####
# 
# # Number of failed recruitment / habitats / Replicates ---
# GG_Rec_1 <- ggplot(data = stats_mean_rep %>%
#                      filter(grepl(paste(c("_failed_diff", "_failed_same"), collapse = "|"), stats)),
#                    aes(x = gen, group = stats)) +
#   geom_line(aes(y = mean_value, color = stats, ), linewidth = 1) +
#   scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) + 
#   facet_grid(rep ~ .) + 
#   ylab(label = "Number of failed recruitment") +
#   xlab(label = "Generation") +
#   scale_color_discrete(name = " ", labels = c("Different","Same")) +
#   guides(fill = "none") +
#   labs(
#     title = "Number of failed recruitment by habitat ~ Replicate.",
#     subtitle = "Number of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from a plot from the same or a different habitat."
#   )
# 
# # Proportion of failed recruitment / habitats ---
# GG_Rec_2 <- ggplot(data = ratio_succes_rep %>%
#                      filter(grepl(paste(c("_diff_hab", "_same_hab"), collapse = "|"), name)),
#                    aes(x = gen, group = name)) +
#   geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
#   facet_grid(rep ~ .) + 
#   scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
#   ylab(label = "Proportion of failed recruitment") +
#   xlab(label = "Generation") +
#   scale_color_discrete(name = " ", labels = c("Different","Same")) +
#   guides(fill = "none") +
#   labs(
#     title = "Proportion of failed recruitment by habitat ~ Replicate.",
#     subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from a plot from the same or a different habitat."
#   )
# 
# # Number of failed recruitment / plots ---
# GG_Rec_3 <- ggplot(data = stats_mean_rep %>%
#                      filter(grepl(paste(c("_failed_within", "_failed_between"), collapse = "|"), stats)),
#                    aes(x = gen, group = stats)) +
#   geom_line(aes(y = mean_value, color = stats), linewidth  = 1) +
#   facet_grid(rep ~ .) +
#   scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
#   ylab(label = "Number of failed recruitment") +
#   xlab(label = "Generation") +
#   scale_color_discrete(name = " ",  labels = c("Different","Same")) +
#   guides(fill = "none") +
#   labs(
#     title = "Number of failed recruitment by plot ~ Replicate.",
#     subtitle = "Number of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from the same or a different plot, regardless of the habitat."
#   )
# 
# # Proportion of failed recruitment / plots ---
# GG_Rec_4 <- ggplot(data = ratio_succes_rep %>%
#                      filter(grepl(paste(c("_within", "_between"), collapse = "|"), name)),
#                    aes(x = gen, group = name)) +
#   geom_line(aes(y = prop_mean, color = name), linewidth  = 1) +
#   facet_grid(rep ~ .) +
#   scale_x_discrete(breaks = pretty(1:Ngen,n = 10)) +
#   ylab(label = "Proportion of failed recruitment") +
#   xlab(label = "Generation") +
#   scale_color_discrete(name = " ", labels = c("Different","Same")) +
#   guides(fill = "none") +
#   labs(
#     title = "Proportion of failed recruitment by plot ~ Replicate.",
#     subtitle = "Proportion of failed recruitments at each generation based on the origin of the individual recruited.\nThis individual could be from the same or a different plot, regardless of the habitat."
#   )
# 
# # Arrange plots ---
# RecArr_Rep <- as.ggplot(ggarrange(
#   GG_Rec_1,GG_Rec_3,common.legend = T, legend = "right")) +
#   labs(title = paste0("Scenario ",Scenario,": Species recruitment data by generation and replicates.")) +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10),
#           hjust = 0.5
#         ))
# 
#     # ----- #


##### ------------ 2.A.c : % Success Distance ------------ #####

  #### . Plots . ####

# First, we summarize the values across the replicates to get the mean and the sd for each distance. For 1 dist, the mean across all the replicates
Rec_Summarize <- Rec_by_Dist_Total %>% group_by(Dist) %>% 
  summarise_at(colnames(Rec_by_Dist_Total[!(colnames(Rec_by_Dist_Total) %in% c("rep","sample","type","Habitat","SampleType","Dist"))]),list(Mean = mean, Sd = sd), na.rm = TRUE)

# Plot the proportion of recruitment in function of the distance
GG_Rec_5 <- ggplot(data = Rec_Summarize, aes(x = Dist, y = Prop_Rec_Mean, group = 1)) +
  geom_line() +
  geom_point(aes(colour = Prop_Rec_Mean, size = Prop_Rec_Mean)) +
  geom_ribbon(aes(
    y = Prop_Rec_Mean,
    ymin = Prop_Rec_Mean - Prop_Rec_Sd,
    ymax = Prop_Rec_Mean + Prop_Rec_Sd,
    alpha = .2)) +
  geom_text_repel(
    data = . %>% filter(Prop_Rec_Mean > 0),
    aes(label = Prop_Rec_Mean),
    fill = rgb(red = 1, green = 1, blue = 1, alpha = 0.75),
    force_pull   = 0, # do not pull toward data points
    nudge_x      = 2,
    direction    = "y",
    xlim  = c(0.2,NA),           # Define a zone where the labels should be present
    ylim  = c(0.1,NA),
    min.segment.length = 0) +  # Always draw line segments
  ylab(label = "Proportion of Recruitment") +
  xlab(label = "Distance") +
  theme(legend.position = "none") +
  labs(
    title = "Proportion of recruitment by distance. ",
    subtitle = "The values represent the proportion of recruitments between plots in function of the distance.\nThe distance is in plot unit between the plot receiving an individual and the plot giving this individual accross all the meta-community. Values of the proportion of recruitment were averaged across the replicates."
  )


# Plot the proportion of successful recruitment in function of the distance
GG_Rec_6 <- ggplot(data = Rec_Summarize, aes(x = Dist, y = Prop_success_Rec_Mean, group = 1)) +
  geom_line() +
  geom_point(aes(colour = Prop_success_Rec_Mean, size = Prop_success_Rec_Mean)) +
  geom_ribbon(aes(
    y = Prop_success_Rec_Mean,
    ymin = Prop_success_Rec_Mean - Prop_success_Rec_Sd,
    ymax = Prop_success_Rec_Mean + Prop_success_Rec_Sd,
    alpha = .2)) +
  geom_text_repel(
    data = . %>% filter(Prop_success_Rec_Mean > 0),
    aes(label = Prop_success_Rec_Mean),
    fill = rgb(red = 1, green = 1, blue = 1, alpha = 0.75),
    force_pull   = 0, # do not pull toward data points
    nudge_x      = 2,
    direction    = "y",
    xlim  = c(0.2,NA),           # Define a zone where the labels should be present
    ylim  = c(0.1,NA),
    min.segment.length = 0) +  # Always draw line segments
  ylab(label = "Proportion of successful Recruitment") +
  xlab(label = "Distance") +
  theme(legend.position = "none") +
  labs(
    title = "Proportion of successful recruitment by distance. ",
    subtitle = "The values represent the proportion of successful recruitments between plots in function of the distance.\nThe distance is in plot unit between the plot receiving an individual and the plot giving this individual accross all the meta-community. Values of the proportion of recruitment were averaged across the replicates.",
  )

# Arrange plots ---
RecArr2 <- ggarrange(
  GG_Rec_5,GG_Rec_6,nrow = 2, ncol = 1) +
  labs(title = paste0("Scenario ",Scenario,": Species recruitment data by distance.")) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))


# Print a message
cli_alert_success("Recruitment plots correctly done !")

    # ----- #

##### ------------ 2.B: ABUNDANCE -------------------- #####

# Print a message
cli_alert_info("Abundance plots ... ")

    #### . Statistics . ####

# Compute total species abundance per replicates.
Abundance_Total_Summed <- Abundance_Total %>%
  group_by(rep) %>%
  dplyr::summarise_all(sum, na.rm = TRUE) %>% # Compute the sum of each species across all samples for each replicates
  mutate(Total_Sum = rowSums(.[,-c(1:2)]),.before = 3) # Compute the total number of individuals for each replicates

Abundance_Total_Summed_Dplyr<- Abundance_Total_Summed %>%
  pivot_longer(!c(rep,sample,Total_Sum), names_to = "sp", values_to = "Abundance")

# Compute species diversity per replicate
Rep_diversity <- diversity(Abundance_Total_Summed[,-c(1:3)],index = "simpson") %>%  # Compute the diversity
  as.data.frame() %>% # Coerce to data frame
  setNames("Diversity") %>%
  rownames_to_column("rep") %>%
  mutate(across(rep,function(x){factor(as.numeric(x))})) # Change into numeric and factor

# Plot the diversity of the replicates
Plot_diversity <- ggplot(Rep_diversity,aes(x = 0, y = Diversity)) +
  geom_boxplot() +
  geom_jitter(aes(color=rep), size= 1.5, alpha=0.9) +
  scale_x_discrete(limits = 0, labels =  c("0"= " ")) +
  labs(title = "Replicate diversity.",
       subtitle = "Simpson diversity of the replicates.",
       x = NULL)+
  guides(color = FALSE) 

# Compute the sum of each species abundance across each plots for each replicates. (Sum by columns. This is not plot abundance but species abundance accross all the plots)
GG_Abundance_Total <- dplyr::select(Abundance_Total,-sample) %>%    # We don't want the column sample
  group_by(rep) %>%
  group_map( ~ sapply(.x,sum), .keep = F)

# Modify the form to be compatible with GGplot (That's very ugly but it works)
for(i in 1:length(GG_Abundance_Total)){
  GG_Abundance_Total[[i]] <- as.data.frame(GG_Abundance_Total[[i]])
  GG_Abundance_Total[[i]] <- cbind(rownames(GG_Abundance_Total[[i]]),GG_Abundance_Total[[i]],rep(i,nrow(GG_Abundance_Total[[i]])))
  colnames(GG_Abundance_Total[[i]]) <- c("sp","abundance","rep") 
}

# Merging all the list elements on the species columns.
GG_Abundance_Total <- rbind.fill(GG_Abundance_Total) %>%
  mutate(across(c(sp,rep),as.factor))

    # ----- #

  #### . Lollipop plot . ####

# /!\ Replicates Number /!\

# # Lollipop species abundance by replicates
# GG_AB_1 <- ggplot(GG_Abundance_Total, aes(x=species,y=abundance)) +
#   geom_point(data = GG_Abundance_Total %>% filter(abundance > 0),
#              aes(color = factor(rep),size = abundance),
#              na.rm=TRUE) +
#   geom_segment(aes(xend=species, y = 0, yend=abundance), alpha = 0.2,na.rm=TRUE) + 
#   geom_hline(yintercept = 0, alpha = 0.2) +
#   guides(size = "none") +
#   labs(color = "Replicate",
#        title = "Species abundance by replicates.",
#        subtitle = "Lollipop plot of the abundance of each species of the simulated community for each replicates.",
#        x = "Species",
#        y = "Abundance (Number of individuals)") + 
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         panel.grid = element_blank())

  # Species are not ordered by replicates, some are abundant in one replicate and rare in another replicate. 
  # It also gives a coarse information on the type of communities: Are a few species very abundant or many species are common ? 

    # ----- #

  #### . Density plot . ####

  # The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 

# Create a vector to be used as tick marks and reused for the quantile plot. 
tick_breaks <- c(pretty(1:max(GG_Abundance_Total$abundance, na.rm = T), n = 3))
# Names them as their values
names(tick_breaks) <- as.character(tick_breaks)

GG_AB_2 <- ggplot(GG_Abundance_Total, aes(x=abundance)) +
  geom_density(na.rm=FALSE, fill = "grey80", alpha = 0.3) + 
  scale_x_continuous(trans = "log10") +
  guides(color = "none") +
  labs(title = "Species abundance density.",
       subtitle = "Density plots of the species abundance accross all replicates.",
       x = "Abundance (Number of individuals)",
       y = "Density")

# # Summary statistics by replicates ---
# GG_AB_Stat <- GG_Abundance_Total %>%
#   # group_by(rep) %>%
#   get_summary_stats(type = "common",) %>%
#   ggsummarytable(x = "Abundance",
#                  y = c("n","min","max","mean","median","sd"),
#                  ggtheme = theme_bw(),
#                  digits = 2,
#                  size = 2.5) +
#   labs(x = NULL) +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"))

# The most of the species are absent from the communities, few are present and very few are abundant. 

# --------------------------------------------------- #

# Lollipop species abundance total
# GG_AB_3 <- ggplot(GG_Abundance_Total, aes(x=species,y=abundance)) +
#   geom_point(color = "grey30",na.rm=TRUE) +
#   geom_segment(aes(xend=species, y = 0, yend=abundance), alpha = 0.2,na.rm=TRUE) + 
#   labs(color = "Replicate",
#        title = "Species abundance (Total)",
#        x = "Species",
#        y = "Abundance (Number of individuals)") + 
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         plot.title = element_text(size=14, face="italic")) ; GG_AB_3

# --------------------------------------------------- #

    # ----- #

  #### . Density plot with quantiles . ####

# The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 

# Compute the density
dens <- density(GG_Abundance_Total$abundance, na.rm = T)
# Create a dataframe witht the density values
df <- data.frame(x=dens$x, y=dens$y)
# Choose the percentiles we want
probs <- c(0.25, 0.5, 0.75, 0.95, 0.99)
# Compute the quantiles
quantiles <- quantile(GG_Abundance_Total$abundance, prob=probs, na.rm = T)
# Add a column of factors for the quantiles
df$quant <- factor(findInterval(df$x,quantiles))
# Modify tick_breaks
tick_breaks <- c(tick_breaks[!(names(tick_breaks) %in% c("0"))],quantiles[c("75%","95%","99%")])
  
# Plot it
GG_AB_4 <- ggplot(df, aes(x,y)) + 
  geom_line(na.rm=TRUE) + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous(breaks = tick_breaks,
                     limits = c(0,NA),
                     guide = guide_axis(n.dodge=2)) + 
  scale_fill_brewer() +
  guides(fill = "none") +
  labs(color = "Replicate",
       title = "Total species abundance",
       subtitle = "Quantiles of the distribution of species abundance for all replicates.",
       x = "Quantiles of species abundance.",
       y = "Density")


# Summary statistics ---
GG_AB_Stat_2 <- GG_Abundance_Total %>%
  get_summary_stats(type = "common") %>%
  ggsummarytable(x = 1,
                 y = c("n","min","max","mean","median","sd"),
                 ggtheme = theme_bw(),
                 digits = 2,
                 size = 4) +
  labs(x = NULL,
       title = "Species abundance statistics.",
       subtitle = "Statistics of species abundance for all replicates combined.") +
  # theme(axis.text.x= element_blank(),
  #       axis.ticks.x = element_blank()) +
  arrange_theme()

# The most of the species are absent from the communities, few are present and very few are abundant.

    # ----- #

  #### . Heatmap . ####

# Prepare the dataframe to be compatible with pheatmap
GG_HM <- GG_Abundance_Total %>%
  pivot_wider(names_from = rep, values_from = abundance) %>% # Modify the dataframe across replicates
  tibble::column_to_rownames('sp') %>% t()

# Create the palette 
paletteLength <- 50 # Choose the length of the palette (the number of breaks)
# Create the palette
myColor <- viridis(paletteLength) 
# Change the first value (the zero) to white 
myColor[1] <- "#EEEEEE" 
# Create a logarithmically spaced vector of breaks. 
myBreaks <- c(0,pracma::logseq(1,max(GG_HM, na.rm = T),paletteLength))

# Create the heatmap
HM_Plot <- pheatmap(GG_HM,
                    color = myColor,
                    breaks = myBreaks,
                    cluster_rows = FALSE, # Do not cluster rows
                    cluster_cols = FALSE, # Do not cluster cols
                    cellwidth = NA,       # Cell_width and Height (if NA, automatic)
                    cellheight = NA,
                    border_color = "grey40", # Color of the border of the cell
                    scale = "none",          # Should the values be centered (rows, columns or none)
                    show_colnames = FALSE,
                    na_col = "grey60",
                    labels_row = rep("",Nrep),
                    silent = T
) %>% 
  # Modify the theme
  as.ggplot() +
  labs(title = "Species abundance by replicates.",
       subtitle = "Heatmap of the species abundance (columns) by replicates (row). 0 values are in light grey, NA values are in dark grey.",
       caption = "Color is on a logarithmic scale.")+
  arrange_theme()

    # ----- #

  #### . Patchworking . ####

# AbArr <- ggarrange(
#   ggarrange(GG_AB_1,Plot_diversity,widths = c(3,1)),HM_Plot, ncol = 1, align = "v",common.legend = F, legend = "right", heights = c(1,1)) +
#   labs(title = paste0("Scenario ",Scenario,": Species abundance data by replicates.")) +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10),
#           hjust = 0.5
#         ))

# Print a message
cli_alert_success("Abundance plots correctly done !")

# ----- #

##### ------------ 2.C: OCCURENCE -------------------- #####

# Print a message
cli_alert_info("Occurence plots ... ")

  #### . Statistics . ####

# Compute the sum of each species abundance across each plots for each replicates. 
GG_Occurrence_Total <-dplyr::select(Occurrence_Total,-sample)  %>%
  group_by(rep) %>%
  group_map( ~ sapply(.x,sum), .keep = F)

# Modify the form to be compatible with GGplot (That's very ugly but it works)
for(i in 1:length(GG_Occurrence_Total)){
  GG_Occurrence_Total[[i]] <- as.data.frame(GG_Occurrence_Total[[i]])
  GG_Occurrence_Total[[i]] <- cbind(rownames(GG_Occurrence_Total[[i]]),GG_Occurrence_Total[[i]],rep(i,nrow(GG_Occurrence_Total[[i]])))
  colnames(GG_Occurrence_Total[[i]]) <- c("species","abundance","rep") 
}

# Merging all the list elements on the species columns.
GG_Occurrence_Total <- rbind.fill(GG_Occurrence_Total)
# Change the wanted columns into factors
GG_Occurrence_Total[c("species","rep")] <- lapply(GG_Occurrence_Total[c("species","rep")],factor)

    # ----- #

  #### . Lollipop plot . ####

# Lollipop species abundance by replicates
# GG_OC_1 <- ggplot(GG_Occurrence_Total, aes(x=species,y=abundance)) +
#   geom_point(data = GG_Occurrence_Total %>% filter(abundance > 0),
#              aes(color = factor(rep), size = abundance), 
#              na.rm = T) +
#   geom_segment(aes(xend=species, y = 0, yend=abundance), alpha = 0.2, na.rm = T) + 
#   geom_hline(yintercept = 0, alpha = 0.2) +
#   guides(size = "none") +
#   labs(color = "Replicate",
#        title = "Species occurence by replicates.",
#        subtitle = "Lollipop plot of the occurence of each species of the simulated community for each replicates.",
#        x = "Species",
#        y = "Occurence") + 
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         panel.grid = element_blank())

# Species are not ordered by replicates, some are abundant in one replicate and rare in another replicate. 
# It also gives a coarse information on the type of communities: Are a few species very abundant or many species are common ? 

    # ----- #

  #### . Density plot . ####

# The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 

# Create a vector to be used as tick marks and reused for the quantile plot. 
tick_breaks <- c(pretty(1:max(GG_Occurrence_Total$abundance, na.rm = T), n = 5))
# Names them as their values
names(tick_breaks) <- as.character(tick_breaks)

GG_OC_2 <- ggplot(GG_Occurrence_Total, aes(x=abundance)) +
  geom_density(na.rm=FALSE, fill = "grey80", alpha = 0.3) + 
  scale_x_continuous(breaks = tick_breaks) +
  guides(color = "none") +
  labs( title = "Species occurence density by replicates.",
        subtitle = "Density plots of the species occurence by replicates.",
        x = "Occurence",
        y = "Species")

# # Summary statistics by replicates ---
# GG_OC_Stat <- GG_Occurrence_Total %>%
#   group_by(rep) %>%
#   get_summary_stats(type = "common") %>%
#   ggsummarytable(x = "rep",
#                  y = c("n","min","max","mean","median","sd"),
#                  ggtheme = theme_bw(),
#                  digits = 2,
#                  size = 2.5) +
#   labs(x = "Replicates") +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"))

# The most of the species are absent from the communities, few are present and very few are abundant. 

# --------------------------------------------------- #

# Lollipop species abundance total
# GG_OC_3 <- ggplot(GG_Occurrence_Total, aes(x=species,y=abundance)) +
#   geom_point(color = "grey30", na.rm = T) +
#   geom_segment(aes(xend=species, y = 0, yend=abundance), alpha = 0.2, na.rm = T) + 
#   labs(color = "Replicate",
#        title = "Species occurrence (Total).",
#        x = "Species",
#        y = "Occurrence") + 
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         plot.title = element_text(size=14, face="italic"))

# Species are not ordered by replicates, some are abundant in one replicate and rare in another replicate. 
# It also gives a coarse information on the type of communities: Are a few species very abundant or many species are common ? 

    # ----- #

#   #### . Density plot with quantiles . ####
# 
#   # The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 
# 
# # Compute the density
# dens <- density(GG_Occurrence_Total$abundance, na.rm = T)
# # Create a dataframe witht the density values
# df <- data.frame(x=dens$x, y=dens$y)
# # Choose the percentiles we want
# probs <- c(0.25, 0.5, 0.75, 0.95, 0.99)
# # Compute the quantiles
# quantiles <- quantile(GG_Abundance_Total$abundance, prob=probs, na.rm = T)
# # Add a column of factors for the quantiles
# df$quant <- factor(findInterval(df$x,quantiles))
# # Modify tick_breaks
# tick_breaks <- c(tick_breaks,quantiles[c("75%","95%","99%")])
# 
# # Plot it
# GG_OC_4 <- ggplot(df, aes(x,y)) + 
#   geom_line(na.rm = T) + 
#   geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
#   scale_x_continuous(breaks = tick_breaks,
#                      limits = c(0,NA),
#                      guide = guide_axis(n.dodge=1)) + 
#   scale_fill_brewer() +
#   guides(fill = "none") +
#   labs(color = "Replicate",
#        title = "Total species occurence",
#        subtitle = "Quantiles of the distribution of species occurence for all replicates.",
#        x = "Quantiles of species occurence.",
#        y = "Species")

# Summary statistics ---
GG_OC_Stat_2 <- GG_Occurrence_Total %>%
  get_summary_stats(type = "common") %>%
  ggsummarytable(x = 1,
                 y = c("n","min","max","mean","median","sd"),
                 ggtheme = theme_bw(),
                 digits = 2,
                 size = 4) +
  labs(x = NULL,
       title = "Species occurence statistics.",
       subtitle = "Statistics of species occurence for all replicates combined.") +
  theme(axis.text.x= element_blank()) +
  #       axis.ticks.x = element_blank()) +
  arrange_theme()

# The most of the species are absent from the communities, few are present and very few are abundant. 

    # ----- #

  #### . Heatmap . ####

# Prepare the dataframe to be compatible with pheatmap
GG_HM <- GG_Occurrence_Total %>%
  pivot_wider(names_from = rep, values_from = abundance) %>% # Modify the dataframe across replicates
  tibble::column_to_rownames('species') %>% t()

# Create the palette 
paletteLength <- 20 # Choose the length of the palette (the number of breaks)
# Create the palette
myColor <- viridis(paletteLength) 
# Change the first value (the zero) to white 
myColor[1] <- "#EEEEEE" 
# Create a logarithmically spaced vector of breaks. 
# myBreaks <- c(0,pracma::logseq(1,max(GG_HM, na.rm = T),paletteLength))

# Create the heatmap
HM_Plot_2 <- pheatmap(GG_HM,
                    color = myColor,
                    # breaks = myBreaks,
                    cluster_rows = FALSE, # Do not cluster rows
                    cluster_cols = FALSE, # Do not cluster cols
                    cellwidth = NA,       # Cell_width and Height (if NA, automatic)
                    cellheight = NA,
                    border_color = "grey40", # Color of the border of the cell
                    scale = "none",          # Should the values be centered (rows, columns or none)
                    na_col = "grey60",
                    show_colnames = FALSE,
                    labels_row = rep("",Nrep),
                    silent = T
) %>% 
  # Modify the theme
  as.ggplot() +
  labs(title = "Species occurence by replicates.",
       subtitle = "Heatmap of the species occurence (column) by replicates (rows). 0 values are in light grey, NA values are in dark grey.")+
  arrange_theme()

    # ----- #

  #### . Patchworking . ####

# AbOcc <- ggarrange(
#   GG_OC_1,HM_Plot, ncol = 1, align = "v",common.legend = T, legend = "right", heights = c(2,2,1)) +
#   labs(title = paste0("Scenario ",Scenario,": Species occurence data by replicates.")) +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10),
#           hjust = 0.5
#         ))

AbArr2 <- ggarrange(
  GG_AB_2,Plot_diversity,GG_AB_4,GG_AB_Stat_2,GG_OC_2,GG_OC_Stat_2, ncol = 2, nrow = 3, align = "v",common.legend = T, legend = "right", widths  = c(2,1))+
  labs(title = paste0("Scenario ",Scenario,": Species abundance data total.")) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))


# Print a message
cli_alert_success("Occurence plots correctly done !")

# ----- #

##### ------------ 2.D: RICHNESS -------------------- #####

# Print a message
cli_alert_info("Richness plots ... ")

# Group by rep.
GG_SR_Total <- group_by(SR_Total,rep,sample)
# Bind the species richness to the sites data
GG_SR_Total <-merge(GG_SR_Total,Sites_Sampled,sort = FALSE)
# Change sample type data from 1-2 to core-ecotone
GG_SR_Total$SampleType[GG_SR_Total$SampleType == 1] <- "Core"
GG_SR_Total$SampleType[GG_SR_Total$SampleType == 2] <- "Borders"
# Change column to factors
GG_SR_Total <-  mutate(GG_SR_Total, across(c(hab,rep,sample,SampleType),as.factor))

    # ----- #

  #### . Density plot: Richness combined Replicates ####

# The goal of this graph is to look how species abundance is distributed. Are a few species very abundant or are all species evenly distributed and common. 
GG_SR_1 <- ggplot(GG_SR_Total, aes(x=SR_Data)) +
  geom_density(na.rm = T, fill = "grey80", alpha = 0.3) + 
  # scale_color_brewer(palette = "Paired") +
  labs(color = "Replicate",
       x = "Species richness",
       subtitle = "Density of the plot richness for each replicates combined.",
       title = "Species richness.")


GG_SR_Stat_1 <- GG_SR_Total %>%
  # group_by(rep) %>%
  get_summary_stats(SR_Data,type = "common") %>%
  ggsummarytable(x = 1,
                 y = c("n","min","max","mean","median","sd"),
                 ggtheme = theme_bw(),
                 digits = 2,
                 size = 4) +
  labs(x = NULL,
       title = "Species richness statistics.",
       subtitle = "Statistics of plot richness for all replicates combined.") +
  theme(axis.text.x= element_blank()) +
  arrange_theme()
  
    # ----- #
 
  #### . Density plot: Richness ~ SampleTypes. ####
 
# The goal of this graph is to look how species richness is distributed between the two possible sampletypes
GG_SR_2 <- ggplot(GG_SR_Total, aes(x=SR_Data)) +
  geom_density(aes(color = SampleType, fill = SampleType, alpha = 0.3), na.rm = T) +
  scale_color_manual(values = c("#CA5310","#FBBA72")) +
  scale_fill_manual(values = c("#CA5310","#FBBA72")) +
  guides(alpha = "none")+
  labs(color = "SampleType", 
       x = "Species richness",
       subtitle = "Density of the plot richness for each sampletypes.  All replicates are merged together.",
       title = "Species richness by sample types.")


GG_SR_Stat_2 <- GG_SR_Total %>%
  group_by(SampleType) %>%
  get_summary_stats(SR_Data,type = "common") %>%
  ggsummarytable(x = "SampleType",
                 y = c("n","min","max","mean","median","sd"),
                 ggtheme = theme_bw()) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))

    # ----- #

  #### . Density plot: Richness ~ Habitat. ####

# The goal of this graph is to look how species richness is distributed between the two possible 
GG_SR_3 <- ggplot(GG_SR_Total, aes(x=SR_Data)) +
  geom_density(aes(color = hab, fill = hab, alpha = 0.3)) +
  scale_color_manual(values = c("#5FAD56","#8B80F9")) +
  scale_fill_manual(values = c("#5FAD56","#8B80F9")) +
  guides(alpha = "none", fill = "none") +
  labs(color = "Habitat",
       x = "Species richness",
       subtitle = "Density of the plot richness for each habitat.  All replicates are merged together.",
       title = "Species richness by habitat.")

GG_SR_Stat_3 <- GG_SR_Total %>%
  group_by(hab) %>%
  get_summary_stats(SR_Data,type = "common") %>%
  ggsummarytable(x = "hab", y = c("n","min","max","mean","median","sd"),
                 ggtheme = theme_bw()) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98")) +
  labs(x = "Habitat")

    # ----- #

  #### . Patchworking . ####

SR_Density <- ggarrange(GG_SR_1,GG_SR_2,GG_SR_3,
                        ncol = 1,nrow = 3, align = "v",common.legend = F, legend = "right") +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))

SR_Stat <- ggarrange(GG_SR_Stat_1,GG_SR_Stat_2,GG_SR_Stat_3, 
                        ncol = 1,nrow = 3, align = "v",common.legend = F, legend = "right") +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))

  # Total
SR_Total_plot <- ggarrange(SR_Density,SR_Stat,widths = c(2,1)) +
  labs(title = paste0("Scenario ",Scenario,": Species richness data.")) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
          hjust = 0.5
        ))

# Print a message
cli_alert_success("Richness plots correctly done !")

    # ----- #

##### ------------ 2.E: TRAITS -------------------- #####

# Print a message
cli_alert_info("Traits plots ... ")

    # ----- #

  #### . Plots . ####  

# Transformation into a ggplot compatible DF --- #
  # extract rep and sp columns. 
Rep_SP <- Traits_Total[,c("rep","sp")]
  # Remove it from the initial database. 
Traits_Total <- subset(Traits_Total, select = -c(rep,sp))
  # Add N times Rep_SP to be the size of the traits
Rep_SP <- Rep_SP[rep(seq_len(nrow(Rep_SP)),ncol(Traits_Total)), ]
  # Add trait column
Rep_SP <- cbind(Rep_SP, "Trait" = rep(1:ncol(Traits_Total),each = nrow(Traits_Total)))
  # Add each trait value
GG_Traits <- cbind(Rep_SP, "Value" = unlist(Traits_Total))
  # Transform wanted columns into factors
GG_Traits[c("rep","Trait","sp")] <- lapply(GG_Traits[c("rep","Trait","sp")],factor)

# # Density plot by Traits and Replicates---
# GG_Traits_1 <- ggplot(GG_Traits, aes(x=Value)) +
#   geom_density(aes(color = factor(Trait))) + 
#   geom_vline(xintercept = HabOpt) +
#   labs(color = "Trait",
#        x = "Trait value") +
#   facet_grid(rep ~ . , labeller = label_both)
# 
# # Summary statistics ---
# GG_Traits_Stat_1 <- GG_Traits %>%
#   group_by(rep,Trait) %>%
#   get_summary_stats(type = "common") %>% 
#   ggsummarytable(x = "Trait",
#                 y = c("n","min","max","mean","sd"),
#                 size = 2.5,
#                 digits = 4,
#                 color = "Trait", 
#                 facet.by = c("rep",NA),
#                 ggtheme = theme_light()) +
#   guides(color = "none") +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"))

# # Density plot by Replicates ---
# GG_Traits_2 <- ggplot(GG_Traits, aes(x=Value)) +
#   geom_density(aes(color = rep)) + 
#   # scale_color_brewer(palette = "Paired")+
#   geom_vline(xintercept = HabOpt) +
#   labs(color = "Replicates",
#        title = "Total species traits by replicates.",
#        x = "Trait value",
#        subtitle = "Distribution of the combined trait values for each replicates. Black lines represent the trait optimums for each habitat.")
# 
# # Summary statistics ---
# GG_Traits_Stat_2 <- GG_Traits %>%
#   group_by(rep) %>%
#   get_summary_stats(type = "common") %>% 
#   ggsummarytable(x = "rep",
#                  y = c("n","min","max","mean","sd"),
#                  digits = 4,
#                  ggtheme = theme_light()) +
#   labs(x = "Replicate") +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"))

# Density plot total ---
GG_Traits_3 <- ggplot(GG_Traits, aes(x=Value)) +
  geom_density() + 
  geom_vline(xintercept = HabOpt) +
  labs(color = "Replicates",
       title = "Total species traits.",
       x = "Trait value",
       subtitle = "Distribution of the combined trait values of combined replicates. Black lines represent the trait optimums for each habitat.")


# Ungroup GG_Traits to work on all the replicates.
GG_Traits <- ungroup(GG_Traits)
# Summary statistics ---
GG_Traits_Stat_3 <- GG_Traits %>%
  get_summary_stats("Value",type = "common") %>%
  ggsummarytable(x = "Value",
                 y = c("n","min","max","mean","sd"),
                 digits = 4,
                 ggtheme = theme_light()) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98")) +
  labs(x = NULL)

    # ----- #

  #### . Patchworking . ####

# STArr_1 <- ggarrange(
#   GG_Traits_1,GG_Traits_Stat_1, ncol = 2, align = "v",common.legend = T, legend = "right",widths = c(3,2)) +
#   labs(title = paste0("Scenario ",Scenario,": Species traits data (1/2)."),
#        subtitle = "Distribution and values of the nine trait for each replicates. Black lines represent the trait optimums for each habitat.") +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10),
#         ),
#         # Customize subtitle appearence
#         plot.subtitle = element_markdown(
#           color = "grey30", 
#           size = 8,
#           lineheight = 1.35,
#           margin = margin(t = 10, b = 20)
#         ))
#   
# STArr_2 <- ggarrange(
#   GG_Traits_2,GG_Traits_Stat_2,GG_Traits_3,GG_Traits_Stat_3,ncol = 2, nrow = 2, align = "v",common.legend = T, legend = "right",widths = c(3,2)) +
#   labs(title = paste0("Scenario ",Scenario,": Species traits data (2/2).")) +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10),
#         ))

# Print a message
cli_alert_success("Traits plots correctly done !")

# ----- #

##### ------------ 2.F: FITNESS -------------------- #####

# Print a message
cli_alert_info("Fitness plots ... ")

# Add the difference of fitness between the two habitats as a new column
Fitness_Total <- Fitness_Total %>%
  mutate(across(c("rep","sp"),factor)) %>%
  mutate(Abs_Difference = abs(Fitness_Total$habitat_1 - Fitness_Total$habitat_2),.before = sp) 

# Create a pivoted version
Fitness_Pivot <- Fitness_Total %>%
  # Remove the unwanted columns
  select(-Abs_Difference) %>%
  # Pivot the Fitness values for the habitats
  pivot_longer(starts_with("habitat"), names_to = "Habitat", values_to = "Fitness")

    # ----- #

  #### . Density plot: Absolute difference . ####

GG_Fitness_1 <- ggplot(Fitness_Total, aes(x = Abs_Difference)) +
  geom_density(alpha = 0.2, fill = "grey80") +
  # scale_fill_brewer(palette = "Paired") +
  labs(subtitle = "Density plot of the absolute difference of fitness between the two habitats for all species for each replicates.",
       x ="Absolute difference of fitness",
       y = "Density")

GG_Fitness_Stat_1 <- Fitness_Total %>%
  # group_by(rep) %>%
  get_summary_stats(Abs_Difference,type = "common") %>% 
  ggsummarytable(x = "rep",
                 y = c("n","min","max","mean","sd"),
                 digits = 4,
                 ggtheme = theme_bw()) +
  xlab("") +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        axis.text.x = element_blank())

    # ----- #

  #### . Density plot: Fitness ~ Habitat . ####

GG_Fitness_2 <- ggplot(Fitness_Pivot, aes(x = Fitness, color = Habitat)) +
  geom_density(alpha = 0.2, fill = "grey80") +
  # scale_fill_brewer(palette = "Paired") +
  # theme(panel.border = element_rect(colour = "#5FAD56", fill=NA, linewidth =3)) +
  labs(subtitle="Species fitness for each habitat.",
       x ="Species fitness",
       y = "Density")

GG_Fitness_Stat_2 <- Fitness_Pivot %>%
  group_by(Habitat) %>%
  get_summary_stats(Fitness,type = "common") %>% 
  ggsummarytable(x = "Habitat",
                 y = c("n","min","max","mean","sd"),
                 digits = 4,
                 ggtheme = theme_bw()) +
  xlab(" ") +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"))

    # ----- #

  #### . Patchworking . ####

SFArr <- ggarrange(
  GG_Fitness_1,GG_Fitness_Stat_1,GG_Fitness_2,GG_Fitness_Stat_2, 
  ncol = 2, nrow = 2, align = "v",common.legend = T, legend = "right", widths = c(2,1)) +
  labs(title = paste0("Scenario ",Scenario,": Species fitness data.")) +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"),
        # Customize title appearence
        plot.title = element_text(
          color = "grey10", 
          size = 12, 
          face = "bold",
          margin = margin(t = 10),
        ))

# Print a message
cli_alert_success("Fitness plots correctly done !")

# ----- #

# ##### ------------ 2.G: TRAITS ~ FITNESS -------------------- #####
# 
# # Print a message
# cli_alert_info("Traits ~ Fitness plots ... ")
# 
  #### . Statistics . ####
   
# Reload Total Species Traits
Traits_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Traits.csv"),row.names = 1)

# Compute total species abundance per replicates.
Abundance_Total_Summed <- Abundance_Total %>%
  # Change into numeric and factor
  mutate(across(c("rep","sample"),as.factor)) %>%
  # Add the sites data to have information about the habitat
  merge(sites, by = "sample") %>%
  # Group the data
  group_by(rep,hab) %>%
  # Compute the sum of each species across all samples for each replicates, splitted between the two possible habitats
  dplyr::summarise(across(starts_with("t"), \(x) sum(x, na.rm = TRUE))) %>%
  # Pivot the data
  pivot_longer(!c(rep,hab), names_to = "sp", values_to = "Abundance") %>%
  # Add the traits values
  merge(Traits_Total,by = c("sp","rep"))
  
  #### . Abundance ~ Traits ~ HabOpt . ####

# Prepare the datase for the Bubble plot
BubblePlot_Data <- Abundance_Total_Summed %>%
  # Select only one replicate
  filter(rep == 1) %>%
  # Pivot the traits
  pivot_longer(starts_with("Trait"), names_to = "Trait", values_to = "Values")
  
# Vector of names
Names <- c("Habitat 1","Habitat 2") ; names(Names) <- c("1","2")

# Create the Bubble Plot
BubblePlot <- BubblePlot_Data %>%
  ggplot(aes(x=Trait, y=Values, size = Abundance, color = hab)) +
  geom_point(alpha=0.3) +
  facet_grid(. ~ hab, labeller = labeller(hab = Names)) + 
  scale_size_area() +
  guides(color = "none")  +
  labs(title = "Species abundance for each traits values.",
       subtitle = "Bubble plot of the species abundance based on thair trait values for the 9 traits. Only one replicate is shown.") +
  arrange_theme()

# 
# # Compute trait distance to optimal values (0.5 or -0.5) for the first three traits, the 4-6 and the 7-9 traits
# Global_data$Dist_POpt_13 = sqrt((Global_data$Trait_1-0.5)^2+(Global_data$Trait_2-0.5)^2+(Global_data$Trait_3-0.5)^2)
# Global_data$Dist_NOpt_13 = sqrt((Global_data$Trait_1-(-0.5))^2+(Global_data$Trait_2-(-0.5))^2+(Global_data$Trait_3-(-0.5))^2)
# Global_data$Dist_POpt_46 = sqrt((Global_data$Trait_4-0.5)^2+(Global_data$Trait_5-0.5)^2+(Global_data$Trait_6-0.5)^2)
# Global_data$Dist_NOpt_46 = sqrt((Global_data$Trait_4-(-0.5))^2+(Global_data$Trait_5-(-0.5))^2+(Global_data$Trait_6-(-0.5))^2)
# Global_data$Dist_POpt_79 = sqrt((Global_data$Trait_7-0.5)^2+(Global_data$Trait_8-0.5)^2+(Global_data$Trait_9-0.5)^2)
# Global_data$Dist_NOpt_79 = sqrt((Global_data$Trait_7-(-0.5))^2+(Global_data$Trait_8-(-0.5))^2+(Global_data$Trait_9-(-0.5))^2)
# 
# # now add total species abundances per species and replicate in Traits data.frame
# Global_data <- merge(Global_data,Abundance_Total_Summed_Dplyr, by = c("rep","sp"))
# 
# # Compute the correlation between fitness of species for each habitat and the distance of the traits with this habitat
# # Plot_correl <- ggpairs(Global_data,
# #                        columns = 12:19,
# #                        lower = list(continuous = "points", size = 0.2),
# #                        ggplot2::aes(colour=rep),
# #                        columnLabels = c("Hab1_Fitness","Hab2_Fitness","Dist_Opt2_13","Dist_Opt1_13","Dist_Opt2_46","Dist_Opt1_46","Dist_Opt2_79","Dist_Opt1_79")) +
# #   labs(title = paste0("Scenario ",Scenario,": Correlogram."),
# #        subtitle = "Correlogram between Habitat fitnesses and distance from the habitat optimums for each traits. Traits for the same process are combined alltogether.") +
# #   theme(  axis.text.x = element_blank(),
# #           axis.text.y = element_blank(),
# #           axis.ticks = element_blank())
# 
# # GGpairs create a ggmatrix object that is not grobabble, therefore, we will do it by hand ...
# 
# # Create a vector of colnames to compare with the habitat fitness
# Dists <- colnames(Global_data)[14:20]
# 
#     # ----- #
# 
#   #### . Scatterplot: Correlation . ####
# 
# # Create a loop for the habitat 1
# Corel_Plot_hab1 <- foreach (i = 1:length(Dists)) %dopar% {
#   
#   # Create the plot
#   Plot <-
#     ggplot(data = Global_data, aes_string(x = "habitat_1", y = Dists[i], color = "rep")) +
#     geom_point(size = 0.25) +
#     if(i == length(Dists)) {
#       labs(x = "Fitness Habitat 1")
#     } else {
#       labs(x = NULL)
#     }
#   }
# 
# # Create a loop for the habitat 2
# Corel_Plot_hab2 <- foreach (i = 1:length(Dists)) %dopar% {
#   
#   # Create the plot
#   Plot <-
#     ggplot(data = Global_data, aes_string(x = "habitat_2", y = Dists[i], color = "rep")) +
#     geom_point(size = 0.25) +
#     if(i == length(Dists)) {
#       labs(x = "Fitness Habitat 2")
#     } else {
#       labs(x = NULL)
#     }
# }
# 
# 
# # Arrange each plots all together
# Corel_Plot_hab1 <- ggarrange(plotlist = Corel_Plot_hab1, nrow=7,ncol = 1,  align = "v", legend = "none")
# Corel_Plot_hab2 <- ggarrange(plotlist = Corel_Plot_hab2, nrow=7,ncol = 1,  align = "v", common.legend = T, legend = "right")
# 
#     # ----- #
# 
#   #### . Patchworking . ####
# 
# # Arrange them together 
# Core_Plot_Total <-  ggarrange(
#   Corel_Plot_hab1,
#   Corel_Plot_hab2) %>% as.ggplot() +
#   labs(title = paste0("Scenario ",Scenario,": Fitness ~ Distance from habitat optimums."),
#        subtitle = "Correlation between the fitness of the species for the habitats and the distance between their traits and the two habitat optimums.") +
#   theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
#         panel.background = element_rect(fill = "grey98", color = "grey98"),
#         # Customize title appearence
#         plot.title = element_text(
#           color = "grey10", 
#           size = 12, 
#           face = "bold",
#           margin = margin(t = 10)),
#           # Customize subtitle appearence
#           plot.subtitle = element_markdown(
#             color = "grey30", 
#             size = 8,
#             lineheight = 1.35,
#             margin = margin(t = 10, b = 20)
#           ))
# 
# 
# 
# # Print a message
# cli_alert_success("Traits ~ Fitness plots correctly done !")
# 
# # ----- #

##### ------------ 2.H: GRID -------------------- #####

# Print a message
cli_alert_info("Grid plots ... ")
  
# Make the grid 
grid <- ggplot(sites,aes(X,Y)) +
  # geom_point(aes(shape = hab, size = 1)) +
  labs(shape="Habitat") +  
  geom_point(aes(color = sample, size = 0.5, shape = hab, stroke = 1)) + 
  scale_color_discrete(type = sample(viridis(max(as.numeric(sites$sample), na.rm = T))), na.value = "grey85") +
  scale_shape_manual(values = c(15,0)) +
  guides(color = "none", size = "none", shape = "none") +
  theme_void() +
  # labs(title = paste0("Scenario ",Scenario,": Community GRID of plots."),
  #     subtitle = "Community grid of plots of the simulation. Matching colored plots represent sub-samples. Core plots are represented by squares, ecotone plots are represented by dots.") +
  theme(plot.background = element_rect(fill = "grey98", color = "grey98"),
        panel.background = element_rect(fill = "grey98", color = "grey98"))

# Print a message
cli_alert_success("Grid plots correctly done !")

# ----- #

# ##### ------------ 2.I: PHYLOGENETIC DISTANCES -------------------- #####
# 
# # Print a message
# cli_alert_info("Phylogenetic distance plots ... ")
# 
#   #### . Statistics . ####
# 
# # Load the phylogenetic trees
# tree <- ape::read.tree(paste0("coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_PhyloTrees.tree"))
# 
# rep <- 1
# # Plot one of the tree
# # Depending on the size of the phylogeny, it could be very useless
# # ggtree(tree[[1]])
# 
# # Compute the phylogenetic distances between each species
# phyd <- cophenetic.phylo(tree[[rep]])
# 
# # Extract the traits for one replicate
# Tr <- Traits_Total[Traits_Total$rep==rep,]
# rownames(Tr) <- Tr$sp
# 
# # Compute the distance between the fitness
# Fit <- Fitness_Total[Fitness_Total$rep==rep,]
# rownames(Fit) <- Fit$sp
# 
# # Compute the distance between each of the traits
# Dist_9 <- as.matrix(dist(Tr[,3:11])) # dist for the 9 traits
# Dist_13 <- as.matrix(dist(Tr[,3:5])) # dist for the 1st set of 3 traits
# Dist_46 <- as.matrix(dist(Tr[,6:8])) # dist for the 2nd set of 3 traits
# Dist_79 <- as.matrix(dist(Tr[,9:11])) # dist for the 3rd set of 3 traits
# 
# # Compute the distance between each of the habitat fitness
# Dist_hab1 <- as.matrix(dist(Fit$habitat_1))
# Dist_hab2 <- as.matrix(dist(Fit$habitat_2))
# # Rename the row and columns names 
# colnames(Dist_hab1) <- Fit$sp
# rownames(Dist_hab1) <- Fit$sp
# colnames(Dist_hab2) <- Fit$sp
# rownames(Dist_hab2) <- Fit$sp
# 
# # Compute the correlation between phylogenetic distance and trait distance
# Cor_9 <- cor(as.dist(Dist_9),as.dist(phyd))
# Cor_13 <- cor(as.dist(Dist_13),as.dist(phyd))
# Cor_46 <- cor(as.dist(Dist_46),as.dist(phyd))
# Cor_79 <- cor(as.dist(Dist_79),as.dist(phyd))
# 
# # Plot the correlation between phylogenetic distance and trait distance
# # Create a data frame compatible with dplyr
# phyd <- as.data.frame(phyd) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# Dist_9 <- as.data.frame(Dist_9) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# Dist_13 <- as.data.frame(Dist_13) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# Dist_46 <- as.data.frame(Dist_46) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# Dist_79 <- as.data.frame(Dist_79) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# 
# # Create a data frame compatible with dplyr
# Dist_hab1 <- as.data.frame(Dist_hab1) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# Dist_hab2 <- as.data.frame(Dist_hab2) %>% rownames_to_column('row') %>% pivot_longer(cols = -row) %>% setNames(c("sp1","sp2","dist"))
# 
# # Merge the data alltogether
# # Verify that the species are ordered equally in the dataframes
# if (AllIdentical(Dist_9[,1:2],Dist_13[,1:2],Dist_46[,1:2],Dist_79[,1:2]) == TRUE){
#   Global_Dist <- cbind(Dist_9,Dist_13[,3],Dist_46[,3],Dist_79[,3]) %>%  # Bind the values of the trait distances
#     merge(phyd, by = c("sp1","sp2")) %>%
#     setNames(c("sp1","sp2","Dist_9","Dist_13","Dist_46","Dist_79","Dist_Phy"))
# }
# 
# if (AllIdentical(Dist_hab1[,1:2],Dist_hab2[,1:2]) == TRUE){
#   Global_Dist_Fit <- cbind(Dist_hab1,Dist_hab2[,3]) %>%  # Bind the values of the trait distances
#     merge(phyd, by = c("sp1","sp2")) %>%
#     setNames(c("sp1","sp2","Fit_hab1","Fit_hab2","Dist_Phy"))
# }
# 
#     # ----- #
# 
#   #### . Scatterplot: Correlation . ####
# 
# # Plot the correlations between the traits and the phylogenetic distances
# Plot_9 <- ggplot(Global_Dist[sample(1:nrow(Global_Dist), size = 10000),], aes(x = Dist_Phy, y = Dist_9)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Dist_9.",
#        subtitle = "Correlation between phylogenetic distances and all trait distances.",
#        x = "Phylogenetic distances",
#        y = "Trait distances") +
#   guides(color = FALSE)
# 
# Plot_13 <- ggplot(Global_Dist[sample(1:nrow(Global_Dist), size = 10000),], aes(x = Dist_Phy, y = Dist_13)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Dist_13.",
#        subtitle = "Correlation between phylogenetic distances and stabilizing trait distances.",
#        x = "Phylogenetic distances",
#        y = "Trait distances") +
#   guides(color = FALSE)
# 
# Plot_46 <- ggplot(Global_Dist[sample(1:nrow(Global_Dist), size = 10000),], aes(x = Dist_Phy, y = Dist_46)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Dist_46.",
#        subtitle = "Correlation between phylogenetic distances and equalizing trait distances.",
#        x = "Phylogenetic distances",
#        y = "Trait distances") +
#   guides(color = FALSE)
# 
# Plot_79 <- ggplot(Global_Dist[sample(1:nrow(Global_Dist), size = 10000),], aes(x = Dist_Phy, y = Dist_79)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Dist_79.",
#        subtitle = "Correlation between phylogenetic distances and stochastic trait distances.",
#        x = "Phylogenetic distances",
#        y = "Trait distances") +
#   guides(color = FALSE)
# 
# Plot_Fit_1 <- ggplot(Global_Dist_Fit[sample(1:nrow(Global_Dist_Fit), size = 10000),], aes(x = Dist_Phy, y = Fit_hab1)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Fitness habitat 1.",
#        subtitle = "Correlation between phylogenetic distances and the fitness distances of each species for the first habitat.",
#        x = "Phylogenetic distances",
#        y = "Fitness distances") +
#   guides(color = FALSE)
# 
# Plot_Fit_2 <- ggplot(Global_Dist_Fit[sample(1:nrow(Global_Dist_Fit), size = 10000),], aes(x = Dist_Phy, y = Fit_hab2)) +
#   geom_point() +
#   labs(title = "Phylogenetic distances ~ Fitness habitat 2.",
#        subtitle = "Correlation between phylogenetic distances and the fitness distances of each species for the second habitat.",
#        x = "Phylogenetic distances",
#        y = "Fitness distances") +
#   guides(color = FALSE)
# 
#     # ----- #
# 
#   #### . Patchworking . ####
# 
# # Arrange the four precedent plots together (using patchwork)
# Dist_Plots <- (Plot_9|Plot_13) / (Plot_46|Plot_79)
# 
# # Arrange the fitness distance plots
# Dist_Fit_Plots <- Plot_Fit_1 | Plot_Fit_2
# 
# # Transform into a ggplot (whithout creating a pdf file)
# pdf(file = NULL)
# 
# Dist_Plots <- as.ggplot(Dist_Plots) +
#   labs(title = "Phylogenetic distances ~ Trait distances.") +
#   theme(
#     plot.background = element_rect(fill = "grey98", color = "grey98"),
#     # Customize subtitle appearence
#     plot.subtitle = element_text(
#       color = "grey30",
#       size = 10,
#       lineheight = 1.35,
#       margin = margin(t = 5, b = 5),
#       hjust = 0.5
#     ),
#     plot.title = element_text(
#       color = "grey10", 
#       size = 12, 
#       face = "bold",
#       margin = margin(t = 10),
#       hjust = 0.5
#     ))
# 
# # Transform into a ggplot
# Dist_Fit_Plots <- as.ggplot(Dist_Fit_Plots) +
#   labs(title = "Phylogenetic distances ~ Fitness distances.") +
#   theme(
#     plot.background = element_rect(fill = "grey98", color = "grey98"),
#     # Customize subtitle appearence
#     plot.subtitle = element_text(
#       color = "grey30",
#       size = 10,
#       lineheight = 1.35,
#       margin = margin(t = 5, b = 5),
#       hjust = 0.5
#     ),
#     plot.title = element_text(
#       color = "grey10", 
#       size = 12, 
#       face = "bold",
#       margin = margin(t = 10),
#       hjust = 0.5
#     ))
# 
# dev.off()
# 
# # Print a message
# cli_alert_success("Phylogenetic distances plots correctly done !")
# 
# # ----- #

##### ______________________________________________________________________________________________________________________________________ #####

# --------------------------- #
##### STEP 3: GRID REPORT #####
# --------------------------- #

# Print a message
cat(rule(left = "LOADING OF THE DATA FOR ABUNDANCE HEATMAPS", line_col = "green", line = "-", col = "br_green"))

# ----- #

#### . Loading files . ####

# Each time, we load the data and we select only one replicate.

# Total Abundance
Abundance_Total  <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Sampled_Communities/Scenario",Scenario,"_Step4_Total_Sampled_Abundance.csv"),row.names = 1) %>%
  filter(rep == 1)

# Total species fitness
Fitness_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Fitness.csv"),row.names = 1) %>%
  filter(rep == 1)

# Total Sites data 
Sites <- read.csv(paste0("coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step2_Sites_Data.csv"), row.names = 1)  

# Total_Community abundance
Whole_Com <- read.csv(paste0("coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Community_Abundance.csv"), row.names = 1) %>%
  filter(rep == 1) %>%
  # In addition, remove the columns of the species not at all presen in that replicate
  select(where(~sum(!is.na(.)) > 0))

# Total Species Traits
Traits_Total <- read.csv(paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario,"/Whole_Communities/Scenario",Scenario,"_Step1_Total_Species_Traits.csv"),row.names = 1) %>%
  filter(rep == 1)

# ----- #

#### . Preparing files . ####

# Add the fitness to the species traits values
Global_data <- merge(Traits_Total,Fitness_Total, by = c("rep","sp")) %>%
  mutate(across(c(rep,sp),as.factor))

# Compute total species abundance per replicates.
Abundance_Total_Summed <- Abundance_Total %>%
  dplyr::summarise_all(sum, na.rm = TRUE) %>% # Compute the sum of each species across all samples for each replicates
  mutate(Total_Sum = rowSums(.[,-c(1:2)]),.before = 3) # Compute the total number of individuals for each replicates

Abundance_Total_Summed_Dplyr<- Abundance_Total_Summed %>%
  pivot_longer(!c(rep,sample,Total_Sum), names_to = "sp", values_to = "Abundance") %>%
  select(-c(rep,sample,Total_Sum))

# now add total species abundances per species and replicate in Traits data.frame
Global_data <- merge(Global_data,Abundance_Total_Summed_Dplyr, by = c("sp"))

### .. HM Data .. ###

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
  slice_head(n = 10)

# Only keep the species that are the most abundant
HM_Data_Lite <- Whole_Com %>%
  # select by the most abundant species
  select(HM_Data_Best$Species) %>%
  # Add the Sites Data
  cbind(.,Sites) %>%
  # Pivot longer the species
  pivot_longer(
    cols = starts_with("t"), names_to = "Species",values_to = "Abundance") %>%
  mutate_at(vars(!Abundance),factor) 

# ----- # 

# Print a message
cat(rule(left = "CREATION OF THE GRID REPORT", line_col = "green", line = "-", col = "br_green"))

# ----- #

# Loop over the Replicates
Total_Grid_Plots <- foreach(Rep = 1) %do% {
  
  ##### ------------ 3.A: Summary statistics ------------ #####
  
  #### . Species frequencies . ####
  
  Abr = Whole_Com[Whole_Com$rep == Rep,] # Select one replicate
  Abr = Abr[,-(1:2)]/Nip # Compute the relative frequencies of species per plot (Number of individuals per plot multiplied by the number of plots combined to create a sample)
  
  # Verify that everything is OK, if sum of relative abundance accross all species is equal to 1
  # if(rowSums(Abr,na.rm =T) != 1) {stop("The sum of the relative abundance of all species is not equal to 1: Please check the file: Scenario",Scenario,"_Step1_Total_Community_Abundance.csv")}
  
  # Sort the colmeans to get the most abundant species first
  Mean_Abr <- sort(colMeans(Abr), decreasing = T)
  
  # Reorder Abr according to this new vector, from the most abundant to the least abundant specie.
  Abr <- Abr[,match(names(Mean_Abr), colnames(Abr))]
  
  # Modify Mean_Abr to be plotted
  Mean_Abr <- as.data.frame(Mean_Abr) %>%
    rownames_to_column("sp")
  
  # Plot it
  Plot_Mean_Abr <-
    ggplot(data = Mean_Abr[1:20,],aes(x = fct_inorder(sp), y = Mean_Abr)) +
    geom_point() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(title = paste0("Scenario ",Scenario,": Species relative abundance."),
         subtitle = "Relative abundance of the ranked 20 most abundant species in the whole community (grid).",
         x = "Species",
         y = "Relative Abundance")
  
  # ----- #
  
  #### . Species fitness . ####
  
  # Species fitness, ranked by species abundances
  Fitnr = Fitness_Total[Fitness_Total$rep == Rep,]
  rownames(Fitnr) = Fitnr$sp
  # Order the fitness data by the species abundance values
  Fitnr = Fitnr[colnames(Abr),]
  
  # Plot the fitness of the 20 most abundant species for both habitats.
  Plot_Abund_fit <- Fitnr[1:20,] %>%
    ggplot(aes(x = fct_inorder(sp), group = 1)) +
    geom_segment(aes(x = fct_inorder(sp), xend = fct_inorder(sp), y=habitat_1, yend=habitat_2), color="grey") +
    geom_point(aes(y = habitat_1, color= "Habitat 1")) +
    geom_point(aes(y = habitat_2, color= "Habitat 2")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_color_manual(values = c("Habitat 1" = rgb(0.2,0.7,0.1,0.5),"Habitat 2" = rgb(0.7,0.2,0.1,0.5))) +
    labs(title = paste0("Scenario ",Scenario,": Species fitness."),
         subtitle = "Relative fitness of the 20 ranked most abundant species in the whole community for the two habitats (grid).",
         x = "Species",
         y = "Fitness",
         color = "Habitat")
  
  # Split data between replicates to check correlation between abundance and traits
  Rep_X = Global_data[Global_data$rep == Rep,]
  
  # Plot the two habitat fitnesses against each other
  Plot_FitVSFit <- ggplot(Rep_X,aes(x = habitat_1, y = habitat_2)) +
    # geom_point() +
    geom_point(aes(color = abs(habitat_1-habitat_2)>0.50)) +
    scale_colour_manual(name = "Difference of fitness > 0.5", values = setNames(c(rgb(0.2,0.7,0.1,0.5),"grey"), c(T, F))) +
    # geom_smooth() +
    labs(title = paste0("Scenario ",Scenario,": Species fitness."),
         subtitle = "Species fitness for the two habitats.",
         x = "Species Habitat 1 fitness",
         y = "Species Habitat 2 fitness")
  
  # Arrange the three precedent plots together (using patchwork)
  Rep_Plots_1 <- Plot_Mean_Abr / (Plot_FitVSFit|Plot_Abund_fit)
  # Transform them into a ggplot object
  Rep_Plots_1 <- as.ggplot(Rep_Plots_1) +
    labs(title = paste0("Scenario ",Scenario,": Replicate: ",Rep,"/",Nrep)) +
    theme(
      plot.background = element_rect(fill = "grey98", color = "grey98"),
      panel.background = element_rect(fill = "grey98", color = "grey98"),
      # Customize title appearence
      plot.title = element_text(
        color = "grey10", 
        size = 12, 
        face = "bold",
        margin = margin(t = 10),
        hjust = 0.5
      ))
  
  ##### ------------ 3.B: Heatmaps -------------------- #####
  
  # Loop over the 10 most abundant species to create the grid of their repartition
  Total_Grid_Plots <- foreach(sp = 1:10) %dopar% {
    
    # POSITIVE VERSION ------ # 
    # Plot <- ggplot(Sites,aes(x = X, y = Y)) +
    #   geom_point(size = 3*Abr[,sp], shape=1, aes(color = ifelse(Abr[,sp] == 0,'transparent','black'))) +
    #   geom_point(shape = 20, size = 0.1, aes(color = ifelse(hab == 1,'transparent','red'))) +
    #   scale_color_identity() +
    #   geom_rect(aes(xmin = 0, xmax = max(X), ymin = 0, ymax = max(Y)),
    #             fill = "transparent", color = "black") +
    #   labs(title = paste0("Scenario ",Scenario,": Replicate: ",Rep,"/",Nrep),
    #        subtitle = paste0("Species ", colnames(Abr)[sp]," relative abundance: ",Mean_Abr[sp,2]),
    #        x = "",
    #        y = "") +
    #   theme_void() +
    #   theme(
    #   # Customize subtitle appearence
    #     plot.subtitle = element_text(
    #       color = "grey30",
    #       size = 10,
    #       lineheight = 1.35,
    #       margin = margin(t = 5, b = 5),
    #       hjust = 0.5
    #     ),
    #     plot.title = element_text(
    #       color = "grey10", 
    #       size = 12, 
    #       face = "bold",
    #       margin = margin(t = 10),
    #       hjust = 0.5
    #     ))
    
    # NEGATIVE VERSION ------ # 
    
    # Select each time one specie to plot (From the most abundant to the least abundant), 
    # Each time for both habitat
    
    # Select the "sp" rank species for both habitats
    SP <- HM_Data_Best %>%
      group_by(hab) %>%
      slice(sp)
    
    # Select the data on which we are gonna work
    HM_Filter <- HM_Data_Lite %>%
      subset(Species %in% SP$Species)
    
    # Create the labels for the plot
    Facet_Names <- c(paste0("Hab_1 / ",SP$Species[1]," / ",SP$Abundance[1]),paste0("Hab_2 / ",SP$Species[2]," / ",SP$Abundance[2]))
    names(Facet_Names) <- c(SP$Species[1:2])
    
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
      # labs(title = paste0(Scenario.List[[y]])) +
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
    
    # return the plot
    return(HM)
    
  } # End of the loop 
  
  #### . Patchworking . ####
  
  # Append the summary plots
  Total_Grid <- append(Total_Grid_Plots,list(Rep_Plots_1),after = 0)

  # Transform the list into grobs
  Total_Grid <- lapply(Total_Grid, as.grob)
  
  # return the grid
  return(Total_Grid)
  
} # End of the loop across all replicates.

# ---- #

# Unlist the graphs
Total_Grid_Plots <- purrr::flatten(Total_Grid_Plots) 

# Print a message
cli_alert_success("Heatmaps of abundance correctly done !")

# ----- #

# ---------------------------- #
##### STEP 3: PATCHWORKING #####
# ---------------------------- #

##### ______________________________________________________________________________________________________________________________________ #####

# Print a message
cat(rule(left = "PATCHWORKING", line_col = "green", line = "-", col = "br_green"))


# Read the config file
Config <- read.csv(paste0("coMet_ComSim_Outputs/",Scenario,"/Scenario",Scenario,"_Parameters.csv"))
# Change the colnames
colnames(Config) <- c("Parameters","Value")
# Grob the parameters 
Config <-tableGrob(Config)
# Create a plot list
pl <- list(Config,grid,RecArr,RecArr2,AbArr2,HM_Plot,HM_Plot_2,SR_Total_plot,SFArr,BubblePlot) %>%
  append(Total_Grid_Plots)

# Transform them into grobs
pdf(file = NULL)
pl <- lapply(pl, as.grob)
dev.off()
# Create the pdf
ml <- marrangeGrob(pl,nrow=1,ncol = 1)

# Save the pdf
ggsave(filename = paste0("Scenario",Scenario,"_CommunityReport.pdf"),
       plot = ml,
       device = "pdf",
       path = paste0(getwd(),"/coMet_ComSim_Outputs/",Scenario),
       width = 35,
       height = 35,
       units = "cm")

cli_alert_success("FINAL PATCHWORKING: DONE!")

##### ______________________________________________________________________________________________________________________________________ ###

cat(rule(left = "SCRIPT COMET_COMSIM_CHECK.R ENDING", line_col = "red", line = "-", col = "br_red"))

