#This is the master script for initial data analysis steps
#It identifies data to use, wrangles it into a usable form, creates additional variables based on this data, and creates plots of certain variables for quality-checking.

library(qualtRics)
library(tidyverse)
library(ddpcr)
library(tidyverse)
library(sigmoid)

##SET MANUALLY
path_to_project_directory <- "~/projects/RVA_3/"
path_to_s22 <- "~/projects/spring_2022_study/"
path_to_s22fu <- "~/projects/s22_follow_up/"
path_to_rva <- "~/projects/RVA/"

ids_to_exclude <- c("6398dc9eee4c333e1712e777") #Ps whose data you don't want to analyze even if it looks good
##############

#clear out results from old analyses
system("mv /Users/dp/projects/RVA_3/analysis_data/*.csv /Users/dp/projects/RVA_3/analysis_data/old_analysis_data") #move any CSV files in this folder to old_analysis_data folder
system("rm -rf /Users/dp/projects/RVA_3/analysis_data/qc_plots/trial_level_vars/*") #clear the folder containing trial-level plots
system("rm -f /Users/dp/projects/RVA_3/analysis_data/qc_plots/sub_level_variables.pdf") #remove old subject level variable plot

#read in necessary functions
source(paste0(path_to_s22,"code/functions/s22_utilities.R"))
source(paste0(path_to_s22fu,"code/functions/s22fu_utilities.R"))
source(paste0(path_to_rva,"code/functions/rva_utilities.R"))
source(paste0(path_to_project_directory,"code/functions/rva3_wrangle_psychopy_data.R"))

csvs_to_analyze <- get_csvs_to_analyze(path_to_project_directory,ids_to_exclude,first_time="2025-03-05") #get the usable CSVs to analyze
all_data <- lapply(csvs_to_analyze, rva3_wrangle_psychopy_data) #reformats each CSV, turning it into a long dataset usable for analysis, and adds all variables of interest that can be created from the raw data alone.
                                                               #Returns a list - one element for each subject - where each element is itself a list containing dfs with the trial-level data and subject-level data

trials_list <- lapply(all_data, function(l) l[[1]]) #get a list of the trial-level dfs only
trial_level_data <- do.call(rbind,trials_list) #stack trial dfs into one big df
sub_list <- lapply(all_data, function(l) l[[2]]) #get a list of the subject-level dfs only
sub_level_data <- do.call(rbind,sub_list) #stack them into one big data frame

#append the Yes/No answer from qualtrics to this data
sub_level_data <- qualtrics_append(sub_level_data,survey="/Users/dp/Downloads/RVA_3_survey.csv",id_col_qualt="participant",join_cols_qualt=c("Q46","Q40","Q41","Q37","Q25"),join_cols_df=c("sig_emot","age","gender","race","ethnicity"))
trial_level_data <- qualtrics_append(trial_level_data,survey="/Users/dp/Downloads/RVA_3_survey.csv",id_col_qualt="participant",join_cols_qualt=c("Q46","Q40","Q41","Q37","Q25"),join_cols_df=c("sig_emot","age","gender","race","ethnicity"))
#write both to CSVs
date_time <- Sys.time() %>% chartr(" ","_",.) %>% chartr(":","_",.) #grab the date and time, reformatting ':' and '-' to  '_' so you can label the files with it
write.csv(trial_level_data,paste0(path_to_project_directory,"analysis_data/trial_level_data_all_subs_",date_time,".csv"),row.names = FALSE)
#write subject-level data
write.csv(sub_level_data,paste0(path_to_project_directory,"analysis_data/sub_level_data_all_subs_",date_time,".csv"),row.names = FALSE)

#Create plot grids of select variables for quality checking. These are saved to the qc_plots folder
tl_hists <- c("val_rat_rt","prob_rat_rt","val_rat","prob_rat") #trial level variables to plot
plot_trial_level_vars(trial_level_data,tl_hists,path_to_project_directory) #create and save plot grids
sl_hists <- c("answers_incorrect","instruct_keypress",
              "feed_check_passed","aff_check_passed","worth_check_passed","att_checks_passed",
              "probrat_skipped_percent","sd_probrat","mean_probrat_rt",
              "valrat_skipped_percent","sd_valrat","mean_valrat_rt",
              "earnings","total_experiment_time") #subject-level variables to plot
plot_sub_level_vars(sub_level_data,sl_hists,path_to_project_directory) #create and save plot grids




