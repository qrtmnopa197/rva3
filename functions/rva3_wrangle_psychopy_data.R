#This function takes a single subject's raw psychopy data and reformats it into a long dataset usable for analysis. 
#It creates all variables of interest that can be created from the raw data alone.
rva3_wrangle_psychopy_data <- function(csv_path){
  if(!is.function(get_csvs_to_analyze)){
    stop("get_csvs_to_analyze() is not loaded into the environment of rva3_wrangle_psychopy_data(). This function is in s22_utilities.R.")
  }
  print(csv_path) #print the csv path for debugging purposes - so you know what subject you're on
  df_full <- read_csv(csv_path, show_col_types = F) #read in subject's data
  #if the subject never pressed the right arrow key on affect probes, the RT columns won't appear. Thus, you need to create fully empty RT columns for them.
  if(!("feed_rate.rt" %in% names(df_full))){
    df_full$feed_rate.rt <- NA
  }
  if(!("worth_rate.rt" %in% names(df_full))){
    df_full$worth_rate.rt <- NA
  }
  #make one df with all the subject-level info, and a second with all the trial-level info
  sub_info <- select(df_full, any_of(c("participant", "date", earnings="earn_str", "total_experiment_time", 
                                       comp_qs = "mouse_3.clicked_name",instruct_keypress = "instruct_key.keys",
                                       feed_check = "pressq_feedrate_2.keys", aff_check = "pressq_feedrate.keys", worth_check = "pressq_feedrate_3.keys"))) 
  trials <- select(df_full, any_of(c("outcome", rew_hid = "reward_reveal", frac_img = "fA_img", val_rat = "val_slider_8.response", 
                                     val_rat_rt = "feed_rate.rt", prob_rat = "worth_slider.response", 
                                     prob_rat_rt = "worth_rate.rt", trial_raw = "trials.thisN", block_raw = "blocks.thisN")),
                   ends_with(".ran"))
  
  #FIRST, REFORMAT THE TRIALS DF...
  
  #Find the rows containing the first and last trials, getting rid of everything above and below those...
  end_instruct_ran <- which(trials$end_instruct_loop.ran == 1) #identify the rows on which end_instruct_loop - the last loop before the main task - ran
  main_task_start <- max(end_instruct_ran) + 1 #the row after the last run of this loop is where the main task starts
  block_runs <- which(trials$blocks.ran == 1) #identify all rows on which the blocks loop ran. 
  main_task_end <- max(block_runs) #the last run of this loop signifies the end of the main task
  trials <- trials[main_task_start:main_task_end,] #keep only rows within the main task
  
  #Create column for block number
  #This is output at the end of each block/repetition loop,
  #so you need to fill out all the rows above each loop end
  #with that value. You can use the custom "fill_vec" function to do this (see s22_functions.R in this folder)
  trials$block_raw <- fill_vec(trials$block_raw, bottom_up = TRUE)
  
  # You want one row per trial, meaning you want to ratchet down one row at the end of every trial. 
  # Fortunately, the end of each trial is marked by a loop end, so psychopy does indeed ratchet down one row
  # at the end of every trial. 
  # However, the block loop sometimes ends in between trials, leading to a redundant
  # ratchet-down. To resolve this issue, you can simply delete all rows on which the blocks loop has run 
  # (since no trial data is collected before the ratchet-down).
  # A second issue is that the aff_rate_loop and worth_rate_loop sometimes end in the middle of trials,
  # so you get a ratchet down before the end of the trial.
  # To address this, you should identify trials on which the aff rate or worth rate loop ran.
  # In these cases, you know that the current row and the row below it represent a single trial’s data, 
  # and that you need to combine them into one row. The simplest way to do this is to copy the trial data
  # on the second row to the first row, thus ensuring that the first row contains the full trial’s data. 
  # Then, delete the second row, which contains only redundant information.
  trials$delete_row <- 0
  for(row in 1:nrow(trials)){
    #On rows where a pre-trial ratchet-down has occurred...
    if (!is.na(trials$blocks.ran[row])){
      trials$delete_row[row] <- 1 #mark row for deletion
    } else if (!(is.na(trials$worth_rate_loop.ran[row])) || !(is.na(trials$aff_rate_loop.ran[row]))){
      #Mid-trial ratchet-down. First, copy all the valuable data from the row below to the row on which the trial started.
      row_below <- row + 1
      trials$trial_raw[row] <- trials$trial_raw[row_below]
      trials$frac_img[row] <- trials$frac_img[row_below]
      #Now, mark the row below (containing only redundant data) for deletion
      trials$delete_row[row_below] <- 1 
    }
  }
  trials <- filter(trials, delete_row == 0) #actually delete the marked rows
  
  #python indexing starts at 0, so this rectifies that
  trials$trial_blk <- trials$trial_raw + 1
  trials$block <- trials$block_raw + 1
  trials <- select(trials,-c(trial_raw,block_raw)) #don't need raw rows anymore
  
  trials$trial_sub <- c(1:nrow(trials)) #get the overall trial numbers for this subject
  
  block_list <- by(trials,trials$block,add_prevrate,rat_col_name="val_rat") #add prev_rate column to the df for each block
  trials <- do.call(rbind,block_list)
  
  trials <- trials %>% mutate(valrat_z = (val_rat - mean(val_rat,na.rm=T))/sd(val_rat,na.rm=T)) #get z-scored valence ratings
  block_list <- by(trials,trials$block,add_prevrate,rat_col_name="valrat_z",pr_col_name="prev_rate_z") #add prev_rate column for valrat_z
  trials <- do.call(rbind,block_list)
  
  #Get fractal numbers
  trials$frac_img <- str_extract(trials$frac_img,"\\d+")
  
  #NOW, THE SUBJECT INFO DF...
  trials_completed <- max(trials$trial_sub)
  
  #get valence rating skipped percent
  valrat_shown <- sum(!is.na(trials$aff_rate_loop.ran))
  valrat_made <- sum(!is.na(trials$val_rat))
  valrat_skipped_percent <- (valrat_shown - valrat_made)/valrat_shown 
  #get probability rating skipped percent
  prob_rat_shown <- sum(!is.na(trials$worth_rate_loop.ran))
  prob_rat_made <- sum(!is.na(trials$prob_rat))
  probrat_skipped_percent <- (prob_rat_shown - prob_rat_made)/prob_rat_shown
  
  trials <- select(trials, -delete_row,-ends_with(".ran")) #get rid of .ran and delete_row columns (don't need anymore)

  earnings <- sub_info$earnings[grepl("\\d+",sub_info$earnings)] #grab earnings from whatever row it's on
  total_experiment_time <- sub_info$total_experiment_time[grepl("\\d+",sub_info$total_experiment_time)] #ditto experiment time
  
  #count the number of each type of attention check was passed
  feed_check_passed <- sum(sub_info$feed_check == "q",na.rm=T)
  aff_check_passed <- sum(sub_info$aff_check == "q",na.rm=T)
  worth_check_passed <- sum(sub_info$worth_check == "q",na.rm=T)
  
  att_checks_passed <- feed_check_passed + aff_check_passed + worth_check_passed #get the total
  
  id <- sub_info$participant[1] #get id from one of the rows
  date <- sub_info$date[1] #get date from one of the rows
  
  #get the standard deviations of the ratings
  sd_valrat <- sd(trials$val_rat,na.rm=TRUE) 
  sd_probrat <- sd(trials$prob_rat,na.rm=TRUE) 
  
  instruct_keypress <- sum(sub_info$instruct_keypress == "right",na.rm=T) #number of manual processions on instructions slides; 
                                                                          #a small number is a red flag they weren't paying attention
  answers_incorrect <- sum(sub_info$comp_qs=="poly_false",na.rm=T)
  
  mean_probrat_rt <- mean(trials$prob_rat_rt, na.rm=TRUE) #grab mean prob rat RT
  mean_valrat_rt <- mean(trials$val_rat_rt, na.rm=TRUE) #get mean RT for valence ratings
  
  valrat_autoprocess <- sum(!is.na(trials$val_rat) & is.na(trials$val_rat_rt)) #number of trials on which a val rating was made but the button wasn't clicked
  probrat_autoprocess <- sum(!is.na(trials$prob_rat) & is.na(trials$prob_rat_rt)) #ditto probability
  
  
  sub_info_final <- data.frame(id,date,
                               answers_incorrect,instruct_keypress,
                               feed_check_passed,worth_check_passed,aff_check_passed,att_checks_passed,
                               probrat_skipped_percent, sd_probrat, mean_probrat_rt, probrat_autoprocess,
                               valrat_skipped_percent,sd_valrat, mean_valrat_rt, valrat_autoprocess,
                               earnings,total_experiment_time,trials_completed,
                               row.names=NULL)
  
  base_data <- cbind(sub_info_final,trials) #merge the two dfs
  
  return(list(base_data,sub_info_final)) #return the trial-level data and subject-level-data
}