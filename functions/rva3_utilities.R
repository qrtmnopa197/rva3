# Returns a list of data for input to stan, given trial-level data
# trials: trial-level data
# n_t: number of trials; must be set manually
stan_data_rva3 <- function(trials,n_t){
  n_s <- length(unique(trials$id)) # Get number of subjects
  n_f <- max(trials$frac_ix)
  
  out_size <- abs(trials$outcome[1]) # Absolute value of outcomes
  
  rew <- sub_by_trial_vec_list(trials,"outcome") # Outcome on each trial
  rew_hid <- sub_by_trial_vec_list(trials,"rew_hid") 
  
  bl_cent <- sub_by_trial_vec_list(trials,"bl_cent") # Block number, mean-centered
  tr_cent <- sub_by_trial_vec_list(trials,"tr_cent") # Trial number, mean-centered
  tr_sub_cent <- sub_by_trial_vec_list(trials,"tr_sub_cent") # Subject-level trial number, mean-centered
  
  
  frac <- sub_by_trial_matrix(trials,"frac_ix") # Fractal index on each trial

  # Valence rating data
  val_rat_num <- sub_by_trial_matrix(trials,"vrat_number")
  n_vrat <- max(trials$vrat_number)
  
  val_rat_trials <- filter(trials,vrat_number != 0)
  val_rat <- val_rat_trials$val_rat
  val_rat_z <- val_rat_trials$valrat_z
  
  # Probability rating data
  prob_rat_num <- sub_by_trial_matrix(trials,"prat_number")
  n_prat <- max(trials$prat_number)
  
  prob_rat_trials <- filter(trials,prat_number != 0)
  prob_rat <- prob_rat_trials$prob_rat
  
  prev_vrat_cent <- sub_by_trial_vec_list(trials,"prev_vrat_cent") # Previous valence rating

  data <- list(
    n_t = n_t,
    n_s = n_s,
    n_f = n_f,
    out_size = out_size,
    rew = rew,
    rew_hid = rew_hid,
    frac = frac,
    val_rat_num = val_rat_num,
    n_vrat = n_vrat,
    val_rat = val_rat,
    val_rat_z = val_rat_z,
    prob_rat_num = prob_rat_num,
    n_prat = n_prat,
    prob_rat = prob_rat,
    bl_cent = bl_cent,
    tr_cent = tr_cent,
    tr_sub_cent = tr_sub_cent,
    prev_vrat_cent = prev_vrat_cent
  )
  return(data)
}