data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=2> n_f; // Number of fractals (i.e., cues)

  array[n_s,n_t] int frac; // The identity of the fractal presented on each trial

  real out_size; // Absolute value of the binary outcomes

  array[n_s] vector[n_t] rew; // The reward on each trial

  array[n_s] vector[n_t] rew_hid; // 1 if the reward was hidden, 0 if shown

  int n_vrat; // Total number of valence ratings made
  array[n_s,n_t] int val_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_vrat] val_rat; // The valence rating on each trial

  int n_prat; // Total number of probability ratings made
  array[n_s,n_t] int prob_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_prat] prob_rat; // The probability rating on each trial

  array[n_s] vector[n_t] bl_cent; // Block number for each trial, mean-centered
  array[n_s] vector[n_t] tr_cent; // Trial number within the block, mean-centered
}

transformed data{
  int<lower=1> n_tb = n_t/2; // number of trials per block

  array[n_s] vector[n_t] cue_pr_prev_cent = rep_array(rep_vector(0,n_t),n_s); // Centered previous cue probability
  array[n_s,n_t] int cue_pr_prev_avail = rep_array(0,n_s,n_t); // 1 if previous cue trial had prob rating

  // 1 if valence observation should be included in likelihood, 0 otherwise
  array[n_vrat] int val_keep = rep_array(0,n_vrat);
  int n_vkeep = 0; // Number of valence observations kept in likelihood
  array[n_vrat] int vkeep_idx = rep_array(0,n_vrat); // Indices of kept valence observations

  for (s in 1:n_s) {
    array[n_f] real last_pr = rep_array(0.0,n_f); // Last probability rating for each cue in block
    array[n_f] int last_was_pr = rep_array(0,n_f); // 1 if most recent cue trial had probability rating

    for (t in 1:n_t) {
      int f = frac[s,t]; // Cue shown on current trial

      // Reset cue history at block boundaries
      if (t == 1 || t == (n_tb + 1)) {
        last_pr = rep_array(0.0,n_f);
        last_was_pr = rep_array(0,n_f);
      }

      cue_pr_prev_avail[s,t] = last_was_pr[f];
      if (last_was_pr[f] == 1) {
        cue_pr_prev_cent[s,t] = last_pr[f] - 0.5; // Center at midpoint
      }

      // Keep this valence observation only if:
      // 1) valence was collected on this trial,
      // 2) outcome was shown on this trial, and
      // 3) previous trial for this cue had a probability rating.
      if (val_rat_num[s,t] != 0 && rew_hid[s,t] == 0 && cue_pr_prev_avail[s,t] == 1) {
        val_keep[val_rat_num[s,t]] = 1;
      }

      // Update cue-specific state based on what happened on this cue trial
      if (prob_rat_num[s,t] != 0) {
        last_was_pr[f] = 1;
        last_pr[f] = prob_rat[prob_rat_num[s,t]];
      } else {
        last_was_pr[f] = 0;
        last_pr[f] = 0;
      }
    }
  }
  
  // Build compact index of valence observations included in likelihood.
  for (i in 1:n_vrat) {
    if (val_keep[i] == 1) {
      n_vkeep += 1;
      vkeep_idx[n_vkeep] = i;
    }
  }
}

parameters {
  // Baseline valence
  real w_0_mu;
  real<lower=0> w_0_sigma;
  vector[n_s] w_0_z;

  // Effect of block number (mean-centered) on valence
  real w_bl_mu;
  real<lower=0> w_bl_sigma;
  vector[n_s] w_bl_z;

  // Effect of trial number (mean-centered) on valence
  real w_tr_mu;
  real<lower=0> w_tr_sigma;
  vector[n_s] w_tr_z;

  // Effect of centered cue-linked probability rating (from previous cue trial)
  real w_prob_mu;
  real<lower=0> w_prob_sigma;
  vector[n_s] w_prob_z;

  // Effect of current-trial outcome (points scale)
  real w_out_mu;
  real<lower=0> w_out_sigma;
  vector[n_s] w_out_z;

  // Effect of interaction: centered cue-linked probability * current-trial outcome
  real w_int_mu;
  real<lower=0> w_int_sigma;
  vector[n_s] w_int_z;
  
  // Effect of previous-trial probability-rating priming (centered)
  real w_pr_mu;
  real<lower=0> w_pr_sigma;
  vector[n_s] w_pr_z;
  
  // Decay rate for probability-rating priming history
  real kappa_mu;
  real<lower=0> kappa_sigma;
  vector[n_s] kappa_z;

  // Effect of shown-outcome history (decayed, excludes current trial)
  real w_rhist_mu;
  real<lower=0> w_rhist_sigma;
  vector[n_s] w_rhist_z;

  // Effect of decayed rew_hid term (as in original predictor 4)
  real w_hid_mu;
  real<lower=0> w_hid_sigma;
  vector[n_s] w_hid_z;

  // Decay rate for history terms
  real gamma_mu;
  real<lower=0> gamma_sigma;
  vector[n_s] gamma_z;

  // Mean valence sigma is implicitly 0 for half-normal prior
  real<lower=0> val_sigma_sigma;
  vector<lower=0>[n_s] val_sigma_z;
}

transformed parameters {
  vector[n_vrat] val_pred = rep_vector(0.5,n_vrat); // Predicted valence rating
  vector[n_vrat] val_sigma_ss = rep_vector(1,n_vrat); // Subject-specific valence residual SD
  {//anonymous scope start

    vector[n_s] gamma = inv_logit(gamma_mu + gamma_sigma*gamma_z); // Subject-level history decay
    vector[n_s] kappa = inv_logit(kappa_mu + kappa_sigma*kappa_z); // Subject-level priming decay
    vector[n_s] val_sigma = val_sigma_sigma*val_sigma_z; // Subject-specific valence residual SD

    vector[n_s] w_0 = w_0_mu + w_0_sigma*w_0_z; // Intercept
    vector[n_s] w_bl = w_bl_mu + w_bl_sigma*w_bl_z; // Block trend
    vector[n_s] w_tr = w_tr_mu + w_tr_sigma*w_tr_z; // Trial trend
    vector[n_s] w_prob = w_prob_mu + w_prob_sigma*w_prob_z; // Previous cue-probability effect
    vector[n_s] w_out = w_out_mu + w_out_sigma*w_out_z; // Current shown outcome effect
    vector[n_s] w_int = w_int_mu + w_int_sigma*w_int_z; // Probability-by-outcome interaction
    vector[n_s] w_pr = w_pr_mu + w_pr_sigma*w_pr_z; // Previous-trial probability priming
    vector[n_s] w_rhist = w_rhist_mu + w_rhist_sigma*w_rhist_z; // Shown-outcome history effect
    vector[n_s] w_hid = w_hid_mu + w_hid_sigma*w_hid_z; // Decayed rew_hid effect

    for (s in 1:n_s) {
      real rew_hist = 0; // Decayed shown-outcome history at t (excludes current outcome)
      real hid_hist = 0; // Decayed rew_hid history including current trial rew_hid
      real pr_prime = 0; // Decayed probability-rating priming history (not cue-specific)

      for (t in 1:n_t) {
        // Reset all history states at block boundaries (no cross-block carryover)
        if (t == 1 || t == (n_tb + 1)) {
          rew_hist = 0;
          hid_hist = 0;
          pr_prime = 0;
        }

        // Update rew_hid history before prediction.
        hid_hist = gamma[s]*hid_hist + rew_hid[s,t];

        // Predict only eligible valence observations; excluded rows keep default placeholders.
        if (val_rat_num[s,t] != 0 && rew_hid[s,t] == 0 && cue_pr_prev_avail[s,t] == 1) {
          val_pred[val_rat_num[s,t]] = inv_logit(
            w_0[s] + w_bl[s]*bl_cent[s,t] + w_tr[s]*tr_cent[s,t]
            + w_prob[s]*cue_pr_prev_cent[s,t]
            + w_out[s]*rew[s,t]
            + w_int[s]*(cue_pr_prev_cent[s,t]*rew[s,t])
            + w_pr[s]*pr_prime
            + w_rhist[s]*rew_hist
            + w_hid[s]*hid_hist
          );
          val_sigma_ss[val_rat_num[s,t]] = val_sigma[s];
        }

        // Update shown-outcome history after prediction (strict history term).
        if (rew_hid[s,t] == 0) {
          rew_hist = gamma[s]*rew_hist + rew[s,t];
        } else {
          // Hidden-outcome trials contribute 0 to shown-outcome history.
          rew_hist = gamma[s]*rew_hist;
        }
        
        // Update global probability-rating priming history for trial t+1.
        // Non-probability trials contribute 0 and only decay the existing history.
        if (prob_rat_num[s,t] != 0) {
          pr_prime = kappa[s]*pr_prime + (prob_rat[prob_rat_num[s,t]] - 0.5);
        } else {
          pr_prime = kappa[s]*pr_prime;
        }
      }
    }
  }//anonymous scope end
}

model {
  // Approx. equal to uniform distribution when passed through inv_logit.
  gamma_mu ~ normal(-.05,1.7);
  gamma_sigma ~ normal(0,4);
  gamma_z ~ std_normal();

  w_0_mu ~ normal(0,2);
  w_0_sigma ~ normal(0,3);
  w_0_z ~ std_normal();

  w_bl_mu ~ normal(0,2);
  w_bl_sigma ~ normal(0,3);
  w_bl_z ~ std_normal();

  w_tr_mu ~ normal(0,0.04);
  w_tr_sigma ~ normal(0,0.06);
  w_tr_z ~ std_normal();

  // Similar prior scale as other probability-scale predictors
  w_prob_mu ~ normal(0,6);
  w_prob_sigma ~ normal(0,10);
  w_prob_z ~ std_normal();

  // Outcome-scale effects
  w_out_mu ~ normal(0,0.03);
  w_out_sigma ~ normal(0,0.05);
  w_out_z ~ std_normal();

  w_int_mu ~ normal(0,0.03);
  w_int_sigma ~ normal(0,0.05);
  w_int_z ~ std_normal();
  
  w_pr_mu ~ normal(0,6);
  w_pr_sigma ~ normal(0,10);
  w_pr_z ~ std_normal();
  
  kappa_mu ~ normal(-.05,1.7); // approx. equal to uniform distribution when passed through inv_logit
  kappa_sigma ~ normal(0,4);
  kappa_z ~ std_normal();

  w_rhist_mu ~ normal(0,0.03);
  w_rhist_sigma ~ normal(0,0.05);
  w_rhist_z ~ std_normal();

  // Decayed rew_hid term prior matches the broad prior used for predictor 4 previously
  w_hid_mu ~ normal(0,4);
  w_hid_sigma ~ normal(0,6);
  w_hid_z ~ std_normal();

  val_sigma_sigma ~ normal(0,2);
  val_sigma_z ~ std_normal();

  // Include only eligible valence observations in likelihood.
  if (n_vkeep > 0) {
    for (k in 1:n_vkeep) {
      int i = vkeep_idx[k];
      val_rat[i] ~ normal(val_pred[i], val_sigma_ss[i]) T[0,1];
    }
  }
}

generated quantities {
  // Log likelihoods for included valence observations only.
  vector[n_vkeep] val_lik;

  if (n_vkeep > 0) {
    for (k in 1:n_vkeep) {
      int i = vkeep_idx[k];
      // Truncated-normal log likelihood to match model block.
      val_lik[k] = normal_lpdf(val_rat[i] | val_pred[i], val_sigma_ss[i]) -
                   log_diff_exp(normal_lcdf(1 | val_pred[i], val_sigma_ss[i]),
                                normal_lcdf(0 | val_pred[i], val_sigma_ss[i]));
    }
  }
}
