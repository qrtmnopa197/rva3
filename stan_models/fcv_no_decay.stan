// Predicts probability ratings with a Rescorla-Wagner algorithm, with a single forgetting parameter
// Predicts valence ratings using reward, cue values, and probability rating residuals,
// along with nuisance parameters. Additionally, this model does not allow the effects of these variables to linger
// and decay across trials. Instead, valence ratings are predicted using only the values of these variables on the
// trial of the rating (as in Studies 1 and 2).

data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=2> n_f; // Number of fractals (i.e., cues)
  
  array[n_s,n_t] int frac; // The identity of the fractal presented on each trial

  real out_size; // Absolute value of the binary outcomes

  array[n_s] vector[n_t] rew; // The reward on each trial
  
  array[n_s] vector[n_t] rew_hid; //1 if the reward was hidden, 0 if shown
  
  int n_vrat; // Total number of valence ratings made
  array[n_s,n_t] int val_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_vrat] val_rat; // The valence rating on each trial
  
  int n_prat; // Total number of probability ratings made
  array[n_s,n_t] int prob_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_prat] prob_rat; // The probability rating on each trial
  
  array[n_s] vector[n_t] bl_cent; // Block number for each trial, mean-centered
  array[n_s] vector[n_t] tr_cent; // Trial number within the block, mean-centered
  
  array[n_s] vector[n_t] prev_vrat_cent; // The previous valence rating on each trial, mean-centered
                                         // Used for autoregressive term
}

parameters {
  // Learning rate
  real alpha_mu;
  real<lower=0> alpha_sigma;
  vector[n_s] alpha_z;
  
  // Forgetting rate when the cue isn't presented (Cue Not Shown)
  real forget_mu;
  real<lower=0> forget_sigma;
  vector[n_s] forget_z;
  
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
  
  // Autoregressive term
  real w_auto_mu; 
  real<lower=0> w_auto_sigma;
  vector[n_s] w_auto_z;
  
  // Weights for valence predictors with temporally extended effects on valence
  vector[6] wx_mu;
  vector<lower=0>[6] wx_sigma;
  matrix[n_s,6] wx_z;

  // Residual SDs
  real<lower=0> prob_sigma;
  
  // Mean valence sigma is implicitly 0 for half-normal prior
  real<lower=0> val_sigma_sigma;
  vector<lower=0>[n_s] val_sigma_z;

}

transformed parameters {
  vector[n_prat] prob_pred; // Predicted probability rating
  vector[n_vrat] val_pred; // Predicted valence rating
  vector[n_vrat] val_sigma_ss; // Subject-specific error SDs for valence ratings
  {//anonymous scope start
  
    array[n_s,n_t+1] vector[n_f] V; // EVs of each fractals
    array[n_s,n_t+1] vector[n_f] V_resid; // residual EV for each fractal on each trial
    
    // Valence predictors that have temporally extended effects on valence
    array[n_s] matrix[n_t,6] x_preds = rep_array(rep_matrix(0,n_t,6),n_s); // 1. reward when shown
                                                                           // 2. V when reward shown
                                                                           // 3. V when reward not shown
                                                                           // 4. rew_hid
                                                                           // 5. V_resid when reward shown
                                                                           // 6. V_resid when reward not shown
    x_preds[,,4] = rew_hid;
  
    // Get subject-specific values based on non-centered parameterization
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z);
    vector[n_s] forget = inv_logit(forget_mu + forget_sigma*forget_z); 
    vector[n_s] w_auto = w_auto_mu + w_auto_sigma*w_auto_z;


    vector[n_s] val_sigma = val_sigma_sigma*val_sigma_z;
    
    vector[n_s] w_0 = w_0_mu + w_0_sigma*w_0_z;
    vector[n_s] w_bl = w_bl_mu + w_bl_sigma*w_bl_z;
    vector[n_s] w_tr = w_tr_mu + w_tr_sigma*w_tr_z;
    
    matrix[n_s,6] wx; 
    for(w in 1:6){
      wx[,w] = wx_mu[w] + wx_sigma[w]*wx_z[,w];
    }
    
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          // On the first trial of each subject, initialize all V-values to 0
          V[s,t] = rep_vector(0,n_f); 
          V_resid[s,t] = rep_vector(0,n_f);
        }
        
        V[s,t+1] = V[s,t]*forget[s]; // Decay toward 0 by defafult
        if(rew_hid[s,t] == 0){
          // If the reward was shown...
          V[s,t+1,frac[s,t]] = V[s,t,frac[s,t]] + alpha[s]*(rew[s,t] - V[s,t,frac[s,t]]); // Update V of presented cue
          x_preds[s,t,1] = rew[s,t]; // Assign the reward to the vector of shown rewards
          x_preds[s,t,2] = V[s,t,frac[s,t]]; // Assign V to the vector of Vs when the reward is shown
          x_preds[s,t,5] = V_resid[s,t,frac[s,t]]; // Assign V_resid to the vector of V_resids when the reward is shown
        } else if (rew_hid[s,t] == 1){
          // If the reward was hidden, don't update V
          x_preds[s,t,3] = V[s,t,frac[s,t]]; // Assign V to the vector of Vs when the reward is not shown
          x_preds[s,t,6] = V_resid[s,t,frac[s,t]]; // Assign V_resid to the vector of V_resids when the reward is not shown
        }

        V_resid[s,t+1] = V_resid[s,t]*forget[s]; // Decay toward 0 by default
        // If the participant made a probability rating...
        if(prob_rat_num[s,t] != 0){
          // Add to the vector of predicted Vs - first transforming V to be on the probability scale
          prob_pred[prob_rat_num[s,t]] = 0.5 + V[s,t+1,frac[s,t]]/(2*out_size); 
          // Set the V residual
          V_resid[s,t+1,frac[s,t]] = (2*out_size*prob_rat[prob_rat_num[s,t]] - out_size) - V[s,t+1,frac[s,t]]; 
        } else{
          V_resid[s,t+1,frac[s,t]] = 0; // If no rating, reset the residual
        }
        
        // If the participant made a valence rating... 
        if(val_rat_num[s,t] != 0){
          // Generate valence prediction
          val_pred[val_rat_num[s,t]] = inv_logit(w_0[s] + w_bl[s]*bl_cent[s,t] + w_tr[s]*tr_cent[s,t] + w_auto[s]*prev_vrat_cent[s,t]
                                                 + dot_product(x_preds[s,t],wx[s]));
          val_sigma_ss[val_rat_num[s,t]] = val_sigma[s]; //the subject-specific error SD for this observation
        }
      }
    }
  }//anonymous scope end
}


model {
  alpha_mu ~ normal(-.05,1.7); //approx. equal to uniform distribution when passed through inv_logit
  alpha_sigma ~ normal(0,4);
  alpha_z ~ std_normal();
  
  forget_mu ~ normal(-.05,1.7); //approx. equal to uniform distribution when passed through inv_logit
  forget_sigma ~ normal(0,4);
  forget_z ~ std_normal();
  
  w_0_mu ~ normal(0,2);
  w_0_sigma ~ normal(0,3);
  w_0_z ~ std_normal();
  
  w_bl_mu ~ normal(0,2);
  w_bl_sigma ~ normal(0,3);
  w_bl_z ~ std_normal();
  
  w_tr_mu ~ normal(0,0.04);
  w_tr_sigma ~ normal(0,0.06);
  w_tr_z ~ std_normal();
  
  w_auto_mu ~ normal(0,2); 
  w_auto_sigma ~ normal(0,3);
  w_auto_z ~ std_normal();
  
  wx_mu[{1, 2, 3, 5, 6}] ~ normal(0,0.03);
  wx_sigma[{1, 2, 3, 5, 6}] ~ normal(0,0.05);
  wx_z[,1] ~ std_normal();
  wx_z[,2] ~ std_normal();
  wx_z[,3] ~ std_normal();
  wx_z[,5] ~ std_normal();
  wx_z[,6] ~ std_normal();
  
  wx_mu[4] ~ normal(0,4);
  wx_sigma[4] ~ normal(0,6);
  wx_z[,4] ~ std_normal();
  
  prob_sigma ~ normal(0,1);
  
  val_sigma_sigma ~ normal(0,2);
  val_sigma_z ~ std_normal();
  
  val_rat ~ normal(val_pred, val_sigma_ss) T[0,1];
  prob_rat ~ normal(prob_pred, prob_sigma) T[0,1];
}

generated quantities{
  vector[n_vrat] val_lik; //log likelihoods of all valence ratings
  vector[n_prat] prob_lik; //log likelihoods of all probability ratings
  
  for(i in 1:n_vrat){
    val_lik[i] = normal_lpdf(val_rat[i] | val_pred[i], val_sigma_ss[i]) - 
                 log_diff_exp(normal_lcdf(1 | val_pred[i], val_sigma_ss[i]), 
                              normal_lcdf(0 | val_pred[i], val_sigma_ss[i]));
  }
  for(i in 1:n_prat){
    prob_lik[i] = normal_lpdf(prob_rat[i] | prob_pred[i], prob_sigma) - 
                  log_diff_exp(normal_lcdf(1 | prob_pred[i], prob_sigma), 
                               normal_lcdf(0 | prob_pred[i], prob_sigma));
  }
}