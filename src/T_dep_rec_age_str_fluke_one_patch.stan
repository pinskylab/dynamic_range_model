functions {
  
  real T_dep(real sst, real Topt, real width){
    return exp(-0.5 * ((sst - Topt)/width)^2); // gussian temperature-dependent function
  }
  
}

data {
  
  //  real spill; //spillover rate between adjacent patches
  
  // int np; // number of patches
  
  int n_ages; // number of ages
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_a_y[n_ages, ny_train];
  
  int proj_init[n_ages];
  
  real sst[ny_train]; // temperature data for training
  
  real sst_proj[ny_proj]; // temperature data for testing 
  
  real z;  // total mortality 
  
  real k;
  
  real loo;
  
  real length_50_sel_guess;
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  vector<lower=0>[n_ages] mean_length_at_age;
  
  matrix[n_ages, n_lbins] length_at_age_key;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  
}

transformed data{
  
  
}

parameters{
  
  real<lower=1e-3> width; // sensitivity to temperature variation
  
  real Topt; //  temp at which recruitment is maximized
  
  real log_sigma_r; // sigma recruits
  
  real<lower = -1, upper = 1> alpha; // autocorrelation term
  
  // real  log_mean_recruits; // log mean recruits per patch, changed to one value for all space/time
  
  real log_mean_recruits;
  
  vector[ny_train-1] raw; // array of raw recruitment deviates, changed to one value per year
  
  //  real<lower = 1e-3> sigma_obs;
  real<lower=0> phi_obs; // I thought it was possible for this parameter to be negtive, but at one point got this error: Exception: neg_binomial_2_lpmf: Precision parameter is -0.317205, but must be > 0!  (in 'model1ee6785925e9_stage_ar_model_adult_dispersal' at line 145)
  
  real<lower = 0> p_length_50_sel; // length at 50% selectivity
  
  real<lower=0, upper=1> theta; // Bernoulli parameter for encounter probability
  
  real theta_threshold; // actual value to trigger positive encounters
}

transformed parameters{
  
  real T_adjust[ny_train]; 
  
  real<lower=0> sigma_r;
  
  real length_50_sel;
  
  real sel_delta;
  
  real mean_recruits;
  
  real n_a_y_hat [n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr
  
  // vector[n_ages] sel_at_age; // vector of selectivity at stage
  
  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  sel_delta = 2;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
  print("mean selectivity at age is ",mean_selectivity_at_age); // check that counts are different for every year
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_a_y_hat[1:n_ages,1] = n_a_y[1:n_ages,1]; // initialize population with "known" values
  
  for(y in 1:ny_train){
    T_adjust[y] = T_dep(sst[y], Topt, width); // calculate temperature-dependence correction factor for each patch and year depending on SST 
  } // close years
  
  for (y in 2:ny_train){
    
    if (y == 2){ 
      rec_dev[y-1]  =  raw[1]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
      
    } // close y==2 case  
    else {
      
      rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y-1]; 
      
      
    } // close ifelse
    
    // density-independent, temperature-dependent recruitment of age 1
    
    n_a_y_hat[1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[y]; //added in T-dependence here
    
    for(a in 2:n_ages){
      
      n_a_y_hat[a,y] = n_a_y_hat[a-1, y-1] * (1-z);
      
    } // close ages for 2+
    
  } // close year 2+ loop
  
  
} // close transformed parameters block

model {
  
  theta ~ uniform(0, 1); // Bernoulli probability of encounter
  
  // RIGHT NOW THIS ISN'T DOING ANYTHING! FIX THAT
  theta_threshold ~ normal(50, 100); // number of total individuals for which encounters are positive
  
  log_mean_recruits ~ normal(1,5);
  
  Topt ~ normal(18, 4);
  
  width ~ normal(4, 4); 
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu); 
  // should we be modeling sigma_obs instead?
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .05);
  
  for(y in 2:ny_train) {
    
    // the n_p_a_y_transform business is here to force Stan to multiply a real (pop size) by a vector (mean sel) to pass to the multinomial
    // it didn't mind this before and I'm not sure what changed!
    vector[n_ages] n_a_y_transform; // should get recreated for every y
    
    for(a in 1:n_ages){
      n_a_y_transform[a] = n_a_y_hat[a,y] .* mean_selectivity_at_age[a];
    }
    
    print("adjusted counts at age in year ",y," are ",n_a_y_transform); // check that counts are different for every year
    
    raw[y-1] ~ normal(0,sigma_r); // prior on raw process error (sorry this is kinda buried in the multinomial stuff)

    // if(sum(n_a_y[sel_100:n_ages,y]) > 0) { // this was still passing some vectors of all zeros to the multinomial, so I replaced  it with:
    // QUESTION: before we were conditioning on the sum true counts being >0. but if the sum *modeled* counts aren't >0, the multinomial doesn't work... so is it OK to condition on those instead?
    
    if(sum(n_a_y_transform[1:n_ages]) > 0) {
      
      // multinomial to estimate relative abundance of stages
      
      n_a_y[1:n_ages,y] ~ multinomial(n_a_y_transform[1:n_ages] / sum(n_a_y_transform[1:n_ages])); // did I parameterize this right?
      
      // negative binomial to estimate absolute abundance (counts) -- was calibrated  to stage 3 before -- now to sum count
      n_a_y[n_ages,y] ~ neg_binomial_2(sum(n_a_y_transform[1:n_ages]) + 1e-3, phi_obs); 
      // FLAG -- is this the right way to recalibrate the magnitude of n_a_y by the absolute scale?
      
      1 ~ bernoulli(theta);
      
      
    } else { // only evaluate length comps if there are length comps to evaluate
    
    0 ~ bernoulli(theta);
    
    } // close else 
    
    
  } // close year loop
  
  
}

generated quantities{
  
  real pp_proj_n_a_y_hat[n_ages, ny_proj];
  
  real tmp_proj[n_ages, ny_proj];
  
  real rec_dev_proj[ny_proj - 1];
  
  real T_adjust_proj[ny_proj];
  
  real pp_theta;
  
  real seen;
  
  for(y in 1:ny_proj){
    T_adjust_proj[y] = T_dep(sst_proj[y], Topt, width); // calculate temperature-dependence correction factor for each patch and year depending on SST 
  }
  
  
  // generate posterior predictives for the testing data aka project recruitment deviates into the future
  
  // left over from when gA/gY was T-dependent
  // for(p in 1:np){
    //   for(y in 1:ny_proj){
      //     gA_proj[p,y] = growth(sst_proj[p, y], Topt, width);
      //     gY_proj[p,y] = growth(sst_proj[p, y], Topt, width);
      //   }
      // }
      
      for(y in 2:ny_proj){
        
        if (y == 2){
          rec_dev_proj[y-1]  = rec_dev[y-1]; // initiate projections with last estimated process errors
        }
        else{
          rec_dev_proj[y-1] = alpha * rec_dev_proj[y-2] + normal_rng(0,sigma_r); // generate autoregressive recruitment deviates
          
        }
      }
      
      
      for(a in 1:n_ages){
        pp_proj_n_a_y_hat[a,1] = proj_init[a]; // initiate projection with fixed observation
        tmp_proj[a, 1] = proj_init[a];
      }
      
      
      // separate loop for projecting pop because we added stage structure
      
      for(y in 2:ny_proj){
        
        // project stage 1 (recruitment)
        tmp_proj[1, y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[y]; // added in T-dependence here 
        //        print("proj recruits in patch: ",p," year: ",y," is ",tmp_proj[p,1,y]);
        
        
        // age 1+ dispersal
        
        for(a in 2:n_ages){
          
          tmp_proj[a, y] = tmp_proj[a-1, y-1] * (1-z);
          
          //  print("tmp_proj at age: ",a," patch: ",p," year: ",y," is ",tmp_proj[p,a,y]);
          
        } // close age loop
        
        // simulate selectivity and sampling error, at the moment only for stage 3 for weird reasons
        // AF: I just pooled these so that both are added for every age but I'm sure it's wrong somehow 
        
        
        pp_theta = theta; // calculate detection probability
        
        seen = bernoulli_rng(pp_theta);
        
        for (a in 1:n_ages){
          
          pp_proj_n_a_y_hat[a, y] = seen * neg_binomial_2_rng(tmp_proj[a,y] + 1e-3, phi_obs) * mean_selectivity_at_age[a]; // is it OK to have added the 0.001 buffer here?
          
        }
        
        // pp_proj_n_p_a_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
        
      } // close year loop
      
} // close generated quantities


