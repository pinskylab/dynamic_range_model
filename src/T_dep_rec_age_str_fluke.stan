functions {
  
  real T_dep(real sst, real Topt, real width){
    return exp(-0.5 * ((sst - Topt)/width)^2); // gussian temperature-dependent function
  }
  
}

data {
  
  //  real spill; //spillover rate between adjacent patches
  
  int np; // number of patches
  
  int n_ages; // number of ages
  
  int ny; // number of years
  
  int n_p_a_y[np,n_ages,ny]; // array of numbers at patch, stage, and year 
  
  int proj_init[np,n_ages]; // array of initial states for the future projections
  
  int ny_proj; // number of years to forecast 
  
  real sst[np, ny]; // temperature data for training
  
  real sst_proj[np, ny_proj]; // temperature data for testing 
  
  // real g; // do we want to estimate these at all? 
  
  real z;  // total mortality 
  
  real k;
  
  real loo;
  
  real length_50_sel_guess;
  
  // PASTING IN FROM SCROOGE:
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  vector<lower=0>[n_ages] mean_length_at_age;
  
  //vector<lower=0>[n_ages] mean_weight_at_age;
  
  //vector<lower=0, upper=1>[n_ages] mean_maturity_at_age;
  
  matrix[n_ages, n_lbins] length_at_age_key;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  //int age_sel; // estimate of the age at first selectivity
  
  //real h; //steepness
  
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
  
  vector[ny-1] raw; // array of raw recruitment deviates, changed to one value per year
  
  //  real<lower = 1e-3> sigma_obs;
  real<lower=0> phi_obs; // I thought it was possible for this parameter to be negtive, but at one point got this error: Exception: neg_binomial_2_lpmf: Precision parameter is -0.317205, but must be > 0!  (in 'model1ee6785925e9_stage_ar_model_adult_dispersal' at line 145)
  
  real seen_intercept;
  
  real seen_sst;
  
  real sel_50;
  
  real<lower = 0> p_length_50_sel; // length at 50% selectivity
  
  
}

transformed parameters{
  
  real T_adjust[np, ny]; 
  
  real<lower=0> sigma_r;
  
  real length_50_sel;
  
  real sel_delta;
  
  real mean_recruits;
  
  real n_p_a_y_hat [np,n_ages,ny]; // array of numbers at patch, stage, and year 
  
  vector[ny-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr
  
  // vector[n_ages] sel_at_age; // vector of selectivity at stage
  
  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  sel_delta = 2;
  
  length_50_sel = loo * p_length_50_sel;
  
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_p_a_y_hat[1:np,1:n_ages,1] = n_p_a_y[1:np,1:n_ages,1]; // initialize population with "known" values
  
  for(p in 1:np){
    for(y in 1:ny){
      T_adjust[p,y] = T_dep(sst[p,y], Topt, width); // calculate temperature-dependence correction factor for each patch and year depending on SST 
    } // close years
  } // close patches
  
  for (y in 2:ny){
    
    if (y == 2){ 
      rec_dev[y-1]  =  raw[1]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
      
      
    } // close y==2 case  
    else {
      
      // rec_dev[y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[y-2] - -pow(sigma_r,2)/2) + raw[y-1]; // why the double minus signs here? -- now not patch-specific
      
      rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y-1]; // why the double minus signs here? -- now not patch-specific
      
      
      
    } // close ifelse
    
    // density-independent, temperature-dependent recruitment of age 1
    for (p in 1:np){
      
      n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y]; //added in T-dependence here
      
    } // close age 1 case 
    
    // NEED TO THINK ABOUT WHICH AGES CAN DISPERSE--RIGHT NOW, NO DISPERSAL
    // 100% TRANSITION EVERY YEAR
    // ALSO NEED TO UNDERSTAND THE SIZE SELECTIVITY FORMULA
    // HOW TO INCORPORATE MAX AGE? 
    
    for(a in 2:n_ages){
      
      for(p in 1:np){
        
        // note there is no growth rate because 100% of individuals transition (unless they die)
        // is z supposed to be logged here? 
        // still need to add dispersal back in
        n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (1-z);
        
      } // close patches
      
    } // close ages for 2+
    
  } // close year 2+ loop
  
  
} // close transformed parameters block

model {
  
  real theta;
  
  // m ~ normal(0.2, 0.1); 
  
  log_mean_recruits ~ normal(1,5);
  
  Topt ~ normal(18, 4);
  
  width ~ normal(4, 4); 
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  //sigma_obs ~ normal(0.75, 0.25);
  
  phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu); 
  
  seen_intercept ~ normal(0,.5);
  
  seen_sst ~ normal(0,.5);
  
  sel_50 ~ normal(1.3,.1);
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .05);
  
  for(y in 2:ny) {
    
    raw[y-1] ~ normal(0,sigma_r); // prior on raw process error
    
    for(p in 1:np){
      
      theta = exp((seen_intercept + seen_sst * sst[p,y])) / (1 +  exp((seen_intercept + seen_sst * sst[p,y]))); // calculate detection probability
      
      // if (n_p_s_y[p,n_ages,y] > 0) { // if any stage three are observed
      
      if(sum(n_p_a_y[p,sel_100:n_ages,y]) > 0) {
        
        // print("hmm",min(n_p_a_y[p,1:ns,y]));
        // // 
        // print("hmm 2 ", sum(to_vector(n_p_a_y_hat[p,1:ns,y])));
        // 
        
        // if (min(n_p_a_y[p,1:ns,y]) > 0) {
          // n_p_a_y[p,1:ns,y] ~ multinomial((to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage + 1) / sum(to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage + 1));
          
          // }
          
          // n_p_a_y[p,1:ns,y] ~ multinomial(to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage);
          
          
          // print("hat stages are ", (to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage) / sum(to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage));
          
          
          // print("observed stages are ",n_p_a_y[p,1:ns,y]);
          
          n_p_a_y[p,n_ages,y] ~ neg_binomial_2(n_p_a_y_hat[p,n_ages,y] * mean_selectivity_at_age[n_ages] + 1e-3, phi_obs); // fit population scale to number of stage three observed
          
          // print("stage 3 seen is ",n_p_a_y_hat[p,ns,y] * sel_at_stage[ns]);
          
          1 ~ bernoulli(theta);
          
      } else { // only evaluate length comps if there are length comps to evaluate
      
      0 ~ bernoulli(theta);
      
      } // close else 
      
    } // close patch loop
    
    
  } // close year loop
  
  
}

generated quantities{
  
  real pp_proj_n_p_a_y_hat[np, n_ages, ny_proj];
  
  real tmp_proj[np, n_ages, ny_proj];
  
  real rec_dev_proj[ny_proj - 1];
  
  real T_adjust_proj[np, ny_proj];
  
  real pp_theta;
  
  real seen;
  
  for(p in 1:np){
    for(y in 1:ny_proj){
      T_adjust_proj[p,y] = T_dep(sst_proj[p,y], Topt, width); // calculate temperature-dependence correction factor for each patch and year depending on SST 
    }
  }
  
  // generate posterior predictives for the testing data aka project recruitment deviates into the future
  
  // left over from when gA/gY was T-dependent
  // for(p in 1:np){
    //   for(y in 1:ny_proj){
      //     gA_proj[p,y] = growth(sst_proj[p, y], Topt, width);
      //     gY_proj[p,y] = growth(sst_proj[p, y], Topt, width);
      //   }
      // }
      
      for(p in 1:np){
        for(y in 2:ny_proj){
          
          if (y == 2){
            rec_dev_proj[y-1]  = rec_dev[y-1]; // initiate projections with last estimated process errors
          }
          else{
            rec_dev_proj[y-1] = alpha * rec_dev_proj[y-2] + normal_rng(0,sigma_r); // generate autoregressive recruitment deviates
            
          }
        }
      }
      
      
      for(p in 1:np){
        for(a in 1:n_ages){
          pp_proj_n_p_a_y_hat[p,a,1] = proj_init[p,a]; // initiate projection with fixed observation
          tmp_proj[p, a, 1] = proj_init[p,a];
        }
      }
      
      // separate loop for projecting pop because we added stage structure
      
      for(y in 2:ny_proj){
        
        for(p in 1:np){
          
          // project stage 1 (recruitment)
          tmp_proj[p, 1, y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[p, y]; // added in T-dependence here 
          //        print("proj recruits in patch: ",p," year: ",y," is ",tmp_proj[p,1,y]);
          
          
          // age 1+ dispersal
          
          for(a in 2:n_ages){
            
            tmp_proj[p, a, y] = tmp_proj[p, a-1, y-1] * (1-z);
            
            //  print("tmp_proj at age: ",a," patch: ",p," year: ",y," is ",tmp_proj[p,a,y]);
            
          } // close age loop
          
        } // close patch loop
        
        // simulate selectivity and sampling error, at the moment only for stage 3 for weird reasons
        // AF: I just pooled these so that both are added for every age but I'm sure it's wrong somehow 
        
        for(p in 1:np){
          
          pp_theta = exp((seen_intercept + seen_sst * sst_proj[p,y])) / (1 +  exp((seen_intercept + seen_sst * sst_proj[p,y]))); // calculate detection probability
          
          seen = bernoulli_rng(pp_theta);
          
          for (a in 1:n_ages){
            
            pp_proj_n_p_a_y_hat[p, a, y] = seen * neg_binomial_2_rng(tmp_proj[p,a,y] + 1e-3, phi_obs) * mean_selectivity_at_age[a]; // is it OK to have added the 0.001 buffer here?
            
          }
          
        } // close patches for pp_proj
        
        // pp_proj_n_p_a_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
        
      } // close year loop
      
} // close generated quantities


