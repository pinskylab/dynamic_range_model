functions {
  
  real T_dep(real sst, real Topt, real width){
    return exp(-0.5 * ((sst - Topt)/width)^2); // gussian temperature-dependent function
  }
  
}

data {
  
  real spill; //spillover rate between adjacent patches
  
  int np; // number of patches
  
  int ns; // number of stages
  
  int ny; // number of years
  
  int n_p_s_y[np,ns,ny]; // array of numbers at patch, stage, and year 
  
  int proj_init[np,ns]; // array of initial states for the future projections
  
  int ny_proj; // number of years to forecast 
  
  real sst[np, ny]; // temperature data for training
  
  real sst_proj[np, ny_proj]; // temperature data for testing 
  
  real g; // do we want to estimate these at all? 
  
  real m;  
  
  vector[ns] mean_length_at_age;
  // vector<lower = 0>[ns] sel_at_stage; // vector of selectivity at stage
  
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
  
}

transformed parameters{
  
  real T_adjust[np, ny]; 
  
  real<lower=0> sigma_r;
  
  // real mean_recruits;
  
  real mean_recruits;

  real n_p_s_y_hat [np,ns,ny]; // array of numbers at patch, stage, and year 
  
  vector[ny-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr
  
  vector[ns] sel_at_stage; // vector of selectivity at stage

  sel_at_stage = 1.0 ./ (1 + exp(-log(19) * ((mean_length_at_age - sel_50) / 1e-1))); // selectivity ogive at age

  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_p_s_y_hat[1:np,1:ns,1] = n_p_s_y[1:np,1:ns,1]; // initialize population with "known" values
  
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
    
    
    for (p in 1:np){
      
      // note life stages are now hard-coded by number, not indexed by s here
      // print("T_adjust is ", T_adjust[p,y]);
      n_p_s_y_hat[p,1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y]; //added in T-dependence here
      
      n_p_s_y_hat[p,2,y] = n_p_s_y_hat[p,1,y-1] * (1 - m) * g + n_p_s_y_hat[p,2,y-1] * (1 - m) * (1 - g); // AF: got rid of exp(-m) here 
      
    } // close patches for stages 1 and 2
    
    // adding adult dispersal - edge cases first (reflecting edges)
    
    // this used to be in a for loop over all p but it caused some syntax issues, so it's hard-coded for the edge patches now 
    
    // p = 1 case 
    n_p_s_y_hat[1,3,y] = n_p_s_y_hat[1,2,y-1] * (1 - m) * g + n_p_s_y_hat[1,3,y-1] * (1 - m) * (1 - spill) + n_p_s_y_hat[2, 3, y-1] * (1 - m) * spill; 
    // young adults from last year that grew
    // adults from last year that survived, minus those that dispersed to one patch only
    // dispersal from patch above
    
    // p = np case
    n_p_s_y_hat[np,3,y] = n_p_s_y_hat[np,2,y-1] * (1 - m) * g + n_p_s_y_hat[np,3,y-1] * (1 - m) * (1 - spill) + n_p_s_y_hat[np-1, 3, y-1] * (1 - m) * spill;
    // young adults from last year that grew
    // adults from last year that survived, minus those that dispersed to one patch only
    // dispersal from patch below
    
    // general case for non-edge patches
    for(p in 2:(np-1)){
      n_p_s_y_hat[p,3,y] = n_p_s_y_hat[p,2,y-1] * (1 - m) * g + n_p_s_y_hat[p,3,y-1] * (1 - m) * (1 - 2*spill) + n_p_s_y_hat[p-1, 3, y-1] * (1 - m) * spill + n_p_s_y_hat[p+1, 3, y-1] * (1 - m) * spill; 
      // young adults from last year that grew
      // adults from last year that survived, minus those that dispersed to two patches
      //dispersal from patch below
      // dispersal from patch above
      
    } // close patches from 2 to np-1
    
  } // close years
  
  
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
  
  //   
  //   print("n_p_s_y is ", n_p_s_y_hat[12,1:ns,2]);
  //   
  //     print("rec_dev is ", rec_dev[12,1]);
  // 
  //   
  // print("mean recruits are ", mean_recruits[12]);
  //  print("Topt is ", Topt);
  //  print("width is ", width);
  
  for(y in 2:ny) {
    
          raw[y-1] ~ normal(0,sigma_r); // prior on raw process error
    
    for(p in 1:np){
      
      theta = exp((seen_intercept + seen_sst * sst[p,y])) / (1 +  exp((seen_intercept + seen_sst * sst[p,y]))); // calculate detection probability
      
      if (n_p_s_y[p,ns,y] > 0) { // if any stage three are observed
        
        // print("hmm",min(n_p_s_y[p,1:ns,y]));
        // // 
        // print("hmm 2 ", sum(to_vector(n_p_s_y_hat[p,1:ns,y])));
        // 
        
        // if (min(n_p_s_y[p,1:ns,y]) > 0) {
        // n_p_s_y[p,1:ns,y] ~ multinomial((to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel_at_stage + 1) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel_at_stage + 1));
        
        // }

      // n_p_s_y[p,1:ns,y] ~ multinomial(to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel_at_stage);

        
              // print("hat stages are ", (to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel_at_stage) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel_at_stage));

       
              // print("observed stages are ",n_p_s_y[p,1:ns,y]);
        
        n_p_s_y[p,ns,y] ~ neg_binomial_2(n_p_s_y_hat[p,ns,y] * sel_at_stage[ns] + 1e-3, phi_obs); // fit population scale to number of stage three observed

      // print("stage 3 seen is ",n_p_s_y_hat[p,ns,y] * sel_at_stage[ns]);

       1 ~ bernoulli(theta);

      } else { // only evaluate length comps if there are length comps to evaluate
      
      0 ~ bernoulli(theta);
      
      } // close else 
      
    } // close patch loop
    
    
  } // close year loop
  
  
}

generated quantities{
  
  real pp_proj_n_p_s_y_hat[np, ns, ny_proj];
  
  real tmp_proj[np, ns, ny_proj];
  
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
            for(s in 1:ns){
              pp_proj_n_p_s_y_hat[p,s,1] = proj_init[p,s]; // initiate projection with fixed observation
              tmp_proj[p, s, 1] = proj_init[p,s];
            }
          }
          
          // separate loop for projecting pop because we added stage structure
          
          for(y in 2:ny_proj){
            
            for(p in 1:np){
              
              // project stage 1 (recruitment)
              tmp_proj[p, 1, y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[p, y]; // added in T-dependence here 
              
              // project other stages
              
              tmp_proj[p,2,y] = tmp_proj[p,1,y-1] * (1 - m) * g + tmp_proj[p,2,y-1] * (1 - m) * (1 - g); 
            } // close patch loop
            
            // adult dispersal
            
            // p = 1
            tmp_proj[1,3,y] = tmp_proj[1,2,y-1] * (1 - m) * g + // young adults from last year that grew
            tmp_proj[1,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
            tmp_proj[1+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
            
            
            // p = np
            tmp_proj[np,3,y] = tmp_proj[np,2,y-1] * (1 - m) * g + // young adults from last year that grew
            tmp_proj[np,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
            tmp_proj[np-1, 3, y-1] * (1 - m) * spill; // dispersal from patch below
            
            
            // general case for non-edge patches
            for(p in 2:(np-1)){
              tmp_proj[p,3,y] = tmp_proj[p,2,y-1] * (1 - m) * g + // young adults from last year that grew
              tmp_proj[p,3,y-1] * (1 - m) * (1 - 2*spill) + // adults from last year that survived, minus those that dispersed to two patches
              tmp_proj[p-1, 3, y-1] * (1 - m) * spill + //dispersal from patch below
              tmp_proj[p+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
            }         // close patches 
            
            
            // simulate selectivity and sampling error, at the moment only for stage 3 for weird reasons
            
            for(p in 1:np){
              
          pp_theta = exp((seen_intercept + seen_sst * sst_proj[p,y])) / (1 +  exp((seen_intercept + seen_sst * sst_proj[p,y]))); // calculate detection probability
              
             seen = bernoulli_rng(pp_theta);
              
              for (s in 1:ns){
                
              pp_proj_n_p_s_y_hat[p, s, y] = seen * tmp_proj[p,s,y] * sel_at_stage[s];
              
              }
              
              pp_proj_n_p_s_y_hat[p,ns,y] = seen * neg_binomial_2_rng(tmp_proj[p,ns,y] + 1e-3, phi_obs); // add in observation error
              
            } // close patches for pp_proj
            
            // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
            
          } // close year loop
          
} // close generated quantities


