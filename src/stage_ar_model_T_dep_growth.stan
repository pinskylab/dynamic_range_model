data {
  real m; // natural mortality
  
  real spill; //spillover rate between adjacent patches
  
  int np; // number of patches
  
  int ns; // number of stages
  
  int ny; // number of years
  
  int n_p_s_y[np,ns,ny]; // array of numbers at patch, stage, and year 
  
  int proj_init[np,ns]; // array of initial states for the future projections
  
  int ny_proj; // number of years to forecast 
  
  matrix[np, ny+ny_proj] gY; // growth rate juveniles -> young adults, for both the training and testing sets
  
  matrix[np, ny+ny_proj] gA; // growth rate young adults -> adults
}

transformed data{
  
  
}

parameters{
  
  real log_sigma_r; // sigma recruits
  
  real<lower = -1, upper = 1> alpha; // autocorrelation term
  
  vector[np] log_mean_recruits; // log mean recruits per patch
  
  matrix[np,ny-1] raw; // array of raw recruitment deviates
  
  //  real<lower = 1e-3> sigma_obs;
  real<lower=0> phi_obs; // I thought it was possible for this parameter to be negative, but at one point got this error: Exception: neg_binomial_2_lpmf: Precision parameter is -0.317205, but must be > 0!  (in 'model1ee6785925e9_stage_ar_model_adult_dispersal' at line 145)
  
}

transformed parameters{
  
  real<lower=0> sigma_r;
  
  vector[np] mean_recruits;
  
  real n_p_s_y_hat [np,ns,ny]; // array of numbers at patch, stage, and year 
  
  matrix[np,ny-1] rec_dev; // array of realized recruitment deviates 
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_p_s_y_hat[1:np,1:ns,1] = n_p_s_y[1:np,1:ns,1]; // initialize population with "known" values
  
  
  for (y in 2:ny){
    
    for (p in 1:np){
      
      if (y == 2){
        
        rec_dev[p,y-1]  =  raw[p,1];
        
      } 
      else {
        
        rec_dev[p,y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[p,y-2] - -pow(sigma_r,2)/2) + raw[p,y-1];
        
        
      } // close ifelse
      
      // note life stages are now hard-coded by number, not indexed by s here
      
      n_p_s_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[p,y-1]);
      
      n_p_s_y_hat[p,2,y] = n_p_s_y_hat[p,1,y-1] * (1 - m) * gY[p,y] + n_p_s_y_hat[p,2,y-1] * (1 - m) * (1 - gA[p,y]); // AF: got rid of exp(-m) here 
      
    } // close patches for stages 1 and 2
    
    // adding adult dispersal - edge cases first (reflecting edges)
    
    // this used to be in a for loop over all p but it caused some syntax issues, so it's hard-coded for the edge patches now 
    
    // p = 1 case 
    n_p_s_y_hat[1,3,y] = n_p_s_y_hat[1,2,y-1] * (1 - m) * gA[1,y] + n_p_s_y_hat[1,3,y-1] * (1 - m) * (1 - spill) + n_p_s_y_hat[2, 3, y-1] * (1 - m) * spill; 
    // young adults from last year that grew
    // adults from last year that survived, minus those that dispersed to one patch only
    // dispersal from patch above
    
    // p = np case
    n_p_s_y_hat[np,3,y] = n_p_s_y_hat[np,2,y-1] * (1 - m) * gA[np,y] + n_p_s_y_hat[np,3,y-1] * (1 - m) * (1 - spill) + n_p_s_y_hat[np-1, 3, y-1] * (1 - m) * spill;
    // young adults from last year that grew
    // adults from last year that survived, minus those that dispersed to one patch only
    // dispersal from patch below
    
    // general case for non-edge patches
    for(p in 2:(np-1)){
      n_p_s_y_hat[p,3,y] = n_p_s_y_hat[p,2,y-1] * (1 - m) * gA[p,y] + n_p_s_y_hat[p,3,y-1] * (1 - m) * (1 - 2*spill) + n_p_s_y_hat[p-1, 3, y-1] * (1 - m) * spill + n_p_s_y_hat[p+1, 3, y-1] * (1 - m) * spill; 
      // young adults from last year that grew
      // adults from last year that survived, minus those that dispersed to two patches
      //dispersal from patch below
      // dispersal from patch above
      
    } // close patches from 2 to np-1
    
  } // close years
  
} // close transformed parameters block

model {
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  //sigma_obs ~ normal(0.75, 0.25);
  
  phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu); 
  
  //   
  //   print("n_p_s_y is ", n_p_s_y_hat[12,1:ns,2]);
  //   
  //     print("rec_dev is ", rec_dev[12,1]);
  // 
  //   
  // print("mean recruits are ", mean_recruits[12]);
  
  
  for(y in 2:ny) {
    
    for(p in 1:np){
      
      raw[p,y-1] ~ normal(0,sigma_r); // prior on raw process error
      
      if (sum(n_p_s_y_hat[p,1:ns,y]) > 0){
        
        // print("hmm",min(n_p_s_y[p,1:ns,y]));
        // // 
        // print("hmm 2 ", sum(to_vector(n_p_s_y_hat[p,1:ns,y])));
        // 
        n_p_s_y[p,1:ns,y] ~ multinomial((to_vector(n_p_s_y_hat[p,1:ns,y])) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]))); 
      } // only evaluate length comps if there are length comps to evaluate
      
      n_p_s_y[p,1,y] ~ neg_binomial_2(n_p_s_y_hat[p,1,y], phi_obs); // fit to mean number of recruits per patch AF: changed from sigma_obs to phi_obs
      
      
    } // close patch loop
    
    
  } // close year loop
  
  
}

generated quantities{
  
  real pp_n_p_s_y_hat[np, ns, ny];
  
  real pp_proj_n_p_s_y_hat[np, ns, ny_proj];
  
  real tmp[np, ns, ny];
  
  real tmp_proj[np, ns, ny_proj];
  
  real rec_dev_proj[np, ny_proj - 1];
  
    int y_sum; // because gA and gY are passed in as one matrix each for both the testing and training periods, we need to index along them as a sum of the testing year plus the whole training period
  
  pp_n_p_s_y_hat[1:np, 1:ns, 1] = n_p_s_y[1:np,1:ns,1]; // initialize posterior predictive for training data with real starting pop
  
  tmp[1:np, 1:ns, 1] = n_p_s_y[1:np,1:ns,1]; // initialize posterior predictive for training data with real starting pop
  
  
  // generate posterior predictives for training data
  
  
  // for(p in 1:np){
    //   for(s in 1:ns){
      //     pp_n_p_s_y_hat[p,s,1] = proj_init[p,s]; // initiate projection with fixed observation
      //     tmp[p, s, 1] = proj_init[p,s];
      //   }
      // }
      
      // separate loop for projecting pop because we added stage structure
      
      for(y in 2:ny){
        
        
        for(p in 1:np){
          
          // project stage 1 (recruitment)
          
          tmp[p, 1, y] = neg_binomial_2_rng(n_p_s_y_hat[p,1,y], phi_obs); // observation error for number of recruits, a little weird since outputting an int into a real, but works. AF: changed from sigma_obs to phi_obs
          
          // derterministically fill in the rest of the numbers at stage
          
          tmp[p,2,y] = tmp[p,1,y-1] * (1 - m) * gY[p,y] + tmp[p,2,y-1] * (1 - m) * (1 - gA[p,y]); 
          
        } // close patch loop for stages 1 and 2
        
        // adult dispersal
        
        // p = 1
        tmp[1,3,y] = tmp[1,2,y-1] * (1 - m) * gA[1,y] + // young adults from last year that grew
        tmp[1,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp[1+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        
        // p = np
        
        tmp[np,3,y] = tmp[np,2,y-1] * (1 - m) * gA[np,y] + // young adults from last year that grew
        tmp[np,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp[np-1, 3, y-1] * (1 - m) * spill; // dispersal from patch below
        
        // general case for non-edge patches
        for(p in 2:(np-1)){
          tmp[p,3,y] = tmp[p,2,y-1] * (1 - m) * gA[p,y] + // young adults from last year that grew
          tmp[p,3,y-1] * (1 - m) * (1 - 2*spill) + // adults from last year that survived, minus those that dispersed to two patches
          tmp[p-1, 3, y-1] * (1 - m) * spill + //dispersal from patch below
          tmp[p+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        }
        
        // fill in pp_n_p_s_y_hat
        for(p in 1:np){
          pp_n_p_s_y_hat[p, 1:ns, y] = tmp[p,1:ns,y];
        } // close pp patches
        
      } // close year loop
      
      // fit multinomial to all stages, not really needed except to add noise from sub-sampling
      
      // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
      
      // generate posterior predictives for the testing data aka project recruitment deviates into the future
      
      for(p in 1:np){
        for(y in 2:ny_proj){
          
          if (y == 2){
            rec_dev_proj[p, y-1]  = rec_dev[p, y-1]; // initiate projections with last estimated process errors
          }
          else{
            rec_dev_proj[p, y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev_proj[p, y-2] - -pow(sigma_r,2)/2) + normal_rng(-pow(sigma_r,2)/2,sigma_r); // generate autoregressive recruitment deviates
            
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
        
        y_sum = y + ny; 
        
        for(p in 1:np){
          
          // project stage 1 (recruitment)
          tmp_proj[p, 1, y] = neg_binomial_2_rng(mean_recruits[p] * exp(rec_dev_proj[p, y-1]),phi_obs); // not sure if this is in the same units as n_p_s_y_hat anymore... double-check this! adding in observation error for number of recruits AF: changed from sigma_obs to phi_obs
          
          // project other stages
          
          
          tmp_proj[p,2,y] = tmp_proj[p,1,y-1] * (1 - m) * gY[p,y_sum] + tmp_proj[p,2,y-1] * (1 - m) * (1 - gA[p,y_sum]); 
        } // close patch loop
        
        // adult dispersal
        
        // p = 1
        tmp_proj[1,3,y] = tmp_proj[1,2,y-1] * (1 - m) * gA[1,y_sum] + // young adults from last year that grew
        tmp_proj[1,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp_proj[1+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        
        
        // p = np
        tmp_proj[np,3,y] = tmp_proj[np,2,y-1] * (1 - m) * gA[np,y_sum] + // young adults from last year that grew
        tmp_proj[np,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp_proj[np-1, 3, y-1] * (1 - m) * spill; // dispersal from patch below
        
        
        // general case for non-edge patches
        for(p in 2:(np-1)){
          tmp_proj[p,3,y] = tmp_proj[p,2,y-1] * (1 - m) * gA[p,y_sum] + // young adults from last year that grew
          tmp_proj[p,3,y-1] * (1 - m) * (1 - 2*spill) + // adults from last year that survived, minus those that dispersed to two patches
          tmp_proj[p-1, 3, y-1] * (1 - m) * spill + //dispersal from patch below
          tmp_proj[p+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        }         // close patches 
        
        
        // fit multinomial to all stages: no need for this except to simulate the slight addition of error from subsampling the numbers at age
        
        for(p in 1:np){
          
          pp_proj_n_p_s_y_hat[p, 1:ns, y] = tmp_proj[p,1:ns,y];
          
        } // close patches for pp_proj
        
        // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
        
      } // close year loop
      
} // close generated quantities

