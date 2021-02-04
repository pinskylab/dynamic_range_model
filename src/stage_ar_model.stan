data {
  real m; // natural mortality
  
  int np; // number of patches
  
  int ns; // number of stages
  
  int ny; // number of years
  
  int n_p_s_y[np,ns,ny]; // array of numbers at patch, stage, and year 
  
  int proj_init[np,ns]; // array of initial states for the future projections
  
  int ny_proj; // number of years to forecast 
  
}

transformed data{
  
  
}

parameters{
  
  real log_sigma_r; // sigma recruits
  
  real<lower = -1, upper = 1> alpha; // autocorrelation term
  
  vector[np] log_mean_recruits; // log mean recruits per patch
  
  matrix[np,ny-1] raw; // array of raw recruitment deviates
  
  real<lower = 1e-3> sigma_obs;
  
}

transformed parameters{
  
  real sigma_r;
  
  vector[np] mean_recruits;
  
  real n_p_s_y_hat [np,ns,ny]; // array of numbers at patch, stage, and year 
  
  matrix[np,ny-1] rec_dev; // array of realized recruitment deviates 
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_p_s_y_hat[1:np,1:ns,1] = n_p_s_y[1:np,1:ns,1]; // initialize population with "known" values
  
  for (p in 1:np){
    
    for (y in 2:ny){
      
      if (y == 2){
        
        rec_dev[p,y-1]  =  raw[p,1];
        
      } 
      else {
        
        rec_dev[p,y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[p,y-2] - -pow(sigma_r,2)/2) + raw[p,y-1];
        
        
      } // close ifelse
      
      n_p_s_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[p,y-1]);
      
      for (s in 2:ns){
        
        n_p_s_y_hat[p,s,y] = n_p_s_y_hat[p,s-1,y-1] * exp(-m); // replace this with real growth models eventually 
        
      } // close stages
      
    } // close years
    
  } // close patches
  
  
}

model {
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  sigma_obs ~ normal(0.75, 0.25);
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
        n_p_s_y[p,1:ns,y] ~ multinomial((to_vector(n_p_s_y_hat[p,1:ns,y])) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]))); // fit to proportions at age by patch AF: this shouldn't have an issue with zeros (?) -- also why not define this as 2:ns? also -- this doesn't have an issue with n_p_s_y_hat being real (not int), but multinomial_rng does... 
      } // only evaluate length comps if there are length comps to evaluate
      
      
      
      n_p_s_y[p,1,y] ~ neg_binomial_2(n_p_s_y_hat[p,1,y], sigma_obs); // fit to mean number of recruits per patch // AF: negative binomial should also be fine with zeros! also, shouldn't sigma_obs be estimated from variance in all the counts, not just smalljuv? there are very few of those and we want a realistic estimate of sigma_obs that we can apply to adults / largejuvs too. maybe we need to move towards estimating sigma_obs by life stage...?
      
      
    } // close patch loop
    
    
  } // close year loop
  
  
}

generated quantities{
  
  real pp_n_p_s_y_hat[np, ns, ny];
  
  real pp_proj_n_p_s_y_hat[np, ns, ny_proj]; 
  
  real tmp[np, ns, ny]; 
  
  real tmp_proj[np, ns, ny_proj]; 

  real rec_dev_proj[np, ny_proj - 1];
  
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
  
  for(p in 1:np){
    for(y in 2:ny){
      
      // project stage 1 (recruitment)
      
      tmp[p, 1, y] = neg_binomial_2_rng(n_p_s_y_hat[p,1,y],sigma_obs); // observation error for number of recruits, a little weird since outputting an int into a real, but works. 
      
      for(s in 2:ns){
        
        // project other stages 
        tmp[p, s, y] = tmp[p,s-1,y-1] * exp(-m); // derterministically fill in the rest of the numbers at stage
        
      }
      

      pp_n_p_s_y_hat[p, 1:ns, y] = tmp[p,1:ns,y];

      // fit multinomial to all stages, not really needed except to add noise from sub-sampling

      // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
      
    } // close proj loop
    
  } // close patch loop
 
  
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
  
  for(p in 1:np){
    for(y in 2:ny_proj){
      
      // project stage 1 (recruitment)
      tmp_proj[p, 1, y] = neg_binomial_2_rng(mean_recruits[p] * exp(rec_dev_proj[p, y-1]),sigma_obs); // not sure if this is in the same units as n_p_s_y_hat anymore... double-check this! adding in observation error for number of recruits
      
      for(s in 2:ns){
        
        // project other stages 
        tmp_proj[p, s, y] = tmp_proj[p,s-1,y-1] * exp(-m);
        
      }
      
      // fit multinomial to all stages: no need for this except to simulate the slight addition of error from subsampling the numbers at age

      pp_proj_n_p_s_y_hat[p, 1:ns, y] = tmp_proj[p,1:ns,y];

      // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
      
    } // close proj loop
    
  } // close patch loop
} // close generated quantities

