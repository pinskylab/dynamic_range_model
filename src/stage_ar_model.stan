data {
  real m; // natural mortality
  
  int np; // number of patches
  
  int ns; // number of stages
  
  int ny; // number of years
  
  int n_p_s_y[np,ns,ny]; // array of numbers at patch, stage, and year 
  
  int proj_init[np,ns,1]; // array of initial states for the future projections
  
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
              n_p_s_y[p,1:ns,y] ~ multinomial((to_vector(n_p_s_y_hat[p,1:ns,y])) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]))); // fit to proportions at age by patch AF: this shouldn't have an issue with zeros (?) -- also why not define this as 2:ns?
      } // only evaluate length comps if there are length comps to evaluate
      

      
      n_p_s_y[p,1,y] ~ neg_binomial_2(n_p_s_y_hat[p,1,y], sigma_obs); // fit to mean number of recruits per patch // AF: negative binomial should also be fine with zeros! also, shouldn't sigma_obs be estimated from variance in all the counts, not just smalljuv? there are very few of those and we want a realistic estimate of sigma_obs that we can apply to adults / largejuvs too. maybe we need to move towards estimating sigma_obs by life stage...?

      
    } // close patch loop
    
    
  } // close year loop
  
  
}

generated quantities{
  
real pp_n_p_s_y_hat[np, ns, ny_proj];

real raw_proj[np, ny_proj - 1]; // why are these one year shorter?

real rec_dev_proj[np, ny_proj - 1];

pp_n_p_s_y_hat[1:np, 1:ns, 1] = n_p_s_y_hat[1:np,1:ns,1]; // initialize with starting pop

  // generate posterior predictives for training data
  
for(p in 1:np){
  for(s in 1:ns){
    for(y in 2:ny){
           pp_n_p_s_y_hat[p,s,y] =  multinomial_rng(to_vector(n_p_s_y_hat[p,1:ns,y]) / sum(to_vector(n_p_s_y_hat[p,1:ns,y])));  // the 1:ns syntax makes it a proportion over all stages
           // getting an error here because multinomial requires integers and we made npsyhat a real array

      }
    }
  }
  
  // generate posterior predictives for the testing data

for(p in 1:np){
  for(s in 1:ns){
    pp_proj_n_p_s_y_hat[p,s,1] = proj_init[p,s,1]; // initiate projection with fixed observation
  }
}

// project recruitment deviates into the future
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

// separate loop for projecting pop because we added stage structure 

for(p in 1:np){
  for(n in 1:ns){
    for(y in 2:ny_proj){
      
      tmp = pp_proj_n_p_s_y_hat[p, s, y] * exp(rec_dev_proj[p, y-1]) + 1; // calculate population in each patch and stage
      
      pp_proj_n_p_s_y_hat[p, s, y] = multinomial_rng(to_vector(tmp[p, s, y]) / sum(to_vector(tmp[p, s, y]))); // is it OK that this is also calculating stage 1?
      
    }
  }
}
      
}

