functions {
  
  real T_dep(real sbt, real Topt, real width){
    return exp(-0.5 * ((sbt - Topt)/width)^2); // gussian temperature-dependent function
  }
  
} // close functions block



data {
  
  // survey data 

  int n_ages; // number of ages
  
  int np; // number of patches
  
  int ny_train; // years for training
  
  // int ny_proj; // number of years to forecast 
  
    int n_lbins; // number of length bins (here just the range of cm values)
  
  int n_p_l_y[np, n_lbins, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
    // int n_p_a_y[np, n_ages, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort

  
  matrix[n_ages, n_lbins] l_at_a_key;
 
  real abund_p_y[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  // matrix[np, ny_train] c_p_y;

  
  // environmental data 
  
  real sbt[np, ny_train]; // temperature data for training
  
  // fish data
  
  real m;  // total mortality 
  
  real k;
  
  real loo;
  
  real t0;
  
  real cv;
  
  real length_50_sel_guess;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  vector[np] patcharea;
  
  int<lower = 0, upper = 1> do_dirichlet;
  
}

transformed data{
  
  // matrix[n_lbins, n_ages] prob_age_at_length;
  
  // int n_p_a_y[np, n_ages, ny_train]; 
  
  // vector[n_ages] age_dist[n_lbins]; // create array of vectors 
  
  // prob_age_at_length = age_at_length_key(
  //   loo=loo,
  //   l0=l0,
  //   k=k,
  //   cv=cv,
  //   n_lbins=n_lbins,
  //   n_ages=n_ages ); 
  
  // for(p in 1:np){
  //   for(l in 1:n_lbins){
  //     for(y in 1:ny_train){
  //       age_dist[l] = prob_age_at_length[l,] * n_p_l_y[p, l, y]; // get vector of probabilities for each length, times counts
  //     }
  //   }
  // }
  
  // for(p in 1:np){
  //   for(a in 1:n_ages){
  //     for(y in 1:ny_train){
  //       n_p_a_y[p,a,y] = sum(age_dist[,a]); // not sure I'm indexing age_dist correctly, and still need to round to int somehow!
  //     }
  //   }
  // }
}

parameters{
  
  
  real<lower = 1e-3> sigma_r;

  real<lower = 1e-3> sigma_obs;

  // real<lower = 1e-3,upper = 10> sigma;
  // 
  // real<lower = 0, upper = 1>  proc_ratio; // sigma recruits
  // 
  vector[np] log_mean_recruits;
  
  vector[ny_train] raw; // array of raw recruitment deviates, changed to one value per year
  
  real p_length_50_sel; // length at 50% selectivity
  
  // real log_scalar;
  
  real log_f;
  
  real <lower = 0> theta_d;
  
}

transformed parameters{
  
  // real sigma_r;
  
  // real sigma_obs;
  
  real length_50_sel;
  
  real sel_delta;
  
  vector[np] mean_recruits;
  
  real n_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  matrix[np, n_lbins] n_p_l_y_hat[ny_train]; // array number of years containing matrices with numbers at patch, length bin, and year 

  real dens_p_y_hat [np, ny_train]; // for tracking sum density 
  
  real n_eff [np, ny_train]; // for tracking sum density 

  real n_raw [np, ny_train]; // for tracking sum density 


  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  real f = exp(log_f);
  
  // real scalar = exp(log_scalar);
  
  real d = 0;
  
  real alpha = 0;
  
  real z = exp(-(f + m));
  
  real c_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year

  real c_p_y_hat [np, ny_train]; // for tracking sum density 

  c_p_a_y_hat = rep_array(0, np, n_ages,ny_train);

  // sigma_r = sigma * proc_ratio;
  
  // sigma_obs = sigma * (1 - proc_ratio);
  
  sel_delta = 1;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  
  mean_recruits = exp(log_mean_recruits);
  
  // fill in year 1 of n_p_a_y_hat, initialized with mean_recruits 
  for(p in 1:np){
    for(a in 1:n_ages){
      if(a==1){
        n_p_a_y_hat[p,a,1] = mean_recruits[p] * exp(raw[1] - pow(sigma_r,2)/2); // initialize age 0 with mean recruitment in every patch
      }
      else{
        n_p_a_y_hat[p,a,1] = n_p_a_y_hat[p,a-1,1] * (z); // initialize population with mean recruitment propogated through age classes with mortality
      }
      
    } // close ages
  } // close patches
  
  // rec_dev[1] = raw[1];
  
  // calculate recruitment deviates every year (not patch-specific)
  for (y in 2:ny_train){
    if (y == 2){
      rec_dev[y-1]  =  raw[y]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific

    } // close y==2 case
    else {
      
      rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y]; 
      
    } // close ifelse
    
    // describe population dynamics
    for(p in 1:np){
      
      n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[y-1] - pow(sigma_r,2)/2);
      
      // pop dy for non-reproductive ages 
      for(a in 2:(age_at_maturity-1)){
        
        n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (z); // these just grow and die in the patch
        c_p_a_y_hat[p,a - 1,y - 1] =  n_p_a_y_hat[p, a-1, y-1] * (z) * (f / (f + m));

        
      } // close ages for 2 to age at maturity

      // pop dy for reproductive adults
      // mortality and dispersal are happening simultaneously here, between generations
      // because neither is patch-specific I don't think the order matters
      
      for(a in age_at_maturity:n_ages){
        
        if (a == n_ages){
                   
            n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (z) + n_p_a_y_hat[p, a, y-1] * z;

        } else {
            
            n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (z);

        }
          
          c_p_a_y_hat[p,a - 1,y - 1] =  n_p_a_y_hat[p, a-1, y-1] * (z) * (f / (f + m));

      }// close ages
      
    } // close patches 
    
    
  } // close year 2+ loop
  
  
  // to lazy to do this better right now
  for (p in 1:np){
    
    for (a in 1:n_ages){
     
      c_p_a_y_hat[p,a,ny_train] =  (n_p_a_y_hat[p, a, ny_train]) * (f / (f + m)); 
      
    }

  }
  

  

    for(p in 1:np){
    for(y in 1:ny_train){
      
      n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_p_a_y_hat[p,1:n_ages,y])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      // n_p_l_y_hat[y,p,1:n_lbins]  =  (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  .* selectivity_at_bin)';
      
      dens_p_y_hat[p,y] = sum((to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
     
      c_p_y_hat[p,y] =  sum((to_vector(c_p_a_y_hat[p,1:n_ages,y])));
      
      n_eff[p,y] = (1 + theta_d * sum(n_p_l_y[p,1:n_lbins, y])) / (1 + theta_d);
      
      n_raw[p,y] = sum(n_p_l_y[p,1:n_lbins, y]);
 
    }
  }
  
} // close transformed parameters block

model {
  
  real n;
  
  vector[n_lbins] prob_hat;
  
  vector[n_lbins] prob;

  real dml_tmp;
  
  real test;
  
  log_mean_recruits ~ normal(10,5);
  
  log_f ~ normal(log(.1),.1);
  
  // sigma ~ normal(1,.1);  // total error prior

  // proc_ratio ~ beta(2,2);
  
  sigma_r ~ normal(.7,.2);
  
  sigma_obs ~ normal(0.1,.2);
  
  // normal(0.5,.1);
  
  theta_d ~ normal(0.5,.1);

  p_length_50_sel ~ normal(length_50_sel_guess/loo, 1);
  
  // log_scalar ~ normal(log(2),1);
  
  raw ~ normal(0, sigma_r);
  
  for(y in 1:ny_train) {
    
    for(p in 1:np){

        if (sum(n_p_l_y[p,1:n_lbins,y]) > 0){
        
        if (do_dirichlet == 1){
        
        prob_hat = (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
        
        prob = (to_vector(n_p_l_y[p,1:n_lbins,y])  / sum(to_vector(n_p_l_y[p,1:n_lbins,y])));

        n = sum(n_p_l_y[p,1:n_lbins,y]);
        
        dml_tmp = lgamma(n + 1) -  sum(lgamma(n * prob + 1)) + lgamma(theta_d * n) - lgamma(n + theta_d * n) + sum(lgamma(n * prob + theta_d * n * prob_hat) - lgamma(theta_d * n * prob_hat)); // see https://github.com/merrillrudd/LIME/blob/9dcfc7f7d5f56f280767c6900972de94dd1fea3b/src/LIME.cpp#L559 for log transformation of dirichlet-multinomial in Thorston et al. 2017

        target += dml_tmp;
        
        } else 
        {
        
        (n_p_l_y[p,1:n_lbins,y]) ~ multinomial((to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins]))));
        } // close dirichlet statement
      
        } // close if any observations
          // print(n_p_l_y[p,1:n_lbins,y]);
          
        if (abund_p_y[p,y] > 0) {

        log(abund_p_y[p,y]) ~ normal(log(dens_p_y_hat[p,y] + 1e-3), sigma_obs); 
        
        }
        // log(c_p_y[p,y]) ~ normal(log(c_p_y_hat[p,y] + 1e-3), .05);
        // 

    } // close patch loop
    
  } // close year loop
  
  
}



