functions {
  
  real T_dep(real sbt, real Topt, real width){
    return exp(-0.5 * ((sbt - Topt)/width)^2); // gussian temperature-dependent function
  }
  
  matrix age_at_length_key(real loo, real l0, real k, real cv, int n_lbins, int n_ages){
    
    vector[n_lbins] mean_age_at_length;
    vector[n_lbins] sigma_age_at_length;
    matrix[n_lbins, n_ages] prob_age_at_length;
    
    // make vector of mean ages for each length 
    for(i in 1:n_lbins){
      mean_age_at_length[i] = log((loo - i) / (loo - l0)) / -k; 
      sigma_age_at_length[i] = mean_age_at_length[i] * cv; 
    }
    
    for(i in 1:n_lbins){
      for(j in 1:n_ages){
        if(j < n_ages){
          prob_age_at_length[i,j] = normal_cdf(j+1, mean_age_at_length[i], sigma_age_at_length[i]) - normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);  // analog of pnorm in R
        }
        else{
          prob_age_at_length[i,j] = normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);
        } // close if/else
      } // close ages
    } // close lengths
    return prob_age_at_length;
  } // close function
  
} // close functions block



data {
  
  // survey data 
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  int ny_train; // years for training
  
  // int ny_proj; // number of years to forecast 
  
    int n_lbins; // number of length bins (here just the range of cm values)
  
  int n_p_l_y[np, n_lbins, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
  matrix[n_ages, n_lbins] l_at_a_key;
 
  real abund_p_y[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
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
  
  real<lower=1e-3> width; // sensitivity to temperature variation
  
  real Topt; //  temp at which recruitment is maximized
  
  real log_sigma_r; // sigma recruits
  
  real<lower = -1, upper = 1> alpha; // autocorrelation term 
  
  // real  log_mean_recruits; // log mean recruits per patch, changed to one value for all space/time
  
  vector[np] log_mean_recruits;
  
  vector[ny_train-1] raw; // array of raw recruitment deviates, changed to one value per year
  
  //  real<lower = 1e-3> sigma_obs;
  //real<lower=0> phi_obs; // I thought it was possible for this parameter to be negtive, but at one point got this error: Exception: neg_binomial_2_lpmf: Precision parameter is -0.317205, but must be > 0!  (in 'model1ee6785925e9_stage_ar_model_adult_dispersal' at line 145)
  real<lower=0> sigma_obs;
  
  real p_length_50_sel; // length at 50% selectivity
  
  real<lower=0, upper=1> theta; // Bernoulli parameter for encounter probability
  
  real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
  real log_f;
  // real log_scalar;
  
}

transformed parameters{
  
  real T_adjust[np, ny_train]; // tuning parameter for sbt suitability in each patch*year
  
  real sigma_r;
  
  real length_50_sel;
  
  real sel_delta;
  
  vector[np] mean_recruits;
  
  real n_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  matrix[np, n_lbins] n_p_l_y_hat[ny_train]; // array number of years containing matrices with numbers at patch, length bin, and year 

  real dens_p_y_hat [np, ny_train]; // for tracking sum density 
  
  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  // vector[n_ages] sel_at_age; // vector of selectivity at stage
  
  // vector[n_ages] mean_selectivity_at_age; // mean selectivity at age
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  // real scalar = exp(log_scalar);
  
  real f = exp(log_f);
  
  real z = exp(-(f + m));
  
  sel_delta = 2;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  // mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
  //  print("mean selectivity at age is ",mean_selectivity_at_age); // check that counts are different for every year
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  //print("mean recruits is ",mean_recruits);
  
  // calculate temperature-dependence correction factor for each patch and year depending on sbt
  for(p in 1:np){
    for(y in 1:ny_train){
      T_adjust[p,y] = T_dep(sbt[p,y], Topt, width);  
    } // close years
  } // close patches
  
  // fill in year 1 of n_p_a_y_hat, initialized with mean_recruits 
  for(p in 1:np){
    for(a in 1:n_ages){
      if(a==1){
        n_p_a_y_hat[p,a,1] = mean_recruits[p] * T_adjust[p,1] * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
      }
      else{
        n_p_a_y_hat[p,a,1] = n_p_a_y_hat[p,a-1,1] * z; // initialize population with mean recruitment propogated through age classes with mortality
      }
      
    } // close ages
  } // close patches
  
  
  // calculate recruitment deviates every year (not patch-specific)
  for (y in 2:ny_train){
    if (y == 2){ 
      rec_dev[y-1]  =  raw[2]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
      // need to fix this awkward burn in
    } // close y==2 case  
    else {
      
      rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y-1]; 
      
    } // close ifelse
    
    // describe population dynamics
    for(p in 1:np){
      
      // density-independent, temperature-dependent recruitment of age 1
      n_p_a_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y-1];
      // why estimate raw and sigma_r? we want to estimate process error
      // if we got rid of sigma_r, we would be saying that raw could be anything
      // that means raw could pick any value, and it would pick deviates to perfectly match abundance index
      // sigma_r scales the amount of process error that we say is possible
      // letting it go to infinity means the model will always fit the data perfectly 
      // raw is the realized process error 
      // exp(rec_dev[y-1] - pow(sigma_r,2)/2) = random variable with mean 0 and SD sigma_r
      // allows us to have a different recruitment deviation every year even though sigma_r is constant 
      
      // pop dy for non-reproductive ages 
      if(age_at_maturity > 1){ // confirm that there are non-reproductive age classes above 1
      for(a in 2:(age_at_maturity-1)){
        
        n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * z; // these just grow and die in the patch
        
      } // close ages for 2 to age at maturity
      } // close if 
      
      // pop dy for reproductive adults
      // mortality and dispersal are happening simultaneously here, between generations
      // because neither is patch-specific I don't think the order matters
      
      for(a in age_at_maturity:n_ages){
        
        // edge cases -- edges are reflecting
        if(p==1){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * z * (1-d) + n_p_a_y_hat[p+1, a-1, y-1] * z * d;
        } // close patch 1 case 
        
        else if(p==np){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * z * (1-d) + n_p_a_y_hat[p-1, a-1, y-1] * z * d;
        } // close highest patch
        
        else{
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * z * (1-2*d) + n_p_a_y_hat[p-1, a-1, y-1] * z * d + n_p_a_y_hat[p+1, a-1, y-1] * z * d;
          
        } // close if/else for all other patches
        
      }// close ages
    } // close patches 
    
    
  } // close year 2+ loop
  
  for(p in 1:np){
    for(y in 1:ny_train){
      
      n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_p_a_y_hat[p,1:n_ages,y])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      // n_p_l_y_hat[y,p,1:n_lbins]  =  (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  .* selectivity_at_bin)';
      
      dens_p_y_hat[p,y] = sum((to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
      
    }
  }
  
  
} // close transformed parameters block

model {
  
  theta ~ uniform(0, 1); // Bernoulli probability of encounter
  
  d ~ normal(0.1, 0.1); // dispersal rate as a proportion of total population size within the patch
  
  log_f ~ normal(log(m / 2),.5);
  
  log_mean_recruits ~ normal(7,5);
  
  Topt ~ normal(18, 4);
  
  width ~ normal(4, 4); 
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  sigma_obs ~ normal(.1, .1); // think more about whether these numbers are reasonable
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .2);
  
  // log_scalar ~ normal(log(2),1);
  
  for(y in 2:ny_train) {
    
    raw[y-1] ~ normal(0,sigma_r); // prior on raw process error (sorry this is kinda buried in the multinomial stuff)
    
    for(p in 1:np){
      if((abund_p_y[p,y]) > 0) {
        
        (n_p_l_y[p,1:n_lbins,y]) ~ multinomial((to_vector(n_p_l_y_hat[y,p,1:n_lbins]) / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins]))));

        log(abund_p_y[p,y]) ~ normal(log( dens_p_y_hat[p,y] + 1e-6), sigma_obs); 
        
        1 ~ bernoulli(theta);
        
        
      } else { // only evaluate length comps if there are length comps to evaluate
      
      0 ~ bernoulli(theta);
      
      } // close else 
    } // close patch loop
    
  } // close year loop
  
  
}



