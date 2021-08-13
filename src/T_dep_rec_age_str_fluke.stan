functions {
  
  real T_dep(real sst, real Topt, real width){
    return exp(-0.5 * ((sst - Topt)/width)^2); // gussian temperature-dependent function
  }
  
}

data {
  
  // survey data 
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_p_a_y[np, n_ages, ny_train]; // SUM number of individuals in each age, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
  real abund_p_y[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  int proj_init[np, n_ages]; // data with which to initialize projection (one year after the end of n_p_a_y)
  
  // environmental data 
  
  real sst[np, ny_train]; // temperature data for training
  
  real sst_proj[np, ny_proj]; // temperature data for testing 
  
  // fish data
  
  real z;  // total mortality 
  
  real k;
  
  real loo;
  
  real length_50_sel_guess;
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  vector<lower=0>[n_ages] mean_length_at_age;
  
  matrix[n_ages, n_lbins] length_at_age_key;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  vector[np] patcharea;
  
  
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
  //real<lower=0> phi_obs; // I thought it was possible for this parameter to be negtive, but at one point got this error: Exception: neg_binomial_2_lpmf: Precision parameter is -0.317205, but must be > 0!  (in 'model1ee6785925e9_stage_ar_model_adult_dispersal' at line 145)
  real<lower=0> sigma_obs;
  
  real<lower = 0> p_length_50_sel; // length at 50% selectivity
  
  real<lower=0, upper=1> theta; // Bernoulli parameter for encounter probability
  
  // CURRENTLY NOT REALLY USING THIS, NEED TO INCORPORATE OR DELETE.
  real theta_threshold; // actual value to trigger positive encounters
  
  real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
}

transformed parameters{
  
  real T_adjust[np, ny_train]; // tuning parameter for SST suitability in each patch*year
  
  real sigma_r;
  
  real length_50_sel;
  
  real sel_delta;
  
  real mean_recruits;
  
  real n_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  real dens_p_y_hat [np, ny_train]; // for tracking sum density 
  
  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  // vector[n_ages] sel_at_age; // vector of selectivity at stage
  
  vector[n_ages] mean_selectivity_at_age; // mean selectivity at age
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  sel_delta = 2;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
  //  print("mean selectivity at age is ",mean_selectivity_at_age); // check that counts are different for every year
  
  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  //print("mean recruits is ",mean_recruits);
  
  // calculate temperature-dependence correction factor for each patch and year depending on SST
  for(p in 1:np){
    for(y in 1:ny_train){
      T_adjust[p,y] = T_dep(sst[p,y], Topt, width);  
    } // close years
  } // close patches
  
  // fill in year 1 of n_p_a_y_hat, initialized with mean_recruits 
  for(p in 1:np){
    for(a in 1:n_ages){
      if(a==1){
        n_p_a_y_hat[p,a,1] = mean_recruits * T_adjust[p,1]; // initialize age 0 with mean recruitment in every patch
      }
      else{
        n_p_a_y_hat[p,a,1] = n_p_a_y_hat[p,a-1,1] * (1-z); // initialize population with mean recruitment propogated through age classes with mortality
      }
      // TRIED BACKING INIT POP SIZE OUT FROM COUNT AT SEL_100, BUT IT CAUSED THIS ERROR:
      // Chain 1: Rejecting initial value:
      //Chain 1:   Error evaluating the log probability at the initial value.
      //Chain 1: Exception: multinomial_lpmf: Probabilities parameter is not a valid simplex. sum(Probabilities parameter) = -nan, but should be 1  (in 'modelb55735c34310_T_dep_rec_age_str_fluke' at line 272)
      // else if(a>=sel_100){
        //   n_p_a_y_hat[p,a,1] = n_p_a_y[p,a,1];
        // }
        // else{
          //   // backing out pre-fully-selected pop sizes to kick off the model
          //   n_p_a_y_hat[p,a,1] = n_p_a_y[p,sel_100,1] / pow(1-z, (sel_100-a)); 
          // }
          //    print("initial pop size for patch ",p," and age ",a," is ",n_p_a_y_hat[p,a,1]);
    } // close ages
  } // close patches
  
  
  // calculate recruitment deviates every year (not patch-specific)
  for (y in 2:ny_train){
    if (y == 2){ 
      rec_dev[y-1]  =  raw[1]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
      
    } // close y==2 case  
    else {
      
      rec_dev[y-1] =  alpha * rec_dev[y-2] + raw[y-1]; 
      
    } // close ifelse
    
    // describe population dynamics
    for(p in 1:np){
      
      // density-independent, temperature-dependent recruitment of age 1
      n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y-1];
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
        
        n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (1-z); // these just grow and die in the patch
        
      } // close ages for 2 to age at maturity
      } // close if 
      
      // pop dy for reproductive adults
      // mortality and dispersal are happening simultaneously here, between generations
      // because neither is patch-specific I don't think the order matters
      
      for(a in age_at_maturity:n_ages){
        
        // edge cases -- edges are reflecting
        if(p==1){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (1-z) * (1-d) + n_p_a_y_hat[p+1, a-1, y-1] * (1-z) * d;
        } // close patch 1 case 
        
        else if(p==np){
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (1-z) * (1-d) + n_p_a_y_hat[p-1, a-1, y-1] * (1-z) * d;
        } // close highest patch
        
        else{
          n_p_a_y_hat[p,a,y] = n_p_a_y_hat[p, a-1, y-1] * (1-z) * (1-2*d) + n_p_a_y_hat[p-1, a-1, y-1] * (1-z) * d + n_p_a_y_hat[p+1, a-1, y-1] * (1-z) * d;
          
        } // close if/else for all other patches
        
      }// close ages
    } // close patches 
    
    
  } // close year 2+ loop
  
  for(p in 1:np){
    for(y in 1:ny_train){
      dens_p_y_hat[p,y] = sum((to_vector(n_p_a_y_hat[p,1:n_ages,y]) .* mean_selectivity_at_age));

    }
  }

  
} // close transformed parameters block

model {
  
  theta ~ uniform(0, 1); // Bernoulli probability of encounter
  
  d ~ normal(0.1, 0.1); // dispersal rate as a proportion of total population size within the patch
  
  // RIGHT NOW THIS ISN'T DOING ANYTHING! FIX THAT
  theta_threshold ~ normal(50, 100); // number of total individuals for which encounters are positive
  
  log_mean_recruits ~ normal(1,5);
  
  Topt ~ normal(18, 4);
  
  width ~ normal(4, 4); 
  
  log_sigma_r ~ normal(log(.5),.1); // process error prior
  
  alpha ~ normal(0,.25); // autocorrelation prior
  
  //  phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu); 
  sigma_obs ~ normal(.1, 2); // think more about whether these numbers are reasonable
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .05);
  
  for(y in 2:ny_train) {
    
    // the n_p_a_y_transform business is here to force Stan to multiply a real (pop size) by a vector (mean sel) to pass to the multinomial
    // it didn't mind this before and I'm not sure what changed!
    //vector[n_ages] n_a_y_transform; // should get recreated for every y
    // matrix[np, n_ages] n_p_a_y_transform; // 
    
    //  for(p in 1:np){
      //    for(a in 1:n_ages){
        //      n_p_a_y_transform[p,a] = n_p_a_y_hat[p,a,y] .* mean_selectivity_at_age[a];
        //    }
        //  }
        
        // print("adjusted counts at age in year ",y," are ",n_a_y_transform); // check that counts are different for every year
        
        raw[y-1] ~ normal(0,sigma_r); // prior on raw process error (sorry this is kinda buried in the multinomial stuff)
        
        // if(sum(n_a_y[sel_100:n_ages,y]) > 0) { // this was still passing some vectors of all zeros to the multinomial, so I replaced  it with:
        // QUESTION: before we were conditioning on the sum true counts being >0. but if the sum *modeled* counts aren't >0, the multinomial doesn't work... so is it OK to condition on those instead?
        
        for(p in 1:np){
          //   if(sum(n_a_y_transform[1:n_ages]) > 0 && sum(n_a_y[sel_100:n_ages,y]) > 0) {
            
            // this should work just with the n_p_a_y clause, but I get this error:
            // Exception: multinomial_lpmf: Probabilities parameter is not a valid simplex. sum(Probabilities parameter) = -nan, but should be 1  (in 'modelea1f73d5b5cc_T_dep_rec_age_str_fluke' at line 261)
            if((abund_p_y[p,y]) > 0) {
              
              // multinomial to estimate relative abundance of stages
              
              /// n_a_y[1:n_ages,y] ~ multinomial(n_a_y_transform[1:n_ages] / sum(n_a_y_transform[1:n_ages])); // did I parameterize this right?
              
              // negative binomial to estimate absolute abundance (counts) -- was calibrated  to stage 3 before -- now to sum count
              //   sum(n_a_y[1:n_ages,y]) ~ neg_binomial_2(sum(n_a_y_transform[1:n_ages]) + 1e-3, phi_obs); 
              // FLAG -- is this the right way to recalibrate the magnitude of n_a_y by the absolute scale?
              // n_p_a_y[p,1:ns,y] ~ multinomial(to_vector(n_p_a_y_hat[p,1:ns,y]) .* sel_at_stage);
              (n_p_a_y[p,1:n_ages,y]) ~ multinomial(softmax(to_vector(n_p_a_y_hat[p,1:n_ages,y]) .*  mean_selectivity_at_age));
              //  sum(n_p_a_y[p,1:n_ages,y]) ~ neg_binomial_2(sum((to_vector(n_p_a_y_hat[p,1:n_ages,y]) .* mean_selectivity_at_age)) + 1e-3, phi_obs); 
              //abund_p_y[p,y] ~ lognormal(sum((to_vector(n_p_a_y_hat[p,1:n_ages,y]) .* mean_selectivity_at_age)), sigma_obs); 
              log(abund_p_y[p,y]) ~ normal(log(dens_p_y_hat[p,y] + 1e-6), sigma_obs); 

              
              1 ~ bernoulli(theta);
              
              
            } else { // only evaluate length comps if there are length comps to evaluate
            
            0 ~ bernoulli(theta);
            
            } // close else 
        } // close patch loop
        
  } // close year loop
  
  
}

generated quantities{
  
  real pp_proj_n_p_a_y_hat[np,n_ages, ny_proj];
  
  real pp_proj_dens_p_y_hat[np, ny_proj];
  
  real tmp_proj[np,n_ages, ny_proj];
  
  real n_p_a_y_proj[np,n_ages, ny_proj]; // adding this as a post-selectivity forecast 
  
  real rec_dev_proj[ny_proj - 1];
  
  real T_adjust_proj[np,ny_proj];
  
  real pp_theta;
  
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
      
      for(y in 2:ny_proj){
        
        if (y == 2){
          rec_dev_proj[y-1]  = rec_dev[y-1]; // initiate projections with last estimated process errors
        }
        else{
          rec_dev_proj[y-1] = alpha * rec_dev_proj[y-2] + normal_rng(0,sigma_r); // generate autoregressive recruitment deviates
          
        }
      }
      
      // project pop dy
      
      for(y in 1:ny_proj){
        
        for(p in 1:np){
          
          if(y==1){
            // initiate projection with fixed observation
            for(a in 1:n_ages){
              //      pp_proj_n_p_a_y_hat[p,a,1] = proj_init[p,a]; 
              tmp_proj[p,a, 1] = fmin(proj_init[p,a] / mean_selectivity_at_age[a], 0); // transform into pre-selectivity units
            }
          } else { // add case for all other years
          
          for(a in 1:n_ages){
            
            // project age 1
            if(a==1){
              tmp_proj[p,1, y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[p,y]; 
            }
            // project non-reproductive ages
            else if(a < age_at_maturity){
              tmp_proj[p,a,y] = tmp_proj[p,a-1,y-1] * (1-z);// these just grow and die in the patch
            }
            // project reproductive ages
            else{
              if(p==1){
                tmp_proj[p,a,y] = tmp_proj[p, a-1, y-1] * (1-z) * (1-d) + tmp_proj[p+1, a-1, y-1] * (1-z) * d;
              } // close patch 1 case 
              
              else if(p==np){
                tmp_proj[p,a,y] = tmp_proj[p, a-1, y-1] * (1-z) * (1-d) + tmp_proj[p-1, a-1, y-1] * (1-z) * d;
              } // close highest patch
              
              else{
                tmp_proj[p,a,y] = tmp_proj[p, a-1, y-1] * (1-z) * 2*(1-d) + tmp_proj[p-1, a-1, y-1] * (1-z) * d + tmp_proj[p+1, a-1, y-1] * (1-z) * d;
                
              } // close if/else for main patches
              
            } // close else for reproductive age group
            
            //     print("tmp_proj for patch ",p,", year ",y,", age ",a,", is: ",tmp_proj[p,a,y]); 
            
          } // close age loop
          } // close else for all years except the initial one 
        } // close year loop
      } // close patch loop -- end of pop dy
      
      // no longer trying to simulate selectivity/sampling process 
      for(p in 1:np){
        for(a in 1:n_ages){
          for(y in 1:ny_proj){
            pp_proj_n_p_a_y_hat[p,a,y] = tmp_proj[p,a,y];
          }
        }
      }
      
  for(p in 1:np){
    for(y in 1:ny_proj){
      pp_proj_dens_p_y_hat[p,y] = sum((to_vector(pp_proj_n_p_a_y_hat[p,1:n_ages,y]) .* mean_selectivity_at_age));
    }
  }
      
      // simulate selectivity and sampling error
      // for(p in 1:np){
        //   for(y in 1:ny_proj){
          //     pp_theta = bernoulli_rng(theta); // not being used for anything right now
          //     
          //     // right now we're just missing the whole negative binomial process because I am not sure how it fits in with the multinomial
          //     // since this is a calculation we can't just set pp_proj_n_p_a_y_hat equal to two things
          // 
          //     if(sum(tmp_proj[p,1:n_ages,y]) > 0){
            //     pp_proj_n_p_a_y_hat[p,1:n_ages,y] = multinomial_rng(to_vector(tmp_proj[p,1:n_ages,y]) / sum(to_vector(tmp_proj[p,1:n_ages,y])), 1000);
            //     
            //     } else {
              //       for(a in 1:n_ages){
                //                 pp_proj_n_p_a_y_hat[p,a,y] = 0;
                //       }
                //     }
                //     // impose selectivity on pp_proj_n_p_a_y_hat so we can compare it to the data 
                //     n_p_a_y_proj[p,1:n_ages,y] = to_array_1d(to_vector(pp_proj_n_p_a_y_hat[p,1:n_ages,y]) .* to_vector(mean_selectivity_at_age));
                //     
                //   } // close year loop
                //   
                // } // close patch loop
                
                
} // close generated quantities


