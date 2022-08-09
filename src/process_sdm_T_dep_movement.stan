functions {
  // Gaussian temperature dependence function
  real T_dep(real sbt, real Topt, real width, real exp_yn){
    if(exp_yn==1){
    return exp(-0.5 * ((sbt - Topt)/width)^2); // gaussian temperature-dependent function
  } else{
        // return (-0.1 * ((sbt - Topt)/width)^2); // gaussian temperature-dependent function, not exponentiated, for temperature-dependent process that are exponentiated later (like movement)
        // I'm taking this out because I can't get it to be non-negative without the exp! 
        
        // return normal_lpdf(sbt | Topt, width);
        return log((1 / sqrt(2 * pi() * width)) * exp(-pow(sbt - Topt,2) / (2 * pow(width,2))));
  }
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
  
  vector colSums(matrix M){
    int ncol; 
    vector[cols(M)] sums; 
    
    ncol = cols(M); 
    
    for(i in 1:ncol){
      sums[i] = sum(M[,i]); //sums[i] = sum(col(M,i)); 
    }
    return(sums);
  }
  

} // close functions block

data {
  
  // survey data 
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  matrix[n_ages, n_lbins] l_at_a_key;
  
  vector[n_ages] wt_at_age;
  
  real abund_p_y[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  // environmental data 
  
  real sbt[np, ny_train]; // temperature data for training
  
  real sbt_proj[np, ny_proj];
  
  // fish data
  
  real m;  // total mortality 
  
  real k;
  
  real loo;
  
  real t0;
  
  real cv;
  
  real f[n_ages, ny_train]; 
  
  real f_proj[n_ages, (ny_proj+1)];
  
  real length_50_sel_guess;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  int<lower = 0, upper = 1> do_dirichlet;
  
  int<lower = 0, upper = 1> T_dep_recruitment;
  
  int<lower = 0, upper = 1> T_dep_mortality;
  
  int<lower = 0, upper = 1> T_dep_movement;
  
  int<lower = 0, upper = 1> eval_l_comps;
  
  int<lower = 0, upper = 1> spawner_recruit_relationship;
  
  int<lower = 0, upper = 1> run_forecast;
  
    int<lower = 0, upper = 1> exp_yn;
  
  int<lower = 0, upper = 1> process_error_toggle;
  
  int n_p_l_y[np, n_lbins, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
}

transformed data{
  
  vector[n_ages] maturity_at_age; // vector of probabilities of being mature at each age, currently binary (0/1) and taken as known
  
  matrix[np, np] adj_m; // adjacency matrix for patches 
  
  matrix[np, np] outer; // outer difference matrix 
  
  //int exp_yn; 
  
  // in R this is just outer(np, np, "-") 
  for(i in 1:np){
    for(j in 1:np){
      outer[i, j] = i - j; // should just create a matrix that counts from 1 to np along rows and columns
      if(fabs(outer[i,j])==1) {
        adj_m[i,j] = 1; // fill with 1s for adjacent patches and 0s otherwise 
      }else{
        adj_m[i,j] = 0; 
      }
    }
  }
  
 // print("the outer matrix is ",outer); 
 // print("the adjacency matrix is ",adj_m); 
  
  for(a in 1:n_ages){
    if(a < age_at_maturity){
      maturity_at_age[a]=0;
    }else{
      maturity_at_age[a]=1;
    }
  }
  
   // exp_yn=1;
  
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
  
  // real<lower = 1e-6> sigma_total;
  
  real<lower = 1e-6> sigma_obs;
  
  real<lower = 1e-6> sigma_r_raw;

  
  real<lower=0.5> width; // sensitivity to temperature variation
  
  real Topt; //  temp at which recruitment is maximized
  
  real<lower = 0, upper = 0.99> alpha; // autocorrelation term
  
  real  log_mean_recruits; // log mean recruits per patch, changed to one value for all space/time
  
  vector[ny_train] raw; // array of raw recruitment deviates, changed to one value per year
  
  real<upper = 0.8> p_length_50_sel; // length at 50% selectivity
  
  real<lower=0> beta_obs; // controls how fast detection goes up with abundance
  
  // real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
  real<lower=0, upper=1> d; // increasing bounds on this for the temperature-dependent movement model 
  
  vector<lower = 0, upper = 1>[np] init_dep;
  
  real<lower = 0> theta_d;
  
  real<lower=0, upper=1> h;
  
  real beta_t; // responsiveness of movement to temperature
  
  real beta_rec; // responsivenses of mean recruits to temperature
  
  real beta_obs_int; // intercept of detection probability
  
  real log_r0;
  
}

transformed parameters{
  
  matrix[np, ny_train] T_adjust; // tuning parameter for sbt suitability in each patch*year
  
  real length_50_sel;
  
  real sel_delta;
  
  real mean_recruits;
  
  matrix[np, ny_train] theta; // Bernoulli probability of encounter  
  
  // real n_p_a_y_hat [np, n_ages,ny_train]; // array of numbers at patch, stage, and year 
  
  matrix[np, n_ages] n_p_a_y_hat[ny_train];
  
  matrix[np, n_lbins] n_p_l_y_hat[ny_train]; // array number of years containing matrices with numbers at patch, length bin, and year 
  
  real dens_p_y_hat [np, ny_train]; // for tracking sum density 
  
  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  real surv[np, n_ages, ny_train];
  
  matrix[np, np] diff_m; // dispersal matrix (setting up for temperature-dependent movement, but doesn't actually need to vary by year)
  
  matrix[np, np] mov_inst_m[ny_train]; // array containing instantaneous movement matrices (dispersal + temperature-dependent taxis)
  
  matrix[np, np] mov_m[ny_train]; // Dan thinks this is the wrong syntax -- need to fix 
  // array[ny_train] matrix[np, np] mov_m; // array containing annualized movement matrices (dispersal + temperature-dependent taxis)
  
  matrix[np, np] T_adjust_m[ny_train]; // array containing temperature preference matrices
  
  matrix[np, np] tax_m[ny_train]; // array containing taxis between patches every year 
  
  real ssb0;
  
  vector[n_ages] unfished;
  
  matrix[np, ny_train] ssb;
  
  vector[n_ages] v_ssb; // placeholder vector for adding up ssb 
  
  real r0;
  
  vector[np] v_in; // vector for matrix multiplication 
  
  vector[np] v_out; // vector for matrix multiplication 
  
  // real sigma_r = sigma_total / 2;
  
  // real sigma_obs = sigma_total / 2;
  
  real sigma_r;
  
  r0 = exp(log_r0);
  
  sigma_r = sigma_r_raw * process_error_toggle;
  
  ssb0 = -999;
  
  // for(a in 1:n_ages){
  //   unfished[a] = 999;
  //   v_ssb[a] = 999;
  // }
  // for(p in 1:np){
  //   for(y in 1:ny_train){
  //     ssb[p,y] = 999;
  //   }
  // }
  // 
  
  unfished[1] = r0;
  if(spawner_recruit_relationship==1){
    for(a in 2:n_ages){

        unfished[a] = unfished[a-1] * exp(-m);
      }

    ssb0 = sum(unfished .* maturity_at_age .* wt_at_age);
  }
  
  
  sel_delta = 1;
  
  length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
  
  selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  // mean_selectivity_at_age = length_at_age_key * selectivity_at_bin; // calculate mean selectivity at age given variance in length at age
  
  //  print("mean selectivity at age is ",mean_selectivity_at_age); // check that counts are different for every year
  
  mean_recruits = exp(log_mean_recruits);
  
  //print("mean recruits is ",mean_recruits);
  
  // calculate temperature-dependence correction factor for each patch and year depending on sbt
  for(p in 1:np){
    for(y in 1:ny_train){
      T_adjust[p,y] =   T_dep(sbt[p,y], Topt, width, exp_yn);  
    } // close years
  } // close patches
  
  
  // calculate total annual mortality from instantaneous natural + fishing mortality data 
  
  for(p in 1:np){
    for(a in 1:n_ages){
      for(y in 1:ny_train){
        
        if(T_dep_mortality==1){
          surv[p,a,y] = exp(-((f[a,y] + m) * T_adjust[p,y]));
        }
        
        if(T_dep_mortality==0){
          surv[p,a,y] = exp(-(f[a,y] + m)) ;
        }
        
      }
    }
  }
  
  //  calculate annual movement matrix if using spatially varying dispersal
  if(T_dep_movement==1){
    
    diff_m = adj_m * d; // set up matrix as adjacency * dispersal
    diff_m = add_diag(diff_m, -1 * colSums(diff_m)); // rescale the columns so they will sum to 1 at the end
    
    // set up array of matrices containing 
    for(y in 1:ny_train){
      for(i in 1:np){
        for(j in 1:np){
          // in R this is just outer(np, np, "-") 
          
            if(exp_yn==1){
                     T_adjust_m[y,i,j] = exp(beta_t * (log(T_adjust[i,y]) - log(T_adjust[j,y]))); 

            } else {
                  
                  T_adjust_m[y,i,j] = fmin(500,exp(beta_t * (T_adjust[i,y] - T_adjust[j,y]))); 
            }
          
         // print("T_adjust in patch ",i," and year ",y," is ",T_adjust[i,y]); 
         //           print("T_adjust in patch ",j," and year ",y," is ",T_adjust[i,y]); 

          // print("T_adjust_m in year ",y," from patch ",i," into patch ",j," is ",T_adjust_m[y,i,j]);
        }
      }
      tax_m[y] = adj_m .* T_adjust_m[y]; 
      
      tax_m[y] = add_diag(tax_m[y], -1 * colSums(tax_m[y])); // fill in the diagonal with within-patch "taxis" so everything sums to 1 

      
      mov_inst_m[y] =  diff_m + tax_m[y]; // movement as a sum of diffusion and taxis (can cancel each other out)
      mov_m[y] = matrix_exp(mov_inst_m[y]); // matrix exponentiate, although see https://discourse.mc-stan.org/t/matrix-exponential-function/9595
      
      // print(colSums(mov_m[y]))

      if ((sum(colSums(mov_m[y])) / np - 1) > .001 ){
        print("Something has gone very wrong, movement matrix columns do not sum to 1")
        print(colSums(mov_m[y]))
        print("width is", width)
        print("Topt is", Topt)
        print(diagonal(mov_inst_m[y]))
        
      }
        // print(colSums(mov_m[y]))

      // print("column sums of the annualized movement matrix in year ",y," is ",colSums(mov_m[y]));
      //  print("the annualized movement matrix in year ",y," is ",mov_m[y]);
    }
  } else {
    for(z in 1:np){
      for(x in 1:np){
    diff_m[z,x] = 999; 
  }
  }
  
    }
  
  // fill in year 1 of n_p_a_y_hat, initialized with mean_recruits 
  for(p in 1:np){
    for(a in 1:n_ages){
      if(a==1){
        // n_p_a_y_hat[p,a,1] = mean_recruits[p] * T_adjust[p,1] * exp(raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
        
        if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
          n_p_a_y_hat[1,p,a] = init_dep[p] * mean_recruits * beta_rec * T_adjust[p,1] * exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
        }
        if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
          n_p_a_y_hat[1,p,a] = init_dep[p] * mean_recruits * exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
        }
        if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
          n_p_a_y_hat[1,p,a] = init_dep[p] * r0 *  exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // scale it down a bit -- historical fishing was still occurring
        }
        if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
          n_p_a_y_hat[1,p,a] = init_dep[p] *r0 *  exp(sigma_r *raw[1] - pow(sigma_r,2) / 2) * T_adjust[p,1] * beta_rec;
        }
      } // close age==1 case
      else{
        n_p_a_y_hat[1,p,a] = n_p_a_y_hat[1,p,a-1] * surv[p,a-1,1]; // initialize population with mean recruitment propogated through age classes with mortality
      }
      
    } // close ages
    
          ssb[p,1] = sum(to_vector(n_p_a_y_hat[1,p,1:n_ages]) .* maturity_at_age .* wt_at_age);

  } // close patches
  
  
  // calculate recruitment deviates every year (not patch-specific)
  for (y in 2:ny_train){
    
      if (y == 2){ 
        rec_dev[y-1]  =  sigma_r * raw[y]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific
        // need to fix this awkward burn in
      } // close y==2 case  
      else {
        
        rec_dev[y-1] =  alpha * rec_dev[y-2] +  sqrt(1 - pow(alpha,2)) *  sigma_r * raw[y]; // why does rec_dev[y-1] use raw[y]?

      } // close ifelse
      
      // print(rec_dev[y-1]);
    
    // describe population dynamics
    for(p in 1:np){
      
      // density-independent, temperature-dependent recruitment of age 1
      
      if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
        n_p_a_y_hat[y,p,1] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y-1] * beta_rec;
      }
      if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
        n_p_a_y_hat[y,p,1] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) ;
      }
      
      if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
        n_p_a_y_hat[y,p,1] = (0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) + ssb[p, y-1] * (h - 0.2));
        
        n_p_a_y_hat[y,p,1] =  n_p_a_y_hat[y,p,1] *  exp(rec_dev[y-1] - pow(sigma_r,2)/2);
        
      }
      if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
        n_p_a_y_hat[y,p,1] = ((0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) +  ssb[p, y-1] * (h - 0.2))) * T_adjust[p,y-1];
        
        n_p_a_y_hat[y,p,1] =  n_p_a_y_hat[y,p,1] *  exp(rec_dev[y-1] - pow(sigma_r,2)/2);
      }
      // 
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
      
       n_p_a_y_hat[y,p,2:(age_at_maturity-1)] = (n_p_a_y_hat[y-1,p, 1:(age_at_maturity-2)]) .* to_row_vector(surv[p,1:(age_at_maturity-2),y-1]);
     
      } // close if 
    } // close patches 
    
    // pop dy for reproductive adults
    
    
    if(T_dep_movement==0){
      
      // mortality and dispersal are happening simultaneously here, between generations
      // because neither is patch-specific I don't think the order matters
      for(p in 1:np){
        for(a in age_at_maturity:n_ages){
          // edge cases -- edges are reflecting
          if(p==1){
            n_p_a_y_hat[y,p,a] = n_p_a_y_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-d) + n_p_a_y_hat[y-1,+1, a-1] * surv[p+1,a-1,y-1] * d;
          } // close patch 1 case 
          
          else if(p==np){
            n_p_a_y_hat[y,p,a] = n_p_a_y_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-d) + n_p_a_y_hat[y-1,p-1, a-1] * surv[p-1,a-1,y-1] * d;
          } // close highest patch
          
          else{
            n_p_a_y_hat[y,p,a] = n_p_a_y_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-2*d) + n_p_a_y_hat[y-1,p-1, a-1] * surv[p-1,a-1,y-1] * d + n_p_a_y_hat[y-1,p+1, a-1] * surv[p+1,a-1,y-1] * d;
            
          } // close if/else for all other patches
          
        }// close ages
      } // close patches 
    } // close T-dep movement if 
    
    // this code block calculates adult population size based on survival and directional movement 
    if(T_dep_movement==1){
      
      
      for (p in 1:np){
        
        n_p_a_y_hat[y,p, age_at_maturity:n_ages] = n_p_a_y_hat[y-1,p,  (age_at_maturity - 1):(n_ages - 1)] .* to_row_vector(surv[p, (age_at_maturity - 1):(n_ages - 1),y-1]);
        
      }
 
      for(a in age_at_maturity:n_ages){
        
        // some acrobatics required here, because Stan won't do matrix multiplication with an array of reals like n_p_a_y_hat
        // instead we do the matrix multiplication with a placeholder vector and then populate n_p_a_y_hat 
        
        // fill in placeholder vector with reproductive ages across patches, and do mortality   
        for(p in 1:np){
          v_in[p] = n_p_a_y_hat[y,p, a]; 
        }
        v_out = mov_m[y] * v_in; // redistribute each age among patches according to the movement matrix 
        
        // fill in n_p_a_y_hat
        for(p in 1:np){
          n_p_a_y_hat[y,p,a] = v_out[p]; 
          
        }
        
      } // close ages
    }// close T-dep movement if 
    
    for (p in 1:np){
      
        ssb[p,y]  =  sum(to_vector(n_p_a_y_hat[y,p,1:n_ages]) .* maturity_at_age .* wt_at_age);

    }
    
    
  } // close year 2+ loop
  
  for(p in 1:np){
    for(y in 1:ny_train){
      
      n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_p_a_y_hat[y,p,1:n_ages])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      // n_p_l_y_hat[y,p,1:n_lbins]  =  (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  .* selectivity_at_bin)';
      
      dens_p_y_hat[p,y] = sum((to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
      
    //  print("for patch ",p," in year ",y," dens_p_y_hat is ",dens_p_y_hat[p, y]);
      
      theta[p,y] = ((1/(1+exp(-(beta_obs_int + beta_obs*log(dens_p_y_hat[p,y] + 1e-6))))));

      // print(theta[p,y])
    }
  }

  
} // close transformed parameters block

model {
  
  real n;
  
  vector[n_lbins] prob_hat;
  
  vector[n_lbins] prob;
  
  real dml_tmp;
  
  real test;
  
  init_dep ~ beta(1.5,3);
  
  beta_obs ~ normal(0.001,0.1); 
  
  beta_obs_int ~ normal(-100,4);

  raw ~ normal(0, 1);

  sigma_r_raw ~ normal(.2,.1);

  sigma_obs ~ normal(0.1,.2);

  log_mean_recruits ~ normal(7,5);
  
  h ~ normal(0.6, 0.25);
  
  log_r0 ~ normal(15,5);

  Topt ~ normal(18, 2);
  
  width ~ normal(4, 2); 
  
  beta_t ~ normal(0,2);
  
  beta_rec ~ normal(0,2);
  
  alpha ~  beta(12,20); // concentrated around 0.4
  
  d ~ normal(0.1, 0.1); // dispersal rate as a proportion of total population size within the patch
  
  // sigma ~ normal(1,.1);  // total error prior
  // 
  // proc_ratio ~ beta(2,2);
  // 
  
  
  // sigma_obs ~ normal(.1, .1); // think more about whether these numbers are reasonable
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .2);
  
  // log_scalar ~ normal(log(2),1);
  
  
  theta_d ~ normal(0.5,.1);
  
  
  for(y in 2:ny_train) {
    
    for(p in 1:np){
      
      if((abund_p_y[p,y]) > 0) {
        
        if(eval_l_comps==1){
          if (sum(n_p_l_y[p,1:n_lbins,y]) > 0) {
            
            if (do_dirichlet == 1){
              
              prob_hat = (to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
              
              prob = (to_vector(n_p_l_y[p,1:n_lbins,y])  / sum(to_vector(n_p_l_y[p,1:n_lbins,y])));
              
              n = sum(n_p_l_y[p,1:n_lbins,y]);
              
              dml_tmp = lgamma(n + 1) -  sum(lgamma(n * prob + 1)) + lgamma(theta_d * n) - lgamma(n + theta_d * n) + sum(lgamma(n * prob + theta_d * n * prob_hat) - lgamma(theta_d * n * prob_hat)); // see https://github.com/merrillrudd/LIME/blob/9dcfc7f7d5f56f280767c6900972de94dd1fea3b/src/LIME.cpp#L559 for log transformation of dirichlet-multinomial in Thorston et al. 2017
              
              // dml_tmp = lgamma(n + 1) - sum(lgamma(n * prob_hat + 1)) + (lgamma(theta_d * n) - lgamma(n + theta_d * n)) * prod(((lgamma(n * prob_hat + theta_d * n * prob))./(lgamma(theta_d * n * prob))));
              
              // test = prod(1:10);
              
              target += (dml_tmp);
              
            } else {
              
              (n_p_l_y[p,1:n_lbins,y]) ~ multinomial((to_vector(n_p_l_y_hat[y,p,1:n_lbins])  / sum(to_vector(n_p_l_y_hat[y,p,1:n_lbins]))));
              
            } // close dirichlet statement
            
          } // close if any positive length comps
          
        } // close eval_length_comps
        
        log(abund_p_y[p,y]) ~ normal(log(dens_p_y_hat[p,y] + 1e-6), sigma_obs); 
        
        1 ~ bernoulli(theta[p,y]);

      
      } else { // only evaluate density if there are length comps to evaluate
      
      0 ~ bernoulli(theta[p,y]);
      
      } // close else 
    } // close patch loop
    
  } // close year loop
  
  
}


generated quantities {
  // real proj_n_p_a_y_hat[np, n_ages, ny_proj+1];
  // real T_adjust_proj[np, ny_proj];
  // vector[ny_proj] rec_dev_proj;
  // vector[ny_proj] raw_proj;
  // real surv_proj[n_ages, (ny_proj+1)];
  // matrix[np, n_lbins] proj_n_p_l_y_hat[ny_proj];
  // real proj_dens_p_y_hat [np, ny_proj];
  // 
  // 
  // if(run_forecast==1){
    //   for(p in 1:np){
      //     for(y in 1:ny_proj){
        //       T_adjust_proj[p,y] = T_dep(sbt_proj[p,y], Topt, width);
        //     } // close years
        //   } // close patches
        //   
        //   for(p in 1:np){
          //     for(a in 1:n_ages){
            //       for(y in 1:(ny_proj+1)){
              //         
              //         if(T_dep_mortality==0){
                //           surv_proj[a,y] = exp(-(f_proj[a,y] + m));
                //         }
                //         if(T_dep_mortality==1){
                  //           surv_proj[a,y] = exp(-(f_proj[a,y] + m))* T_adjust_proj[p,y];
                  //           
                  //         }
                  //       }
                  //     }
                  //   }
                  //   
                  //   // initialize with final year of our model
                  //   proj_n_p_a_y_hat[,,1] = n_p_a_y_hat[,,ny_train];
                  //   rec_dev_proj[1] = rec_dev[ny_train-1];
                  //   raw_proj[1] = raw[ny_train];
                  //   
                  //   // project pop dy
                  //   for(y in 2:ny_proj){
                    //     raw_proj[y] = normal_rng(0, sigma_r);
                    //     //  print("raw_proj in year ",y," is ",raw_proj[y]);
                    //     rec_dev_proj[y] = alpha * rec_dev_proj[y-1] + raw_proj[y];
                    //     //  print("rec_dev_proj in year ",y," is ",rec_dev_proj[y]);
                    //     
                    //   }
                    //   
                    //   for(y in 2:(ny_proj+1)){
                      //     for(p in 1:np){
                        //       
                        //       if(T_dep_recruitment==1){
                          //         proj_n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[p,y-1];
                          //       }
                          //       if(T_dep_recruitment==0){
                            //         proj_n_p_a_y_hat[p,1,y] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2);
                            //       }
                            //       
                            //       if(age_at_maturity > 1){
                              //         for(a in 2:(age_at_maturity-1)){
                                //           proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1];
                                //         } // close ages for 2 to age at maturity
                                //       } // close if
                                //       
                                //       for(a in age_at_maturity:n_ages){
                                  //         if(p==1){
                                    //           proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-d) + proj_n_p_a_y_hat[p+1, a-1, y-1] * surv_proj[a-1,y-1] * d;
                                    //         } // close patch 1 case
                                    //         
                                    //         else if(p==np){
                                      //           proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-d) + proj_n_p_a_y_hat[p-1, a-1, y-1] * surv_proj[a-1,y-1] * d;
                                      //         } // close highest patch
                                      //         
                                      //         else{
                                        //           proj_n_p_a_y_hat[p,a,y] = proj_n_p_a_y_hat[p, a-1, y-1] * surv_proj[a-1,y-1] * (1-2*d) + proj_n_p_a_y_hat[p-1, a-1, y-1] * surv_proj[a-1,y-1] * d + proj_n_p_a_y_hat[p+1, a-1, y-1] * surv_proj[a-1,y-1] * d;
                                        //           
                                        //         } // close if/else for all other patches
                                        //         
                                        //       }// close ages
                                        //     } // close patches
                                        //     
                                        //     
                                        //   } // close year 2+ loop
                                        //   
                                        //   for(p in 1:np){
                                          //     for(y in 1:(ny_proj)){
                                            //       
                                            //       proj_n_p_l_y_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(proj_n_p_a_y_hat[p,1:n_ages,y])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
                                            //       proj_dens_p_y_hat[p,y] = sum((to_vector(proj_n_p_l_y_hat[y,p,1:n_lbins])));
                                            //       
                                            //     }
                                            //   }
                                            // }
                                            
}



