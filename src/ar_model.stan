data{
int<lower=1> len_t; // the number of time points

int<lower=1> len_t_proj; // the number of time points


int<lower=1> len_i; // number of patches 

int<lower=0> y[len_i, len_t]; // defining y as an array of integers with patches as rows and years as columns

// data inputs
vector<lower=0>[len_i] z0; // vector of starting pop. values, one per patch 

vector<lower=0>[len_i] proj_init; // vector of starting pop. values, one per patch 


} 

//transformed data{
//}

parameters{ 

real<lower=0> sigma_r; // recruitment deviates CV 

real<lower = 1e-3> phi_obs; // observation error 

matrix[len_i,len_t-1] raw; // array of raw recruitment deviates

real<lower = -.999, upper = .999> alpha; // autocorrelation parameter

} 

transformed parameters {

matrix[len_i, len_t] y_hat; // estimated numbers in patch i and timestep t

matrix[len_i,len_t-1] rec_dev; // array of realized recruitment deviates

y_hat[1:len_i,1] = z0;

for (t in 2:len_t){
  
  if (t == 2){
  rec_dev[1:len_i,t-1]  =  raw[1:len_i,1];
  } else {
  
  // rec_dev[1:len_i,t-1] = alpha * rec_dev[1:len_i,t-2] + sqrt(1 - pow(alpha,2)) * raw[1:len_i,t-1]; // recruitment deviates as an autocorrelated process

  rec_dev[1:len_i,t-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[1:len_i,t-2] - -pow(sigma_r,2)/2) + raw[1:len_i,t-1];
  

  }
  
  y_hat[1:len_i,t] = y_hat[1:len_i,t-1] .* exp(rec_dev[1:len_i,t-1]); // calculate population in each patch

} // close time loop

} // close transformed parameters

model{

phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu)

sigma_r ~ normal(.5,.25); // process error prior

// alpha ~ beta(.5,2);

alpha ~ normal(0,.25);

// observation model
for(t in 2:len_t) {
  
    for(i in 1:len_i){
      
      raw[i,t-1] ~ normal(-pow(sigma_r,2)/2,sigma_r); // prior on raw process error
      
      
      y[i,t] ~ neg_binomial_2(y_hat[i,t], phi_obs); // this version of neg binom has a more familiar form

    }

  
}

}

generated quantities{
  
matrix[len_i, len_t] pp_y_hat; // pp estimated numbers in patch i and timestep t

matrix[len_i, len_t_proj] pp_proj_y_hat; // pp projected estimated numbers in patch i and timestep t

matrix[len_i, len_t_proj - 1] raw_proj; // raw projected deviates

matrix[len_i, len_t_proj - 1] rec_dev_proj; //pp projected rec devs

real tmp;

pp_y_hat[1:len_i,1] = z0;

// generate posterior predictives for training data
for (i in 1:len_i){
  
  for (t in 2:len_t){
    
   pp_y_hat[i,t] =  neg_binomial_2_rng(y_hat[i,t] + 1, phi_obs);
    
  }
  
}

// generate posterior predictives for the testing data


pp_proj_y_hat[1:len_i,1] = proj_init; // initiate projection with fixed observation

for (i in 1:len_i){
  
  for (t in 2:len_t_proj){
    
      if (t == 2){
        
        rec_dev_proj[i,t-1]  = normal_rng(-pow(sigma_r,2)/2,sigma_r); // generate initial recruitment deviate
        
        } else {
  
  rec_dev_proj[i,t-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev_proj[i,t-2] - -pow(sigma_r,2)/2) + normal_rng(-pow(sigma_r,2)/2,sigma_r); // generate autoregressive recruitment deviates
  

  }
  
  tmp = pp_proj_y_hat[i,t-1] * exp(rec_dev_proj[i,t-1]) + 1; // calculate population in each patch
    
  pp_proj_y_hat[i,t] =  neg_binomial_2_rng(tmp, phi_obs);
  
    // pp_proj_y_hat[i,t] = tmp;
  
  } // close projection time loop
  
} // close projection space loop

  
} // close generated quantities
