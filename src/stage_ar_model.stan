data {
real m; // natural mortality

int np; // number of patches

int ns; // number of stages

int ny; // number of years

int n_p_s_y[np,ns,ny]; // array of numbers at patch, stage, and year 
 
}

transformed data{
  
  
}

parameters{
  
  real log_sigma_r; // sigma recruits
  
  real<lower = -1, upper = 1> alpha; // autocorrelation term
  
  vector[np] log_mean_recruits; // log mean recruits per patch
  
  matrix[np,ny-1] raw; // array of raw recruitment deviates
  
  real<lower = 1e-6> sigma_obs;

}

transformed parameters{
  
  real sigma_r;
  
  vector[np] mean_recruits;
  
  real n_p_s_y_hat [np,ns,ny]; // array of numbers at patch, stage, and year 

  matrix[np,ny-1] rec_dev; // array of realized recruitment deviates

  sigma_r = exp(log_sigma_r);
  
  mean_recruits = exp(log_mean_recruits);
  
  n_p_s_y_hat[1:np,1:ns,1] = n_p_s_y[1:np,1:ns,1]; // initialze population with "known" values
  
  for (p in 1:np){
  
    for (y in 2:ny){
    
    if (y == 2){
      
      rec_dev[p,y-1]  =  raw[p,1];
     
     } 
     else {
  
      rec_dev[p,y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[p,y-2] - -pow(sigma_r,2)/2) + raw[p,y-1];
  

      } // close ifesle
  
    n_p_s_y_hat[p,1,y] = mean_recruits[p] * exp(rec_dev[p,y-1]);
      
      for (s in 2:ns){
    
        n_p_s_y_hat[p,s,y] = n_p_s_y_hat[p,s-1,y-1] * exp(-m);

      } // close stages
    
    } // close years
  
  } // close patches
  

}

model {
  
log_sigma_r ~ normal(log(.5),.1); // process error prior

alpha ~ normal(0,.25); // autocorrelation prior

sigma_obs ~ normal(0.75, 0.25);

for(y in 2:ny) {
  
    for(p in 1:np){
      
      raw[p,y-1] ~ normal(0,sigma_r); // prior on raw process error
      
      n_p_s_y[p,1:ns,y] ~ multinomial(to_vector(n_p_s_y_hat[p,1:ns,y]) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]))); // fit to proportions at age by patch

     n_p_s_y[p,1,y] ~ neg_binomial_2(n_p_s_y_hat[p,1,y], sigma_obs); // git to mean number of recruits per patch
    
    
    } // close patch loop

  
} // close year loop

  
}

generated quantities{
  
  
  
}

