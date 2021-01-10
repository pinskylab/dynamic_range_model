data{
int<lower=1> len_t; // the number of time points

int<lower=0> len_i; // number of patches 

int<lower=0> y[len_i, len_t]; // defining y as an array of integers with patches as rows and years as columns

// data inputs
vector<lower=0>[len_i] z0; // vector of starting pop. values, one per patch 

} 

//transformed data{
//}

parameters{ 
row_vector[len_i] log_mean_rec; // average number of recruits per year

real<lower=0> sigma_r; // recruitment deviates CV 

real<lower=0> phi_obs; // observation error 

matrix[len_i,len_t-1] raw; // array of raw recruitment deviates

real<lower = -.999, upper = .999> alpha; // autocorrelation parameter

real log_m; // average mortality, this is really not right... eventually would need to fix M and estimate F

} 

transformed parameters {

matrix[len_i, len_t] y_hat; // estimated numbers in patch i and timestep t

matrix[len_i,len_t-1] rec_dev; // array of realized recruitment deviates

real m;

row_vector[len_i] mean_rec;

mean_rec = exp(log_mean_rec);

m = exp(log_m);

y_hat[1:len_i,1] = z0;

for (t in 2:len_t){
  
  if (t == 2){
  rec_dev[1:len_i,t-1]  =  raw[1:len_i,1];
  } else {
  
  // rec_dev[1:len_i,t-1] = alpha * rec_dev[1:len_i,t-2] + sqrt(1 - pow(alpha,2)) * raw[1:len_i,t-1]; // recruitment deviates as an autocorrelated process

  rec_dev[1:len_i,t-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev[1:len_i,t-2] - -pow(sigma_r,2)/2) + raw[1:len_i,t-1];
  

  }
  y_hat[1:len_i,t] = y_hat[1:len_i,t-1] * exp(-m) + mean_rec * exp(rec_dev[1:len_i,t-1]); // calculate population in each patch

} // close time loop

}

model{

phi_obs ~ normal(0.75, 0.25); // from https://mc-stan.org/docs/2_20/functions-reference/nbalt.html  phi = mu^2 / (sigma^2-mu)

log_m ~ normal(log(0.2),.5); // natural mortality prior

sigma_r ~ normal(.75,.1); // process error prior

log_mean_rec ~ normal(log(10),.1); // prior on mean number of recruits per patch

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
