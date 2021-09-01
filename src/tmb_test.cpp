#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_ARRAY(n_p_l_y);
  
  DATA_MATRIX(a_p_y);
  
  DATA_MATRIX(l_at_a_key);
  
  DATA_SCALAR(m);
  
  DATA_SCALAR(loo);
  
  DATA_SCALAR(sel_delta);
  
  DATA_INTEGER(np);
  
  DATA_INTEGER(ny);
  
  DATA_INTEGER(na);
  
  
  DATA_INTEGER(nl);
  
  DATA_VECTOR(bin_mids);
  
  PARAMETER(log_sigma_obs);
  
  PARAMETER(log_sigma_r);
  
  PARAMETER_VECTOR(log_mean_recruits);
  
  PARAMETER_VECTOR(raw);
  
  PARAMETER(log_p_length_50_sel);
  
  PARAMETER(log_f);
  
  // model
  
  Type sigma_obs = exp(log_sigma_obs);
  
  Type sigma_r = exp(log_sigma_r);
  
  vector<Type> mean_recruits(np);
  
  vector<Type> selectivity_at_bin(nl);
  
  mean_recruits = exp(log_mean_recruits);
  
  Type f = exp(log_f);
    
  Type length_50_sel = loo * exp(log_p_length_50_sel); // Dan made a note to change this sometime
  
  selectivity_at_bin = 1.0 / (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
  
  array<Type> n_p_a_y_hat(np, na, ny);
  
  array<Type> n_p_l_y_hat(np, nl, ny);
  
  Type z = exp(-(f + m));
  
  for(int p = 0; p < np; p++){
    for(int a = 0; a < na; a++){
      if(a==1){
        n_p_a_y_hat(p,a,0) = mean_recruits(p) * exp(raw(0) - pow(sigma_r,2)/2); // initialize age 0 with mean recruitment in every patch
      }
      else{
        n_p_a_y_hat(p,a,0) = n_p_a_y_hat(p,a-1,0) * (z); // initialize population with mean recruitment propogated through age classes with mortality
      }
      
    } // close ages
  } // close patches
  

  // calculate recruitment deviates every year (not patch-specific)
  for (int y = 1; y < ny; y++){

    // describe population dynamics
    for(int p = 0; p < np; p++){
      
      n_p_a_y_hat(p,0,y) = mean_recruits(p) * exp(raw(y-1) - pow(sigma_r,2)/2);
    
    for(int a = 1; a < na; a++){
      
      n_p_a_y_hat(p,a,y) = n_p_a_y_hat(p, a-1, y-1) * (z); // these just grow and die in the patch
      
      
      }
    
    } // close patches 
    
    
  } // close year 2+ loop
  
  
  for(int p = 0; p < np; p++){
    for(int y = 0; y < ny; y++){
      for (l = 0; l < nl; l++){
        for (a = 0; a < na; a++){
          
          
          
          
          n_p_l_y_hat(p,l,y) = ((l_at_a_key.transpose() %*% (n_p_a_y_hat(p,0:n_ages,y)) * selectivity_at_bin); 
                                          
                                          dens_p_y_hat[p,y] = sum((to_vector(n_p_l_y_hat[y,p,1:n_lbins])));
        
        }
        
      }
      

                                       
    }
  }
  
  
  
  return nll;
}
