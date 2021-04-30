generated quantities{
  
  real pp_n_p_s_y_hat[np, ns, ny];
  
  real pp_proj_n_p_s_y_hat[np, ns, ny_proj];
  
  real tmp[np, ns, ny];
  
  real tmp_proj[np, ns, ny_proj];
  
  real rec_dev_proj[ny_proj - 1];
  
  real T_adjust_proj[np, ny_proj];
  
  for(p in 1:np){
    for(y in 1:ny_proj){
      T_adjust_proj[p,y] = T_dep(sst_proj[p,y], Topt, width); // calculate temperature-dependence correction factor for each patch and year depending on SST 
    }
  }
  
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
      
      for(y in 2:ny){
        
        
        for(p in 1:np){
          
          // project stage 1 (recruitment)
          
          // generate both distribution and scale.. ah so just generate from multinomial based on proportions, then scale up by total numbers, which I think should do it? though also need error in the scale. sigh. 
          
          // multinomial_rng((to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel) / sum(to_vector(n_p_s_y_hat[p,1:ns,y]) .* sel)); 
          // 
          tmp[p, 1, y] = neg_binomial_2_rng(n_p_s_y_hat[p,1,y], phi_obs); // observation error for number of recruits, a little weird since outputting an int into a real, but works. AF: changed from sigma_obs to phi_obs
          
          // derterministically fill in the rest of the numbers at stage
          
          tmp[p,2,y] = tmp[p,1,y-1] * (1 - m) * g + tmp[p,2,y-1] * (1 - m) * (1 - g); 
          
        } // close patch loop for stages 1 and 2
        
        // adult dispersal
        
        // p = 1
        tmp[1,3,y] = tmp[1,2,y-1] * (1 - m) * g + // young adults from last year that grew
        tmp[1,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp[1+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        
        // p = np
        
        tmp[np,3,y] = tmp[np,2,y-1] * (1 - m) * g + // young adults from last year that grew
        tmp[np,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
        tmp[np-1, 3, y-1] * (1 - m) * spill; // dispersal from patch below
        
        // general case for non-edge patches
        for(p in 2:(np-1)){
          tmp[p,3,y] = tmp[p,2,y-1] * (1 - m) * g + // young adults from last year that grew
          tmp[p,3,y-1] * (1 - m) * (1 - 2*spill) + // adults from last year that survived, minus those that dispersed to two patches
          tmp[p-1, 3, y-1] * (1 - m) * spill + //dispersal from patch below
          tmp[p+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
        }
        
        // fill in pp_n_p_s_y_hat
        for(p in 1:np){
          pp_n_p_s_y_hat[p, 1:ns, y] = tmp[p,1:ns,y];
        } // close pp patches
        
      } // close year loop
      
      // fit multinomial to all stages, not really needed except to add noise from sub-sampling
      
      // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
      
      // generate posterior predictives for the testing data aka project recruitment deviates into the future
      
      // left over from when gA/gY was T-dependent
      // for(p in 1:np){
        //   for(y in 1:ny_proj){
          //     gA_proj[p,y] = growth(sst_proj[p, y], Topt, width);
          //     gY_proj[p,y] = growth(sst_proj[p, y], Topt, width);
          //   }
          // }
          
          for(p in 1:np){
            for(y in 2:ny_proj){
              
              if (y == 2){
                rec_dev_proj[y-1]  = rec_dev[y-1]; // initiate projections with last estimated process errors
              }
              else{
                rec_dev_proj[y-1] = -pow(sigma_r,2)/2 + alpha *  (rec_dev_proj[y-2] - -pow(sigma_r,2)/2) + normal_rng(-pow(sigma_r,2)/2,sigma_r); // generate autoregressive recruitment deviates
                
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
          
          for(y in 2:ny_proj){
            
            for(p in 1:np){
              
              // project stage 1 (recruitment)
              tmp_proj[p, 1, y] = neg_binomial_2_rng(mean_recruits * exp(rec_dev_proj[y-1]) * T_adjust_proj[p, y],phi_obs); // added in T-dependence here 
              
              // project other stages
              
              
              tmp_proj[p,2,y] = tmp_proj[p,1,y-1] * (1 - m) * g + tmp_proj[p,2,y-1] * (1 - m) * (1 - g); 
            } // close patch loop
            
            // adult dispersal
            
            // p = 1
            tmp_proj[1,3,y] = tmp_proj[1,2,y-1] * (1 - m) * g + // young adults from last year that grew
            tmp_proj[1,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
            tmp_proj[1+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
            
            
            // p = np
            tmp_proj[np,3,y] = tmp_proj[np,2,y-1] * (1 - m) * g + // young adults from last year that grew
            tmp_proj[np,3,y-1] * (1 - m) * (1 - spill) + // adults from last year that survived, minus those that dispersed to one patch only
            tmp_proj[np-1, 3, y-1] * (1 - m) * spill; // dispersal from patch below
            
            
            // general case for non-edge patches
            for(p in 2:(np-1)){
              tmp_proj[p,3,y] = tmp_proj[p,2,y-1] * (1 - m) * g + // young adults from last year that grew
              tmp_proj[p,3,y-1] * (1 - m) * (1 - 2*spill) + // adults from last year that survived, minus those that dispersed to two patches
              tmp_proj[p-1, 3, y-1] * (1 - m) * spill + //dispersal from patch below
              tmp_proj[p+1, 3, y-1] * (1 - m) * spill; // dispersal from patch above
            }         // close patches 
            
            
            // fit multinomial to all stages: no need for this except to simulate the slight addition of error from subsampling the numbers at age
            
            for(p in 1:np){
              
              pp_proj_n_p_s_y_hat[p, 1:ns, y] = tmp_proj[p,1:ns,y];
              
            } // close patches for pp_proj
            
            // pp_proj_n_p_s_y_hat[p, 1:ns, y] = multinomial_rng(to_vector(tmp[p, 1:ns, y]) / sum(to_vector(tmp[p, 1:ns, y])), 100);
            
          } // close year loop
          
} // close generated quantities
