# Alternate optimization ####
# Cleans up signal-free regional curve standardization algorithm and allows expansion to N dimensions

standardize_alternate <- function (tra, model=c("Time", "Age"), link="log", cor_threshold=0.999999)
{
  
  # Create storage for the estimated effects
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Dummy starting effects
  for (i in model){
    dim_i <- nlevels(tra[[i]])
    
    if (link=="log")
    {
      effects[[i]] <-  rep.int(1,  dim_i)
    } else
    {
      effects[[i]] <-  rep.int(0,  dim_i)
    }
    
    names(effects[[i]]) <- levels(tra[[i]])
  }
  
  # Don't attempt to fit null models
  if (length(model)==0)
  {
    converged <- TRUE 
  } else 
  {
    converged <- FALSE
  }
  
  # Loop controls
  iteration <- 0
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    
    # Determine effect order from order in which model is given  
    for (id in model){
      
      # Estimate the effects across each dimension
      est_j <- est_effect(working_tra, id, link)
      
      # Remove them from the signal-free data
      working_tra <- remove_effect (working_tra, est_j, id, link)
      
      # Combine them with previously determined effects for that dimension
      if (link == "identity")
      {
        effects[[id]] <- effects[[id]] + est_j
      } else
      {
        effects[[id]] <- effects[[id]] * est_j
      }
    }
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal    
    if (link=="identity"){
      conv_cor <- cor(working_tra$Growth, last_tra$Growth)       
    } else {
        conv_cor <- cor(log(working_tra$Growth), log(last_tra$Growth))
    }
        
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  return (effects)
  
}