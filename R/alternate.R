# Alternate optimization ####
# Cleans up signal-free regional curve standardization algorithm and allows expansion to N dimensions

standardize_alternate <- function (tra, model=c("Time", "Age"), group_by=NA, link="log", dep_var="Growth", cor_threshold=0.999999, max_iter=25)
{
  
  # Construct skeleton effects to scaffold algorithm
  effects <- make_skeleton_effects(tra, model, group_by)
  
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
  
  while (!converged & iteration < max_iter){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    if (iteration == max_iter)
    {
      print ("Maximum number of iterations reached.")
    }
    
    # Determine effect order from order in which model is given  
    for (id in model){
      if (id %in% group_by)
      {
        groups <- levels(tra[[paste(id, "Group", sep="_")]])
        for (group in groups)
        {
          # Estimate the effects across each dimension
          current_est <- est_effect(working_tra, id, link, dep_var, group)
          
          # Remove the effect
          working_tra <- remove_effect(working_tra, current_est, id, link, dep_var, group)
          
          # Combine them with previously determined effects for that dimension
          if (link == "identity")
          {
            effects[[id]][[group]] <- effects[[id]][[group]] + current_est
          } else
          {
            effects[[id]][[group]] <- effects[[id]][[group]] * current_est
          }
        }
      } else 
      {
        # Estimate the effects across each dimension
        current_est <- est_effect(working_tra, id, link, dep_var)
        
        # Remove them from the signal-free data
        working_tra <- remove_effect (working_tra, current_est, id, link, dep_var)
        
        # Combine them with previously determined effects for that dimension
        if (link == "identity")
        {
          effects[[id]] <- effects[[id]] + current_est
        } else
        {
          effects[[id]] <- effects[[id]] * current_est
        }
      }
    }
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal    
    if (link=="identity"){
      conv_cor <- cor(working_tra[[dep_var]], last_tra[[dep_var]])       
    } else {
        conv_cor <- cor(log(working_tra[[dep_var]]), log(last_tra[[dep_var]]))
    }
        
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  return (effects)
  
}