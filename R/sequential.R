# Sequential standardization ####
# Used in traditional regional curve standardization or flat detrending
# effect_order: the order in which effects are sequentially estimated
standardize_sequential <- function(tra, model=c("Age", "Time"), link="log")
{
  
  # Convert information about link function to type of mean and form used
  if(link=="log"){
    error <- "lnorm"
    form <- "multiplicative"
  } else if (link=="identity"){
    error <- "norm"
    form <- "additive"    
  }
  
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- model
  
  # Create storage for the estimated effects
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  for (id in model){
    # Estimate an effect    
    effects[[id]] <- est_effect(working_tra, id, mean_type)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[id]], effect, link)
  }
  
  return (effects)   
}
