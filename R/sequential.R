# Sequential standardization ####
# Used in traditional regional curve standardization or flat detrending
# effect_order: the order in which effects are sequentially estimated
standardize_sequential <- function(tra, model=c("Age", "Time"), form="multiplicative", error="lnorm")
{
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- model
  
  # Make a dummy list of effects
  effects <- list("Tree"=NULL, "Time"=NULL, "Age"=NULL)
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  for (id in model){
    # Estimate an effect    
    effects[[id]] <- est_effect(working_tra, id, mean_type)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[id]], effect, form)
  }
  
  return (effects)   
}
