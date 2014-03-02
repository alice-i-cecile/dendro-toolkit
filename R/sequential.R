# Sequential standardization ####
# Used in traditional regional curve standardization or flat detrending
# effect_order: the order in which effects are sequentially estimated
standardize_sequential <- function(tra, model=c("Age", "Time"), link="log", dep_var="Growth")
{
  
  # Create storage for the estimated effects
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  # Effects are removed in the order they are listed in the model
  for (id in model){
    # Estimate an effect    
    effects[[id]] <- est_effect(working_tra, id, link, dep_var)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[id]], id, link, dep_var)
  }
  
  return (effects)   
}



