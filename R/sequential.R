# Sequential standardization ####
# Used in traditional regional curve standardization or flat detrending
# effect_order: the order in which effects are sequentially estimated
standardize_sequential <- function(tra, model=c("Age", "Time"), group_by=NA, link="log", dep_var="Growth")
{
  
  # Create storage for the estimated effects
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  # Effects are removed in the order they are listed in the model
  for (id in model)
  {
    if (id %in% group_by)
    {
      groups <- levels(tra[[paste(id, "Group", sep="_")]])
      for (group in groups)
      {
        # Estimate an effect    
        effects[[id]][[group]] <- est_effect(working_tra, id, link, dep_var, group)
        
        # Remove the effect
        working_tra <- remove_effect(working_tra, effects[[id]][[g]], id, link, dep_var, group)
      }
    } else 
    {
      # Estimate an effect    
      effects[[id]] <- est_effect(working_tra, id, link, dep_var)
      
      # Remove the effect
      working_tra <- remove_effect(working_tra, effects[[id]], id, link, dep_var)
    }
  }
  
  return (effects)   
}
