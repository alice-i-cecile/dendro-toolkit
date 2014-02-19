# Regional curve standardization ####
# effect_order: the order in which effects are sequentially estimated
standardize_rcs <- function(tra, model=list(I=FALSE, A=TRUE, T=TRUE), form="multiplicative", error="lnorm")
{
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- names(model[unlist(model)==TRUE])
  
  # Make a dummy list of effects
  effects <- list(I=NULL, T=NULL, A=NULL)
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  for (id in inc_effects){
    # Estimate an effect    
    effects[[id]] <- est_effect(working_tra, id, mean_type)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[id]], effect, form)
  }
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, method="rcs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)   
}
