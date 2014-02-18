# Regional curve standardization ####
# effect_order: the order in which effects are sequentially estimated
standardize_rcs <- function(tra, model=list(I=FALSE, A=TRUE, T=TRUE), form="multiplicative", error="lnorm", sparse=TRUE)
{
  
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse)
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- names(model[unlist(model)==TRUE])
  name_dim_dict <- c(I="I", T="T", A="A")
  effect_order <- match(inc_effects, name_dim_dict)
  
  # Make a dummy list of effects
  effects <- list(I=NULL, T=NULL, A=NULL)
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  for (effect in effect_order){
    # Estimate an effect    
    effects[[effect]] <- est_effect(working_tra, effect, mean_type, sparse)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[effect]], effect, form, sparse)
  }
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error, sparse)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="rcs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)   
}
