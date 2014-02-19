# Maximum likelihood fixed effects standardization ####
# Fits models using maximum likelihood
# Searches for solutions with simulated annealing

# Flattening effects to a vector for use with optim()
flatten_effects <- function(effects, model)
{
  flat_effects <- vector()
  
  if (model$I)
  {
    effects_I <- dummy_effects$I
    names(effects_I) <- paste("I", names(dummy_effects$I), sep="@")
    flat_effects <- c(flat_effects, effects_I)
  }
  if (model$T)
  {
    effects_T <- dummy_effects$T
    names(effects_T) <- paste("T", names(dummy_effects$T), sep="@")
    flat_effects <- c(flat_effects, effects_T)
  }
  if (model$A)
  {
    effects_A <- dummy_effects$A
    names(effects_A) <- paste("A", names(dummy_effects$A), sep="@")
    flat_effects <- c(flat_effects, effects_A)
  }
  return(flat_effects)
}

unflatten_effects <- function (flat_effects)
{
  effects <- list()
  
  id_code <- substr(names(flat_effects), 1, 1)
  names(flat_effects) <- substr(names(flat_effects), 3, length(names(flat_effects)))
  
  effects$I <- flat_effects[id_code=="I"]
  effects$T <- flat_effects[id_code=="T"]
  effects$A <- flat_effects[id_code=="A"]
  
  return(effects)
}

standardize_mle <- function(tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", ...)
{
  
  # Determine starting estimates for the effects
  dummy_effects <- pad_effects(list(), tra, form)
  flat_effects0 <- flatten_effects(dummy_effects, model)
  
  # Likelihood function, optimize this!
  fes_likelihood <- function(flat_effects)
  {
    # Unflatten effects
    effects <- unflatten_effects(flat_effects)
    
    # Pad the effects
    effects <- pad_effects(effects, tra, form)
    
    # Find the predicted values
    predicted <- predicted_tra(effects, tra, form)
    
    # Find the residuals
    residuals <- residuals_tra(tra, predicted, error)
    
    # Find the likelihood
    llh <- llh_tra(residuals, error)
    
    # Return the negative likelihood (optim is a minimizer)
    return(-llh)
  }
  
  # Optimize the model
  mle_solution <- optim(flat_effects0, fes_likelihood, ...)
  
  # Report the optimizer's output
  print(mle_solution[2:5])
  
  # Extract the effects
  effects <- unflatten_effects(mle_solution$par)
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, method="likelihood", ...)
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return(out)
}