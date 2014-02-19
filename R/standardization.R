# Wrapper for standardizing tree ring data
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# form: are the effects added together or multiplied together
# error: is the data drawn from a normal additive PDF or a multiplicative log-normal PDF
# optim: which algorithm should we use to standardize the data?

standardize_tra <- function(tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", optim="sfs", post_hoc=TRUE, ...)
{
  
  # Exception handling ####
  if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra<=0], na.rm=TRUE) > 0))
  {
    # Raise a warning if negative values found for multiplicative models
    if (form == "multiplicative")
    {
      warning("Zero or negative values found in data when using a multiplicative model. Estimated effects will not be stable.")
    }
    
    # Raise a warning if negative or 0 values found for lognormal errors
    if (error == "lnorm")
    {
      warning("Zero or negative values found in data when using a lognormal error form. Model fitting will fail, as non-positive values can not be generated using a log-normal distribution.")
    }
  }  
  
  # Raise a warning if error family and form do not match
  if (
      (form=="additive" & error=="lnorm")
      |
      (form=="multiplicative" & error=="norm")
  )
  {
    warning("Model form and error distribution are an unrealistic match. Be sure to check the model residuals. Model fitting problems may arise.")
  }

  # Model fitting
  if (optim=="likelihood")
  {
    effects <- standardize_likelihood(tra, model, form, error, ...)
  }
  else if(optim == "least_squares")
  {
    effects <- standardize_least_squares(tra, model, form, error, ...)
  }
  else if(optim == "sfs")
  {
    effects <- standardize_sfs(tra, model, form, error, ...)
  }
  else if(optim == "rcs")
  {
    effects <- standardize_rcs(tra, model, form, error, ...)
  }
  else if(optim == "gam")
  {
    effects <- standardize_gam(tra, model, form, error, ...)
  }
  
  # Check for 3 effect model
  if (sum(unlist(model))==3)
  {
    if (post_hoc)
    {
      effects <- post_hoc_intercession(out$effects, out$tra, form)
      warning("Three effect model selected. Post-hoc effect selection was used to stabilize parameter estimates.")
    } else {
      warning("Three effect model selected. Parameter estimates are wildly unreliable. Consider using post-hoc effect selection.")
    }
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
  settings <- list(model=model, form=form, error=error, optim=optim
                   
  # Compile and output all relevant information               
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
}