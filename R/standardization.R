# Wrapper for standardizing tree ring data
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# link: combines information about the structure and error term in the model as in generalized linear models
# optim: which algorithm should we use to standardize the data?
# post_hoc: apply a post-hoc correction to reduce instability in 3-effect models

standardize_tra <- function(tra, model=c("Age", "Time"), link="log", optim="sfs", post_hoc=TRUE, ...)
{
  
  # Exception handling ####
  if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra<=0], na.rm=TRUE) > 0))
  {
    # Raise a warning if negative values found for multiplicative models
    if (link == "log")
    {
      stop("Zero or negative values cannot be use. Estimated effects will not be stable.")
    }
  }  

  # Model fitting
  if (optim=="likelihood")
  {
    effects <- standardize_likelihood(tra, model, link, ...)
  }
  # TODO
  #else if(optim == "glm")
  #{
  #  effects <- standardize_glm(tra, model, link, ...)
  #}
  else if(optim == "alternate")
  {
    effects <- standardize_alternate(tra, model, link, ...)
  }
  else if(optim == "rcs")
  {
    effects <- standardize_sequential(tra, model, link, ...)
  }
  else if(optim == "gam")
  {
    effects <- standardize_gam(tra, model, link, ...)
  }
  
  # Check for 3 effect model
  if (("Tree" %in% model) & ("Time" %in% model) & ("Age" %in% model))
  {
    if (post_hoc)
    {
      effects <- post_hoc_intercession(out$effects, out$tra, link)
      warning("Tree-time-age model selected. Post-hoc effect selection was used to stabilize parameter estimates.")
    } else {
      warning("Tree-time-age model selected. Parameter estimates are wildly unreliable. Consider using post-hoc effect selection.")
    }
  }
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, link, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, link)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, link)
  
  # Record model fitting settings
  settings <- list(model=model, link=link, optim=optim
                   
  # Compile and output all relevant information               
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
}