# Wrapper for standardizing tree ring data
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# link: combines information about the structure and error term in the model as in generalized linear models
# optim: which algorithm should we use to standardize the data?
# post_hoc: apply a post-hoc correction to reduce instability in 3-effect models
# ...: further arguments to control model fitting in optimization algorithm

standardize_tra <- function(tra, model=c("Age", "Time"), link="log", dep_var="Growth", optim="glm", post_hoc=TRUE, return_data=FALSE, make_plots=TRUE, ...)
{
  
  # Exception handling
  if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra<=0], na.rm=TRUE) > 0))
  {
    # Raise a warning if negative values found for multiplicative models
    if (link == "log")
    {
      stop("Zero or negative values cannot be use. Estimated effects will not be stable.")
    }
  }  
  
  # Fitting the model
  if(optim == "alternate")
  {
    effects <- standardize_alternate(tra, model, link, dep_var, ...)
  }
  else if(optim == "sequential")
  {
    effects <- standardize_sequential(tra, model, link, dep_var, ...)
  }
  else if(optim == "glm")
  {
   effects <- standardize_glm(tra, model, link, dep_var, ...)
  }
 
  else if(optim == "gam")
  {
    effects <- standardize_gam(tra, model, link, dep_var, ...)
  }
  
  # Check for 3 effect model
  if (("Tree" %in% model) & ("Time" %in% model) & ("Age" %in% model))
  {
    if (post_hoc)
    {
      effects <- post_hoc_intercession(effects, tra, link)
      warning("Tree-time-age model selected. Post-hoc effect selection was used to stabilize parameter estimates.")
    } else {
      warning("Tree-time-age model selected. Parameter estimates will be unreliable. Consider using post-hoc effect selection.")
    }
  }
  
  # Standardization complete
  print("Standardization complete")
  
  # Make sure elements of effects are in the right order
  effects <- sort_effects(effects, tra, dep_var)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, link, dep_var)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, link,  dep_var)
  print("Model fit computed")  
  
  # Seperate predicted and residuals from fit data
  data <- list(original=tra, predicted=fit$predicted, residuals=fit$residuals)
  fit <- fit[-which(names(fit)=="predicted" | names(fit)=="residuals")]
  
  # Make plots if needed
  if (make_plots){
      plots <- make_standardization_plots(effects, data, link, dep_var)
  
      # Display the plots immediately
      print("Plots constructed")        
      print(plots)
  }
  
  # Record model fitting settings
  settings <- list(model=model, link=link, optim=optim,  dep_var=dep_var)
                   
  # Compile and output all relevant information
  out <- list(effects=effects, fit=fit, settings=settings)
  
  if (return_data){
    out <- c(out, list(data=data))
  }
  
  if(make_plots){
    out <- c(out, list(plots=plots))    
  }
  
  return (out)
}