# Wrapper for standardizing tree ring data
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# link: combines information about the structure and error term in the model as in generalized linear models
# optim: which algorithm should we use to standardize the data?
# post_hoc: apply a post-hoc correction to reduce instability in 3-effect models
# ...: further arguments to control model fitting in optimization algorithm

standardize_tra <- function(tra, model=c("Age", "Time"), group_by=NA, link="log", dep_var="Growth", optim="glm", ci_size=0.95, post_hoc=TRUE, return_data=FALSE, make_plots=TRUE, show_plots=TRUE, ...)
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
  
  # Data cleaning
  # Ensure groups and ids are formatted as factors
  for (id in model){
    tra[[id]] <- as.factor(tra[[id]])
  }
  if (!is.na(group_by)){
    for (g in group_by){
      cname <- paste(g, "Group", sep="_")
      tra[[cname]] <- as.factor(tra[[cname]])
    }
  }

  # Fitting the model
  if(optim == "alternate")
  {
    effects <- standardize_alternate(tra, model, group_by, link, dep_var, ...)
  }
  else if(optim == "sequential")
  {
    effects <- standardize_sequential(tra, model, group_by, link, dep_var, ...)
  }
  else if(optim == "glm")
  {
   effects <- standardize_glm(tra, model, group_by, link, dep_var, ...)
  }
 
  else if(optim == "gam")
  {
    results <- standardize_gam(tra, model, group_by, link, dep_var, ...)
    effects <- results$effects
    k <- results$k
  }
  
  # Standardization complete
  print("Standardization complete")
  
  # Make sure elements of effects are in the right order
  effects <- sort_effects(effects, tra, group_by)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, link, group_by)
  
  # Compute model fit statistics
  if (optim=="gam"){
    fit <- model_fit_tra (effects, tra, model, group_by, link, dep_var, optim, k)
  } else {
    fit <- model_fit_tra (effects, tra, model, group_by, link, dep_var, optim)
  }
  print("Model fit computed")  
  
  # Compute standard errors of estimated coefficients
  se <- est_se(fit$residuals, model, group_by, link, dep_var)
  print("Standard errors of estimates computed")
  
  # Seperate predicted and residuals from fit data
  dat <- list(original=tra, predicted=fit$predicted, residuals=fit$residuals)
  fit <- fit[-which(names(fit)=="predicted" | names(fit)=="residuals")]
  
  # Make plots if needed
  if (make_plots){
      plots <- make_standardization_plots(effects, se, dat, group_by, link, dep_var, ci_size)
  
      # Display the plots immediately
      print("Plots constructed")        
      if (show_plots){
        print(plots)
      }
  }
  
  # Record model fitting settings
  settings <- list(model=model, group_by=group_by, link=link, optim=optim,  dep_var=dep_var)
                   
  # Compile and output all relevant information
  out <- list(effects=effects, se=se, fit=fit, settings=settings)
  
  if (return_data){
    out <- c(out, list(dat=dat))
  }
  
  if(make_plots){
    out <- c(out, list(plots=plots))    
  }
  
  return (out)
}