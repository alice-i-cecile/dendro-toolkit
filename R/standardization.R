# Wrapper for standardizing tree ring data
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# link: combines information about the structure and error term in the model as in generalized linear models
# optim: which algorithm should we use to standardize the data?
# post_hoc: apply a post-hoc correction to reduce instability in 3-effect models
# ...: further arguments to control model fitting in optimization algorithm

standardize_tra <- function(tra, model=c("Age", "Time"), split=NA, link="log", dep_var="Growth", optim="alternate", auto_cluster=FALSE, cluster_criteria="BIC", distance="euclidean", clust="kmeans", ci_size=0.95, return_data=FALSE, make_plots=TRUE, show_plots=TRUE, ...)
{
  
  # Exception handling
  if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra<=0], na.rm=TRUE) > 0))
  {
    # Raise a warning if negative values found for multiplicative models
    if (link == "log")
    {
      stop("Zero or negative values cannot be used. Estimated effects will not be stable.")
    }
  }
  
  # First run of automatic clustering should be done with no splits
  if (auto_cluster){
    original_split <- split
    split <- NA
  }
  
  # Fitting the model
  if(optim == "alternate")
  {
    effects <- standardize_alternate(tra, model, split, link, dep_var, ...)
  }
  else if(optim == "sequential")
  {
    effects <- standardize_sequential(tra, model, split, link, dep_var, ...)
  }
  else if(optim == "glm")
  {
   effects <- standardize_glm(tra, model, split, link, dep_var, ...)
  }
 
  else if(optim == "gam")
  {
    results <- standardize_gam(tra, model, split, link, dep_var, ...)
    effects <- results$effects
    k <- results$k
  }
  
  else if(optim == "direct_search")
  {
    standardize_direct_search(tra, model, split, link, dep_var, ...)
  }
  
  # Standardization complete
  print("Standardization complete")
  
  # Make sure elements of effects are in the right order
  effects <- sort_effects(effects, tra, split)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, link, split)
  
  # Compute model fit statistics
  if (optim=="gam"){
    fit <- model_fit_tra (effects, tra, model, split, link, dep_var, optim, k)
  } else {
    fit <- model_fit_tra (effects, tra, model, split, link, dep_var, optim)
  }
  print("Model fit computed")  
  
  # Compute standard errors of estimated coefficients
  se <- est_se(fit$residuals, model, split, link, dep_var)
  print("Standard errors of estimates computed")
  
  # Seperate predicted and residuals from fit data
  dat <- list(original=tra, predicted=fit$predicted, residuals=fit$residuals)
  fit <- fit[-which(names(fit)=="predicted" | names(fit)=="residuals")]
  
  # Automatic clustering
  if (auto_cluster){
    results <- auto_cluster_tra(tra, resids, fit, model, original_split, link, dep_var, optim, cluster_criteria, ...)
    
    # Use clustered results instead of initial run
    effects <- results$effects
    se <- results$se
    fit <- results$fit
    dat <- results$dat
    
    # Reset split parameter
    split <- original_split
  }
  
  # Make plots if needed
  if (make_plots){
      plots <- make_standardization_plots(effects, se, dat, split, link, dep_var, ci_size)
  
      # Display the plots immediately
      print("Plots constructed")        
      if (show_plots){
        print(plots)
      }
  }
  
  # Record model fitting settings
  settings <- list(model=model, split=split, link=link, optim=optim,  dep_var=dep_var, auto_cluster=auto_cluster, cluster_criteria=cluster_criteria, distance=distance, clust=clust)
                   
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