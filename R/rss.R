# Maximum likelihood fixed effects standardization ####
# Fits models using maximum likelihood
# Searches for solutions with various optimizing algorithms

standardize_rss <- function(tra, model=c("Time", "Age"), link="log", method="CG", report_optim=FALSE, ...)
{
  
  # Create storage for the estimated effects
  effects_skeleton <- vector(mode="list", length=length(model))
  names(effects_skeleton) <- model
  
  # Dummy starting effects
  for (i in model){
    dim_i <- nlevels(tra[[i]])
    
    if (link=="log")
    {
      effects_skeleton[[i]] <-  rep.int(1,  dim_i)
    } else
    {
      effects_skeleton[[i]] <-  rep.int(0,  dim_i)
    }
    
    names(effects_skeleton[[i]]) <- levels(tra[[i]])
  }
  
  # Likelihood function, optimize this!
  fes_rss <- function(flat_effects)
  {
    effects <- relist(flat_effects, skeleton=effects_skeleton)
    
    # Find the predicted values
    predicted <- predicted_tra(effects, tra, model, link)
    
    # Find the residuals
    residuals <- residuals_tra(tra, predicted, link)
    
    # Find the likelihood
    rss <- rss_tra(residuals, link)
    
    # Use the negative likelihood because optim() is a minimizer
    return(rss)
  }
  
  # Optimize the model
  mle_solution <- optim(unlist(effects_skeleton), fes_rss, method, ...)
  
  # Report the optimizer's output
  if (report_optim){
    print(mle_solution[2:5])
  }  
  # Extract and return effects
  effects <- relist(mle_solution$par, effects_skeleton)
  
  return(effects)
}