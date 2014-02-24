# Maximum likelihood fixed effects standardization ####
# Fits models using maximum likelihood
# Searches for solutions with various optimizing algorithms

standardize_mle <- function(tra, model=c("Time", "Age"), link="log", method="CG", ...)
{
  
  # Create storage for the estimated effects
  effects0 <- vector(mode="list", length=length(model))
  names(effects0) <- model
  
  # Dummy starting effects
  for (i in model){
    dim_i <- nlevels(tra[[i]])
    
    if (link=="log")
    {
      effects[[i]] <-  rep.int(1,  dim_i)
    } else
    {
      effects[[i]] <-  rep.int(0,  dim_i)
    }
    
    names(effects[[i]]) <- levels(tra[[i]])
  }
  
  # Likelihood function, optimize this!
  fes_likelihood <- function(effects)
  {
    # Find the predicted values
    predicted <- predicted_tra(effects, tra, link)
    
    # Find the residuals
    residuals <- residuals_tra(tra, predicted, link)
    
    # Find the likelihood
    llh <- llh_tra(residuals, error)
    
    return(llh)
  }
  
  # Optimize the model
  # Use the negative likelihood because optim() is a minimizer
  mle_solution <- optim(effects0, -fes_likelihood, method, ...)
  
  # Report the optimizer's output
  print(mle_solution[2:5])
    
  return(effects)
}