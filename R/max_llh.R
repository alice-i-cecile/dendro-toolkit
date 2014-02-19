# Maximum likelihood fixed effects standardization ####
# Fits models using maximum likelihood
# Searches for solutions with various optimizing algorithms

standardize_mle <- function(tra, model=c("Time", "Age"), form="multiplicative", error="lnorm", method="CG" ...)
{
  
  # Create storage for the estimated effects
  effects <-vector(mode="list", length=3)
  names (effects) <- c("Tree","Time","Age")
  
  # Determine starting estimates for the effects
  for (i in names(tra)[2:ncol(tra)]){
    dim_i <- nlevels(tra[[i]])
    
    if (form=="additive")
    {
      effects[[i]] <-  rep.int(0,  dim_i)
    } else
    {
      effects[[i]] <-  rep.int(1,  dim_i)
    }
    
    names(effects[[i]]) <- levels(tra[[i]])
  }
  
  # Likelihood function, optimize this!
  fes_likelihood <- function(effects)
  {
    # Find the predicted values
    predicted <- predicted_tra(effects, tra, form)
    
    # Find the residuals
    residuals <- residuals_tra(tra, predicted, error)
    
    # Find the likelihood
    llh <- llh_tra(residuals, error)
    
    return(llh)
  }
  
  # Optimize the model
  # Use the negative likelihood because optim() is a minimizer
  mle_solution <- optim(effects0, -fes_likelihood, ...)
  
  # Report the optimizer's output
  print(mle_solution[2:5])
    
  return(effects)
}