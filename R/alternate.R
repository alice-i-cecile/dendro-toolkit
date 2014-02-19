# Ssignal-free standardization ####
# Cleans up signal-free regional curve standardization algorithm and allows expansion to N dimensions

standardize_sfs <- function (tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", cor_threshold=0.999999)
{
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Create storage for the estimated effects
  effects <-vector(mode="list", length=3)
  names (effects) <- c("I","T","A")
  
  # Dummy starting effects
  for (i in 1:3){
    dim_i <- nlevels(tra[[i+1]])
    
    if (form=="additive")
    {
      effects[[i]] <-  rep.int(0,  dim_i)
    } else
    {
      effects[[i]] <-  rep.int(1,  dim_i)
    }

    names(effects[[i]]) <- levels(tra[[i+1]])
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- names(model[unlist(model)==TRUE])
  name_dim_dict <- c(I="I", T="T", A="A")
  effect_order <- match(inc_effects, name_dim_dict)
  
  # Loop controls
  
  # Don't attempt to fit null models
  if (length(inc_effects)==0)
  {
    converged <- TRUE 
  } else 
  {
    converged <- FALSE
  }
  iteration <- 0
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    
    for (j in effect_order){
      
      # Estimate the effects across each dimension
      est_j <- est_effect(working_tra, j, mean_type)
      
      # Remove them from the signal-free data
      working_tra <- remove_effect (working_tra, est_j, j, form)
      
      # Combine them with previously determined effects for that dimension
      if (form == "additive")
      {
        effects[[j]] <- effects[[j]] + est_j
      } else
      {
        effects[[j]] <- effects[[j]] * est_j
      }
    }
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal    
    if (error=="norm"){
      conv_cor <- cor(working_tra$G, last_tra$G)       
    } else {
        conv_cor <- cor(log(working_tra$G), log(last_tra$G))
    }
    
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, method="sfs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
}