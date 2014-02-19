# Wrapper for standardizing tree ring data ####
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# form: are the effects added together or multiplied together
# error: is the data drawn from a normal additive PDF or a multiplicative log-normal PDF
# method: which algorithm should we use to standardize the data?

standardize_tra <- function(tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", post_hoc=TRUE, ...)
{
  
  # Exception handling
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

 
  if (method=="likelihood")
  {
    out <- standardize_likelihood(tra, model, form, error, ...)
  }
  else if(method == "least_squares")
  {
    out <- standardize_least_squares(tra, model, form, error, ...)
  }
  else if(method == "sfs")
  {
    #out <- standardize_sfs(tra, model, form, error, ...)
    out <- standardize_tsfs(tra, model, form, error, ...)
  }
  else if(method == "rcs")
  {
    out <- standardize_rcs(tra, model, form, error, ...)
  }
  else if(method == "gam")
  {
    out <- standardize_gam(tra, model, form, error, ...)
  }
  
  # Check for 3 effect model
  if (sum(unlist(model))==3)
  {
    if (post_hoc)
    {
      out$effects <- post_hoc_intercession(out$effects, out$tra, form)
      warning("Three effect model selected. Post-hoc selection was used to stabilize parameter estimates.")
    } else {
      warning("Three effect model selected. Parameter estimates are wildly unreliable. Consider using post-hoc selection.")
    }
  }
  
  
  return (out)
}

# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, factor.dim, mean_type="arithmetic")
{  
  id <- factor.dim + 1
  
  estimate_effect <- function(id_i)
  {
    data <- tra[tra[[id]]==id_i, "G"]
    if (mean_type=="geometric"){
      est_effect <- geomMean(data) 
    } else {
      est_effect <- mean(data, na.rm=TRUE)
    }
    return(est_effect)
  }
  
  est_effect <- sapply(levels(tra[[id]]), estimate_effect)
    
  return (est_effect)
}

# Remove an effect from a tree ring array
remove_effect <- function (tra, effect, factor.dim, form="multiplicative")
{
  
  id <- factor.dim + 1
      
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    relevant_rows <- tra[[id]]==effect_id
    if (form=="additive")
    {
      removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] - effect[effect_id]
    }
    else
    {
      removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] / effect[effect_id]
    }
  }
  
  return (removed_tra)
  
}

# Add dummy effect vectors if some are missing
pad_effects <- function(effects, tra, form="multiplicative")
{
  # Set the value to fill dummy coefficients with
  if (form=="multiplicative"){
    na.value <- 1
  } else {
    na.value <- 0
  }
  
  # Initialize dummy effects lists
  new_effects <- list(I=NA, T=NA, A=NA)
  
  # Fill empty values  
  for(i in c("I", "T", "A")){
    
    j <- c(I="i", T="t", A="a")[i]
    
    if (length(effects[[i]] > 0)){
      new_effects[[i]] <-  effects[[i]]
      
      # Fill in dummy levels
      if(length(effects[[i]]) < nlevels(tra[[j]]))
      {
        missing_names <- levels(tra[[j]])[!(levels(tra[[j]]) %in% names(effects[[i]]))]
        missing_effects <- rep(na.value, times=length(missing_names))
        names (missing_effects) <- missing_names
        new_effects[[i]] <- c(effects[[i]], missing_effects)
      }
      
    } else {
      tra_dim <- which(c("I", "T", "A")==i)
      
      new_effects[[i]] <- rep.int(na.value, nlevels(tra[[tra_dim + 1]]))
      effect_names <- levels(tra[[tra_dim + 1]])
      if(i=="T" | i=="A")
      {
        effect_names <- as.numeric(as.character(effect_names))
      }
      names(new_effects[[i]]) <- sort(effect_names)
    }
  }
  
  return(new_effects)
}

# Correctly order effect vectors
sort_effects <- function(effects, tra)
{
  # Only ascending sorting as no "standard" order is retained
  sorted_effects <- list()
  
  # Sort I
  effect_names <- names(effects$I)
  sorted_effects$I <- effects$I[sort(effect_names)]
  
  # Sort T and A
  for (j in c("T", "A"))
  {
    effect_names <- as.numeric(names(effects[[j]]))
    sorted_effects[[j]] <- effects[[j]][as.character(sort(effect_names))]
  }
  
  return(sorted_effects)
}

# Rescale effect vectors to canonical form
rescale_effects <- function (effects, form="multiplicative")
{
  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (form=="multiplicative")
  {
    effects <- lapply(effects, log)
  }
  
  mean_effects <- lapply(effects, mean, na.rm=T)  
  
  # Scale I and T to mean of 0
  # Scale A so sum of effects stays the same
  
  # If A is missing, leave effects uncscaled
  if(!sum(!is.null(effects[[3]])))
  {
    return (effects)
  }
  
  # I
  if (sum (!is.null(effects[[1]]))){
    effects[[1]] <- effects[[1]]-mean_effects[[1]]
    effects[[3]] <- effects[[3]]+mean_effects[[1]]
  }
  
  # T
  if (sum (!is.null(effects[[2]]))){
    effects[[2]] <- effects[[2]]-mean_effects[[2]]
    effects[[3]] <- effects[[3]]+mean_effects[[2]]
  }
  
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
  }
  
  return (effects)
} 

# Model fit statistics ####

# Compute all the relevant model fit statistics for fixed-effects standardization
model_fit_tra <- function(effects, tra, model, form, error, method="sfs", k=NA)
{
  
  fit <- list()
  
  # Need to compute appropriate model fit stats for null models
  if (sum(unlist(model))==0)
  {
    # Predictions of the null model are the null value
    fit$predicted <- tra
    
    if (form=="additive")
    {
      fit$predicted$G <- 0
    } else
    {
      fit$predicted$G <- 1
    }
    
    # Residuals of the null model are the observed data
    fit$residuals <- tra
    
  } else {
    fit$predicted <- predicted_tra(effects, tra, form)
    fit$residuals <- residuals_tra (tra, fit$predicted, error)
  }
  
  fit$n <- n_tra(tra)
  if (method=="gam")
  {
    fit$k <- k
  } else{
    fit$k <- k_tra(tra, model)
  }
  fit$sigma <- sigma_tra(fit$residuals, error)
  
  fit$tss <- tss_tra(tra, error)
  fit$rss <- rss_tra(fit$residuals, error)
  fit$llh <- llh_tra(fit$residuals, error)
  
  fit$Rsq <- Rsq_tra(fit$rss, fit$tss)
  fit$adj.Rsq <- adj.Rsq_tra(fit$rss, fit$tss, fit$n, fit$k)
  
  fit$AIC <- AIC_tra(fit$llh, fit$k)
  fit$AICc <- AICc_tra(fit$llh, fit$k, fit$n)
  fit$BIC <- BIC_tra(fit$llh, fit$k, fit$n)
  
  return(fit)
}

# Predicted values
predicted_tra <- function (effects, tra, form)
{
    predicted <- tra
    
    for (r in 1:nrow(predicted))
    {
      i <- as.character(tra[r, "i"])
      t <- as.character(tra[r, "t"])
      a <- as.character(tra[r, "a"])
      
      if (form=="additive")
      {
        predicted[r, "G"] <- effects$I[i] + effects$T[t] + effects$A[a]
      }
      else
      {
        predicted[r, "G"] <- effects$I[i] * effects$T[t] * effects$A[a]
      }
    }
     
  return(predicted)
}

# Residuals
residuals_tra <- function (tra, predicted, error)
{
  residuals <- tra
  
  if (error=="norm")
  {
    residuals$G <- tra$G - predicted$G
  }
  else
  {
    residuals$G <- tra$G / predicted$G
  }
  
  return (residuals)
}

# Number of data points in the model
n_tra <- function (tra)
{

  n <- nrow(tra)
  
  return(n)
}

# Number of parameters estimated (k)
k_tra <- function (tra, model)
{
 
  # One parameter automatically for estimate of error
  k <- 1
  
  # One parameter is estimated for each index of the effect vectors
  if (model$I)
  {
    k <- k + nlevels(tra$i)
  }
  if (model$T)
  {
    k <- k + nlevels(tra$t)
  }
  if (model$A)
  {
    k <- k + nlevels(tra$a)
  }
  
  # Information about some parameters is lost due to rescaling (dummy variable trap)
  num_effects <- sum(unlist(model))
  
  k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  
  
  return(k)
}

# Calculate sigma, the level of dispersal in the noise PDF as the RMSE (the standard deviation of the residuals)
sigma_tra <- function (residuals, error)
{
 
  val <- residuals$G
  if (error=="lnorm")
  {
    val <- log(val)
  }
  
  # sigma <- sqrt(mean(val^2))
  sigma <- sd(val)
  
  return(sigma)
}

# Total sum of squares
tss_tra <- function (tra, error)
{
  val <- tra$G
  if (error=="lnorm")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val)
  tss <- sum((val-mean_val)^2)
  
  return(tss)
}

# Residual sum of squares
rss_tra <- function (residuals, error)
{
  val <- residuals$G
  if (error=="lnorm")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val)
  rss <-sum((val-mean_val)^2)
    
  return(rss)
}

# Log-likelihood
llh_tra <- function (residuals, error)
{
 
  sigma <- sigma_tra(residuals, error)
  
  val <- residuals$G
  
  # Likelihood is proportional to the probability of observing the data, given the parameters
  if(error=="norm")
  {
    llh <- sum(dnorm(val, sd=sigma, log=TRUE))
  }
  else
  {
    llh <- sum(dlnorm(val, sd=sigma, log=TRUE))
  }
  
  return(llh)
  
}

# R^2
Rsq_tra <- function (rss, tss)
{
  Rsq <- 1-rss/tss
  return(Rsq)
}

# Adjusted R^2
adj.Rsq_tra <- function (rss, tss, n, k)
{
  # p doesn't include estimates of the variance
  p <- k - 1
  adj.Rsq <- 1-(rss/tss)*(n-1)/(n-p-1)
  return(adj.Rsq)
}

# AIC
AIC_tra <- function (llh, k)
{
  AIC <- -2*llh + 2*k
  return(AIC)
}

# AICc
AICc_tra <- function (llh, k, n)
{
  AICc <- -2*llh + 2*k + 2*k*(k+1)/(n-k-1)
  return(AICc)
}

# BIC
BIC_tra <- function (llh, k, n)
{
  BIC <- -2*llh + k*log(n)
  return(BIC)
}