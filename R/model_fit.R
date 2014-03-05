
# Master function ####

# Compute all the relevant model fit statistics for fixed-effects standardization
model_fit_tra <- function(effects, tra, model, group_by=NA, link, dep_var, optim, k=NA)
{
  
  fit <- list()
  
  # Need to compute appropriate model fit stats for null models
  if (length(model)==0)
  {
    # Predictions of the null model are the null value
    fit$predicted <- tra
    
    if (link=="identity")
    {
      fit$predicted[[dep_var]] <- 0
    } else
    {
      fit$predicted[[dep_var]] <- 1
    }
    
    # Residuals of the null model are the observed data
    fit$residuals <- tra
    
  } else {
    fit$predicted <- predicted_tra(effects, tra, model, group_by, link, dep_var)
    fit$residuals <- residuals_tra (tra, fit$predicted, link, dep_var)
  }
  
  fit$n <- n_tra(tra)
  if (optim=="gam")
  {
    fit$k <- k
  } else{
    fit$k <- k_tra(tra, model, group_by)
  }
  fit$sigma_sq <- sigma_sq_tra(fit$residuals, link, dep_var)
  
  fit$tss <- tss_tra(tra, link, dep_var)
  fit$rss <- rss_tra(fit$residuals, link, dep_var)
  fit$llh <- llh_tra(fit$residuals, link, dep_var)
  
  fit$Rsq <- Rsq_tra(fit$rss, fit$tss)
  fit$Adj.Rsq <- adj.Rsq_tra(fit$rss, fit$tss, fit$n, fit$k)
  
  fit$AIC <- AIC_tra(fit$llh, fit$k)
  fit$AICc <- AICc_tra(fit$llh, fit$k, fit$n)
  fit$BIC <- BIC_tra(fit$llh, fit$k, fit$n)
  
  return(fit)
}

# Predicted and residuals ####

# Predicted values
predicted_tra <- function (effects, tra, model, group_by, link, dep_var)
{
  
  # Set starting values to null value
  predicted <- tra
  predicted[[dep_var]] <- 0
  
  # Transform effects according to link
  if(link=="log"){
    
    for (e in model){
      if (e %in% group_by)
      {
        effects[[e]] <- lapply(effects[[e]], log)
      } else {
        effects[[e]] <- log(effects[[e]])
      }
    }
  }
  
  # Add effects one at a time
  for (r in 1:nrow(predicted))
  {
    for (e in model)
    {
      if (e %in% group_by)
      {
        cname <- paste(e, "Group", sep="_")
        group <- tra[r, cname]
        i <- as.character(tra[r, e])
        
        predicted[r, dep_var] <- predicted[r, dep_var] + effects[[e]][[group]][i]
      } else {
        i <- as.character(tra[r, e])
        
        predicted[r, dep_var] <- predicted[r, dep_var] + effects[[e]][i]
      }
    }
  }
  
  # Untransform data
  if(link=="log"){
    predicted[[dep_var]] <- exp(predicted[[dep_var]])
  }
  
  return(predicted)
}

# Residuals
residuals_tra <- function (tra, predicted, link, dep_var)
{
  residuals <- tra
  
  if (link=="identity")
  {
    residuals[[dep_var]] <- tra[[dep_var]] - predicted[[dep_var]]
  }
  else
  {
    residuals[[dep_var]] <- tra[[dep_var]] / predicted[[dep_var]]
  }
  
  return (residuals)
}

# Data points and parameters ####

# Number of data points in the model
n_tra <- function (tra)
{
  
  n <- nrow(tra)
  
  return(n)
}

# Number of parameters estimated (k)
k_tra <- function (tra, model, group_by=NA)
{
  
  # One parameter automatically for estimate of error term
  k <- 1
  
  # One parameter is estimated for each index of the effect vectors
  for (E in model)
  {
    if (E %in% group_by)
    {
      cname <- paste(E, "Group", sep="_")
      groups <- levels(tra[[cname]])
      
      # Each group adds unique parameters
      # But only count parameters that are actually in that group
      for (group in groups){
        k <- k + length(unique(tra[tra[[cname]]==group, E]))
      }
    } else 
    {
      k <- k + nlevels (tra[[E]])  
    }
  }
  
  # Information about some parameters is lost due to rescaling (dummy variable trap)
  num_effects <- length(model)
  
  k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  
  
  return(k)
}

# Noise and R^2 ####

# Calculate sigma_sq, the level of dispersal in the noise PDF as the RMSE (the standard deviation of the residuals)
sigma_sq_tra <- function (residuals, link, dep_var)
{
  
  val <- residuals[[dep_var]]
  if (link=="log")
  {
    val <- log(val)
  }
  
  # sigma_sq <- sqrt(mean(val^2))
  sigma_sq <- sd(val, na.rm=T)
  
  return(sigma_sq)
}

# Total sum of squares
tss_tra <- function (tra, link, dep_var)
{
  val <- tra[[dep_var]]
  if (link=="log")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val)
  tss <- sum((val-mean_val)^2)
  
  return(tss)
}

# Residual sum of squares
rss_tra <- function (residuals, link, dep_var)
{
  val <- residuals[[dep_var]]
  if (link=="log")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val, na.rm=T)
  rss <-sum((val-mean_val)^2, na.rm=T)
  
  return(rss)
}

# Likelihood and *IC ####

# Log-likelihood
llh_tra <- function (residuals, link, dep_var)
{
  
  sigma_sq <- sigma_sq_tra(residuals, link, dep_var)
  
  if (link=="log"){
    residuals[[dep_var]] <- log(residuals[[dep_var]])
  }
  
  # Likelihood is proportional to the probability of observing the data, given the parameters
  llh <- sum(exp(-residuals[[dep_var]]^2/(2*sigma_sq))/(2*pi*sigma_sq), na.rm=T)
    
  return(llh)
  
}

# R^2
# Follow Nakagawa, S., Schielzeth, H. (2013), A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4: 133–142.
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