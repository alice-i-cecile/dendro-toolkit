
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