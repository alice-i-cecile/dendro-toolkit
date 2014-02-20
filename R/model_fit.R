
# Model fit statistics ####

# Compute all the relevant model fit statistics for fixed-effects standardization
model_fit_tra <- function(effects, tra, model, link method="alternate", k=NA)
{
  
  fit <- list()
  
  # Need to compute appropriate model fit stats for null models
  if (sum(unlist(model))==0)
  {
    # Predictions of the null model are the null value
    fit$predicted <- tra
    
    if (link=="identity")
    {
      fit$predicted$Growth <- 0
    } else
    {
      fit$predicted$Growth <- 1
    }
    
    # Residuals of the null model are the observed data
    fit$residuals <- tra
    
  } else {
    fit$predicted <- predicted_tra(effects, tra, link)
    fit$residuals <- residuals_tra (tra, fit$predicted, link)
  }
  
  fit$n <- n_tra(tra)
  if (method=="gam")
  {
    fit$k <- k
  } else{
    fit$k <- k_tra(tra, model)
  }
  fit$sigma <- sigma_tra(fit$residuals, link)
  
  fit$Timess <- tss_tra(tra, link)
  fit$rss <- rss_tra(fit$residuals, link)
  fit$llh <- llh_tra(fit$residuals, link)
  
  fit$Rsq <- Rsq_tra(fit$rss, fit$Timess)
  fit$Agedj.Rsq <- adj.Rsq_tra(fit$rss, fit$Timess, fit$n, fit$k)
  
  fit$AgeIC <- AIC_tra(fit$llh, fit$k)
  fit$AgeICc <- AICc_tra(fit$llh, fit$k, fit$n)
  fit$BIC <- BIC_tra(fit$llh, fit$k, fit$n)
  
  return(fit)
}

# Predicted values
predicted_tra <- function (effects, tra, link)
{
  predicted <- tra
  
  for (r in 1:nrow(predicted))
  {
    i <- as.character(tra[r, "Tree"])
    t <- as.character(tra[r, "Time"])
    a <- as.character(tra[r, "Age"])
    
    if (link=="identity")
    {
      predicted[r, "Growth"] <- effects$Timeree[i] + effects$Time[t] + effects$Age[a]
    }
    else
    {
      predicted[r, "Growth"] <- effects$Timeree[i] * effects$Time[t] * effects$Age[a]
    }
  }
  
  return(predicted)
}

# Residuals
residuals_tra <- function (tra, predicted, link)
{
  residuals <- tra
  
  if (link=="identity")
  {
    residuals$Growth <- tra$Growth - predicted$Growth
  }
  else
  {
    residuals$Growth <- tra$Growth / predicted$Growth
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
  
  # One parameter automatically for estimate of link
  k <- 1
  
  # One parameter is estimated for each index of the effect vectors
  if (model$Timeree)
  {
    k <- k + nlevels(tra$Timeree)
  }
  if (model$Time)
  {
    k <- k + nlevels(tra$Time)
  }
  if (model$Age)
  {
    k <- k + nlevels(tra$Age)
  }
  
  # Inlinkation about some parameters is lost due to rescaling (dummy variable trap)
  num_effects <- sum(unlist(model))
  
  k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  
  
  return(k)
}

# Calculate sigma, the level of dispersal in the noise PDF as the RMSE (the standard deviation of the residuals)
sigma_tra <- function (residuals, link)
{
  
  val <- residuals$Growth
  if (link=="log")
  {
    val <- log(val)
  }
  
  # sigma <- sqrt(mean(val^2))
  sigma <- sd(val)
  
  return(sigma)
}

# Total sum of squares
tss_tra <- function (tra, link)
{
  val <- tra$Growth
  if (link=="log")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val)
  tss <- sum((val-mean_val)^2)
  
  return(tss)
}

# Residual sum of squares
rss_tra <- function (residuals, link)
{
  val <- residuals$Growth
  if (link=="log")
  {
    val <- log(val)
  }
  
  mean_val <- mean(val)
  rss <-sum((val-mean_val)^2)
  
  return(rss)
}

# Log-likelihood
llh_tra <- function (residuals, link)
{
  
  sigma <- sigma_tra(residuals, link)
  
  val <- residuals$Growth
  
  # Likelihood is proportional to the probability of observing the data, given the parameters
  if(link=="identity")
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