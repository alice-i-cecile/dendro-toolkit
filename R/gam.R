# GAM fixed effects standardization ####

# Main gam function
standardize_gam <- function (tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", sparse=TRUE, age_k=10, ...)
{
  # Confirm that form and error match, otherwise GAMs can't be used
  if (
    (form=="additive" & error=="lnorm")
    |
      (form=="multiplicative" & error=="norm")
  )
  {
    simpleError("Model form and error distribution are an unrealistic match. Generalized additive models cannot be used for this configuration.")
  }
  
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse) 
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    # Regression uses table form data
    simpleError("GAM standardization can only handle sparse tree-ring arrays.")
  }
  
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model))
  
  # max_k determines the maximum flexibility of the spline fitting the age effect
  # Increased max_k greatly increases computation time
  # Absolute maximum flexibility is:
  # max_k <- nlevels(tra$a)
  
  print ("Model initialized.")
  print (growth_formula)
  
  # If form is additive and error is additive and normal:
  # Use GAM regression with a log-link
  if (form == "additive")
  {
    family <- gaussian(link="identity")
  }
  
  # If form is multiplicative and error is multiplicative log-normal:  
  # Use GAM regression with a log-link
  if (form == "multiplicative")
  {
    # Extremely slow, but more correct
    #family <- Gamma(link="log")
    # Log transformed response variable hack
    # Breaks gam-estimated AIC
    tra$G <- log(tra$G)
    family <- gaussian(link="identity")
  }
  
  # Estimate the growth model
  growth_model <- gam (growth_formula, family=family, data=tra, ...)
  print ("Generalized additive model used to standardize data.")
  
  # Extract estimates of the effect
  # Correct way
  #effects <- extract_effects_gam(growth_model, model, form, tra)
  
  # Log-transformation hack
  effects <- extract_effects_gam(growth_model, model, form="additive", tra)
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
    tra$G <- exp(tra$G)
  }
  print ("Effects extracted.")
  
  return (effects)
  
}

# Formula construction for GAM standardization
make_gam_formula <- function (model){
  dep.str <- "G"
  ind.str <- "0"
  
  if(model$I){
    ind.str <- paste(ind.str, "i", sep="+")
  }
  if(model$T){
    ind.str <- paste(ind.str, "t", sep="+")
  }
  if(model$A){
    # Change smoothing terms here
    ind.str <- paste(ind.str, "s(as.numeric(as.character(a)), k=age_k, ...)", sep="+")
  }
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

# Extracting effects for gam models
extract_effects_gam <- function(growth_model, model, form, tra)
{
  
  # Extract coefficients for I and T directly
  raw_effects <- growth_model$coefficients
  
  effects <- list()
  effects$I <- raw_effects[substr(names(raw_effects), 1, 1)=="i"]
  effects$T <- raw_effects[substr(names(raw_effects), 1, 1)=="t"]
  
  names(effects$I) <- substr(names(effects$I), 2, length(names(effects$I)))
  names(effects$T) <- substr(names(effects$T), 2, length(names(effects$T)))
  
  # Find A by process of elimination
  # Generate predicted values
  dummy_data <- data.frame(i=levels(tra$i)[2], t=levels(tra$t)[2], a=levels(tra$a))
  predicted_by_age <- predict(growth_model, dummy_data)
  
  # Remove the known effects of time and individuals
  # Convoluted partial matching code to handle variable name truncation
  if (form=="additive")
  {
    base_line <- 0
    if (model$I)
    {
      base_line <- base_line + effects$I[which(!is.na(pmatch(names(effects$I), levels(tra$i)[2])))]
    }
    if (model$T)
    {
      base_line <- base_line + effects$T[which(!is.na(pmatch(names(effects$T), levels(tra$t)[2])))]
    }
    effects$A <- predicted_by_age - base_line
  } else 
  {
    base_line <- 1
    if (model$I)
    {
      base_line <- base_line * effects$I[which(!is.na(pmatch(names(effects$I), levels(tra$i)[2])))]
    }
    if (model$T)
    {
      base_line <- base_line * effects$T[which(!is.na(pmatch(names(effects$T), levels(tra$t)[2])))]
    }
    effects$A <- predicted_by_age / base_line
  }
  
  names(effects$A) <- levels(tra$a)
  
  return(effects)
}