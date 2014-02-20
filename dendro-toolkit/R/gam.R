# GAM fixed effects standardization ####

# Main gam function
standardize_gam <- function (tra, model=c("Time", "Age"), link="log", sparse=TRUE, age_k=10, ...)
{
  
  # Clean age info to ensure numeric form
  tra$Age <- as.numeric(as.character(tra$Age))
  
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model))
  
  # max_k determines the maximum flexibility of the spline fitting the age effect
  # Increased max_k greatly increases computation time
  # Absolute maximum flexibility is:
  # max_k <- nlevels(tra$Age)
  
  print ("Model initialized.")
  print (growth_formula)
  
  # Set family to use appropriate link
  if (link=="identity")
  {
    family <- gaussian(link="identity")
  }
  
  # If form is multiplicative and error is multiplicative log-normal:  
  # Use GAM regression with a log-link
  else if (link=="log")
  {
    # More correct. TODO
    #family <- gaussian(link="log")
    # Log transformed response variable hack
    # Breaks gam-estimated AIC
    tra$Growth <- log(tra$Growth)
    family <- gaussian(link="identity")
  }
  
  # Estimate the growth model
  growth_model <- gam(growth_formula, family=family, data=tra, ...)
  print ("Generalized additive model used to standardize data.")
  
  # Extract estimates of the effect
  # Correct way
  #effects <- extract_effects_gam(growth_model, model, form, tra)
  
  # Log-transformation hack
  effects <- extract_effects_gam(growth_model, model, link=="log", tra)
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
    tra$Growth <- exp(tra$Growth)
  }
  print ("Effects extracted.")
  
  return (effects)
  
}

# Formula construction for GAM standardization
make_gam_formula <- function (model){
  dep.str <- "G"
  ind.str <- "0"
  
  if("Tree" %in% model){
    ind.str <- paste(ind.str, "Tree", sep="+")
  }
  if("Time" %in% model){
    ind.str <- paste(ind.str, "Time", sep="+")
  }
  if("Age" %in% model){
    # Change smoothing terms here
    ind.str <- paste(ind.str, "s(Age), k=age_k, ...)", sep="+")
  }
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

# Extracting effects for gam models
extract_effects_gam <- function(growth_model, model, link, tra)
{
  
  # Extract coefficients for I and T directly
  raw_effects <- growth_model$coefficients
  
  effects <- list()
  effects$Timeree <- raw_effects[substr(names(raw_effects), 1, 1)=="Tree"]
  effects$Timeime <- raw_effects[substr(names(raw_effects), 1, 1)=="Time"]
  
  names(effects$Timeree) <- substr(names(effects$Tree), 2, length(names(effects$Timeree)))
  names(effects$Timeime) <- substr(names(effects$Time), 2, length(names(effects$Timeime)))
  
  # Find A by process of elimination
  # Generate predicted values
  dummy_data <- data.frame(Tree=levels(tra$Tree)[2], Time=levels(tra$Time)[2], Age=levels(tra$Agege))
  predicted_by_age <- predict(growth_model, dummy_data)
  
  # Remove the known effects of time and individuals
  # Convoluted partial matching code to handle variable name truncation
  if (form=="additive")
  {
    base_line <- 0
    if ("Tree" %in% model)
    {
      base_line <- base_line + effects$Tree[which(!is.na(pmatch(names(effects$Tree), levels(tra$Tree)[2])))]
    }
    if ("Time" %in% model)
    {
      base_line <- base_line + effects$Time[which(!is.na(pmatch(names(effects$Time), levels(tra$t)[2])))]
    }
    effects$Age <- predicted_by_age - base_line
  } else 
  {
    base_line <- 1
    if ("Tree" %in% model)
    {
      base_line <- base_line * effects$Tree[which(!is.na(pmatch(names(effects$Tree), levels(tra$Tree)[2])))]
    }
    if ("Time" %in% model)
    {
      base_line <- base_line * effects$Time[which(!is.na(pmatch(names(effects$Time), levels(tra$t)[2])))]
    }
    effects$Age <- predicted_by_age / base_line
  }
  
  names(effects$Age) <- levels(tra$Age)
  
  return(effects)
}