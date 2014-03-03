# GAM fixed effects standardization ####

# Main GLM function
standardize_gam <- function (tra, model=c("Time", "Age"), group_by=NA, link="log", dep_var="Growth", ...)
{
  
  # Clean age info to ensure numeric form
  tra$Age <- as.numeric(as.character(tra$Age))
  
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model, group_by, dep_var))
  
  print ("Using a generalized additive model used to standardize data")
  print (paste("Gaussian family, link is set to", link))
  print (growth_formula)
  
  # Set family to use appropriate link function
  family <- gaussian(link=link)
  
  # Estimate the growth model
  growth_model <- gam(growth_formula, family=family, data=tra, ...)
  
  # Extract estimates of the effect
  effects <- extract_effects_gam(growth_model, model, group_by, link, tra)
  
  # Extract the number of degrees of freedom
  # Custom because of smoothing
  k <- extract_k_gam(growth_model, tra, model, group_by)
  
  return (list(effects=effects, k=k))
  
}

# Formula construction for GAM standardization
make_gam_formula <- function (model, group_by, dep_var)
{
  dep_str <- dep_var
  
  ind_effects <- c("0", model)
  ind_effects <- str_replace(ind_effects, "Age", "s(Age, ...)")
  ind_str <- Reduce(function(...){paste(..., sep="+")}, ind_effects)
  
  # Combine the two sides of the formula
  formula_str <- paste(dep_str, ind_str, sep="~")
  
  return (formula_str)
}


# Extracting effects for glm models
extract_effects_gam <- function(growth_model, model, group_by, link, tra)
{
  # Reset age index to factor
  tra$Age <- as.factor(tra$Age)
  
  # Skeleton effects for relisting coefficients
  skele <- make_skeleton_effects(tra, model, group_by, link)
  
  # Grab the coefficients from the regression model
  effect_coef <- coef(growth_model)
  
  # Fix the names
  for (key_phrase in model){
    names(effect_coef) <- str_replace(names(effect_coef), key_phrase, paste0(key_phrase, "."))
  }
  
  # Put the coefficients into a list of the appropriate size
  boneless_effects <- unlist(skele)
  matching_effect_names <- intersect(names(effect_coef), names(boneless_effects))
  
  for (n in matching_effect_names){
    boneless_effects[n] <- effect_coef[n]
  }
  
  # Missing effects are set to base levels
  # Returned invisibly by regression
  boneless_effects[!(names(boneless_effects) %in% matching_effect_names)] <- 0
  
  # Fix structure of effects
  effects <- relist(boneless_effects, skele)
  
  # Add in effects for age
  if ("Age" %in% model)
  {
    # Find estimates of effect by comparing predictions to predictions from other variables
    Tree_index <- unique(tra$Tree)[1]
    Time_index <- unique(tra$Time)[1]
    Age_levels <- as.numeric(unique(tra$Age))
    
    dummy_data <- data.frame(Tree=Tree_index, Time=Time_index, Age=Age_levels)
    predicted_by_age <- predict(growth_model, dummy_data)
    names(predicted_by_age) <- Age_levels
    
    # Predictions - other effects = age effect
    baseline <- 0
    if ("Tree" %in% model){
      base_Tree <- effect_coef[paste("Tree", Tree_index, sep=".")]
      baseline <- baseline + base_Tree    
    }    
    if ("Time" %in% model){
      base_Time <- effect_coef[paste("Time", Time_index, sep=".")]
      baseline <- baseline + base_Time
    }      
    
    effects$Age <- predicted_by_age - baseline
    
  }
  
  # Account for link
  if (link=="log"){
    effects <- lapply(effects, exp)
  }
  
  return(effects)
  
}

# Number of parameters estimated (k)
# Custom for GAM models
# Borrows from k_tra()
extract_k_gam <- function (growth_model, tra, model, group_by)
{
  
  # One parameter automatically for estimate of link
  k <- 1
  
  # One parameter is estimated for each index of the effect vectors
  if ("Tree" %in% model)
  {
    k <- k + nlevels(tra$Tree)
  }
  if ("Time" %in% model)
  {
    k <- k + nlevels(tra$Time)
  }
  if ("Age" %in% model)
  {
    k <- k + summary(growth_model)$edf
  }
  
  # Inlinkation about some parameters is lost due to rescaling (dummy variable trap)
  num_effects <- length(model)
  
  k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  
  
  return(k)
}
