# GAM fixed effects standardization ####

# Main GLM function
standardize_gam <- function (tra, model=c("Time", "Age"), link="log", ...)
{
  
  # Clean age info to ensure numeric form
  tra$Age <- as.numeric(as.character(tra$Age))
  
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model))
  
  print ("Using a generalized additive model used to standardize data")
  print (paste("Gaussian family, link is set to", link))
  print (growth_formula)
  
  # Set family to use appropriate link function
  family <- gaussian(link=link)
  
  # Estimate the growth model
  growth_model <- gam(growth_formula, family=family, data=tra, ...)
  
  # Extract estimates of the effect
  # Correct way
  effects <- extract_effects_gam(growth_model, model, link, tra)
  
  return (effects)
  
}

# Formula construction for GAM standardization
make_gam_formula <- function (model)
{
  dep.str <- "Growth"
  
  ind_effects <- c("0", model)
  ind_effects <- str_replace(ind_effects, "Age", "s(Age, ...)")
  ind.str <- Reduce(function(...){paste(..., sep="+")}, ind_effects)
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}


# Extracting effects for glm models
extract_effects_gam <- function(growth_model, model, link, tra)
{
  # Reset age index to factor
  tra$Age <- as.factor(tra$Age)
  
  # Skeleton effects for relisting coefficients
  skeleton_effects <- vector(mode="list", length=length(model))
  names(skeleton_effects) <- model
  
  for (i in model){
    dim_i <- nlevels(tra[[i]])
    skeleton_effects[[i]] <-  rep.int(NA,  dim_i)
    names(skeleton_effects[[i]]) <- levels(tra[[i]])
  }
  
  # Grab the coefficients from the regression model
  effect_coef <- coef(growth_model)
  
  # Fix the names
  for (key_phrase in model){
    names(effect_coef) <- str_replace(names(effect_coef), key_phrase, paste0(key_phrase, "."))
  }
  
  # Put the coefficients into a list of the appropriate size
  boneless_effects <- unlist(skeleton_effects)
  matching_effect_names <- intersect(names(effect_coef), names(boneless_effects))
  
  for (n in matching_effect_names){
    boneless_effects[n] <- effect_coef[n]
  }
  
  # Missing effects are set to base levels
  # Returned invisibly by regression
  boneless_effects[is.na(boneless_effects)] <- 0
  
  # Fix structure of effects
  effects <- relist(boneless_effects, skeleton_effects)
  
  # Add in effects for age
  if ("Age" %in% model)
  {
    # Find estimates of effect by comparing predictions to predictions from other variables
    Tree_index <- levels(tra$Tree)[1]
    Time_index <- levels(tra$Time)[1]
    Age_levels <- as.numeric(levels(tra$Age))
    
    dummy_data <- data.frame(Tree=Tree_index, Time=Time_index, Age=Age_levels)
    predicted_by_age <- predict(growth_model, dummy_data)
    names(predicted_by_age) <- Age_levels
    
    # Predictions - other effects = age effect
    baseline <- 0
    if ("Tree" %in% model){
      base_Tree <- effect_coef[paste("Tree", Tree_index, sep=".")]
      baseline <- baseline + base_Tree    }    
    if ("Time" %in% model){
      base_Time <- effect_coef[paste("Time", Time_index, sep=".")]
      baseline <- baseline + base_Time
    }      
    
    effects$Age <- predicted_by_age - baseline
    
  }
  
  return(effects)
  
}