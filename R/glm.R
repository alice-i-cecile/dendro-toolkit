# GLM fixed effects standardization ####

# Main GLM function
standardize_glm <- function (tra, model=c("Time", "Age"), link="log", ...)
{
  
  # Construct formula for regression
  growth_formula <- as.formula(make_glm_formula(model))
  
  print ("Using a generalized linear model used to standardize data")
  print (paste("Gaussian family, link is set to", link))
  print (growth_formula)
  
  # Set family to use appropriate link function
  family <- gaussian(link=link)
    
  # Estimate the growth model
  growth_model <- glm(growth_formula, family=family, data=tra, ...)
    
  # Extract estimates of the effect
  # Correct way
  effects <- extract_effects_glm(growth_model, model, link, tra)
  
  return (effects)
  
}

# Formula construction for GAM standardization
make_glm_formula <- function (model)
{
  dep.str <- "Growth"
  ind.str <- Reduce(function(...){paste(..., sep="+")}, c("0", model))
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}


# Extracting effects for glm models
extract_effects_glm <- function(growth_model, model, link, tra)
{
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
  
  return(effects)
    
}