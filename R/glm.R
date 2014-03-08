# GLM fixed effects standardization ####

# Main GLM function
standardize_glm <- function (tra, model=c("Time", "Age"), split=NA, link="log", dep_var="Growth", ...)
{
  
  # Construct formula for regression
  growth_formula <- as.formula(make_glm_formula(model, split, dep_var))
  
  print ("Using a generalized linear model to standardize data")
  print (paste("Gaussian family, link is set to", link))
  print (growth_formula)
  
  # Set family to use appropriate link function
  family <- gaussian(link=link)
    
  # Estimate the growth model
  growth_model <- glm(growth_formula, family=family, data=tra, ...)
    
  # Extract estimates of the effect
  # Correct way
  effects <- extract_effects_glm(growth_model, model, split, link, tra)
  
  return (effects)
  
}

# Formula construction for GAM standardization
make_glm_formula <- function (model, split, dep_var)
{
  dep_str <- dep_var
  
  if (any(model %in% split)){
    # Code grouping effects as an interaction with no main effects
    m_terms <- model
    for (e in model){
      if (e %in% split){
        m_terms[which(m_terms==e)] <- paste0(e, "*", e, "_Split")
      }
    }
    ind_str <- Reduce(function(...){paste(..., sep="+")}, c("0", m_terms))
  } else {
    ind_str <- Reduce(function(...){paste(..., sep="+")}, c("0", model))
  }
  
  # Combine the two sides of the formula
  formula_str <- paste(dep_str, ind_str, sep="~")
  
  return (formula_str)
}


# Extracting effects for glm models
extract_effects_glm <- function(growth_model, model, split, link, tra)
{
  # Skeleton effects for relisting coefficients
  skele <- make_skeleton_effects(tra, model, split, link)
    
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
  boneless_effects[is.na(boneless_effects)] <- 0
  
  # Fix structure of effects
  effects <- relist(boneless_effects, skele)
  
  # Account for link
  if (link=="log"){
    effects <- lapply(effects, exp)
  }
  
  return(effects)
    
}