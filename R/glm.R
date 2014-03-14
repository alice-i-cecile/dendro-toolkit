# GLM fixed effects standardization ####

# Main GLM function
standardize_glm <- function (tra, model=c("Time", "Age"), split=NA, link="log", dep_var="Growth", ...)
{
  
  # Construct formula for regression
  growth_formula <- as.formula(make_glm_formula(model, split, dep_var))
  
  print ("Using a generalized linear model to standardize data")
  print (paste("Gaussian family with a", link, "link"))
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

# Formula construction for GLM standardization
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


# Extracting effects for GLM models
extract_effects_glm <- function(growth_model, model, split, link, tra)
{
  # Skeleton effects for relisting coefficients
  skele <- make_skeleton_effects(tra, model, split, link)
    
  # Grab the coefficients from the regression model
  effect_coef <- coef(growth_model)
  
  # Fix the names for split effects
  for (e in model[model %in% split]){
    split_name <- paste(e, "Split", sep="_")
    
    # Grab all names within the relevant effect
    all_e <- str_match_all(names(effect_coef), paste0(e,".*"))
    
    # Find included groups
    raw_groups <- lapply(all_e, str_match, pattern=paste0(split_name, ".*"))
    raw_group_names <- unlist(lapply(raw_groups, str_replace, pattern=split_name, replacement=""))
    inc_groups <- unique(raw_group_names[!is.na(raw_group_names)])    
    
    # Find all groups for the relevant effect
    all_groups <- names(skele[[e]])
    
    # Determine base group by process of elimination
    base_group <- all_groups[!(all_groups %in% inc_groups)]
    
    # Assign group label to base coefficients
    base_id <- sapply(names(effect_coef), str_detect, pattern=e) & !sapply(names(effect_coef), str_detect, pattern=split_name)
    
    names(effect_coef)[base_id] <- paste0(names(effect_coef[base_id]), ":", split_name, base_group)
    
    # Fix labels for dummy levels of original variable 
    all_id <- names(skele[[e]][[base_group]])
    
    inc_id <- names(effect_coef)[base_id]
    inc_id <- sapply(inc_id, str_replace, pattern=e, replacement="")
    inc_id <- sapply(inc_id, str_replace, pattern=":.*", replacement="")
    
    dummy_id <- all_id[!(all_id %in% inc_id)]
    
    dummy_group_id <- !sapply(names(effect_coef), str_detect, pattern=":") & sapply(names(effect_coef), str_detect, pattern=e)
    
    names(effect_coef)[dummy_group_id] <- paste0(e, dummy_id, ":", names(effect_coef)[dummy_group_id])
    
    # Reorder parts of names to proper nested list structure
    str_raw_group <- unlist(str_match(names(effect_coef), paste0(split_name, ".*")))
    str_group <-  str_replace(str_raw_group, split_name, "")
    
    str_raw_id <- str_match(names(effect_coef), paste0(e, ".*", ":"))
    str_id <- str_replace_all(str_raw_id, paste0(e, "|:"), "")
    
    e_id <- str_detect(names(effect_coef), e)
    
    names(effect_coef)[e_id] <- paste(e, str_group[e_id], str_id[e_id], sep=".")
    
  }
  
  # Fix the names of unsplit effects
  for (e in model[!(model %in% split)]){
    e_id <- str_detect(names(effect_coef), e)

    names(effect_coef)[e_id] <- str_replace(names(effect_coef)[e_id], e, paste0(e, "."))
  }
  
  # Put the coefficients into a list of the appropriate size
  boneless_effects <- unlist(skele)
  boneless_effects <- sapply(boneless_effects, function(x){NA})
  
  matching_effect_names <- intersect(names(effect_coef), names(boneless_effects))

  for (n in matching_effect_names){
    boneless_effects[n] <- effect_coef[n]
  }
  
  # Missing effects are set to base levels
  # Returned invisibly by regression
  boneless_effects[is.na(boneless_effects)] <- 0
  
  # Fix structure of effects
  effects <- relist(boneless_effects, skele)
  
  # Correctly combine base values and interaction term for split effects
  for (e in model[model %in% split]){
    groups <- names(effects[[e]])
    base_group <- groups[1]
    effect_b <- effects[[e]][[base_group]]
    
    for (group in groups){
      if (group != base_group){
        base <-  effect_b[names(effects[[e]][[group]])]
        effects[[e]][[group]] <- base + effects[[e]][[group]]
      }
    }
  }
  
  # Account for link
  if (link=="log"){
    effects <- lapply(names(effects), function(x){
      if (x %in% split){
        return(lapply(effects[[x]], exp))
      } else {
        return (exp(effects[[x]]))
      }
    })
    names(effects) <- names(skele)
  }
  
  
  
  return(effects)
    
}