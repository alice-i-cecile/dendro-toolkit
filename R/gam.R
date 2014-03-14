# GAM fixed effects standardization ####

# Main GLM function
standardize_gam <- function (tra, model=c("Time", "Age"), split=NA, link="log", dep_var="Growth", ...)
{
  
  # Clean age info to ensure numeric form
  if ("Age" %in% model)
  {
    tra$Age <- as.numeric(as.character(tra$Age))
  }
  
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model, split, dep_var))
  
  print ("Using a generalized additive model to standardize data")
  print (paste("Gaussian family with a", link, "link"))
  print (growth_formula)
  
  # Set family to use appropriate link function
  family <- gaussian(link=link)
  
  # Estimate the growth model
  growth_model <- gam(growth_formula, family=family, data=tra, ...)
  
  # Extract estimates of the effect
  effects <- extract_effects_gam(growth_model, model, split, link, tra)
  
  # Extract the number of degrees of freedom
  # Custom because of smoothing
  k <- extract_k_gam(growth_model, tra, model, split)
  
  return (list(effects=effects, k=k))
  
}

# Formula construction for GAM standardization
# Smoothing choices
# No split: s(Age)
# Simple factor smooth: s(Age, by=Age_Split)
# Factor smooth basis: s(Age, Age_Split, bs="fs")

make_gam_formula <- function (model, split, dep_var)
{
  dep_str <- dep_var
  
  
  if (any(model %in% split)){
    # Code grouping effects as an interaction with no main effects
    m_terms <- model
    for (e in model){
      if (e %in% split){
        if (e == "Age"){
          m_terms[which(m_terms=="Age")] <- "s(Age, by=Age_Split,...)"
        } else {
          m_terms[which(m_terms==e)] <- paste0(e, "*", e, "_Split")
        }
      }
    }
    ind_str <- Reduce(function(...){paste(..., sep="+")}, c("0", m_terms))
  } else {
    ind_effects <- c("0", model)
    ind_effects <- str_replace(ind_effects, "Age", "s(Age, ...)")
    ind_str <- Reduce(function(...){paste(..., sep="+")}, ind_effects)
  }
  
  # Combine the two sides of the formula
  formula_str <- paste(dep_str, ind_str, sep="~")
  
  return (formula_str)
}

# Extracting effects for GAM models
extract_effects_gam <- function(growth_model, model, split, link, tra)
{
  # Reduced model excludes age, useful for iteration
  reduced_model <- model[model!="Age"]
  
  # Skeleton effects for relisting coefficients
  skele <- make_skeleton_effects(tra, model, split, link)
  
  # Grab the coefficients from the regression model
  effect_coef <- coef(growth_model)
  
  # Smoothed coefficients are not what we're looking for, discard them
  smooth_coef <- sapply(names(effect_coef), str_detect, pattern="s(.*)")
  effect_coef <- effect_coef[!smooth_coef]
  
  # Fix the names for split effects
  for (e in reduced_model[reduced_model %in% split]){
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
  for (e in reduced_model[!(reduced_model %in% split)]){
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
  for (e in reduced_model[reduced_model %in% split]){
    groups <- names(effects[[e]])
    
    # Identify base group
    # Grab all names within the relevant effect
    all_e <- str_match_all(names(coef(growth_model)), paste0(e,".*"))
    
    # Find included groups
    raw_groups <- lapply(all_e, str_match, pattern=paste0(split_name, ".*"))
    raw_group_names <- unlist(lapply(raw_groups, str_replace, pattern=split_name, replacement=""))
    inc_groups <- unique(raw_group_names[!is.na(raw_group_names)])    
    
    # Find all groups for the relevant effect
    all_groups <- names(skele[[e]])
    
    # Determine base group by process of elimination
    base_group <- all_groups[!(all_groups %in% inc_groups)]
    
    # Add base effect to interactions
    effect_b <- effects[[e]][[base_group]]
    
    
    for (group in groups){
      if (group != base_group){
        base <-  effect_b[names(effects[[e]][[group]])]
        effects[[e]][[group]] <- base + effects[[e]][[group]]
      }
    }
  }
  
  # Add in effects for age
  if ("Age" %in% model)
  {
    # Find estimates of effect by comparing predictions to predictions from other variables
    base_indices <- sapply(reduced_model, function(x){unique(tra[[x]][1])})
    
    
    if (!("Age" %in% split))
    {
      Age_levels <- as.numeric(unique(tra$Age))
      
      dummy_data <- matrix(NA, nrow=length(Age_levels), ncol=length(model))
      colnames(dummy_data) <- c(reduced_model, "Age")
      dummy_data[,1:length(reduced_model)] <- base_indices
      dummy_data <- as.data.frame(dummy_data)
      dummy_data$Age <- Age_levels
      
      predicted_by_age <- predict(growth_model, dummy_data)
      names(predicted_by_age) <- dummy_data$Age
      
      # Predictions - other effects = age effect
      baseline <- 0
      for (e in reduced_model){
        baseline <- baseline + effects[[e]][base_indices[e]]
      }     
      
      effects$Age <- predicted_by_age - baseline
    } else {
      groups <- names(effects$Age)
      
      # Generate predictions for each group / age combo
      for (group in groups)
      {
        Age_levels <- as.numeric(unique(tra[tra$Age_Split==group, "Age"]))
        
        dummy_data <- matrix(NA, nrow=length(Age_levels), ncol=length(model)+1)
        colnames(dummy_data) <- c(reduced_model, "Age", "Age_Split")
        dummy_data[,1:length(reduced_model)] <- base_indices
        dummy_data <- as.data.frame(dummy_data)
        dummy_data$Age <- Age_levels
        dummy_data$Age_Split <- group
        
        predicted_by_age <- predict(growth_model, dummy_data)
        names(predicted_by_age) <- dummy_data$Age
        
        # Predictions - other effects = age effect
        baseline <- 0
        for (e in reduced_model){
          baseline <- baseline + effects[[e]][base_indices[e]]
        }     
        
        effects$Age[[group]] <- predicted_by_age - baseline
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

# Number of parameters estimated (k)
# Custom for GAM models
# Borrows from k_tra()
extract_k_gam <- function (growth_model, tra, model, split)
{
  
  # One parameter automatically for estimate of error term
  k <- 1
  
  # One parameter is estimated for each index of the effect vectors
  for (E in model[model!="Age"])
  {
    if (E %in% split)
    {
      cname <- paste(E, "Split", sep="_")
      groups <- levels(tra[[cname]])
      
      # Each group adds unique parameters
      # But only count parameters that are actually in that group
      for (group in groups){
        k <- k + length(unique(tra[tra[[cname]]==group, E]))
      }
    } else 
    {
      k <- k + length(unique(tra[[E]]))
    }
  }
    
  if ("Age" %in% model)
  {
    k <- k + sum(summary(growth_model)$edf)
  }
  
  # Information about some parameters is lost due to rescaling (dummy variable trap)
  # Null models have only one parameter
  num_effects <- length(model)
  
  k <- ifelse (num_effects > 0, k - (num_effects-1), 1)
  
  
  return(k)
}

