# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, id, link, dep_var="Growth")
{  

  # Estimate effect using averages (crude method of moments)
  estimate_effect <- function(id_i)
  {
    data <- tra[tra[[id]]==id_i, dep_var]
    if (link=="log"){
      est_effect <- geomMean(data)
    } else {
      est_effect <- mean(data, na.rm=TRUE)
    }
    return(est_effect)
  }
  
  est_effect <- sapply(levels(tra[[id]]), estimate_effect)
  
  return (est_effect)
}

# Remove an effect from a tree ring array
remove_effect <- function (tra, effect, id, link="log", dep_var="Growth")
{
    
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    relevant_rows <- tra[[id]]==effect_id
    if (link=="identity")
    {
      removed_tra[relevant_rows,"Growth"] <- removed_tra[relevant_rows, dep_var] - effect[effect_id]
    }
    else
    {
      removed_tra[relevant_rows,"Growth"] <- removed_tra[relevant_rows, dep_var] / effect[effect_id]
    }
  }
  
  return (removed_tra)
  
}

# Cleaning up effects ####

# Correctly sort elements of the effect vectors
sort_effects <- function(effects, tra)
{
  sorted_effects <- list()
  
  for (i in names(effects)){
    effect_names <- names(effects[[i]])
    
    # Sort time and age by ascending time
    # Other variables are just sorted alphabetically
    if (i == "Time" | i == "Age"){
      ordered_names <- as.character(sort(as.numeric(as.character(effect_names))))
    } else {
      ordered_names <- sort(effect_names)
    }
    
    sorted_effects[[i]] <- effects[[i]][ordered_names]
    
  }
      
  return(sorted_effects)
}

# Rescale effect vectors to canonical form
rescale_effects <- function (effects, link="log")
{
  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (link=="log")
  {
    effects <- lapply(effects, log)
  }
  
  mean_effects <- lapply(effects, mean, na.rm=T)  
  
  # Scale Tree effect and Time effect to mean of 0
  # Scale Age effect so sum of effects stays the same
  
  # If Age effect is missing, try to add magnitude information to I
  # Otherwise return raw effects
  if(!("Age" %in% names(effects)))
  {
    if (all(c("Tree", "Time") %in% names(effects))){
      effects$Tree <- effects$Tree-mean_effects$Tree
      effects$Time <- effects$Time+mean_effects$Tree
    } else{
      return (effects)
    }
  }
  
  # I
  if ("Tree" %in% names(effects)){
    effects$Tree <- effects$Tree-mean_effects$Tree
    effects$Age <- effects$Age+mean_effects$Tree
  }
  
  # T
  if ("Time" %in% names(effects)){
    effects$Time <- effects$Time-mean_effects$Time
    effects$Age <- effects$Age+mean_effects$Time
  }
  
  if (link=="log")
  {
    effects <- lapply(effects, exp)
  }
  
  return (effects)
}

# Crude estimation of standard errors for each effect ####
est_se <- function(residuals, model, link="log", dep_var){
  
  # Method of moments confidence intervals
  # Compute classic standard errors for each effect
  mom_se <- function(x){
    se <- sd(x, na.rm=TRUE)/length(x[!is.na(z)])
    return(se)
  }
  
  grab_res <- function(e, E){
    residuals[residuals[[E]]==e, dep_var]
  }
  
  effect_se <- function(E)
  {
    e_list <- levels(residuals[[E]])
    res_list <- lapply(e_list, grab_res)
    se_list <- lapply(res_list, mom_se)
    return(se)
  }
  
  all_se <- lapply(model, function(E){sapply(E, effect_se)})
  
  # Sort so names line up
  all_se <- sort_effects(all_se, residuals)
  
  return(all_se)
}