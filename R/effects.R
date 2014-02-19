# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, factor.dim, mean_type="arithmetic")
{  
  id <- factor.dim + 1
  
  estimate_effect <- function(id_i)
  {
    data <- tra[tra[[id]]==id_i, "G"]
    if (mean_type=="geometric"){
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
remove_effect <- function (tra, effect, factor.dim, form="multiplicative")
{
  
  id <- factor.dim + 1
  
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    relevant_rows <- tra[[id]]==effect_id
    if (form=="additive")
    {
      removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] - effect[effect_id]
    }
    else
    {
      removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] / effect[effect_id]
    }
  }
  
  return (removed_tra)
  
}

# Add dummy effect vectors if some are missing
pad_effects <- function(effects, tra, form="multiplicative")
{
  # Set the value to fill dummy coefficients with
  if (form=="multiplicative"){
    na.value <- 1
  } else {
    na.value <- 0
  }
  
  # Initialize dummy effects lists
  new_effects <- list(I=NA, T=NA, A=NA)
  
  # Fill empty values  
  for(i in c("I", "T", "A")){
    
    j <- c(I="i", T="t", A="a")[i]
    
    if (length(effects[[i]] > 0)){
      new_effects[[i]] <-  effects[[i]]
      
      # Fill in dummy levels
      if(length(effects[[i]]) < nlevels(tra[[j]]))
      {
        missing_names <- levels(tra[[j]])[!(levels(tra[[j]]) %in% names(effects[[i]]))]
        missing_effects <- rep(na.value, times=length(missing_names))
        names (missing_effects) <- missing_names
        new_effects[[i]] <- c(effects[[i]], missing_effects)
      }
      
    } else {
      tra_dim <- which(c("I", "T", "A")==i)
      
      new_effects[[i]] <- rep.int(na.value, nlevels(tra[[tra_dim + 1]]))
      effect_names <- levels(tra[[tra_dim + 1]])
      if(i=="T" | i=="A")
      {
        effect_names <- as.numeric(as.character(effect_names))
      }
      names(new_effects[[i]]) <- sort(effect_names)
    }
  }
  
  return(new_effects)
}

# Correctly order effect vectors
sort_effects <- function(effects, tra)
{
  # Only ascending sorting as no "standard" order is retained
  sorted_effects <- list()
  
  # Sort I
  effect_names <- names(effects$I)
  sorted_effects$I <- effects$I[sort(effect_names)]
  
  # Sort T and A
  for (j in c("T", "A"))
  {
    effect_names <- as.numeric(names(effects[[j]]))
    sorted_effects[[j]] <- effects[[j]][as.character(sort(effect_names))]
  }
  
  return(sorted_effects)
}

# Rescale effect vectors to canonical form
rescale_effects <- function (effects, form="multiplicative")
{
  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (form=="multiplicative")
  {
    effects <- lapply(effects, log)
  }
  
  mean_effects <- lapply(effects, mean, na.rm=T)  
  
  # Scale I and T to mean of 0
  # Scale A so sum of effects stays the same
  
  # If A is missing, leave effects uncscaled
  if(!sum(!is.null(effects[[3]])))
  {
    return (effects)
  }
  
  # I
  if (sum (!is.null(effects[[1]]))){
    effects[[1]] <- effects[[1]]-mean_effects[[1]]
    effects[[3]] <- effects[[3]]+mean_effects[[1]]
  }
  
  # T
  if (sum (!is.null(effects[[2]]))){
    effects[[2]] <- effects[[2]]-mean_effects[[2]]
    effects[[3]] <- effects[[3]]+mean_effects[[2]]
  }
  
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
  }
  
  return (effects)
} 
