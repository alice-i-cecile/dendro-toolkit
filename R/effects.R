# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, id, mean_type="arithmetic")
{  

  # Estimate effect using averages (crude method of moments)
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
remove_effect <- function (tra, effect, id, link="log")
{
    
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    relevant_rows <- tra[[id]]==effect_id
    if (link=="identity")
    {
      removed_tra[relevant_rows,"Growth"] <- removed_tra[relevant_rows,"Growth"] - effect[effect_id]
    }
    else
    {
      removed_tra[relevant_rows,"Growth"] <- removed_tra[relevant_rows,"Growth"] / effect[effect_id]
    }
  }
  
  return (removed_tra)
  
}

# Add dummy effect vectors if some are missing
pad_effects <- function(effects, tra, link="log")
{
  # Set the value to fill dummy coefficients with
  if (link=="log"){
    na.value <- 1
  } else {
    na.value <- 0
  }
  
  # Initialize dummy effects lists
  new_effects <- list(I=NA, T=NA, A=NA)
  
  # Fill empty values  
  for(i in c("Tree", "Time", "Age")){
        
    if (length(effects[[i]] > 0)){
      new_effects[[i]] <-  effects[[i]]
      
      # Fill in dummy levels
      if(length(effects[[i]]) < nlevels(tra[[i]]))
      {
        missing_names <- levels(tra[[i]])[!(levels(tra[[i]]) %in% names(effects[[i]]))]
        missing_effects <- rep(na.value, times=length(missing_names))
        names (missing_effects) <- missing_names
        new_effects[[i]] <- c(effects[[i]], missing_effects)
      }
      
    } else {
      tra_dim <- which(c("Tree", "Time", "Age")==i)
      
      new_effects[[i]] <- rep.int(na.value, nlevels(tra[[tra_dim + 1]]))
      effect_names <- levels(tra[[tra_dim + 1]])
      if(i=="Time" | i=="Age")
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
  effect_names <- names(effects$Tree)
  sorted_effects$Tree <- effects$Time[sort(effect_names)]
  
  # Sort T and A
  for (j in c("Time", "Age"))
  {
    effect_names <- as.numeric(names(effects[[j]]))
    sorted_effects[[j]] <- effects[[j]][as.character(sort(effect_names))]
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
  
  # Scale I and T to mean of 0
  # Scale A so sum of effects stays the same
  
  # If A is missing, leave effects uncscaled
  if(!sum(!is.null(effects$Age)))
  {
    return (effects)
  }
  
  # I
  if (sum (!is.null(effects$Tree))){
    effects$Tree <- effects$Tree-mean_effects$Tree
    effects$Age <- effects$Age+mean_effects$Tree
  }
  
  # T
  if (sum (!is.null(effects$Time))){
    effects$Time <- effects$Time-mean_effects$Time
    effects$Age <- effects$Age+mean_effects$Time
  }
  
  if (link=="log")
  {
    effects <- lapply(effects, exp)
  }
  
  return (effects)
} 
