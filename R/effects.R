# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, id, link, dep_var="Growth", group=NA)
{  

  # Estimate effect using averages (crude method of moments)
  estimate_effect <- function(id_i)
  {
    # Estimate using rows from relevant group
    if (!is.na(group)){
      cname <- paste(id, "Split", sep="_")
      data <- tra[tra[[id]]==id_i & tra[[cname]]==group, dep_var]
    } else {
      data <- tra[tra[[id]]==id_i, dep_var]
    }
    
    # Geometric mean is log equivalent of the arithmetic mean
    if (link=="log"){
      estim_effect <- geomMean(data)
    } else {
      estim_effect <- mean(data, na.rm=TRUE)
    }
    
    names (estim_effect) <- id_i
    
    return(estim_effect)
  }
  
  
  if (!is.na(group)){
    cname <- paste(id, "Split", sep="_")
    id_levels <- unique(tra[tra[[cname]]==group,][[id]])
  } else {
    id_levels <- unique(tra[[id]])
  }
  estim_effect <- sapply(id_levels, estimate_effect)
  
  # Ensure correct ordering
  list_effect <- list(estim_effect)
  names(list_effect) <- id
  estim_effect <- sort_effects(list_effect, tra, split=NA)[[id]]
  
  return (estim_effect)
}

# Remove an effect from a tree ring array
remove_effect <- function (tra, effect, id, link="log", dep_var="Growth", group=NA)
{
    
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    # Only affect rows from relevant group
    if (!is.na(group)){
      cname <- paste(id, "Split", sep="_")
      relevant_rows <- tra[[id]]==effect_id & tra[[cname]]==group
    } else {
      relevant_rows <- tra[[id]]==effect_id
    }
    
    if (link=="identity")
    {
      removed_tra[relevant_rows, dep_var] <- removed_tra[relevant_rows, dep_var] - effect[effect_id]
    }
    else{
      removed_tra[relevant_rows, dep_var] <- removed_tra[relevant_rows, dep_var] / effect[effect_id]
    }
  }
  
  return (removed_tra)
  
}

# Cleaning up effects ####

# Correctly sort elements of the effect vectors
sort_effects <- function(effects, tra, split=NA)
{
  sorted_effects <- list()
  
  for (i in names(effects)){
    
    if (i %in% split)
    {
      for (j in names(effects[[i]])){
        effect_names <- names(effects[[i]][[j]])
        # Sort time and age by ascending time
        # Other variables are just sorted alphabetically
        if (i == "Time" | i == "Age"){
          ordered_names <- as.character(sort(as.numeric(as.character(effect_names))))
        } else {
          ordered_names <- sort(effect_names)
        }
        
        sorted_effects[[i]][[j]] <- effects[[i]][[j]][ordered_names]
        
      }
    } else 
    {
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
  }
      
  return(sorted_effects)
}

# Rescale effect vectors to canonical form
rescale_effects <- function (effects, link="log", split=NA)
{
  # Set aside effects so we can grab names later
  rescaled_effects <- effects

  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (link=="log")
  {
      rescaled_effects <- lapply(names(rescaled_effects), function(x){
        if (x %in% split){ 
          return(lapply(rescaled_effects[[x]], log))
        } else {
          return (log(rescaled_effects[[x]]))
        }
      })
      names(rescaled_effects) <- names(effects)
  }
      
  mean_effects <- lapply(names(rescaled_effects), function(x){
    if (x %in% split){
      return(mean(unlist(rescaled_effects[[x]]), na.rm=T))
    } else {
      return (mean(rescaled_effects[[x]], na.rm=T))
    }
  }
  )
  names(mean_effects) <- names(effects)
  
  # Scale Tree effect and Time effect to mean of 0
  # Scale Age effect so sum of effects stays the same
  
  # If Age effect is missing, try to add magnitude information to I
  # Otherwise return raw effects
  if(!("Age" %in% names(effects)))
  {
    if (all(c("Tree", "Time") %in% names(effects))){
      if ("Tree" %in% split){
        rescaled_effects$Tree <- lapply(rescaled_effects$Tree, function(x){x+mean_effects$Time})
      } else {
        rescaled_effects$Tree <- rescaled_effects$Tree+mean_effects$Time
      }
      
      if ("Time" %in% split){
        rescaled_effects$Time <- lapply(rescaled_effects$Time, function(x){x-mean_effects$Time})
      } else {
        rescaled_effects$Time <- rescaled_effects$Time-mean_effects$Time
      }
    }
  } else {
    # I
    if ("Tree" %in% names(effects)){
      if ("Tree" %in% split){
        rescaled_effects$Tree <- lapply(rescaled_effects$Tree, function(x){x-mean_effects$Tree})
      } else {
        rescaled_effects$Tree <- rescaled_effects$Tree-mean_effects$Tree
      }
      
      if ("Age" %in% split){
        rescaled_effects$Age <- lapply(rescaled_effects$Age, function(x){x+mean_effects$Tree})
      } else {
        rescaled_effects$Age <- rescaled_effects$Age+mean_effects$Tree
      }
    }
    
    # T
    if ("Time" %in% names(rescaled_effects)){
      if ("Time" %in% split){
        rescaled_effects$Time <- lapply(rescaled_effects$Time, function(x){x-mean_effects$Time})
      } else {
        rescaled_effects$Time <- rescaled_effects$Time-mean_effects$Time
      }
      
      if ("Age" %in% split){
        rescaled_effects$Age <- lapply(rescaled_effects$Age, function(x){x+mean_effects$Time})
      } else {
        rescaled_effects$Age <- rescaled_effects$Age+mean_effects$Time
      }
    }
  }
  
  if (link=="log")
  {
    rescaled_effects <- lapply(names(rescaled_effects), function(x){
      if (x %in% split){
        return(lapply(rescaled_effects[[x]], exp))
      } else {
        return (exp(rescaled_effects[[x]]))
      }
    })
  }
  
  # Reset first level names  
  names(rescaled_effects) <- names(effects)
  
  return (rescaled_effects)

}

# Create skeleton effects ####

make_skeleton_effects <- function(tra, model, split, link)
{
  # Create initial list
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Dummy starting effects
  for (i in model){
    if (i %in% split)
    {
      # Each group is a sub-list
      cname <- paste(i, "Split", sep="_")
      effects[[i]] <- vector(mode="list", length=length(unique(tra[[cname]])))
      names(effects[[i]]) <- unique(tra[[cname]])
      
      for (j in unique(tra[[cname]]))
      {
        dim_j <- length(unique(tra[tra[[cname]]==j,][[i]]))
        
        if (link=="log")
        {
          effects[[i]][[j]] <-  rep.int(1,  dim_j)
        } else
        {
          effects[[i]][[j]] <-  rep.int(0,  dim_j)
        }
        
        names(effects[[i]][[j]]) <- unique(tra[tra[[cname]]==j,][[i]])
        
      }
    } else {
      dim_i <- length(unique(tra[[i]]))
      
      if (link=="log")
      {
        effects[[i]] <-  rep.int(1,  dim_i)
      } else
      {
        effects[[i]] <-  rep.int(0,  dim_i)
      }
      
      names(effects[[i]]) <- unique(tra[[i]])
    }
  }
  
  effects <- sort_effects(effects, tra, split)
  
  return(effects)
  
}

# Crude estimation of standard errors for each effect ####
est_se <- function(resids, model, split=NA, link="log", dep_var="Growth"){
  
  # E: name of effect
  # e: index
  
  # Method of moments confidence intervals
  # Compute classic standard errors for each effect
  mom_se <- function(x){
    se <- sd(x, na.rm=TRUE)/length(x[!is.na(x)])
    return(se)
  }
  
  grab_res <- function(e, E, group=NA){
    if (!is.na(group)){
      cname <- paste(E, "Split", sep="_")
      return(resids[resids[[E]]==e & resids[[cname]]==group, dep_var])
    } else {
      return(resids[resids[[E]]==e, dep_var])
    }
  }
  
  effect_se <- function(E)
  {
    if (E %in% split){
      skele <- make_skeleton_effects(resids,model,split,link)
      groups <- names(skele[[E]])
      cname <- paste(E, "Split", sep="_")
      
      e_list <- lapply(groups, function(group){unique(resids[resids[[cname]]==group,][[E]])})
      names(e_list) <- groups
      
      res_list <- lapply(groups, function(group)
      {
        sel_resid <- sapply(e_list[[group]], grab_res, E=E, group=group)
        names(sel_resid) <- e_list[[group]]
        return(sel_resid)
      })
      names(res_list) <- groups
      
      se_list <- lapply(groups, function(group)
      {
        sapply(res_list[[group]], mom_se) 
      })
      names(se_list) <- groups
      
    } else {
      e_list <- unique(resids[[E]])
      res_list <- lapply(e_list, grab_res, E=E)
      se_list <- sapply(res_list, mom_se)
      names(se_list) <- e_list
    }
    return(se_list)
  }
  
  all_se <- lapply(model, effect_se)
  names(all_se) <- model
  
  # Sort so names line up
  all_se <- sort_effects(all_se, resids, split)
  
  return(all_se)
}