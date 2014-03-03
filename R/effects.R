# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, id, link, dep_var="Growth", group=NA)
{  

  # Estimate effect using averages (crude method of moments)
  estimate_effect <- function(id_i)
  {
    # Estimate using rows from relevant group
    if (!is.na(group)){
      cname <- paste(id, "Group", sep="_")
      data <- tra[tra[[id]]==id_i & tra[[cname]]==group, dep_var]
    } else {
      data <- tra[tra[[id]]==id_i, dep_var]
    }
    
    # Geometric mean is log equivalent of the arithmetic mean
    if (link=="log"){
      est_effect <- geomMean(data)
    } else {
      est_effect <- mean(data, na.rm=TRUE)
    }
    
    names (est_effect) <- id_i
    
    return(est_effect)
  }
  
  
  if (!is.na(group)){
    cname <- paste(id, "Group", sep="_")
    id_levels <- unique(tra[tra[[cname]]==group,][[id]])
  } else {
    id_levels <- unique(tra[[id]])
  }
  est_effect <- sapply(id_levels, estimate_effect)
  
  return (est_effect)
}

# Remove an effect from a tree ring array
remove_effect <- function (tra, effect, id, link="log", dep_var="Growth", group=NA)
{
    
  removed_tra <- tra
  
  for (effect_id in names(effect))
  {
    # Only affect rows from relevant group
    if (!is.na(group)){
      cname <- paste(id, "Group", sep="_")
      relevant_rows <- tra[[id]]==effect_id & tra[[cname]]==group
    } else {
      relevant_rows <- tra[[id]]==effect_id
    }
    
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
sort_effects <- function(effects, tra, group_by=NA)
{
  sorted_effects <- list()
  
  for (i in names(effects)){
    
    if (i %in% group_by)
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
rescale_effects <- function (effects, link="log", group_by=NA)
{
  if (i %in% group_by)
  {
    # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
    # So log transform, set mean to 0, then unlog
    if (link=="log")
    {
      effects <- lapply(names(effects), function(x){
        if (x %in% group_by){
          return(lapply(effects[[x]], log))
        } else {
          return (log(effects[[x]]))
        }
      }
      )
    }
    
    mean_effects <- lapply(effects, function(x){lapply(x,mean, na.rm=T)})  
    
    mean_effects <- lapply(names(effects), function(x){
      if (x %in% group_by){
        return(mean(unlist(effects[[x]]), na.rm=T))
      } else {
        return (mean(effects[[x]], na.rm=T))
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
        if ("Tree" %in% group_by){
          effects$Tree <- lapply(effects$Tree, function(x){x+mean_effects$Time})
        } else {
          effects$Tree <- effects$Tree+mean_effects$Time
        }
        
        if ("Time" %in% group_by){
          effects$Time <- lapply(effects$Time, function(x){x-mean_effects$Time})
        } else {
          effects$Time <- effects$Time+mean_effects$Tree
        }
      } else{
        return (effects)
      }
    }
    
    # I
    if ("Tree" %in% names(effects)){
      if ("Tree" %in% group_by){
        effects$Tree <- lapply(effects$Tree, function(x){x-mean_effects$Tree})
      } else {
        effects$Tree <- effects$Tree-mean_effects$Tree
      }
      
      if ("Age" %in% group_by){
        effects$Age <- lapply(effects$Age, function(x){x+mean_effects$Tree})
      } else {
        effects$Age <- effects$Age+mean_effects$Tree
      }
    }
    
    # T
    if ("Time" %in% names(effects)){
      if ("Time" %in% group_by){
        effects$Time <- lapply(effects$Time, function(x){x-mean_effects$Time})
      } else {
        effects$Time <- effects$Time-mean_effects$Time
      }
      
      if ("Age" %in% group_by){
        effects$Age <- lapply(effects$Age, function(x){x+mean_effects$Time})
      } else {
        effects$Age <- effects$Age+mean_effects$Time
      }
    }
    
    if (link=="log")
    {
      effects <- lapply(names(effects), function(x){
        if (x %in% group_by){
          return(lapply(effects[[x]], exp))
        } else {
          return (exp(effects[[x]]))
        }
      }
      )
    }
  } else {
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
  }
  
  
  return (effects)
}

# Create skeleton effects ####

make_skeleton_effects <- function(tra, model, group_by, link)
{
  # Create initial list
  effects <- vector(mode="list", length=length(model))
  names(effects) <- model
  
  # Dummy starting effects
  for (i in model){
    if (i %in% group_by)
    {
      # Each group is a sub-list
      cname <- paste(i, "Group", sep="_")
      effects[[i]] <- vector(mode="list", length=nlevels(tra[[cname]]))
      names(effects[[i]]) <- levels(tra[[cname]])
      
      for (j in levels(tra[[cname]]))
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
      dim_i <- nlevels(tra[[i]])
      
      if (link=="log")
      {
        effects[[i]] <-  rep.int(1,  dim_i)
      } else
      {
        effects[[i]] <-  rep.int(0,  dim_i)
      }
      
      names(effects[[i]]) <- levels(tra[[i]])
    }
  }
  
  return(effects)
  
}

# Crude estimation of standard errors for each effect ####
est_se <- function(resids, model, group_by=NA, link="log", dep_var="Growth"){
  
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
      cname <- paste(E, "Group", sep="_")
      return(resids[resids[[E]]==e & resids[[cname]]==group, dep_var])
    } else {
      return(resids[resids[[E]]==e, dep_var])
    }
  }
  
  effect_se <- function(E)
  {
    if (E %in% group_by){
      skele <- make_skeleton_effects(resids,model,group_by,link)
      groups <- names(skele[[E]])
      cname <- paste(E, "Group", sep="_")
      
      e_list <- lapply(groups, function(group){unique(tra[tra[[cname]]==group,][[E]])})
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
      e_list <- levels(resids[[E]])
      res_list <- lapply(e_list, grab_res, E=E)
      se_list <- sapply(res_list, mom_se)
      names(se_list) <- e_list
    }
    return(se_list)
  }
  
  all_se <- lapply(model, effect_se)
  names(all_se) <- model
  
  # Sort so names line up
  all_se <- sort_effects(all_se, resids, group_by)
  
  return(all_se)
}