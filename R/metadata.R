# Grabbing ages of trees ####
grab_age <- function(i, tra)
{
  age <- max(as.numeric(as.character(tra[tra$Tree==i,"Age"])))
  return(as.numeric(age))
}


# Grabbing birth years of trees ####
grab_birth_year <- function(i, tra)
{
  slice <- tra[tra$Tree==i,][1,]
  birth_year <- as.numeric(as.character(slice$Time)) -  as.numeric(as.character(slice$Age)) + 1 
  return(birth_year)
}

get_birth_years <- function(tra)
{
  birth_years <- sapply(unique(tra$Tree), grab_birth_year, tra=tra)
  return (birth_years)
}

get_birth_index <- function(birth_year, tra)
{
  years <- sort(as.numeric(as.character(unique(tra$Time))))
  
  year_index <- which(birth_year==years)
  
  # Account for out of range years
  if (length(year_index)==0){
    base_index <- which(years[1]==years)
    year_index <- base_index - (years[1] - birth_year)
  }
  age_index <- 1
  birth_index <- year_index - age_index
  
  return (birth_index)
}

# Get sample size for a characteristic ####
sample_depth_tra <- function(tra, id="Time", split=NA){
  
  ids <- unique(tra[[id]])
  
  if(id %in% split){
    cname <- paste(id, "Split", sep="_")
    groups <- unique(tra[[cname]])
    
    sample_depth <- vector(mode="list", length=length(groups))
    names(sample_depth) <- groups
    
    for (group in groups){
      sample_depth[[group]] <- sapply(ids, function(x){sum(tra[[id]]==x & tra[[cname]]==group)})
      names(sample_depth[[group]]) <- ids
      
      if(id=="Time" | id=="Age")
      {
        ordering <- sort(as.numeric(names(sample_depth[[group]])))
      }  else
      {
        ordering <- sort(names(sample_depth[[group]]))
      }
      
      sample_depth[[group]] <- sample_depth[[group]][sapply(ordering, as.character)]
    }
  } else {
    sample_depth <- sapply(ids, function(x){sum(tra[[id]]==x)})
    names(sample_depth) <- ids
    
    if(id=="Time" | id=="Age")
    {
      ordering <- sort(as.numeric(names(sample_depth)))
    }  else
    {
      ordering <- sort(names(sample_depth))
    }
    
    sample_depth <- sample_depth[sapply(ordering, as.character)]
  }
  
  return (sample_depth)
}

# Find mean series length by time ####
# Age, tree not at all useful

series_length_tra <- function(tra, split){
  ids <- unique(tra$Tree)
  series_lengths <- sapply(ids, function(id){sum(tra$Tree==id)})
  names(series_lengths) <- ids
  years <- unique(tra$Time)
  
  if ("Time" %in% split){
    groups <- unique(tra$Time_Split)
    
    mean_series_length <- vector(mode="list", length=length(groups))
    names(mean_series_length) <- groups
    
    for (group in groups){
      inc_trees <- lapply(years, function(year){
        as.character(unique(tra[tra$Time==year & tra$Time_Split==group, "Tree"]))
      })
      
      mean_series_length[[group]] <- sapply(inc_trees, function(x){
        mean(series_lengths[x])
      })
      
      names(mean_series_length) <- years
    }
  } else {    
    inc_trees <- lapply(years, function(year){
      as.character(unique(tra[tra$Time==year, "Tree"]))
    })
    
    mean_series_length <- sapply(inc_trees, function(x){
      mean(series_lengths[x])
    })
    
    names(mean_series_length) <- years
    
  }
  
  # Sorting to correct order
  msl_list <- list(mean_series_length)
  names(msl_list) <- "Time"
  mean_series_length <- sort_effects(msl_list, split)$Time
  
  return(mean_series_length)
  
}
