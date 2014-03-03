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
  birth_years <- sapply(unique(tra$Tree), grab_birth_year, tra=tra, sparse)
  return (birth_years)
}

get_birth_index <- function(birth_year, tra)
{
  years <- sort(as.numeric(unique(tra$Time)))
  
  year_index <- which(birth_year==years)
  age_index <- 1
  birth_index <- year_index - age_index
  
  return (birth_index)
}

# Get sample size for a characteristic ####
sample_depth_tra <- function(tra, id="Time", group_by=NA){
  
  ids <- unique(tra[[id]])
  
  if(id %in% group_by){
    cname <- paste(id, "Group", sep="_")
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
