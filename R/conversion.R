# utility ####

# Get endpoints of a rwl series
getEndpoints <- function (series, side="start"){
  #side: "start" or "stop" of series
  rings <- subset(series, !is.na(series))
  if (side=="start"){
    return (head(rownames(rings),1))
  }
  if (side=="stop"){
    return (tail(rownames(rings),1))
  }      
}

# RWL to TRA ####

# Go directly from a rwl to tra
# Avoid unnecessary memory bottlenecks
rwl_to_tra <- function (rwl, birth_years=NULL, dep_var="Growth")
{
  if (is.null(birth_years)){
    # Determine birth for each tree
    birth_years <- foreach(i=colnames(rwl)) %do% {getEndpoints(rwl[i], side="start")}
    names (birth_years) <- names(rwl)
    
  }  else {
    # If rwl file is too short, add empty rows to it
    birth_limit <- min(sapply(birth_years, as.numeric))
    rwl_limit <- min(sapply(rownames(rwl), as.numeric))
    
    if (birth_limit < rwl_limit){
      num_missing <- rwl_limit-birth_limit
      
      empty_rows <- rbind(rwl[0,], matrix(NA, num_missing, ncol(rwl)))
      
      colnames(empty_rows) <- colnames(rwl)
      rownames(empty_rows) <- birth_limit:(rwl_limit-1)
      
      rwl <- rbind(empty_rows, rwl)
    }
  }
  
  # Find the indices for each year
  birth_index <- vector ()
  
  for (i in 1:length(birth_years)){
    birth_index[i] <- which(rownames(rwl)==toString(birth_years[i]))
  }
  
  names (birth_index) <- names (birth_years)
  
  # Find the full cells
  filled <- which(!is.na(rwl), arr.ind=TRUE)
  n_data_points <- nrow(filled)
  
  
  add_growth_data <- function (year, tree)
  {
    birth <- birth_index[tree]
    # Extract growth data
    G <- rwl [year, tree]
    
    # Find tree of growth data
    i <- names(rwl)[tree]
    
    # Find year of growth data
    t <- rownames(rwl)[year]
    
    # Find age of ring data
    a <- year - birth + 1
    
    out <- list(G, i, t, a)
    return(out)
  }
  
  # Find the growth data at each filled cell and record its position
  raw_tra <- mapply(add_growth_data, year=filled[,1], tree=filled[,2], SIMPLIFY=FALSE)
  
  # Final formatting
  tra <- data.frame(matrix(unlist(raw_tra),ncol=4, byrow=TRUE))
  names(tra) <- c(dep_var, "Tree", "Time", "Age")
  tra[[dep_var]] <- as.numeric(tra[[dep_var]])
  
  return(tra)
  
}

# TRA to RWL ####

# Often loses information about age, sometimes metadata
tra_to_rwl <- function(tra, dep_var="Growth")
{
  # Create the base table
  tree_names <- unique(tra$Tree)

  years <- unique(tra$Year)
  years <- sort(as.numeric(as.character(year)))
  
  num_trees <- length(tree_names)
  num_years <- length(years)
  
  rwl <- as.data.frame(matrix(data=NA, nrow=num_year, ncol=num_trees, dimnames=list(years, trees)))
  
  # Fill in data
  for (i in 1:nrows(tra))
  {
    
    # Extract year and tree coordinates
    tree <- tra[i, "Tree"]
    year <- tra[i, "Year"]
    
    tree_index <- which(names(rwl)==tree)
    year_index <- which(rownames(rwl)==year)
      
    # Convert growth data
    rwl[year_index, tree_index] <- tra[i, dep_var]
  }
  
  return(rwl)
  
}
