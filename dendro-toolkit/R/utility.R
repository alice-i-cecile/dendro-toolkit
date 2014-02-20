# Geometric mean utility function ####
geomMean <- function(x){
  if (length(x)==0){
    return(NA)
  }
  val <- x[!is.na(x)]
  l_val <- log(val)
  out <- exp(mean(l_val))
  return (out)
}

# Get sample size for a characteristic ####
sample_depth_tra <- function(tra, id="Time"){
    
  positions <- tra[[id]]
  ids <- levels (positions)
  
  sample_depth <- sapply(ids, function(x){sum(positions==x)})
  
  if(id=="Time" | id=="Age")
  {
    ordering <- sort(as.numeric(names(sample_depth)))
  }  else
  {
    ordering <- sort(names(sample_depth))
  }
 
  sample_depth <- sample_depth[sapply(ordering, as.character)]
  
  return (sample_depth)
}
