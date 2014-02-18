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

# Get sample size along a dimension ####
sample_depth_tra <- function(tra, factor.dim=2, sparse=FALSE){ #1 is tree, 2 is time, 3 is age
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse)
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  if(sparse)
  {
    positions <- tra[[factor.dim+1]]
    ids <- levels (positions)
    
    sample_depth <- sapply(ids, function(x){sum(positions==x)})
    
    if(factor.dim==1)
    {
      ordering <- sort(names(sample_depth))
    }  else
    {
      ordering <- sort(as.numeric(names(sample_depth)))
    }
   
    sample_depth <- sample_depth[sapply(ordering, as.character)]
  } else 
  {
    filled_cells <- !is.na(tra)
    
    sample_depth <- apply(filled_cells, factor.dim, sum)
  }
  
  return (sample_depth)
}
