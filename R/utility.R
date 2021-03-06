# Truncate tree ring array to ensure minimum sample depth ####
truncate_tra <- function(tra, min_depth=1, id="Time", split=NA)
{
  depth <- sample_depth_tra(tra, id, split)
  
  if (id %in% split){
    cname <- paste(id, "Split", sep="_")
    valid_indices <- lapply(depth, function(x){names(x[x>=min_depth])})
    
    subset_tra <- vector(mode="list", length=length(depth))
    names(subset_tra) <- names(depth)
    for (group in names(depth)){
      subset_tra[[group]] <- tra[tra[[id]] %in% valid_indices[[group]] & tra[[cname]]==group,]
    }
    
    trunc_tra <- Reduce(rbind, subset_tra)
    
  } else {
    valid_indices <- names(depth[depth>=min_depth])
    
    trunc_tra <- tra[tra[[id]] %in% valid_indices,]
  }

  return(trunc_tra)
}
  
  
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