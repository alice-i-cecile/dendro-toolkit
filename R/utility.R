# Truncate tree ring array to ensure minimum sample depth ####
truncate_tra <- function(tra, min_depth=1, id="Time", group_by=NA)
{
  depth <- sample_depth_tra(tra, id, group_by)
    
  valid_indices <- names(depth[depth>=min_depth])
  
  trunc_tra <- tra[tra[[id]] %in% valid_indices,]
  
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