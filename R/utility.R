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