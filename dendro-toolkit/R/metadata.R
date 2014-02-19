# Grabbing ages of trees ####
grab_age <- function(i, tra, sparse=TRUE)
{
  if (sparse)
  {
    age <- max(as.numeric(as.character(tra[tra$i==i,"a"])))
  } else{
    age_row <- min(which(!is.na(tra[i,,]), arr.ind=T)[[2]])
    age <- dimnames(tra)[[3]][age_row]
  }
  return(as.numeric(age))
}


# Grabbing birth years of trees ####
grab_birth_year <- function(i, tra, sparse=TRUE)
{
  if (sparse)
  {
    slice <- tra[tra$i==i,][1,]
    birth_year <- as.numeric(as.character(slice$t)) -  as.numeric(as.character(slice$a)) + 1
  } else{
    slice <- which(!is.na(tra[i,,]), arr.ind=T)[1,]
    birth_index <- slice[1]-slice[2]+1
    birth_year <- as.numeric(dimnames(tra)[[2]][birth_index])
  }
  
  return(birth_year)
}


get_birth_years <- function(tra, sparse=TRUE)
{
  
  if (sparse){
    birth_years <- sapply(levels(tra$i), grab_birth_year, tra=tra, sparse)
  } else
  {
    birth_years <- sapply(1:dim(tra)[1], grab_birth_year, tra=tra, sparse)
    names(birth_years) <- dimnames(tra)[[1]]
  }
  
  return (birth_years)
}


get_birth_index <- function(birth_year, tra, sparse)
{
  if (sparse){
    years <- sort(as.numeric(levels(tra$t)))
  } else {
    years <- sort(as.numeric(dimnames(tra)[[2]]))
  }
  
  year_index <- which(birth_year==years)
  age_index <- 1
  birth_index <- year_index - age_index
  
  return (birth_index)
}