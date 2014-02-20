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
  birth_years <- sapply(levels(tra$Tree), grab_birth_year, tra=tra, sparse)
  return (birth_years)
}


get_birth_index <- function(birth_year, tra)
{
  years <- sort(as.numeric(levels(tra$Time)))
  
  year_index <- which(birth_year==years)
  age_index <- 1
  birth_index <- year_index - age_index
  
  return (birth_index)
}