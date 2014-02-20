# Demonic ritual ####
demonic_ritual <- function(effects, tra=TRUE){
  birth_years <- get_birth_years(tra)
  birth_index <- sapply(birth_years, get_birth_index, tra=tra)
  k <- -birth_index
  
 years <- sort(as.numeric(levels(tra$Time)))
  
  year_index <- match(names(effects$Time), as.character(years))
  names(year_index) <- names(effects$Time)
  
  ages <- sort(as.numeric(levels(tra$Age)))
  
  age_index <- match(names(effects$Age), as.character(ages))
  names(age_index) <- names(effects$Age)
  
  k <- k[!is.na(k)]
  year_index <- year_index[!is.na(year_index)]
  age_index <- age_index[!is.na(age_index)]
  
  ridge_ids <- list(I=k, T=year_index, A=age_index)
  
  return(ridge_ids)
  
}


# Add a demonic intrusion ####
add_demons <- function (m=1, effects, tra=TRUE, link="log"){
  ridge_ids <- demonic_ritual(effects, tra)
  
  corrupted_effects <- effects
  
  if(link=="log")
  {
    corrupted_effects <- lapply(corrupted_effects, log)
  }
  
  
  # I = I +m*k
  corrupted_effects$Timeree <- corrupted_effects$Timeree + m*ridge_ids$Timeree
  
  # T = T + m*C
  corrupted_effects$Time <- corrupted_effects$Time + m*ridge_ids$Time
  
  # A = A - m*R
  corrupted_effects$Age <- corrupted_effects$Age - m*ridge_ids$Age
  
  if(link=="log")
  {
    corrupted_effects <- lapply(corrupted_effects, exp)
  }
  
  corrupted_effects <- rescale_effects(corrupted_effects, link)
  return(corrupted_effects)
}

# Estimate a demonic intrusion ####
estimate_demons <- function (corrupted_effects, tra=TRUE, link="log"){
  ridge_ids <- demonic_ritual(corrupted_effects, tra)
  
  
  if(link=="log")
  {
    corrupted_effects <- lapply(corrupted_effects, log)
  }
  
  # Assuming normal random variation around the demonic line
  demonic_llh <- function (purified_effects)
  {
    sd_I <- sd(purified_effects$Timeree)
    sd_T <- sd(purified_effects$Time)
    sd_A <- sd(purified_effects$Age)
    
    llh_I <- sum(dnorm(purified_effects$Timeree, sd=sd_I, log=TRUE))
    llh_T <- sum(dnorm(purified_effects$Time, sd=sd_T, log=TRUE))
    llh_A <- sum(dnorm(purified_effects$Age, mean=mean(purified_effects$Age), sd=sd_A, log=TRUE))
    
    llh <- sum(llh_I, llh_T, llh_A)
    
    return(llh)
  } 

  purify_probe <- function(m){
    purified_effects <- remove_demons(m, corrupted_effects, tra, link="additive")
    llh <- demonic_llh(purified_effects)
    return(-llh)
  }
  
  optimal_m <- optimize(f=purify_probe, interval=c(-1, 1))$minimum
  
  return(optimal_m)
}


# Remove a demonic intrusion ####
remove_demons <- function (m=1, corrupted_effects, tra=TRUE, link="log"){
  ridge_ids <- demonic_ritual(corrupted_effects, tra)
  
  effects <- corrupted_effects
  
  if(link=="log")
  {
   effects <- lapply(effects, log)
  }
  
  # I = I + m*k
  effects$Timeree <- effects$Timeree - m*ridge_ids$Timeree
  
  # T = T + m*C
  effects$Time <- effects$Time - m*ridge_ids$Time
  
  # A = A - m*R
  effects$Age <- effects$Age + m*ridge_ids$Age
  
  if(link=="log")
  {
    effects <- lapply(effects, exp)
  }
  
  effects <- rescale_effects(effects, link)
  return(effects)
}

# Full post-hoc correction of demonic intrusion ####
post_hoc_intercession <- function(corrupted_effects, tra=TRUE, link=="log")
{
  optimal_m <- estimate_demons(corrupted_effects, tra, link)
  
  # Ignore correction in the case of suspiciously high trends
  if (abs(optimal_m) > 0.99){
    return(corrupted_effects)
  }
  
  print(paste("m of", optimal_m, "selected in choosing solution on likelihood ridge."))
  purified <- remove_demons(optimal_m, corrupted_effects, tra, link)
  
  
  return(purified)
}

# Testing ####
#effects <- ta_2$effects
#effects <- ita_2$effects
#tra <- ta_2$Timera

#corrupted_effects <- add_demons(m=1, effects, tra, link="additive")
#purified <- remove_demons(m=1, corrupted_effects, tra, link="additive")

#corrupted_effects <- add_demons(m=1, effects, tra, link="log")
#purified <- remove_demons(m=1, corrupted_effects, tra, link="log")

#purified <-  remove_demons(m=-0.00633, effects, tra, link="log")
#purified <-  post_hoc_intercession(effects, tra, link="log")

#purified <- post_hoc_intercession(corrupted_effects, tra, link="log")


#plot(ta_2$effects$Age)
#plot(purified$Time)
#plot(purified$Age)