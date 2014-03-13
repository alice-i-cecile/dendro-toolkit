standardize_direct_search <- function (tra, model=c("Time", "Age"), split=NA, link="log", dep_var="Growth", criteria="llh", ...)
{
  
  # Skeleton effects for starting and relisting
  skele <- make_skeleton_effects(tra, model, split, link)
  
  # Function to optimize
  # Finds model fit given effects
  search_wrapper <- function(boneless_effects){
    effects <- relist(boneless_effects, skele)
    
    predicted <- predicted_tra(effects, tra, model, split, link, dep_var)
    
    residuals <- residuals_tra(tra, predicted, link, dep_var)
    
    if (criteria=="llh"){
      llh <- llh_tra(residuals, link, dep_var)
      print(llh)
      return(-llh)
    } else if(criteria=="rss"){
      rss <- rss_tra(residuals, link, dep_var)
      print(rss)      
      return(rss)
    }
  }
  
  # Search the solution space
  search_solution <- optim(unlist(skele), search_wrapper)  
    
  # Reshape and return the optimal effects
  effects <- relist(search_solution$par, skele)
  
  return(effects)
}