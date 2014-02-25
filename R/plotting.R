# Automation and plots in bulk ####
# Make all of the possible plots at once!

make_standardization_plots <- function(effects, data, link="log"){
  
  plots <- list()
  
  # Sample depth
  sample_depth_time_plot <- make_sample_depth_plot(data$original, "Time")
  sample_depth_age_plot <- make_sample_depth_plot(data$original, "Age")
  
  plots <- c(plots, list(sample_depth_time_plot=sample_depth_time_plot, sample_depth_age_plot=sample_depth_age_plot))
  
  # Tree effect
  if ("Tree" %in% names(effects))
  {
    tree_effect_plot <- make_effect_plot(effects, "Tree")
    tree_effect_density_plot <- make_effect_density_plot(effects, "Tree", link)
    
    tree_effect_age_plot <- make_tree_effect_age_plot (effects, data$original)
    tree_effect_year_plot <- make_tree_effect_year_plot (effects, data$original)
    
    plots <- c(plots, list(tree_effect_plot=tree_effect_plot, tree_effect_density_plot=tree_effect_density_plot))
  }
  
  # Time effect
  if ("Time" %in% names(effects))
  {
    time_effect_plot <- make_effect_plot(effects, "Time", temporal=TRUE)
    time_effect_density_plot <- make_effect_density_plot(effects, "Time", link)
    
    plots <- c(plots, list(time_effect_plot=time_effect_plot, time_effect_density_plot=time_effect_density_plot))
  }
  
  # Age effect
  if ("Age" %in% names(effects))
  {
    age_effect_plot <- make_effect_plot(effects, "Age", temporal=TRUE)
    age_effect_density_plot <- make_effect_density_plot(effects, "Age", link)
    
    plots <- c(plots, list(age_effect_plot=age_effect_plot, age_effect_density_plot=age_effect_density_plot))
  }
  
  # Residuals
  residual_density_plot <- make_residual_density_plot(data$residuals, link)
  plots <- c(plots, list(residual_density_plot=residual_density_plot))
  
  # Return all plots as a list
  return (plots)
  
}

# Sample depth plotting ####
# Only defined / useful for time and age
make_sample_depth_plot <- function(tra, id="Time"){
  sample_depth <- sample_depth_tra(tra, id)
  
  df <- data.frame(sample_depth=sample_depth, id=levels(tra[[id]]))
  
  # Coerce id to numeric value
  df$id <- as.numeric(as.character(df$id))
  
  # Make the plot
  my_plot <- ggplot(df, aes(x=id, y=sample_depth)) + geom_area() + xlab(id) + ylab("Sample depth") + theme_bw()

  return(my_plot)
}

# Effects plotting ####

# Simple plotting of effects
make_effect_plot <- function(effects, effect_name, temporal=FALSE){
  
  # Load data into a data frame
  df <- data.frame(effect=effects[[effect_name]], id=names(effects[[effect_name]]))
  
  # Coerce temporal data ids to numeric form
  if (temporal){
    df$id <- as.numeric(as.character(df$id))
  }
  
  # Make the base graphic
  my_plot <- ggplot(df, aes(x=id,y=effect)) + theme_bw() + ylab(paste(effect_name, "effect")) + xlab(effect_name)
  
  # Nontemporal data is unordered and works best as a bar graph or similar
  # Temporal data is ordered and should be a line graph
  if (temporal){
    my_plot <- my_plot + geom_line()
  } else {
    my_plot <- my_plot + geom_bar(stat="identity")
  }  
  return(my_plot)
}

# Looking at the frequency of each magnitude of effects
make_effect_density_plot <- function(effects, effect_name, link){
  
  # Load data into a data frame
  my_effect <- effects[[effect_name]]
  df <- data.frame(effect=my_effect)
  
  # Make the base graphic
  my_plot <- ggplot(df, aes(x=effect)) + theme_bw() + ylab("Density estimate") + xlab(paste(effect_name, "effect")) + geom_density(colour="red", fill="red", alpha=0.5)
    
  # Estimate idealized PDF
  x_min <- min(min(my_effect), 0)
  x_max <- max(my_effect)
  
  x_ticks <- seq(from=x_min, to=x_max, length.out=100)
  
  # Generate PDF
  if (link=="identity"){
    pdf <- dnorm(x_ticks, mean=mean(my_effect),sd=sd(my_effect))    
  } else if (link=="log"){
    pdf <- dlnorm(x_ticks, meanlog=mean(log(my_effect)), sdlog=sd(log(my_effect)))
  }
  
  pdf_data <- data.frame(effect=x_ticks, density=pdf)
  
  # Add pdf information
  my_plot <- my_plot + geom_area(data=pdf_data, aes(x=effect, y=density), colour="black", fill="black", alpha=0.5)
  
  return(my_plot)
}

# Plotting tree effect vs age of tree at sampling to check for ecological effects of productivity on survival
make_tree_effect_age_plot <- function(effects, tra){
  
  ages <- sapply(names(effects$Tree), grab_age, tra=tra)
  
  df <- data.frame(tree=effects$Tree, age=ages)
  
  my_plot <- ggplot(df, aes(x=age, y=tree)) + geom_point() + geom_smooth() + theme_bw() + ylab("Tree effect") + xlab("Age at sampling") + geom_hline(y=1)
  
  return (my_plot)
}

# Plotting tree effect vs year of birth at sampling to check for modern sample bias
make_tree_effect_year_plot <- function(effects, tra){
  
  birth_years <- sapply(names(effects$Tree), grab_birth_year, tra=tra)
  
  df <- data.frame(tree=effects$Tree, birth_year=birth_years)
  
  my_plot <- ggplot(df, aes(x=birth_year, y=tree)) + geom_point() + geom_smooth() + theme_bw() + ylab("Tree effect") + xlab("Year of birth") + geom_hline(y=1)
  
  return (my_plot)
}

# Residuals plotting ####

# Do the residuals conform to our expectations of their density?
make_residual_density_plot <- function(residuals, link){
  
  # Load data into a data frame
  resid <- residuals$Growth[!is.na(residuals$Growth)]
  
  # Make the base graphic
  my_plot <- ggplot(data.frame(residuals=resid), aes(x=residuals)) + theme_bw() + ylab("Density estimate") + xlab("Residuals") + geom_density(colour="red", fill="red", alpha=0.5)
  
  # Estimate idealized PDF
  x_min <- min(min(resid), 0)
  x_max <- max(resid)
  
  x_ticks <- seq(from=x_min, to=x_max, length.out=100)
  
  # Generate PDF
  if (link=="identity"){
    pdf <- dnorm(x_ticks, sd(resid))    
  } else if (link=="log"){
    pdf <- dlnorm(x_ticks, sdlog=sd(log(resid)))
  }
  
  pdf_data <- data.frame(effect=x_ticks, density=pdf)
  
  # Add pdf information
  my_plot <- my_plot + geom_area(data=pdf_data, aes(x=effect, y=density), colour="black", fill="black", alpha=0.5)

  return(my_plot)
}