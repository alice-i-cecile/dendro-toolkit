# Automation and plots in bulk ####
# Make all of the possible plots at once!

make_standardization_plots <- function(effects, se=NULL, dat, group_by=NA,link="log", dep_var="Growth", ci_size=0.95){
  
  plots <- list()
  
  # Sample depth
  sample_depth_time_plot <- make_sample_depth_plot(dat$original, "Time", group_by)
  sample_depth_age_plot <- make_sample_depth_plot(dat$original, "Age", group_by)
  
  plots <- c(plots, list(sample_depth_time_plot=sample_depth_time_plot, sample_depth_age_plot=sample_depth_age_plot))
  
  # Tree effect
  if ("Tree" %in% names(effects))
  {
    tree_effect_plot <- make_effect_plot(effects, se, "Tree", FALSE, link, ci_size, group_by)
    tree_effect_density_plot <- make_effect_density_plot(effects, "Tree", link, group_by)
    
    tree_effect_age_plot <- make_tree_effect_age_plot(effects, se, dat$original, link, ci_size)
    tree_effect_year_plot <- make_tree_effect_year_plot(effects, se, dat$original, link, ci_size)
    
    plots <- c(plots, list(tree_effect_plot=tree_effect_plot, tree_effect_age_plot=tree_effect_age_plot, tree_effect_year_plot=tree_effect_year_plot, tree_effect_density_plot=tree_effect_density_plot))
  }
  
  # Time effect
  if ("Time" %in% names(effects))
  {
    time_effect_plot <- make_effect_plot(effects, se, "Time", temporal=TRUE, link,  ci_size, group_by)
    time_effect_density_plot <- make_effect_density_plot(effects, "Time", link, group_by)
    
    plots <- c(plots, list(time_effect_plot=time_effect_plot, time_effect_density_plot=time_effect_density_plot))
  }
  
  # Age effect
  if ("Age" %in% names(effects))
  {
    age_effect_plot <- make_effect_plot(effects, se, "Age", temporal=TRUE, link, ci_size, group_by)
    age_effect_density_plot <- make_effect_density_plot(effects, "Age", link, group_by)
    
    plots <- c(plots, list(age_effect_plot=age_effect_plot, age_effect_density_plot=age_effect_density_plot))
  }
  
  # Residuals
  residual_density_plot <- make_residual_density_plot(dat$residuals, link, dep_var)
  plots <- c(plots, list(residual_density_plot=residual_density_plot))
  
  # Return all plots as a list
  return (plots)
  
}

# Sample depth plotting ####
# Only defined / useful for time and age
make_sample_depth_plot <- function(tra, id="Time", group_by=NA){
  
  # Extract and format sample depth
  sample_depth <- sample_depth_tra(tra, id, group_by)
  
  if(id %in% group_by)
  {
    groups <- names (sample_depth)
    dat <- vector(mode="list", length=length(groups))
    names(dat) <- groups
  
    for (group in groups){
      dat[[group]] <- data.frame(sample_depth=sample_depth[[group]], id=names(sample_depth[[group]]), group=group)
    }
    
    dat <- Reduce(rbind, dat)
  
    } else {
    dat <- data.frame(sample_depth=sample_depth, id=names(sample_depth))
  }

  # Coerce id to numeric value
  dat$id <- as.numeric(as.character(dat$id))
  
  # Make the plot
  my_plot <- ggplot(dat, aes(x=id, y=sample_depth)) + geom_area() + xlab(id) + ylab("Sample depth") + theme_bw()
  
  # Facet different groups
  if(id %in% group_by)
  {
    my_plot <- my_plot + facet_grid(group~.)
  }

  return(my_plot)
}

# Effects plotting ####

# Simple plotting of effects
make_effect_plot <- function(effects, se=NULL, id, temporal=FALSE, link="log", ci_size=0.95, group_by=NA){
  
  # Load data into a data frame
  if(id %in% group_by)
  {
    dat <- vector(mode="list", length=length(effects[[id]]))
    groups <- names (effects[[id]])
    names(dat) <- groups
    
    for (group in groups){
      dat[[group]] <- data.frame(effect=effects[[id]][[group]], id=names(effects[[id]][[group]]), group=group)
      
      if(!is.null(se)){
        dat[[group]]$se <- se[[id]][[group]]
      }
    }
    
    dat <- Reduce(rbind, dat)
    
  } else {
    dat <- data.frame(effect=effects[[id]], id=names(effects[[id]]))
    
    if(!is.null(se)){
      dat$se <- se[[id]]
    }
  }
  
  # Coerce temporal data ids to numeric form
  if (temporal){
    dat$id <- as.numeric(as.character(dat$id))
  }
  
  # Transform standard errors for coefficients into confidence intervals
  if(!is.null(se)){
 
    # qnorm requires % of values that are included in each tail
    # Therefore half should be missing on each side
    ci_Z <- 1 - (1-ci_size)/2
    
    
    if (link=="identity"){
      dat$upper <- dat$effect + dat$se*qnorm(ci_Z)
      dat$lower <- dat$effect - dat$se*qnorm(ci_Z)
    } else if(link=="log"){
      dat$upper <- exp(log(dat$effect) + dat$se*qnorm(ci_Z))
      dat$lower <- exp(log(dat$effect) - dat$se*qnorm(ci_Z))
    }
  }
  
  # Make the base graphic
  if(!is.null(se)){
    my_plot <- ggplot(dat, aes(x=id,y=effect, ymax=upper, ymin=lower))
  } else {
    my_plot <- ggplot(dat, aes(x=id,y=effect))
  }
  my_plot <- my_plot + theme_bw() + ylab(paste(id, "effect")) + xlab(id)
  
  # Nontemporal data is unordered and works best as a bar graph or similar
  # Temporal data is ordered and should be a line graph
  if (temporal){
    my_plot <- my_plot + geom_line()
    
    # Add ribbon to show CI
    if(!is.null(se)){
      my_plot <- my_plot + geom_ribbon(alpha=0.5)
    }
  } else {
    my_plot <- my_plot + geom_bar(stat="identity")
    
    # Add linerange to show CI
    if(!is.null(se)){
      my_plot <- my_plot + geom_linerange()
    }
  }  
  
  # Facet plot by group
  if(id %in% group_by)
  {
    my_plot <- my_plot + facet_grid(group~.)
  }  
  
  return(my_plot)
}

# Looking at the frequency of each magnitude of effects
make_effect_density_plot <- function(effects, id, link, group_by=NA){
  
  # Load data into a data frame
  if(id %in% group_by)
  {
    dat <- vector(mode="list", length=length(effects[[id]]))
    groups <- names (effects[[id]])
    names(dat) <- groups
    
    for (group in groups){
      dat[[group]] <- data.frame(effect=effects[[id]][[group]], group=group)
    }
    
    dat <- Reduce(rbind, dat)
    
  } else {
    dat <- data.frame(effect=effects[[id]])
  }
  
  # Make the base graphic
  my_plot <- ggplot(dat, aes(x=effect)) + theme_bw() + ylab("Density estimate") + xlab(paste(id, "effect")) + geom_density(colour="red", fill="red", alpha=0.5)
    
  if (id %in% group_by){
    x_min <- min(min(unlist(effects[[id]])), 0)
    x_max <- max(unlist(effects[[id]]))
    
    x_ticks <- seq(from=x_min, to=x_max, length.out=100)
    
    pdf_data <- vector(mode="list", length=length(groups))
    names(pdf_data) <- groups
    
    for (group in groups){
      # Estimate theoretical PDF
      if (link=="identity"){
        pdf <- dnorm(x_ticks, mean=mean(effects[[id]][[group]]),sd=sd(effects[[id]][[group]]))    
      } else if (link=="log"){
        pdf <- dlnorm(x_ticks, meanlog=mean(log(effects[[id]][[group]])), sdlog=sd(log(effects[[id]][[group]])))
      }
      
      pdf_data[[group]] <- data.frame(effect=x_ticks, density=pdf, group=group)
    }
    
    pdf_data <- Reduce(rbind, pdf_data)
    
  } else {
    x_min <- min(min(effects[[id]]), 0)
    x_max <- max(effects[[id]])
    
    x_ticks <- seq(from=x_min, to=x_max, length.out=100)
    
    # Estimate theoretical PDF
    if (link=="identity"){
      pdf <- dnorm(x_ticks, mean=mean(effects[[id]]),sd=sd(effects[[id]]))    
    } else if (link=="log"){
      pdf <- dlnorm(x_ticks, meanlog=mean(log(effects[[id]])), sdlog=sd(log(effects[[id]])))
    }
    
    pdf_data <- data.frame(effect=x_ticks, density=pdf)
  }

  
  # Add PDF
  my_plot <- my_plot + geom_area(data=pdf_data, aes(x=effect, y=density), colour="black", fill="black", alpha=0.5)
  
  # Facet plot by group
  if(id %in% group_by)
  {
    my_plot <- my_plot + facet_grid(group~.)
  }
  
  return(my_plot)
}

# Plotting tree effect vs age of tree at sampling to check for ecological effects of productivity on survival
make_tree_effect_age_plot <- function(effects, se=NULL, tra, link, ci_size=0.95){
  
  ages <- sapply(names(effects$Tree), grab_age, tra=tra)
  
  dat <- data.frame(tree=effects$Tree, age=ages)
  
  
  # Transform standard errors for coefficients into confidence intervals
  if(!is.null(se)){
    dat$se <- se$Tree
    
    # qnorm requires % of values that are included in each tail
    # Therefore half should be missing on each side
    ci_Z <- 1 - (1-ci_size)/2
    
    
    if (link=="identity"){
      dat$upper <- dat$tree + dat$se*qnorm(ci_Z)
      dat$lower <- dat$tree - dat$se*qnorm(ci_Z)
    } else if(link=="log"){
      dat$upper <- exp(log(dat$tree) + dat$se*qnorm(ci_Z))
      dat$lower <- exp(log(dat$tree) - dat$se*qnorm(ci_Z))
    }
  }
  
  # Base plot
  if (!is.null(se)){
    my_plot <- ggplot(dat, aes(x=age, y=tree,  ymax=upper, ymin=lower))
  } else {
    my_plot <- ggplot(dat, aes(x=age, y=tree))
  }
   my_plot <- my_plot + geom_point() + geom_smooth() + theme_bw() + ylab("Tree effect") + xlab("Age at sampling")
  
  # Add linerange to show CI
  if(!is.null(se)){
    my_plot <- my_plot + geom_linerange()
  }
  
  # Horizontal line to show typical value
  if (link=="log"){
    my_plot <- my_plot + geom_hline(y=1)
  } else {
    my_plot <- my_plot + geom_hline(y=0)
  }
  
  
  return (my_plot)
}

# Plotting tree effect vs year of birth at sampling to check for modern sample bias
make_tree_effect_year_plot <- function(effects,  se=NULL, tra, link="log", ci_size=0.95){
  
  birth_years <- sapply(names(effects$Tree), grab_birth_year, tra=tra)
  
  dat <- data.frame(tree=effects$Tree, birth_year=birth_years)
  
  # Transform standard errors for coefficients into confidence intervals
  if(!is.null(se)){
    dat$se <- se$Tree
    
    # qnorm requires % of values that are included in each tail
    # Therefore half should be missing on each side
    ci_Z <- 1 - (1-ci_size)/2
    
    
    if (link=="identity"){
      dat$upper <- dat$tree + dat$se*qnorm(ci_Z)
      dat$lower <- dat$tree - dat$se*qnorm(ci_Z)
    } else if(link=="log"){
      dat$upper <- exp(log(dat$tree) + dat$se*qnorm(ci_Z))
      dat$lower <- exp(log(dat$tree) - dat$se*qnorm(ci_Z))
    }
  }
  
  # Base plot
  if (!is.null(se)){
    my_plot <- ggplot(dat, aes(x=birth_year, y=tree,  ymax=upper, ymin=lower))
  } else {
    my_plot <- ggplot(dat, aes(x=birth_year, y=tree))
  }
  
  my_plot <-my_plot + geom_point() + geom_smooth() + theme_bw() + ylab("Tree effect") + xlab("Year of birth")
  
  # Add linerange to show CI
  if(!is.null(se)){
    my_plot <- my_plot + geom_linerange()
  }
  
  # Horizontal line to show typical value
  if (link=="log"){
    my_plot <- my_plot + geom_hline(y=1)
  } else {
    my_plot <- my_plot + geom_hline(y=0)
  }
  
  return (my_plot)
}

# Residuals plotting ####

# Do the residuals conform to our expectations of their density?
make_residual_density_plot <- function(residuals, link="log", dep_var="Growth"){
  
  # Load data into a data frame
  resid <- residuals[[dep_var]][!is.na(residuals[[dep_var]])]
  
  # Make the base graphic
  my_plot <- ggplot(data.frame(residuals=resid), aes(x=residuals)) + theme_bw() + ylab("Density estimate") + xlab("Residuals") + geom_density(colour="red", fill="red", alpha=0.5)
  
  # Estimate idealized pdf
  x_min <- min(min(resid), 0)
  x_max <- max(resid)
  
  x_ticks <- seq(from=x_min, to=x_max, length.out=100)
  
  # Generate pdf
  if (link=="identity"){
    pdf <- dnorm(x_ticks, sd=sd(resid))    
  } else if (link=="log"){
    pdf <- dlnorm(x_ticks, sdlog=sd(log(resid)))
  }
  
  pdf_data <- data.frame(effect=x_ticks, density=pdf)
  
  # Add pdf information
  my_plot <- my_plot + geom_area(data=pdf_data, aes(x=effect, y=density), colour="black", fill="black", alpha=0.5)

  # Add line showing null value
  if (link=="log"){
    my_plot <- my_plot + geom_vline(x=1)
  } else {
    my_plot <- my_plot + geom_vline(x=0)
  }
  
  return(my_plot)
}