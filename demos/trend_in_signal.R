# Setting up effects ####
set.seed(42)
n_tree <- 20
n_time <- 20
n_age <- n_time
noise <- 0.1

constBAI_trend <- function (x, k, p0=0){
  # c is the constant basal area increment
  # p0 is the pith offset
  out <- sqrt(k*(x)/pi+p0^2) - sqrt(k*(x-1)/pi+p0^2)
  return (out)
}

# Generate effects
true_tree <- rlnorm(n_tree, sdlog=0)
names(true_tree) <- paste0("T", 1:n_tree)

true_time <- rlnorm(n_time, sdlog=0.2) * sin((1:n_time)/30)
names(true_time) <- (2001-n_time):(2000)

true_age <-  constBAI_trend(1:n_age, 1)
names(true_age) <- 1:n_age

# Combine effects
true_effects <- list(Tree=true_tree, Time=true_time, Age=true_age)

true_effects <- rescale_effects(true_effects, link="log")

# Make complete data set, balanced design ####

effect_names <- lapply(true_effects, function(x){names(x)})

complete_tra <- expand.grid(effect_names)
  
complete_tra$Growth <- NA

# Combine effects
for (i in 1:nrow(complete_tra)){
  Tree_i <- true_tree[complete_tra[i, "Tree"]]
  Time_i <- true_time[complete_tra[i, "Time"]]
  Age_i  <- true_age[complete_tra[i, "Age"]]
  complete_tra[i, "Growth"] <- Tree_i * Time_i * Age_i
}

# Add noise
complete_tra$Growth <- complete_tra$Growth * rlnorm(nrow(complete_tra),sdlog=noise)

# Make incomplete, randomly unbalanced design ####
# Typical trees would live for half of the chronology length
# Therefore # of points remaining should be around |t|/2*|i|

num_retained <- round(0.5*length(true_time)*length(true_tree))

rows_retained <- sample.int(nrow(complete_tra), num_retained)

random_incomplete_tra <- complete_tra[rows_retained,]

# Make realistic incomplete unbalanced design ####
# Trees were born at a random point in the chronology
# Live until time of sampling

date_of_birth <- sample(x=as.numeric(names(true_time)), size=length(true_tree), replace=TRUE)
names(date_of_birth) <- names(true_tree)

year_of_sampling <- max(as.numeric(names(true_time)))

get_valid_time_age_combos<- function(tree){
  valid_time_age_combos <- data.frame(Tree=tree, Time=date_of_birth[tree]:year_of_sampling, Age=1:(year_of_sampling-date_of_birth[tree]+1))
  return(valid_time_age_combos)
}

valid_df <- Reduce(rbind, lapply(names(true_tree), get_valid_time_age_combos))

# Select rows on the basis of ID
valid_df$ID <- paste0(valid_df$Tree, valid_df$Time, valid_df$Age)
complete_tra$ID <- paste0(complete_tra$Tree, complete_tra$Time, complete_tra$Age)

# Extract the final_values
realistic_incomplete_tra <- complete_tra[complete_tra$ID %in% valid_df$ID, ]
realistic_incomplete_tra <- realistic_incomplete_tra[-which(names(realistic_incomplete_tra)=="ID")]

# Saving datasets ####
trend_in_signal_dataset <- list(true_effects=true_effects, complete_tra=complete_tra, random_incomplete_tra=random_incomplete_tra,  realistic_incomplete_tra=realistic_incomplete_tra)

save(trend_in_signal_dataset, file="./Data/trend_in_signal_dataset.RData")

# Testing out various standardization options ####

# Loading in relevant dataset
# load(file="./Data/trend_in_signal_dataset.RData")

data(trend_in_signal_dataset)

# Check optimizers on the complete data
seq_1 <- standardize_tra(complete_tra, optim="sequential")
alt_1 <- standardize_tra(complete_tra, optim="alternate")
glm_1 <- standardize_tra(complete_tra, optim="glm")
gam_1 <- standardize_tra(complete_tra, optim="gam")

# Check optimizers on the random incomplete data
seq_2 <- standardize_tra(random_incomplete_tra, optim="sequential")
alt_2 <- standardize_tra(random_incomplete_tra, optim="alternate")
glm_2 <- standardize_tra(random_incomplete_tra, optim="glm")
gam_2 <- standardize_tra(random_incomplete_tra, optim="gam")

# Check optimizers on the random realistic data
seq_3 <- standardize_tra(realistic_incomplete_tra, optim="sequential")
alt_3 <- standardize_tra(realistic_incomplete_tra, optim="alternate")
mle_3 <- standardize_tra(realistic_incomplete_tra, optim="mle")
rss_3 <- standardize_tra(realistic_incomplete_tra, optim="rss")
glm_3 <- standardize_tra(realistic_incomplete_tra, optim="glm")
gam_3 <- standardize_tra(realistic_incomplete_tra, optim="gam")

# Timing ####

# Check optimizers on the complete data
system.time(standardize_tra(complete_tra, optim="sequential", make_plot=FALSE))
system.time(standardize_tra(complete_tra, optim="alternate", make_plot=FALSE))
system.time(standardize_tra(complete_tra, optim="glm", make_plot=FALSE))
system.time(standardize_tra(complete_tra, optim="gam", make_plot=FALSE))

# Check optimizers on the random incomplete data
system.time(standardize_tra(random_incomplete_tra, optim="sequential", make_plot=FALSE))
system.time(standardize_tra(random_incomplete_tra, optim="alternate", make_plot=FALSE))
system.time(standardize_tra(random_incomplete_tra, optim="glm", make_plot=FALSE))
system.time(standardize_tra(random_incomplete_tra, optim="gam", make_plot=FALSE))

# Check optimizers on the random realistic data
system.time(standardize_tra(realistic_incomplete_tra, optim="sequential", make_plot=FALSE))
system.time(standardize_tra(realistic_incomplete_tra, optim="alternate", make_plot=FALSE))
system.time(standardize_tra(realistic_incomplete_tra, optim="glm", make_plot=FALSE))
system.time(standardize_tra(realistic_incomplete_tra, optim="gam", make_plot=FALSE))

