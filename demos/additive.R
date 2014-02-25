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
true_tree <- rnorm(n_tree, sd=0)
names(true_tree) <- paste0("T", 1:n_tree)

true_time <- rnorm(n_time, sd=0.2) + sin((1:n_time)/30)
names(true_time) <- (2001-n_time):(2000)

true_age <-  constBAI_trend(1:n_age, 1)
names(true_age) <- 1:n_age

# Combine effects
true_effects <- list(Tree=true_tree, Time=true_time, Age=true_age)

true_effects <- rescale_effects(true_effects, link="identity")

# Make complete data set, balanced design ####

effect_names <- lapply(true_effects, function(x){names(x)})

complete_tra <- expand.grid(effect_names)

complete_tra$Growth <- NA

# Combine effects
for (i in 1:nrow(complete_tra)){
  Tree_i <- true_tree[complete_tra[i, "Tree"]]
  Time_i <- true_time[complete_tra[i, "Time"]]
  Age_i  <- true_age[complete_tra[i, "Age"]]
  complete_tra[i, "Growth"] <- Tree_i + Time_i + Age_i
}

# Add noise
complete_tra$Growth <- complete_tra$Growth + rnorm(nrow(complete_tra),sd=noise)

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
additive_dataset <- list(true_effects=true_effects, complete_tra=complete_tra, random_incomplete_tra=random_incomplete_tra,  realistic_incomplete_tra=realistic_incomplete_tra)

save(additive_dataset, file="./Data/additive_dataset.RData")

# Testing out various standardization options ####

# Loading in relevant dataset
# load(file="./Data/additive_dataset.RData")

data(additive_dataset)

# Check optimizers on the complete data
seq_1 <- standardize_tra(complete_tra, optim="sequential", link="identity")
alt_1 <- standardize_tra(complete_tra, optim="alternate", link="identity")
glm_1 <- standardize_tra(complete_tra, optim="glm", link="identity")
gam_1 <- standardize_tra(complete_tra, optim="gam", link="identity")

# Check optimizers on the random incomplete data
seq_2 <- standardize_tra(random_incomplete_tra, optim="sequential", link="identity")
alt_2 <- standardize_tra(random_incomplete_tra, optim="alternate", link="identity")
glm_2 <- standardize_tra(random_incomplete_tra, optim="glm", link="identity")
gam_2 <- standardize_tra(random_incomplete_tra, optim="gam", link="identity")

# Check optimizers on the random realistic data
seq_3 <- standardize_tra(realistic_incomplete_tra, optim="sequential", link="identity")
alt_3 <- standardize_tra(realistic_incomplete_tra, optim="alternate", link="identity")
glm_3 <- standardize_tra(realistic_incomplete_tra, optim="glm", link="identity")
gam_3 <- standardize_tra(realistic_incomplete_tra, optim="gam", link="identity")


# Model fit ####
model_fit_1 <- data.frame(seq=unlist(seq_1$fit),
                          alt=unlist(alt_1$fit),
                          glm=unlist(glm_1$fit),
                          gam=unlist(gam_1$fit)
                )

model_fit_2 <- data.frame(seq=unlist(seq_2$fit),
                          alt=unlist(alt_2$fit),
                          glm=unlist(glm_2$fit),
                          gam=unlist(gam_2$fit)
)

model_fit_3 <- data.frame(seq=unlist(seq_3$fit),
                          alt=unlist(alt_3$fit),
                          glm=unlist(glm_3$fit),
                          gam=unlist(gam_3$fit)
)

print(model_fit_1)
print(model_fit_2)
print(model_fit_3)

# Collating data ####

# Setup
model_list <- c("True", "Seq", "Alt", "GLM", "GAM")

# Tree_model <- rep(model_list, each=length(true_effects$Tree))
Time_model <- rep(model_list, each=length(true_effects$Time))
Age_model <- rep(model_list, each=length(true_effects$Age))

# Tree_id <- rep(true_effects$Tree, times=length(model_list))
Time_id <- rep(names(true_effects$Time), times=length(model_list))
Age_id <- rep(names(true_effects$Age), times=length(model_list))


# Case 1
# Tree_1 <- data.frame(effect=c(true_effects$Tree, seq_1$effects$Tree, alt_1$effects$Tree, glm_1$effects$Tree, gam_1$effects$Tree), model=Tree_model, id=Tree_id, case=1)

Time_1 <- data.frame(effect=c(true_effects$Time, seq_1$effects$Time, alt_1$effects$Time, glm_1$effects$Time, gam_1$effects$Time), model=Time_model, id=Time_id, case=1)

Age_1 <- data.frame(effect=c(true_effects$Age, seq_1$effects$Age, alt_1$effects$Age, glm_1$effects$Age, gam_1$effects$Age), model=Age_model, id=Age_id, case=1)
                    
# Case 2
# Tree_2 <- data.frame(effect=c(true_effects$Tree, seq_2$effects$Tree, alt_2$effects$Tree, glm_2$effects$Tree, gam_2$effects$Tree), model=Tree_model, id=Tree_id, case=2)

Time_2 <- data.frame(effect=c(true_effects$Time, seq_2$effects$Time, alt_2$effects$Time, glm_2$effects$Time, gam_2$effects$Time), model=Time_model, id=Time_id, case=2)

Age_2 <- data.frame(effect=c(true_effects$Age, seq_2$effects$Age, alt_2$effects$Age, glm_2$effects$Age, gam_2$effects$Age), model=Age_model, id=Age_id, case=2)

# Case 3
# Tree_3 <- data.frame(effect=c(true_effects$Tree, seq_3$effects$Tree, alt_3$effects$Tree, glm_3$effects$Tree, gam_3$effects$Tree), model=Tree_model, id=Tree_id, case=3)

Time_3 <- data.frame(effect=c(true_effects$Time, seq_3$effects$Time, alt_3$effects$Time, glm_3$effects$Time, gam_3$effects$Time), model=Time_model, id=Time_id, case=3)

Age_3 <- data.frame(effect=c(true_effects$Age, seq_3$effects$Age, alt_3$effects$Age, glm_3$effects$Age, gam_3$effects$Age), model=Age_model, id=Age_id, case=3)

# Collating cases
# Tree_df <- rbind(Tree_1, Tree_2, Tree_3)
Time_df <- rbind(Time_1, Time_2, Time_3)
Age_df <- rbind(Age_1, Age_2, Age_3)

Time_df$id <- as.numeric(as.character(Time_df$id))
Age_df$id <- as.numeric(as.character(Age_df$id))

# Tree_df_nt <- Tree_df[-Tree_df$model=="True",]
Time_df_diff <- Time_df[-which(Time_df$model=="True"),]
Age_df_diff <- Age_df[-which(Age_df$model=="True"),]

Time_df_diff$effect <- Time_df_diff$effect - true_effects$Time
Age_df_diff$effect <- Age_df_diff$effect - true_effects$Age

# Summary plots ####

# Raw effects
# Tree_raw_plot <- ggplot(Tree_df, aes(x=id, y=effect)) + geom_bar(stat="identity")  + facet_grid(case~model) + theme_bw() + xlab("Tree") + ylab("Effect") + geom_h_line(x=0)

Time_raw_plot <- ggplot(Time_df, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Year") + ylab("Effect") + geom_hline(y=0)

Age_raw_plot <- ggplot(Age_df, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Age") + ylab("Effect")

# Relative to true
# Tree_scatter_plot <- ggplot(Tree_df_nt, aes(x=effect, y=true_effects$Tree)) + geom_point()  + facet_grid(case~model) + theme_bw() + xlab("Estimated Tree Effect") + ylab("True Tree Effect") + geom_smooth() + geom_abline(intercept=0, slope=1)

Time_diff_plot <- ggplot(Time_df_diff, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Year") + ylab("Difference between estimated and true effect") + geom_hline(y=0)

Age_diff_plot <- ggplot(Age_df_diff, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Age") + ylab("Difference between estimated and true effect") + geom_hline(y=0)

# Printing
# print(Tree_raw_plot)
print(Time_raw_plot)
print(Age_raw_plot)

# print(Tree_scatter_plot)
print(Time_diff_plot)
print(Age_diff_plot)

# Performance ####

# Check optimizers on the complete data
system.time(standardize_tra(complete_tra, optim="sequential", make_plot=FALSE, link="identity"))
system.time(standardize_tra(complete_tra, optim="alternate", make_plot=FALSE, link="identity"))
system.time(standardize_tra(complete_tra, optim="glm", make_plot=FALSE, link="identity"))
system.time(standardize_tra(complete_tra, optim="gam", make_plot=FALSE, link="identity"))

# Check optimizers on the random incomplete data
system.time(standardize_tra(random_incomplete_tra, optim="sequential", make_plot=FALSE, link="identity"))
system.time(standardize_tra(random_incomplete_tra, optim="alternate", make_plot=FALSE, link="identity"))
system.time(standardize_tra(random_incomplete_tra, optim="glm", make_plot=FALSE, link="identity"))
system.time(standardize_tra(random_incomplete_tra, optim="gam", make_plot=FALSE, link="identity"))

# Check optimizers on the random realistic data
system.time(standardize_tra(realistic_incomplete_tra, optim="sequential", make_plot=FALSE, link="identity"))
system.time(standardize_tra(realistic_incomplete_tra, optim="alternate", make_plot=FALSE, link="identity"))
system.time(standardize_tra(realistic_incomplete_tra, optim="glm", make_plot=FALSE, link="identity"))
system.time(standardize_tra(realistic_incomplete_tra, optim="gam", make_plot=FALSE, link="identity"))
