# Setting up effects ####
set.seed(42)

# Trees
n_tree <- 50

# Time
n_time <- 1000
start_year <- 1000
end_year <- start_year + n_time - 1

# Age
n_age <- n_time
typical_lifespan <- 300
constBAI_trend <- function (x, k=1, p0=0){
  # c is the constant basal area increment
  # p0 is the pith offset
  out <- sqrt(k*(x)/pi+p0^2) - sqrt(k*(x-1)/pi+p0^2)
  return (out)
}

logistic_trend <- function(x, a=1, k=1/20, t=0.1){
  # a is scaling
  #t is asymptote
  
  out <- a / (1+exp(1)^(k*x)) + t
  return(out)
}

# Sampling
sampling_weights <- 1 / (1+exp(1)^(-((-n_time/2 - 500):(n_time/2-1-500))/300))
names(sampling_weights) <- start_year:end_year
# plot(sampling_weights)

# Noise
noise <- 0.1

# Generate effects
true_tree <- rlnorm(n_tree, sdlog=0)
names(true_tree) <- paste0("T", 1:n_tree)

true_time <- exp(rnorm(n_time, sd=0.2) + 0.15*sin((1:n_time)/30) + 1*sin((1:n_time)/500))
names(true_time) <- start_year:end_year
# plot(true_time)

true_age <-  list(A=constBAI_trend(1:n_age), B=logistic_trend(1:n_age))
true_age <- lapply(true_age, function(x){names(x) <- 1:n_age; x})

# Assign trees to approriate age group
age_split_df <- data.frame(Tree=names(true_tree), Age_Split=sample(c("A", "B"), size=n_tree, replace=T))

# Combine effects
raw_effects <- list(Tree=true_tree, Time=true_time, Age=true_age)

# Setting up demographics and structure ####

# When were the trees born
date_of_birth <- sample(x=start_year:end_year, size=n_tree, replace=TRUE, prob=sampling_weights)
names(date_of_birth) <- names(true_tree)

# How long did they live for
lifespan <- rpois(n=n_tree, lambda=typical_lifespan)
names(lifespan) <- names(true_tree)

# Translate this information into a structure
generate_demographic_structure <- function(tree){
  timespan <- date_of_birth[tree]:(date_of_birth[tree] + lifespan[tree] - 1)
  
  valid_time_age_combos <- data.frame(Tree=tree, Time=timespan, Age=1:lifespan[tree])
  return(valid_time_age_combos)
}

tra_structure <- Reduce(rbind, lapply(names(true_tree), generate_demographic_structure))

# Exclude data points from after the year of sampling
tra_structure <- tra_structure[!(tra_structure$Time > end_year), ]

# Finalizing dataset ####

# Adding information on age split
tra_structure <- merge(tra_structure, age_split_df, by="Tree")

# Setup final tree ring array
tra <- tra_structure
tra$Growth <- NA

# Combine effects
for (i in 1:nrow(tra)){
  Tree_i <- true_tree[tra[i, "Tree"]]
  Time_i <- true_time[tra[i, "Time"]]
  Age_i  <- true_age[[tra[i, "Age_Split"]]][tra[i, "Age"]]
  tra[i, "Growth"] <- Tree_i * Time_i * Age_i
}

# Add noise
tra$Growth <- tra$Growth * rlnorm(nrow(tra),sdlog=noise)

# Truncate true effects to match structure
skele_effects <- make_skeleton_effects(tra, model=c("Tree", "Time", "Age"), split="Age", link="log")

true_effects <- rescale_effects(synchronize_effects(raw_effects, skele_effects, "Age"), link="log", split="Age")

# Saving dataset ####
segment_length_curse_dataset <- list(true_effects=true_effects, tra=tra)

save(segment_length_curse_dataset, file="./Data/segment_length_curse_dataset.RData")

# Loading in relevant dataset ####
load(file="./Data/segment_length_curse_dataset.RData")
# data(segment_length_curse_dataset)

true_effects <- segment_length_curse_dataset$true_effects

tra <- segment_length_curse_dataset$tra

# Cleaning dataset ####

# Truncate by time
tra <- truncate_tra(tra, min_depth=3, id="Time", split="Age")

# Truncate by age
tra <- truncate_tra(tra, min_depth=3, id="Age", split="Age")

# Clean up true effects to match
trunc_skele_effects <- make_skeleton_effects(tra, model=c("Tree", "Time", "Age"), split="Age", link="log")

true_effects <- rescale_effects(synchronize_effects(true_effects, trunc_skele_effects, "Age"), link="log", split="Age")

# Testing out various standardization options ####

# No splitting
seq_1 <- standardize_tra(tra, optim="sequential")
alt_1 <- standardize_tra(tra, optim="alternate")
glm_1 <- standardize_tra(tra, optim="glm")
gam_1 <- standardize_tra(tra, optim="gam")

# True splitting
seq_2 <- standardize_tra(tra, optim="sequential", split="Age")
alt_2 <- standardize_tra(tra, optim="alternate", split="Age")
glm_2 <- standardize_tra(tra, optim="glm", split="Age")
gam_2 <- standardize_tra(tra, optim="gam", split="Age")

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

# Time_3 <- data.frame(effect=c(true_effects$Time, seq_3$effects$Time, alt_3$effects$Time, glm_3$effects$Time, gam_3$effects$Time), model=Time_model, id=Time_id, case=3)
# 
# Age_3 <- data.frame(effect=c(true_effects$Age, seq_3$effects$Age, alt_3$effects$Age, glm_3$effects$Age, gam_3$effects$Age), model=Age_model, id=Age_id, case=3)

# Collating cases
# Tree_df <- rbind(Tree_1, Tree_2, Tree_3)
Time_df <- rbind(Time_1, Time_2)#, Time_3)
Age_df <- rbind(Age_1, Age_2)#, Age_3)

Time_df$id <- as.numeric(as.character(Time_df$id))
Age_df$id <- as.numeric(as.character(Age_df$id))

# Tree_df_nt <- Tree_df[-Tree_df$model=="True",]
Time_df_ratio <- Time_df[-which(Time_df$model=="True"),]
Age_df_ratio <- Age_df[-which(Age_df$model=="True"),]

Time_df_ratio$effect <- Time_df_ratio$effect / true_effects$Time
Age_df_ratio$effect <- Age_df_ratio$effect / true_effects$Age

# Summary plots ####

# Raw effects
# Tree_raw_plot <- ggplot(Tree_df, aes(x=id, y=effect)) + geom_bar(stat="identity")  + facet_grid(case~model) + theme_bw() + xlab("Tree") + ylab("Effect")

Time_raw_plot <- ggplot(Time_df, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Year") + ylab("Effect") + geom_hline(y=1)

Age_raw_plot <- ggplot(Age_df, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Age") + ylab("Effect")

# Relative to true
# Tree_scatter_plot <- ggplot(Tree_df_nt, aes(x=log(effect), y=log(true_effects$Tree))) + geom_point()  + facet_grid(case~model) + theme_bw() + xlab("Log-Estimated Tree Effect") + ylab("Log-True Tree Effect") + geom_smooth() + geom_abline(intercept=0, slope=1)

Time_ratio_plot <- ggplot(Time_df_ratio, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Year") + ylab("Ratio between estimated and true effect") + geom_hline(y=1)

Age_ratio_plot <- ggplot(Age_df_ratio, aes(x=id, y=effect)) + geom_line()  + facet_grid(case~model) + theme_bw() + xlab("Age") + ylab("Ratio between estimated and true effect") + geom_hline(y=1)

# Printing
# print(Tree_raw_plot)
print(Time_raw_plot)
print(Age_raw_plot)

# print(Tree_scatter_plot)
print(Time_ratio_plot)
print(Age_ratio_plot)
