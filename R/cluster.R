# Relevant to paper!
# Specifically, the Mclust( ) function in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models.

# Automated clustering ####

auto_cluster_tra <- function(tra, resids, fit, model=c("Age", "Time"), split="Age", link="log", dep_var="Growth", optim="alternate", cluster_criteria="BIC", distance="euclidean", clust="kmeans", ...){
  
  # Start number of clusters at sqrt(n/2)
  # Common rule of thumb
  n_trees <- length(unique(tra$Tree))
  k <- round(sqrt(n_trees/2))
  
  # Loop controls
  solution_found <- FALSE
  iteration <- 0
  
  # Keeping track of progress  
  cluster_scores <- data.frame(k=1, criteria=fit[[cluster_criteria]])
  
  
  while (!solution_found)
  {
    # Updates
    iteration <- iteration + 1
    print(paste("Cluster iteration", iteration))
    print(paste("Testing model given", k, "clusters"))
    
    # Find clusters 
    clusters <- cluster_tra(resids, k, split, link, dep_var, distance, clust)#, ...)
    
    # Assign clusters to data
    new_tra <- merge(tra, clusters, by="Tree")
    cname <- paste(split, "Split", sep="_")
    new_tra[[cname]] <- new_tra$clusters
    
    # Standardize using clusters
    results <- standardize_tra(new_tra, model, split, link, dep_var, optim, auto_cluster=F, return_data=T, make_plots=F)#, ...)
    
    # Record results
    score <- results$fit[[cluster_criteria]]
    cluster_scores <- rbind(cluster_scores, data.frame(k=k, criteria=score))
    
    # Decide on next direction to search
    
    # Try (in vain) to ensure that each age group has more than 1 tree
    # If age trends are unsmoothed
    # Avoid singularities
    lower_bound <- 1
    upper_bound <- ifelse(optim=="gam", n_trees, floor(n_trees/2))
    
    next_smallest <- max(cluster_scores$k[cluster_scores$k < k], lower_bound)    
    next_largest  <- min(cluster_scores$k[cluster_scores$k > k], upper_bound)
    
    next_smallest_score <- cluster_scores[cluster_scores$k==next_smallest, "criteria"]
    gradient <- score - next_smallest_score
    
    # Adj.Rsq should be maximized
    # All others should be minimized
    if (cluster_criteria == "Adj.Rsq"){
      gradient <- -gradient
    }
    
    # Assumes strictly concave up reward function
    # This is actually a very good assumption
    # LLH should always increase with clusters
    # But number of parameters will increase as well
    # Thus holds for AIC / BIC
    
    direction <- ifelse (gradient > 0, "less", "more")
    
    # Binary search for naive speed / coverage
    if (direction == "less"){
      next_k <- ceiling((k + next_smallest)/2)
    } else if (direction == "more") {
      next_k <- floor((k + next_largest)/2)
    }
    
    # Once step size < 1 the local space is exhaustively searched
    if (next_k == k){
      solution_found <- TRUE
    }
    
    # Change k for next iteration
    k <- next_k
    
  }
  
  # Find the final solution given the optimal number of clusters
  if (cluster_criteria == "Adj.Rsq"){
    optimal_k <- cluster_scores[which.max(cluster_scores$criteria), "k"]
  } else {
    optimal_k <- cluster_scores[which.min(cluster_scores$criteria), "k"]
  }
  clusters <- cluster_tra(resids, k, split, link, dep_var, distance, clust)#, ...)
  
  optimal_tra <- merge(tra, clusters, by="Tree")
  cname <- paste(split, "Split", sep="_")
  optimal_tra[[cname]] <- new_tra$clusters
  optimal_tra <- optimal_tra[,-which(names(optimal_tra)=="clusters")]
  
  print (paste("The optimal number of clusters was found to be", optimal_k))
  print ("Computing final solution")
  
  solution <- standardize_tra(optimal_tra, model, split, link, dep_var, optim, auto_cluster=F, return_data=T, make_plots=F)#,...)  
  
  
  # Report full list of clusters / scores searched
  names(cluster_scores)[2] <- cluster_criteria
  print(cluster_scores)
  
  return(solution)
  
}

#  Wrappers ####

# Find clusters given k
cluster_tra <- function(tra,  num_groups=2, group_by="Age", link="log", dep_var="Growth",  distance="euclidean", clust="kmeans", ...){
  
  # Compute distance matrix
  dist_matrix <- find_dist_tra(tra, group_by, link, distance, dep_var)
  
  # Identify clusters
  clusters <- find_clusters(dist_matrix, num_groups, clust, ...)
  
  return(clusters)
  
}

# Clustering algorithms ####

# Wrapper for clustering
find_clusters <- function(dist_matrix, num_groups=2, clust="kmeans", ...){
  
  # Find clusters
  
  # kmeans (classical)
  if (clust=="kmeans"){
    clust_model <- kmeans(dist_matrix, num_groups, ...)
    clusters <- clust_model$cluster
  }
  
  # pam (robust k-means)
  else if(clust=="pam"){
    clust_model <- pam(dist_matrix, num_groups, diss=TRUE, ...)
    clusters <- clust_model$clustering
  }
  
  # fanny (fuzzy)
  else if(clust=="pam"){
    clust_model <- fanny(dist_matrix, num_groups, diss=TRUE, ...)
    clusters <- clust_model$clustering
  }
  
  # agnes (hierarchical, combines)
  else if(clust=="agnes"){
    clust_model <- agnes(dist_matrix, diss=TRUE, ...)    
    clusters <- cutree(clust_model, num_groups)
  }
  # diana (hierarchical, splits)
  else if(clust=="diana"){
    clust_model <- diana(dist_matrix, num_groups, diss=TRUE, ...)    
    clusters <- cutree(clust_model, num_groups)
  }
  
  # Cleanup
  names(clusters) <- colnames(dist_matrix)
  
  clust_df <- data.frame(Tree=names(clusters), clusters)
  
  clust_df$clusters <- as.factor(clust_df$clusters)
    
  return(clust_df)

}

# Distance functions ####
# In general, more useful when used on residuals than raw data

# Distance wrapper
find_dist_tra <- function(tra, group_by="Age", link="log", distance="euclidean", dep_var="Growth") {
  
  # Compute distance between trees
  if (distance=="euclidean"){
    dist_matrix <- dist_tra_euclidean(tra, group_by, link, dep_var)
  }
  
  # Impute missing values
  if (any(is.na(dist_matrix))){
    dist_matrix <- impute_indirect_path(dist_matrix)
    
    print("Some pairs of trees had no years in common. Missing similarity distances were imputed.")
  }
  
  return(dist_matrix)
}


# Mean Euclidean distance between pairs
dist_tra_euclidean <- function(tra, group_by="Age", link="log", dep_var="Growth"){

  trees <- unique(tra$Tree)
  
  # Transform data according to link
  d_tra <- tra
  if (link=="log"){
    d_tra[[dep_var]] <- log(tra[[dep_var]])
  }
  
  # Pairwise distance function
  pairwise_euclidean <- function(tree_a, tree_b){
    ids_a <- unique(tra[tra$Tree==tree_a, group_by])
    ids_b <- unique(tra[tra$Tree==tree_b, group_by])
    
    common_ids <- intersect(ids_a, ids_b)
    
    common_a <- tra[tra$Tree==tree_a & tra[[group_by]] %in% common_ids, ]
    common_b <- tra[tra$Tree==tree_b & tra[[group_by]] %in% common_ids, ]
    
    # Sort to ensure they line up
    common_a <- common_a[order(common_a[[group_by]]),]
    common_b <- common_b[order(common_b[[group_by]]),]
    
    
    # Compute total Euclidean distance
    tot_dist <- sqrt(sum((common_a[[dep_var]] - common_b[[dep_var]])^2))
    
    # Use mean distance to control for varying overlap
    mean_dist <- tot_dist/length(common_ids)
    
    return(mean_dist)
  }
  
  # Construct the base distance matrix
  d_matrix <- matrix(NA, length(trees), length(trees), dimnames=list(trees, trees))
  
  for (a in trees){
    for (b in trees){
      
      pos_a <- which(trees==a)
      pos_b <- which(trees==b)
      
      # Only fill upper triangle by hand
      if (pos_a < pos_b){
        d_matrix[a,b] <- pairwise_euclidean(a,b)
      } else if (pos_a == pos_b){
        # Diagonals are 0
        d_matrix[a,b] <- 0
      } else{
        # Lower triangle is the same as upper triangle
        d_matrix[a,b] <- d_matrix[b,a]
      }
    }
  }
  
  # NaN values are given when no overlap exists
  # More reasonable to code as NA
  d_matrix[!is.finite(d_matrix)] <- NA
  
  return(d_matrix)
  
} 


# Likelihood distance
# Likelihood of observing A if B is the true value
# Use sigma from whole dataset
# Smaller values mean less similar so... invert? 


# Imputing missing values in a distance matrix ####
# Uses high-dimensional trilateration
# Start with: http://stackoverflow.com/questions/10963054/finding-the-coordinates-of-points-from-distance-matrix
# See also: multidimensional scaling
# See also: distance geometry
# Projection onto subspace


# Find maximum number of dimensions
# Generate position matrix for that
# Place incomplete objects on the basis of minimizing strain (maybe happens automatically?)
# Find missing distances

impute_trilateration <- function(dist_matrix){
  
  o_names <- rownames(dist_matrix)
  new_dist_matrix <- dist_matrix
  
  # Complete objects have no NA values
  complete <- o_names[all(apply(dist_matrix, !is.na))]
  incomplete <- o_names[!(o_names %in% complete)]
  
  # Build position matrix for all complete objects
  # Columns are dimension, rows are points
  # Rownames are object names
  pos_matrix <- NA
  new_pos_matrix <- pos_matrix
  
  # For each incomplete object
  for (i in incomplete){
    
    # Reduce complete distance matrix to projection
    # Such that valid pairwise distances form the axes
    r_pos_matrix <- NA
    
    # Find coordinates of incomplete object in this embedding
    r_pos_i <- NA
    r_pos_matrix_i <- rbind(r_pos_matrix, r_pos_i)
    
    # Find originally missing distances
    missing_partners <- intersect(o_names(is.na(dist_matrix[i, ])), complete)
    
    for (m in missing_partners){
      new_dist_matrix[i,m] <- dist(r_pos_i, r_pos_matrix[m,])
      # Symmetry!
      new_dist_matrix[i,m] <- new_dist_matrix[m,i]
    }
    
    # Transform back into original coordinate system
    pos_i <- NA
    new_pos_matrix <- rbind(new_pos_matrix, pos_i)
    
  }
  
  # Find missing pairwise distances for incomplete objects
  # Use imputed coordinates
  i_i_pairs <- which(is.na(new_dist_matrix), arr.ind=T)
  
  for (j in 1:nrow(i_i_pairs)){
#     row_j <- o_names[i_i_pairs[j, "row"]]
#     col_j <- o_namesi_i_pairs[j, "col"]]

    new_dist_matrix[row_j, col_j] <- dist(new_pos_matrix[row_j,], new_pos_matrix[col_j,])
  }

  return(new_dist_matrix)

}

# Impute distances based on length of indirect path
impute_indirect_path <- function(dist_matrix){  
  
  # Construct base graph representing links between trees
  N <- nrow(dist_matrix)
  dist_graph <- graph.full(N)
  
  # Add weights
  E(dist_graph)$weight <- dist_matrix[lower.tri(dist_matrix)]
  
  # Delete missing edges
  missing_pairs <- which(is.na(dist_matrix), arr.ind=T)
  
  for (r in 1:nrow(missing_pairs)){
    i <- missing_pairs[r,1]
    j <- missing_pairs[r,2]
    dist_graph[i,j] <- FALSE
  }
  
  
  # Find the shortest path between each missing node
  for (r in 1:nrow(missing_pairs)){
    i <- missing_pairs[r,1]
    j <- missing_pairs[r,2]
    
    # Search for shortest weighted path
    path <- get.shortest.paths(dist_graph, i, j)
    
    # Grab path lengths
    # Using only first of shortest paths (ties should be very rare)
    path_lengths <- sapply(1:(length(path$vpath[[1]])-1),function(x){
      dist_graph[path$vpath[[1]][x], path$vpath[[1]][x+1]]
    })
    
    # Lower bound, assumes perfect doubling back
#     path_dist <- 0
#     for (k in path_lengths){
#       if (path_dist<=0){
#         path_dist <- path_dist + k
#       } else {
#         path_dist <- path_dist - k
#       }
#     }
#     path_dist <- abs(path_dist)
    
    
    # Reasonable estimate, assumes right angles every time
    path_dist <- sum(sapply(1:(length(path_lengths)-1), function(x){
      sqrt(path_lengths[x]^2 + path_lengths[x+1]^2)
    }))
    
    # Upper bound, assumes straight line between two missing values 
    # Equivalent to shortest.paths(dist_graph, i, j)
#     path_dist <- sum(path_lengths)
    
    # Use the path distance to impute missing distance values
    dist_matrix[i,j] <- path_dist
    
  }
  
  return(dist_matrix)
}

# Impute distances based on bounds set by indirect paths
impute_bounded_path <- function(dist_matrix, max_paths=100){  
  
  
  # Return upper and lower bounds of all missing values
  dist_matrix_lower <- dist_matrix
  dist_matrix_upper <- dist_matrix
  
  # Return original distance matrix if none are missing
  if (all(!is.na(dist_matrix))){
    return(list(lower=dist_matrix_lower, upper=dist_matrix_upper))
  }

  # Construct base graph representing links between trees
  N <- nrow(dist_matrix)
  dist_graph <- graph.full(N)
  
  # Add weights
  E(dist_graph)$weight <- dist_matrix[lower.tri(dist_matrix)]
  
  # Delete missing edges
  missing_pairs <- which(is.na(dist_matrix), arr.ind=T)
  
  for (r in 1:nrow(missing_pairs)){
    i <- missing_pairs[r,1]
    j <- missing_pairs[r,2]
    dist_graph[i,j] <- FALSE
  }
  
  # Function to find all unique, non-cycle paths between to nodes
  # Algorithm from:
  # http://mathoverflow.net/questions/18603/finding-all-paths-on-undirected-graph
  find_paths <- function(i,j,max_paths){
    
    # Check if we're in a dead end
    stuck <- function(x, seen){
      
      # Not stuck if we've reached the target
      if (x==j){
        return(list(stuck=FALSE, seen=seen))
      }
      
      neighbours <- (1:N)[dist_graph[x,]!=0]
      for (n in neighbours){
        # Don't double back
        if (!(n %in% seen)){
          seen <- c(seen, n)
          s <- stuck(n, seen)
          seen <- union(seen, s$seen)
          if (!s$stuck){
            return(list(stuck=FALSE, seen=seen))
          }
        }
      }
      return(list(stuck=FALSE, seen=seen)) 
    }
    
    # Recursive search function
    search <- function(x, path, seen){      
      
      found <- 0
      paths <- list()
      
      # Break out of path if we can't reach the target
      if (stuck(x, seen)$stuck){
        return(list(paths=list(path), seen=seen, found=0))
      } 
      
      # Our path is valid if it has reached its destination
      if (x==j){
        return(list(paths=list(path), seen=seen, found=1))
      }
      
      
      # Randomize which neighbour to search
      # Avoid weird results when subsampling
      neighbours <- (1:N)[dist_graph[x,]!=0]
      neighbours <- sample(neighbours, length(neighbours))
      
      for (n in neighbours){
        
        # Stop looking if you've found enough
        if (found >= max_paths){
          break
        }
        
        # Don't double back
        if (!(n %in% path)){
          # Add next step to path
          path <- c(path,n)
          
          # Add all valid paths down this branch
          results <- search(n, path, seen)
          paths <- c(paths,  results$paths)            
          found <- found + results$found
          
          # Remove it from the path to continue searching from that node
          path <- path[-which(path==n)]
        }
      }
      
      return(list(paths=paths, seen=seen, found=found))
    }
    
    # Start your search at i
    all_paths <- search(i,i,i)$paths
    
    return(all_paths)
  }
  
  # Use only unique paths
  # Order doesn't matter
  unique_paths <- function(paths){
    
    new_paths <- list(paths[[1]])
    for (i in 2:length(paths)){
      
      matches <- sum(sapply(new_paths, function(x){
        setequal(x, paths[[i]]) 
      }))
      
      # Add unique paths to new dataset
      if (matches == 0){
        new_paths <- c(new_paths, list(paths[[i]]))
      } 
    }
    return(new_paths)
  }
  
  # Lower bound on distance, assumes perfect doubling back
  lower_bound <- function(path){
    path_lengths <- sapply(1:(length(path)-1),function(x){
      dist_graph[path[x], path[x+1]]
    })
    
    path_dist <- 0
    for (k in path_lengths){
      if (path_dist<=0){
        path_dist <- path_dist + k
      } else {
        path_dist <- path_dist - k
      }
    }
    path_dist <- abs(path_dist)
    
    return(path_dist)
  }
  
  # Finding lower / upper bounds for each missing value
  for (r in 1:nrow(missing_pairs)){
    
    i <- missing_pairs[r,1]
    j <- missing_pairs[r,2]
    
    # Don't repeat calculations if mirror has already been done
    if (is.na(dist_matrix_upper[i,j])){
      # Search for shortest weighted path
      # Smallest upper bound, assumes straight line between two missing values 
      upper <- shortest.paths(dist_graph, i, j)
      
      # Grab all the paths between i and j
      paths <- unique_paths(find_paths(i,j,max_paths))
      
      lowers <- sapply(paths, lower_bound)
      lower <- max(lowers)
      
      # Use the path distance to bound missing distance values
      dist_matrix_upper[i,j] <- upper
      dist_matrix_upper[j,i] <- upper
      
      dist_matrix_lower[i,j] <- lower
      dist_matrix_lower[j,i] <- lower
    }    
  }
  
  return(list(lower=dist_matrix_lower, upper=dist_matrix_upper))
}

test_imputation <- function(dist_matrix){
  df <- data.frame(i=NA, j=NA, Pyth=NA, True=NA)[0,]
  
  for (i in 1:ncol(dist_matrix)){
    for(j in 1:ncol(dist_matrix)){
      if (i!=j){
        foo <- dist_matrix
        foo[i,j] <- NA
        foo[j,i] <- NA
        
        bat <- impute_indirect_path(foo)
        df <- rbind(df, data.frame(i=i, j=j, 
              Imputed=bat[i,j], 
              True=dist_matrix[i,j]))
      }
    }
  }
  return(df)
}
