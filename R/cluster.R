# Clustering based on residuals ####

generate_clusters <- function(tra, num_groups=2, group_by="Age", link="log", distance="euclidean", clust="kmeans", dep_var="Growth", ...){
  
  # Compute distance between trees
  dist_matrix <- dist_tra_euclidean(tra, group_by, link, dep_var)
  
  # Impute missing values? Use principal components to understand variance structure?
  # Or maybe impute all values in the distance matrix
  # Actually, that's impossible because requires estimates of time or age outside of range
  
  
  # Find clusters
  
  # kmeans (classical, no support for NA)
  
  # pam (robust k-means)
  
  # fanny (fuzzy)
  
  # agnes (hierarchical, combines)
  agnes(dist_matrix, num_groups, diss=TRUE)
  # diana (hierarchical, splits)
  
  # use cutree to split into groups
  
  # clusplot to view clusters
  
  # Relevant to paper!
  # Specifically, the Mclust( ) function in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models.
  
  # Return cluster object and group_id column
}

# Distance functions ####
# In general, more useful when used on residuals than raw data


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

