#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Utility functions

library(deldir)

# function to compute log-sum-exp
logSumExp <- function(x) {
  x_max = max(x)
  return(x_max + log( sum(exp(x - x_max)) ))
}

# function to obtain indices of n smallest values in a vector
WhichNMin <- function(x, n = 1) {
  n = min(n, length(x))
  threshold = sort(x, partial = n)[n]
  which_n_min = which(x <= threshold)
  # deal with ties
  if(length(which_n_min) > n) {
    x = x[which_n_min]
    which_n_min = which_n_min[order(x)[1:n]]
  }
  return(which_n_min)
}

# function to generate grid points
genGrid <- function(nx, ny, xmin = 0, ymin = 0, xmax = 1, ymax = 1) {
  x_grid = seq(xmin, xmax, length.out = nx + 2)
  y_grid = seq(ymin, ymax, length.out = ny + 2)
  x_grid = x_grid[-c(1, nx + 2)]
  y_grid = y_grid[-c(1, ny + 2)]
  coords_grid = expand.grid(x_grid, y_grid)
  colnames(coords_grid) = c("lon", "lat")
  return(coords_grid)
}

# function to obtain cell index of a grid (column major)
getGridIdx <- function(coords, grid_x, grid_y) {
  n_x = length(grid_x) - 1
  i = sapply(coords[, 1], FUN = function(v) Rfast::binary_search(grid_x, v, T))
  j = sapply(coords[, 2], FUN = function(v) Rfast::binary_search(grid_y, v, T))
  return((i - 2) * n_x + (j - 1))
}

# function to create Delaunay triangulation graph
dentrigraph <- function(coords, threshold = 1000){
  triangulation <- deldir(coords[,1], coords[,2])
  distances <- abs(triangulation$delsgs$x1 - triangulation$delsgs$x2) +
    abs(triangulation$delsgs$y1 - triangulation$delsgs$y2)
  #remove edges that are greater than .001
  edge.list <- as.matrix(triangulation$delsgs[distances < threshold,5:6])
  
  sorted <- sort(c(edge.list), index.return = TRUE)
  run.length <- rle(sorted$x)
  indices <- rep(1:length(run.length$lengths),times=run.length$lengths)
  edge.list.reformatted <- edge.list
  edge.list.reformatted[sorted$ix] <- indices
  edge.list.reformatted <- ord.mat(edge.list.reformatted);
  #create graph from list of edges
  distances <- sqrt((coords[edge.list.reformatted[,1],1] - coords[edge.list.reformatted[,2],1])^2 +
                      (coords[edge.list.reformatted[,1],2] - coords[edge.list.reformatted[,2],2])^2)
  graph0 <- graph_from_edgelist(edge.list.reformatted, directed = F)
  E(graph0)$weight <- distances
  return(graph0)
}

ord.mat <- function(M, decr = F, cols = NULL){
  if(is.null(cols))
    cols = 1:ncol(M)
  out = do.call("order", as.data.frame(M[,cols]))
  if (decr)
    out = rev(out)
  return(M[out, ])
} 
