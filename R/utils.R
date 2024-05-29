#' Calculate Volumes of Neighborhoods of Given Sizes
#' @description Calculate volumes in pixels for different neighborhood radii
#' @param L     List: named vectors of radii of neighborhoods in pixels in X, Y, Z
#' @param dims  Vector: dimensions of overall image in pixels in X, Y, and Z
#' @return      Numeric: area of neighborhood in pixels
#' @export
calc_vols <- function(L, dims) {
  return(unlist(lapply(L, calc_vol, dims=dims)))                                # Pass each radius to calc_vol
}

#' Calculate Volume of a Neighborhood of a Given Size
#' @description Calculate the volume in pixels for a single neighborhood radius
#' @param L     Named Vector: radius of neighborhood in pixels in X, Y, and Z
#' @param dims  Vector: dimensions of overall image in pixels in X, Y, and Z
#' @return      Vector: volumes of neighborhoods in pixels
#' @export
calc_vol <- function(L, dims) {
  D <- pmin((2*L + 1), dims)                                                    # Largest possible diam in each dim
  bound_box <- expand.grid(c(-floor((D[1]-1)/2):ceiling((D[1]-1)/2)),           # All pix coors in bounding box of
                           c(-floor((D[2]-1)/2):ceiling((D[2]-1)/2)),           # biggest possible neighborhood
                           c(-floor((D[3]-1)/2):ceiling((D[3]-1)/2)))           # Then express coors as fractions
  bound_box <- t(t(bound_box) / L)                                              # of dim-specific radii
  return(sum(sqrt(rowSums((bound_box)^2)) <= 1))                                # Volume of sphere in # of whole pix
}

#' Partition a Number into Addends as Evenly as Possible
#' @description Partition a number of variables into addends as evenly as possible
#' @param total  Numeric: whole number to be partitioned into additive components
#' @param comps  Numeric: number of additive components to split the total into
#' @return      Vector: addends into which the total decomposes
partition_evenly <- function(total, comps) {
  out <- c(rep(floor(total/comps), comps-1), ceiling(total/comps))
  if (sum(out) < total) {
    out[1:(total-sum(out))] <- out[1:(total-sum(out))] + 1
  }
  out <- out[out > 0]
  return(out)
}

#' Find Maximum Total Number of Compositional Bins
#' @description Find the total number of possible joint bins for a joint distribution of compositional variables
#' @param dimension     Numeric: total number of compositional variables in the joint distribution
#' @param bins_per_var  Numeric: number of bins to split each compositional variable into
#' @param min_per_var   Vector: minimum value of each compositional variable
#' @param max_per_var   Vector: maximum value of each compositional variable
#' @return      Numeric: the maximum number of possible compositional bins
total_comp_bins <- function(dimension, bins_per_var, min_per_var, max_per_var) {
  if ((missing(min_per_var) && missing(max_per_var)) ||  # If mins and maxes per variable are absent or all 0 and 100
      (sum(min_per_var == 0) == length(min_per_var) && sum(max_per_var == 100) == length(max_per_var))) {
    if (dimension == 1) {                                # Dimension 1: total bins = bins per cell type
      return(bins_per_var)
    } else if (dimension == 2) {                         # Dimension 2: total bins = sum of 1 through bins per cell type
      return(sum(1:bins_per_var))
    } else {                                             # Dimension > 2: invoke a formula
      v <- 1:bins_per_var                                # Start w/ vector of 1 through bins per cell type
      nv <- v                                            # Make a "new" copy of this vector
      for (i in 1:(dimension - 2)) {                     # For each dimension > 2, multiply v by (itself + (dimension-2))
        nv <- nv*(v+i)
      }
      return(sum(nv)/factorial(dimension-1))             # Then sum up v and divide by the factorial of (dimension-1)
    }
  } else {                                               # Given mins & maxes per variable, this might affect possible bins
    combos <- combn(length(min_per_var), dimension)      # All combinations of variables possible
    out <- 1                                             # Maximum number of bins possible - initialize at 1
    bin_vals <- vector("list", dimension)                # List to hold bin values for each variable in the current ensemble
    for (i in 1:ncol(combos)) {                          # For each combination of variables
      for (j in 1:dimension) {                           # List possible bin values for each variable
        bin_vals[[j]] <- seq(from = min_per_var[combos[j,i]], to = max_per_var[combos[j,i]], length.out = bins_per_var)
      }
      all_combos <- expand.grid(bin_vals)                # All combos of bin values for the variables in this combination
      poss_bins <- sum(rowSums(all_combos) <= 100)       # Permissible bins = those that total <= 100
      if (poss_bins > out) {                             # Get highest count of permissible bins across all variable combos
        out <- poss_bins
      }
    }
    return(out)
  }
}

#' Round Values of a Variable
#' @description Round the values of a variable found in a column of a census
#' @param focal_col  Numeric: a specific column of a census to round
#' @param col_min    Numeric: minimum value to round to
#' @param col_max    Numeric: maximum value to round to
#' @param bin_num    Numeric: number of bins to round into
#' @param bin_val    Logical: whether to include bin values or just give them integer IDs
#' @return      Vector: the rounded values of the census column
round_column <- function(focal_col, col_min, col_max, bin_num, bin_val) {
  if (col_min != col_max) {
    round_bins <- seq(col_min, col_max, length.out = bin_num)
    if (bin_val) {
      out <- round_bins[FNN::knnx.index(round_bins, focal_col, k=1)]
    } else {
      out <- as.vector(FNN::knnx.index(round_bins, focal_col, k=1))
    }
  } else {
    out <- focal_col
  }
  return(out)
}

#' Create a Logical Array
#' @description From a list of True coordinates, create a logical array
#' @param xy_coors    Matrix: xy coordinates of points in an array
#' @param matrix_dims Vector: dimensions of a matrix
#' @return      Array: an array of the specified dimensions where the specified coordinates are T and all else F
unwhich2d <- function (xy_coors, matrix_dims) {
  ar <- array(logical(prod(matrix_dims)), matrix_dims)
  xy_coors <- cbind(xy_coors, rep(1, nrow(xy_coors)))
  ar[xy_coors] <- TRUE
  return(ar)
}

#' Calculate the Histogram Entropy Estimate of a Distribution
#' @description Calculate entropy or joint entropy given for an ensemble of variables
#' @param cen   Array 2d: matrix of neighborhoods as rows and variables as columns, giving bin ID for each entry
#' @param bin   Integer: number of bins per variable
#' @return      Numeric: histogram entropy estimate for the ensemble
entropy <- function(cen, bin) {
  jd <- rep(0, nrow(cen))                                                       # joint distribution
  for (i in 1:ncol(cen)) {                                                      # turn multi-dimensional bins
    jd <- jd + (bin ^ (i-1))*(cen[,i] - 1)                                      # into one-dimensional bins
  }
  jd <- unname(table(jd))                                                       # omits 0 counts naturally
  jd <- jd / sum(jd)                                                            # histogram becomes prob dist
  out <- -sum(jd*log2(jd))                                                      # calculate & return entropy
  return(out)
}

#' Find Connected Patches in a 3D Mask
#' @description  Find connected patches and their sizes in a 3d mask. Based on buchsbaum's neuroim2::connComp3D
#' @param mask   Array 3d: logical array of pixels of a particular object or positive for a particular scalar
#' @return list of the 3d array of patch indices and the 3d array of patch sizes
#' @export
patch_3D <- function(mask) {
  if (length(dim(mask)) != 3 || !is.logical(mask)) {
    stop("Error: mask must be a 3d logical array.")
  }
  nodes <- numeric(length(mask)/9)
  labels <- array(0, dim(mask))
  DIM <- dim(mask)
  local.mask <- as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1), z=c(-1,0,1)))   # 27 coors of nbr vox plus self
  dimnames(local.mask) <- NULL
  local.mask <- local.mask[-ceiling(nrow(local.mask)/2),,drop=FALSE]            # Remove self vox
  neighbors <- function(vox) {                                                  # Function to find nbrs of focal vox
    vox.hood <- t(t(local.mask) + vox)                                          # Coors of vox's 26 neighbors
    if (any(vox == 1) || any(vox == DIM)) {                                     # If vox is at edge of mask
      vox.hood <- vox.hood[apply(vox.hood, 1, function(coords) {                # remove nbrs outside of mask
        all(coords >= 1 & coords <= DIM)
      }),,drop=FALSE]
    }
    vox.hood[labels[vox.hood] != 0,,drop=F]                                     # Remove nbrs labeled 0
  }
  find <- function(i) {                                                         # Function to ?
    while (nodes[i] != i) {
      i <- nodes[i]
    }
    nodes[i]
  }
  nextlabel <- 1
  grid <- arrayInd(which(mask), dim(mask))                                      # Coors of all vox of interest
  for (i in 1:NROW(grid)) {                                                     # For each vox of interest
    vox <- grid[i,]                                                             # Get 3d coors
    nbrs <- neighbors(vox)
    if (nrow(nbrs) == 0) {
      nodes[nextlabel] <- nextlabel
      labels[vox[1],vox[2],vox[3]] <- nextlabel
    } else {
      L <- labels[nbrs]
      ML <- min(L)
      labels[vox[1],vox[2], vox[3]] <- ML
      nodes[nextlabel] <- ML
      for (lab in L) {
        rootx <- find(lab)
        nodes[rootx] <- find(ML)
      }
    }
    nextlabel <- nextlabel + 1
  }
  for (k in 1:DIM[3]) {
    for (j in 1:DIM[2]) {
      for (i in 1:DIM[1]) {
        if (labels[i,j,k] > 0) {
          labels[i,j,k] <- find(labels[i,j,k])
        }
      }
    }
  }
  labs <- labels[labels!=0]
  forelabs <- labels > 0
  clusters <- sort(table(labs), decreasing=TRUE)
  SVol <- array(0, dim(mask))
  SVol[forelabs] <- clusters[as.character(labs)]
  indices <- 1:length(clusters)
  names(indices) <- names(clusters)
  IVol <- array(0, dim(mask))
  IVol[forelabs] <- indices[as.character(labs)]
  return(list(index=IVol, size=SVol))
}
