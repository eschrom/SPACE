#' Build a Joint Distribution
#'
#' After the variables in a census have been rounded into bins, assemble a joint distribution for an ensemble.
#'
#' @description Assemble a joint distribution for an ensemble from a rounded census.
#' @param census       Data frame: rounded percentages of each variable across neighborhoods
#' @param vars         Vector: variables for which to build the distribution
#' @param focal_vars   Vector: all variables in a broader analysis
#' @return      Data frame: all possible quantitative combinations of the variables and their counts
build_dist <- function(census, vars, focal_vars="all") {
  if (focal_vars == "all" || (focal_vars != "all" && sum(!is.na(match(vars, focal_vars))) > 0)) {
    cols <- which(!is.na(match(colnames(census), vars)))  # If >1 var is focal, counts of combo frequencies
    joint_dist <- plyr::count(census[,cols])
    colnames(joint_dist)[1:length(cols)] <- colnames(census)[cols]
  } else {                                                # If 0 vars are a focal var, don't do anything
    joint_dist <- NA
  }
  return(joint_dist)
}

#' Regularize a Joint Distribution
#'
#' Regularize a joint distribution to assign small probabilities to all bins, even those with measured value 0.
#'
#' @description Regularize a joint distribution to insure that all possible bins have positive probability
#' @param joint_dist   Data frame: all possible quantitative combinations of variables and their counts
#' @param bin_num      Numeric: number of rounding bins per variable; can be a vector of different numbers for each variable
#' @param min_max      Data frame: minimum and maximum value for each variable in the joint distribution
#' @param full_dist    Logical: whether to return the full distribution or just the vector of frequencies
#' @return      Data frame: regularized frequencies of all possible quantitative combinations of the variables
smooth_dist <- function(joint_dist, bin_num, min_max, full_dist) {
  num_obs <- sum(joint_dist$freq) + 1                     # Total # of observations, # of singletons, and # of doubletons,
  num_sing <- sum(joint_dist$freq[joint_dist$freq == 1])  # Number of singleton counts
  if (num_sing == 0) {                                    # If there are 0 singletons, assume 0.5 singletons
    num_sing <- 0.5                                       # i.e. a sample size twice as large would produce 1 singleton
  }
  num_doub <- sum(joint_dist$freq[joint_dist$freq == 2])  # Number of doubleton counts
  missing_prob <- num_sing / num_obs                      # Original Chao-Shen "missing probability" estimate due to under-sampling
  f1 <- (num_obs - 1) * num_sing                          # Two extra terms to refine "missing probability" to the Chao-Jost version
  f2 <- num_doub * 2
  if ((f1 != 0) || (f2 != 0)) {                           # As long as the two extra terms are not zero
    missing_prob <- missing_prob * (f1 / (f1+f2))         # Chao-Jost estimate of "missing probability"
  }
  num_vars <- ncol(joint_dist) - 1                        # Get number & names of variables, and set up smoothed joint distribution
  var_names <- colnames(joint_dist)[1:num_vars]
  smooth_joint_dist <- expand.grid(mapply(seq, min_max[1,], min_max[2,], length.out=bin_num, SIMPLIFY=F))
  if (num_vars == 1) {                                    # If the distribution only has 1 var, must coerce it to a data frame
    smooth_joint_dist <- as.data.frame(smooth_joint_dist)
  }
  colnames(smooth_joint_dist) <- var_names                # For smoothed distribution, start with original frequencies of all combos
  smooth_joint_dist <- plyr::join(smooth_joint_dist, joint_dist, by=colnames(smooth_joint_dist))
  smooth_joint_dist$freq[is.na(smooth_joint_dist$freq)] <- 0 #plyr::join enforces the row order from expand.grid
  if (num_vars > 1) {                                     # If the distribution includes more than 1 variable
    i <- 1                                                # Look for groups of 2+ vars from the same object map
    while (sum(stringr::str_count(var_names, paste("O", i, sep=""))) > 0) {
      cols <- which(stringr::str_count(var_names, paste("O", i, sep="")) == 1)
      if (length(cols) > 1) {                             # If such groups are found, remove percentage combos that sum > 100
        smooth_joint_dist$freq[rowSums(smooth_joint_dist[,cols]) > 100] <- NA
      }
      i <- i + 1
    }
  }
  if (sum(smooth_joint_dist$freq == 0, na.rm=T) == 0) {   # If there are no zero counts, no smoothing necessary
    smooth_joint_dist$freq <- smooth_joint_dist$freq / sum(smooth_joint_dist$freq, na.rm=T)
    if (full_dist) {
      return(smooth_joint_dist)                           # Return either the full data frame or just the freq column
    } else {
      return(smooth_joint_dist$freq)
    }
  }                                                       # If there are zero counts, perform smoothing
  smooth_joint_dist$freq <- (smooth_joint_dist$freq / sum(smooth_joint_dist$freq, na.rm=T)) * (1 - missing_prob)
  smooth_joint_dist$freq[smooth_joint_dist$freq==0] <- missing_prob / sum(smooth_joint_dist$freq == 0, na.rm=T)
  if (full_dist) {                                        # Normalize non-0s to sum to coverage, & spread missing prob across all 0s
    return(smooth_joint_dist)                             # Return either the full data frame or just the freq column
  } else {
    return(smooth_joint_dist$freq)
  }
}
