#' Calculate Cis Mutual Information
#'
#' Quantify the non-randomness of patterning for each ensemble that is not explained by patterning among any
#' subset of those variables, for a single biological specimen. This is statistically compared to null expectations
#' in the absence of patterning, via simulation of randomized censuses.
#'
#' @description  Calculate cis mutual information for many ensembles and radii from a census
#' @param census       Data frame: census from the biological specimen
#' @param patch_list   List: patch list from the biological specimen
#' @param depth        Numeric: maximum number of variables per ensemble
#' @param radii        Vector: radii to examine, defaulting to all radii present in the provided census
#' @param bootstraps   Numeric: number of randomized censuses to generate for estimating confidence intervals
#' @param all          Vector: variables which all must be included in every maximum-depth ensemble
#' @param alo          Vector: variables of which at least one must be included in every maximum-depth ensemble
#' @param not          Vector: Variables which should not be included in any maximum-depth ensemble
#' @param max_bins     Numeric: maximum number of bins per variable for rounding
#' @param cores        Numeric: number of cores to use for parallel calculations
#' @return       List: data frames of cisMI for all ensembles at each radius
#' @export
measure_cisMI <- function(census, patch_list, depth, radii=NULL, bootstraps=100,
                          all=NULL, alo=NULL, not=NULL, max_bins = 100, cores=NULL) {
  if ((length(all) + !is.null(alo)) > depth) {
    stop("Error: the number of required variables exceeds ensemble depth.")
  }
  if (is.null(radii)) {                                                         # If using all radii
    radii <- unique(census$Radius)                                              # get all numeric radii
  }
  census <- census[census$Radius %in% radii, ]                                  # Drop nbhds of irrelevant radii
  var_cols <- colnames(census)[which(stringr::str_count(colnames(census), "O|S") >= 1)]
  cols_to_drop <- rep(F, length(var_cols))
  if (!is.null(not)) {                                                          # If variables are explicitly excluded
    cols_to_drop[var_cols %in% not] <- T                                        # drop them
  }                                                                             # If variables are implicitly excluded
  if ((length(all) + !is.null(alo)) == depth) {                                 # drop them
    cols_to_drop[!(var_cols %in% all) & !(var_cols %in% alo)] <- T
  }                                                                             # Drop coor cols but not radius col
  census <- census[ , !c(cols_to_drop,T,T,T,F)]                                 # Infer object-scalar pairs of interest
  os_pairs <- colnames(census)[which(stringr::str_count(colnames(census), "O|S") == 2)]
  if (length(os_pairs) > 0) {
    obj_maps <- unique(gsub("\\..*","", os_pairs))                              # Object maps with attached scalars
    OS_pairs <- vector("list", length=length(obj_maps))                         # List for pairing matrices
    names(OS_pairs) <- obj_maps                                                 # Name matrices by their object map
    for (i in 1:length(OS_pairs)) {                                             # For each object map represented
      sub_os_pairs <- os_pairs[which(stringr::str_count(os_pairs, obj_maps[i]) == 1)]
      sub_scalars <- sub(".*_", "", sub_os_pairs)                               # Get scalar IDs and object IDs
      sub_objs <- sub("_.*", "", sub_os_pairs)
      OS_pairs[[i]] <- matrix(0, nrow=length(unique(sub_objs)), ncol=length(unique(sub_scalars)))
      rownames(OS_pairs[[i]]) <- unique(sub_objs)                               # Create and fill in a pairing matrix
      colnames(OS_pairs[[i]]) <- unique(sub_scalars)
      for (j in 1:length(sub_os_pairs)) {
        OS_pairs[[i]][rownames(OS_pairs[[i]]) == sub_objs[j], colnames(OS_pairs[[i]]) == sub_scalars[j]] <- 1
      }
    }
  } else {
    OS_pairs <- NULL
  }
  ss <- min(table(census$Radius))
  b <- min(floor(((ss+1)/10) ^ (1/depth)), max_bins)
  print(paste("Based on a minimum sample size of ", ss, ", there will be ", b, " bins per variable.", sep=""))
  all_col_mins <- matrix(NA, nrow=length(radii), ncol=sum(!cols_to_drop))       # Prep to record min & max of each variable
  all_col_maxs <- matrix(NA, nrow=length(radii), ncol=sum(!cols_to_drop))       # at each neighborhood radius
  out <- vector(mode="list", length=length(radii))                              # Output meta-list of ensembles by radius, depth
  ent <- vector(mode="list", length=length(radii))                              # Output meta-list of entropy by radius, depth
  mi <- vector(mode="list", length=length(radii))                               # Output meta-list of mutual info by radius, depth
  for (i in 1:length(radii)) {                                                  # Address 1 radius at a time
    out[[i]] <- vector(mode="list", length=depth)
    ent[[i]] <- vector(mode="list", length=depth)
    mi[[i]] <- vector(mode="list", length=depth)
    print(paste("Beginning radius ", radii[i], " microns.", sep=""))
    sub_cen <- census[census$Radius == radii[i], -ncol(census)]                 # Narrow census to this radius & drop radius column
    sub_pls <- patch_list[[which(as.numeric(names(patch_list)) == radii[i])]]   # Narrow patch list to correct radius
    print(paste("Generating ", bootstraps, " bootstrap random censuses.", sep=""))
    gc()
    num_cores <- min(parallel::detectCores() - 1, cores, floor(as.numeric(memuse::Sys.meminfo()$freeram)/
                                                               as.numeric(object.size(sub_pls)))-1)
    print(paste("     Preparing ", num_cores, " nodes to simulate random censuses in parallel.", sep=""))
    make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = make_cluster)            # For each local census, excise its bounding box
    `%dopar%` <- foreach::`%dopar%`                              # Attempt to help R find what %dopar% means
    print(paste("     Simulating ", bootstraps, " random censuses.", sep=""))
    all_cen <- foreach::foreach (j = c(1:bootstraps), .packages=c('dplyr','tidyr'), .export=c('random_census', 'summarize_patches')) %dopar% {
      random_census(sub_pls, OS_pairs)                           # Simulate a random census
    }
    parallel::stopCluster(cl = make_cluster)
    all_cen <- lapply(all_cen, function(x) {x[ , !cols_to_drop, drop=F]})       # Drop columns for unused variables
    all_cen <- c(list(sub_cen), all_cen)                                        # Put true census first in list
    all_col_mins[i,] <- apply(dplyr::bind_rows(lapply(all_cen, function(x) {apply(x,2,min)})), 2, min) # Overall min
    all_col_maxs[i,] <- apply(dplyr::bind_rows(lapply(all_cen, function(x) {apply(x,2,max)})), 2, max) # & max per var
    print(paste("Rounding censuses into ", b, " bins per variable.", sep=""))   # Round each variable in each census into bins
    all_cen <- lapply(all_cen, function(x) {mapply(round_column, x[,,drop=F], all_col_mins[i,], all_col_maxs[i,], b, F, SIMPLIFY=T)})
    all_cen <- abind::abind(all_cen, along=3)                                   # Reshape list of data frames into 3d array
    var_names <- colnames(sub_cen)
    for (j in 1:depth) {                                                        # For each depth, collate ensembles & calculate MI
      print(paste("Collating all requested ", j, "-ensembles.", sep=""))
      out[[i]][[j]] <- as.data.frame(t(combn(colnames(sub_cen), j)))
      colnames(out[[i]][[j]]) <- paste("V", LETTERS[1:j], sep="")
      if (j == depth) {                                                         # If this is the highest depth
        if (any(all != "none")) {                                               # Make sure all ensembles have all "all" variables
          out[[i]][[j]] <- out[[i]][[j]][apply(out[[i]][[j]], 1, function(x,y){sum(x %in% y) == length(y)}, y=all), , drop=F]
        }
        if (any(alo != "none")) {                                               # Make sure all ensembles have >= 1 "alo" variables
          out[[i]][[j]] <- out[[i]][[j]][apply(out[[i]][[j]], 1, function(x,y){any(x %in% y)}, y=alo), , drop=F]
        }                                                                       # Below the highest depth, ensembles that don't
      }                                                                         # match demands will still be used in subtraction
      print(paste("Calculating entropy for all ", j, "-ensembles.", sep=""))    # Calculate entropy for every ensemble
      if (nrow(out[[i]][[j]]) < min(parallel::detectCores() - 1, cores)) {      # If more ensembles than cores, parallelize
        ent[[i]][[j]] <- array(NA, dim=c(nrow(out[[i]][[j]]), bootstraps+1))
        for (e in 1:nrow(out[[i]][[j]])) {
          v <- var_names %in% out[[i]][[j]][e, ]
          ent[[i]][[j]][e,] <- apply(all_cen[ , v, ,drop=F], 3, entropy, bin = b)
        }
      } else {
        num_cores <- min(parallel::detectCores() - 1, cores)
        print(paste("     and parallelizing over ", num_cores, " nodes.", sep=""))
        make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
        doParallel::registerDoParallel(cl = make_cluster)
        `%dopar%` <- foreach::`%dopar%`                                         # Tell R what %dopar% means
        ent[[i]][[j]] <- foreach::foreach(e = 1:nrow(out[[i]][[j]]), .packages=c('foreach','parallel','doParallel'),
                                          .export=c('entropy')) %dopar% {
                                            v <- var_names %in% out[[i]][[j]][e, ]
                                            apply(all_cen[ ,v, ,drop=F], 3, entropy, bin = b)
                                          }
        parallel::stopCluster(cl = make_cluster)
        ent[[i]][[j]] <- aperm(abind::abind(ent[[i]][[j]], along=2), c(2,1))
      }
      if (j == 1) {                                                             # At depth 1, all ensembles
        B <- rep(b, nrow(out[[i]][[j]]))                                        # have b bins
      } else {                                                                  # At higher depths, ensembles
        B <- rep(0, nrow(out[[i]][[j]]))                                        # may have different bin #s
        for (k in 1:nrow(out[[i]][[j]])) {                                      # For each ensemble separately
          vars <- as.character(out[[i]][[j]][k,1:j])                            # Get the variables involved
          scls <- vars[which(stringr::str_count(vars, "S") == 1)]               # Split into scalars (ind or linked)
          objs <- vars[which(stringr::str_count(vars, "S") == 0)]               # and objects
          act_B <- 1                                                            # Start with 1 total bin
          if (length(scls) > 0) {                                               # If any scalars, there are
            act_B <- act_B * b^length(scls)                                     # b bins for each
          }
          if (length(objs) > 0 ) {                                              # If any objects,
            obj_maps <- unique(gsub("\\..*", "",objs))                          # get unique object maps
            for (p in 1:length(obj_maps)) {                                     # For each obj map in this ensemble
              obj_vars <- objs[stringr::str_count(objs, obj_maps[p]) == 1]      # Only objs in this map
              wch_obj_vars <- which(colnames(sub_cen) %in% obj_vars)            # which census cols they are in
              act_B <- act_B * total_comp_bins(length(obj_vars), b, all_col_mins[i,wch_obj_vars],
                                               all_col_maxs[i, wch_obj_vars])   # Update actual bin number
            }
          }
        }
      }
      ent[[i]][[j]] <- ent[[i]][[j]] + array(B / (1.5 * nrow(sub_cen)), dim=dim(ent[[i]][[j]]))
      if (j > 1) {                                                              # If depth > 1, correct for compositionality where necessary
        mi[[i]][[j]] <- array(0, dim=dim(ent[[i]][[j]]))
      } else {
        mi[[i]][[j]] <- ent[[i]][[j]]
      }
    }
    if (depth > 1) {                                                            # MI only necessary if depth > 1
      print("Combining entropies to calculate mutual information.")
      for (j in 2:depth) {                                                      # For depths 2 and greater
        for (k in 1:nrow(mi[[i]][[j]])) {                                       # For each ensemble
          for (p in 1:j) {                                                      # For each depth, find sub-ensembles
            rows <- which(apply(out[[i]][[p]], 1, function(x,y) {sum(x %in% y)}, y=out[[i]][[j]][k,]) == p)
            mi[[i]][[j]][k,] <- mi[[i]][[j]][k,] + ((-1)^(p-1))*colSums(ent[[i]][[p]][rows,,drop=F]) # +/- entropies
          }
        }
      }
    }
    print("Testing mutual information for significance.")
    for (j in 1:depth) {                                                        # CisMI = true MI - mean(random MI)
      out[[i]][[j]]$CisMI <- mi[[i]][[j]][,1,drop=F] - rowMeans(mi[[i]][[j]][,2:(bootstraps+1),drop=F])
      out[[i]][[j]]$Zscore <- out[[i]][[j]]$CisMI / apply(mi[[i]][[j]][,2:(bootstraps+1),drop=F], 1, sd)
      out[[i]][[j]]$Zscore[is.nan(out[[i]][[j]]$Zscore)] <- 0                   # Missing variable -> sd=0 -> Z=nan
      out[[i]][[j]]$Pvalue <- pnorm(abs(out[[i]][[j]]$Zscore), lower.tail=F)
      out[[i]][[j]]$Pvalue <- out[[i]][[j]]$Pvalue * 2                          # Correct for two-sided testing
    }
    out[[i]] <- dplyr::bind_rows(out[[i]])
    vcols <- stringr::str_count(colnames(out[[i]]), "V") == 1
    out[[i]] <- cbind(out[[i]][,vcols,drop=F], out[[i]][,!vcols])
    out[[i]]$Padjust <- p.adjust(out[[i]]$Pvalue, "BH")
  }
  names(out) <- as.character(radii)                                             # Idx = nbhd radius
  return(out)
}

#' Generate a Random Census
#'
#' Generate a randomized census from a census and a patch list at a single radius.
#'
#' @description Generate a random census to simulate expectations in the absence of patterning
#' @param pls   List: for each object map or independent scalar, patches of contiguous pixels, at just one radius
#' @param osp   List: for each object map with attached scalars, a matrix giving the combinations of interest
#' @return      Data frame: randomized census simulated by scrambling contiguous patches
random_census <- function(pls, osp) {
  for (i in 1:length(pls)) {                                                    # For each object map & independent scalar
    num_nbhds <- max(pls[[i]]$Nbhd)                                             # Get number of neighborhoods
    elig_pix <- dplyr::summarise(pls[[i]], .by=Nbhd, Area = sum(Area))$Area     # Get # non-bkgd pix from each nbhd
    cum_elig_pix_all <- cumsum(elig_pix)                                        # Cumulative sums of non-bkgd pix
    cum_elig_pix <- unique(cum_elig_pix_all)
    pls[[i]] <- pls[[i]][pls[[i]]$Area > 0, ]                                   # Remove empty nbhds
    pls[[i]] <- pls[[i]][sample(nrow(pls[[i]])), ]                              # Randomly reorder patches
    scl_cols <- grepl("S", colnames(pls[[i]]))                                  # Cols w/ scalars, attached or independent
    if (any(scl_cols)) {                                                        # If there are scalars
      pls[[i]][,scl_cols] <- pls[[i]][,scl_cols] / pls[[i]]$Area                # Convert total expression to avg per pix
    }                                                                           # This makes splitting patches easy
    cum_ptch_size <- cumsum(pls[[i]]$Area)                                      # Cumulative sums of patch sizes
    cum_elig_pix_naln <- cum_elig_pix[!(cum_elig_pix %in% cum_ptch_size)]       # New patch sizes, w/ those straddling nbhds broken
    new_ptch_size <- c(min(cum_ptch_size[1], cum_elig_pix[1]), diff(sort(c(cum_ptch_size, cum_elig_pix_naln))))
    ptch_brks <- which(cumsum(new_ptch_size) %in% cum_elig_pix_naln)            # Which cum patch sizes are new do to breaking
    ptch_brks <- plyr::count(ptch_brks - c(0:(length(ptch_brks)-1)))            # Which of orig patches were broken, & how many times
    patch_IDs <- rep(1, length(pls[[i]]$Area))                                  # List all original patch IDs
    patch_IDs[ptch_brks$x] <- patch_IDs[ptch_brks$x] + ptch_brks$freq           # Broken patches have ID repeated for each break
    pls[[i]] <- pls[[i]][rep(1:nrow(pls[[i]]), times=patch_IDs), ]              # Repeat each patch which must be broken
    pls[[i]]$Area <- new_ptch_size                                              # New patch sizes same as old except for broken patches
    ends <- which(cumsum(new_ptch_size) %in% cum_elig_pix_all)                  # Index of end of each non-empty neighborhood
    begs <- c(1, (ends[1:(length(ends)-1)] + 1))                                # Index of beginning of each non-empty neighborhood
    new_nbhds <- mapply(function(s,e){s:e}, begs, ends)                         # Which new patches belong to each non-empty nbhd
    new_nbhds <- unname(unlist(lapply(new_nbhds, length)))                      # Number of new patches in each new neighborhood
    new_nbhds <- rep(1:length(new_nbhds), times=new_nbhds)                      # New nbhd assignment for each patch
    pls[[i]]$Nbhd <- new_nbhds                                                  # Add new nbhd assignment to data frame
    if (any(scl_cols)) {                                                        # If there are scalars
      pls[[i]][,scl_cols] <- pls[[i]][,scl_cols] * pls[[i]]$Area                # Convert avg per pix expr back to total
    }                                                                           # Now that patches have been split
    if (num_nbhds - max(new_nbhds) > 0) {                                       # If there were empty nbhds with no objects
      empty_nbhds <- as.data.frame(matrix(0,nrow=(num_nbhds - max(new_nbhds)), ncol=ncol(pls[[i]])))
      colnames(empty_nbhds) <- colnames(pls[[i]])                               # Make these empty nbhds with obj = 0
      empty_nbhds$Nbhd <- c((max(new_nbhds)+1):num_nbhds)                       # Assign them the remaining nbhd IDs
      pls[[i]] <- dplyr::bind_rows(list(pls[[i]], empty_nbhds))                 # Append them to the frame of all nbhds
    }
  }
  out <- summarize_patches(pls, osp)
  return(out)
}

#' Calculate Trans Mutual Information
#'
#' Quantify the difference in patterning for every ensemble, across groups of multiple biological specimens.
#' This is statistically compared to null expectations in which sample groups are defined at random.
#'
#' @description  Calculate trans mutual information for every ensemble given a grouping structure of specimens
#' @param censuses     List: censuses for multiple biological specimens
#' @param groups       Frame: data frame with grouping factors as columns and specimens as rows
#' @param depth        Numeric: maximum number of variables per ensemble
#' @param radii        Numeric: radii at which to examine spatial patterns
#' @param bootstraps   Numeric: number of randomized transMI networks to simulate for estimating P values
#' @param all          Vector: variables which all must be included in every maximum-depth ensemble
#' @param alo          Vector: variables of which at least one must be included in every maximum-depth ensemble
#' @param not          Vector: Variables which should not be included in any maximum-depth ensemble
#' @param max_bins     Numeric: maximum number of bins per variable for rounding
#' @param cores        Numeric: number of cores to use for parallel calculations
#' @return       List: data frames of transMI for all ensembles at each radius
#' @export
measure_transMI <- function(censuses, groups, depth, radii, bootstraps=100,
                            all=NULL, alo=NULL, not=NULL, max_bins=100, cores=NULL) {
  if ((length(all) + !is.null(alo)) > depth) {
    stop("Error: the number of required variables exceeds ensemble depth.")
  }
  if (is.character(radii)) {                                                    # If radii = "all"
    radii <- lapply(censuses, function(x) {unique(x$Radius)})                   # Use radii present in all censuses
    radii <- Reduce(intersect, radii)
  }
  if (length(colnames(censuses[[1]])) < length(Reduce(union, lapply(censuses, colnames)))) {
    stop("Error: the censuses do not all contain the same variables.")
  }
  var_cols <- colnames(censuses[[1]])[which(stringr::str_count(colnames(censuses[[1]]), "O|S") >= 1)]
  cols_to_drop <- rep(F, length(var_cols))
  if (!is.null(not)) {                                                          # If variables are explicitly excluded
    cols_to_drop[var_cols %in% not] <- T                                        # drop them
  }                                                                             # If variables are implicitly excluded
  if ((length(all) + !is.null(alo)) == depth) {                                 # drop them
    cols_to_drop <- cols_to_drop[!(var_cols %in% all) & !(var_cols %in% alo)] <- T
  }                                                                             # Drop coor cols but not radius col
  for (i in 1:length(censuses)) {
    censuses[[i]] <- censuses[[i]][ , !c(cols_to_drop,T,T,T,F)]
  }
  os_pairs <- colnames(censuses[[1]])[which(stringr::str_count(colnames(censuses[[1]]), "O|S") == 2)]
  if (length(os_pairs) > 0) {
    obj_maps <- unique(gsub("\\..*","", os_pairs))                              # Object maps with attached scalars
    OS_pairs <- vector("list", length=length(obj_maps))                         # List for pairing matrices
    names(OS_pairs) <- obj_maps                                                 # Name matrices by their obj map
    for (i in 1:length(OS_pairs)) {                                             # For each object map represented
      sub_os_pairs <- os_pairs[which(stringr::str_count(os_pairs, obj_maps[i]) == 1)]
      sub_sclrs <- sub(".*_", "", sub_os_pairs)                                 # Get scalar IDs and obj IDs
      sub_objs <- sub("_.*", "", sub_os_pairs)
      OS_pairs[[i]] <- matrix(0, nrow=length(unique(sub_objs)), ncol=length(unique(sub_sclrs)))
      rownames(OS_pairs[[i]]) <- unique(sub_objs)                               # Create and fill in a pairing matrix
      colnames(OS_pairs[[i]]) <- unique(sub_sclrs)
      for (j in 1:length(sub_os_pairs)) {
        OS_pairs[[i]][rownames(OS_pairs[[i]]) == sub_objs[j], colnames(OS_pairs[[i]]) == sub_sclrs[j]] <- 1
      }
    }
  }
  ss <- min(unlist(lapply(censuses, function(x) {min(table(x$Radius))})))       # Smallest sample size across all radii and object tables
  b <- floor(((ss+1)/10) ^ (1/depth))
  print(paste("Based on a minimum sample size of ", ss, ", there will be ", b, " bins per variable.", sep=""))
  out <- vector(mode="list", length=length(radii))                              # Nested list by radius, depth
  for (i in 1:length(radii)) {                                                  # Address 1 radius at a time
    print(paste("Beginning radius ", radii[i], " pixels.", sep=""))
    sub_min <- vector("list", length(censuses))
    sub_max <- vector("list", length(censuses))
    sub_cen <- vector("list", length(censuses))                                 # For all censuses, narrow down to
    for (j in 1:length(censuses)) {                                             # correct radius, and drop radius col
      sub_cen[[j]] <- censuses[[j]][censuses[[j]]$Radius == radii[i], -ncol(censuses[[j]])]
      sub_min[[j]] <- as.numeric(apply(sub_cen[[j]], 2, min))
      sub_max[[j]] <- as.numeric(apply(sub_cen[[j]], 2, max))
      sub_cen[[j]] <- mapply(round_column, sub_cen[[j]][,,drop=F], sub_min[[j]], sub_max[[j]], b, T)
    }
    census_combos <- t(combn(length(censuses), 2))                              # All pairs of images
    var_names <- colnames(sub_cen[[1]])
    for (j in 1:depth) {                                                        # For each depth, collate ensembles & calculate MI
      print(paste("Collating all requested ", j, "-ensembles.", sep=""))
      out[[i]][[j]] <- as.data.frame(t(combn(var_names, j)))
      colnames(out[[i]][[j]]) <- paste("V", LETTERS[1:j], sep="")
      if (j == depth) {                                                         # If this is the highest depth
        if (any(all != "none")) {                                               # Make sure all ensembles have all "all" variables
          out[[i]][[j]] <- out[[i]][[j]][apply(out[[i]][[j]], 1, function(x,y){sum(x %in% y) == length(y)}, y=all), , drop=F]
        }
        if (any(alo != "none")) {                                               # Make sure all ensembles have >= 1 "alo" variables
          out[[i]][[j]] <- out[[i]][[j]][apply(out[[i]][[j]], 1, function(x,y){any(x %in% y)}, y=alo), , drop=F]
        }                                                                       # Below the highest depth, ensembles that don't
      }
      out[[i]][[j]][,(j + 1):(j + nrow(census_combos)*2)] <-
        matrix(NA, nrow=nrow(out[[i]][[j]]), ncol=nrow(census_combos)*2)
      colnames(out[[i]][[j]])[(j + 1):(j + nrow(census_combos)*2)] <-
        paste(rep(paste("MI_Im", census_combos[,1], "_Im", census_combos[,2], sep=""), times=2),
              "_", rep(c("A","B"), each=nrow(census_combos)), sep="")
      joint_dists <- array(NA, dim=c(nrow(out[[i]][[j]]), b^j, length(censuses)))
      if (nrow(out[[i]][[j]]) > 1000) {                                         # If >1000 ensembles, find joint dists in parallel
        num_cores <- min(parallel::detectCores() - 1, cores)
        print(paste("     and parallelizing over ", num_cores, " nodes.", sep=""))
        make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
        doParallel::registerDoParallel(cl = make_cluster)
        `%:%` <- foreach::`%:%`                                                 # Attempt to help R find what %:% and %dopar% means
        `%dopar%` <- foreach::`%dopar%`
        joint_dists <- foreach::foreach(k = 1:length(censuses), .packages=c('plyr','stringr','foreach','parallel','doParallel'),
                                        .export=c('build_dist','smooth_dist')) %:% # For each census and each ensemble
          foreach::foreach(p = 1:nrow(out[[i]][[j]]), .combine=rbind) %dopar% { # Get col mins and maxes for this ensemble
            mins_and_maxs <- as.data.frame(rbind(sub_min[[k]], sub_max[[j]]))[,match(out[[i]][[j]][p, 1:j], var_names), drop=F]
            smooth_dist(build_dist(sub_cen[[k]], out[[i]][[j]][p,1:j], "all"), b, mins_and_maxs, full_dist=F)
          }
        parallel::stopCluster(cl = make_cluster)
        joint_dists <- abind::abind(joint_dists, along=3)                       # Turn list of joint dist matrices into 3d array
      } else {                                                                  # If <= 1000 ensembles, find joint dists serially
        for (p in 1:nrow(out[[i]][[j]])) {                                      # For each ensemble, get col mins and maxes
          for (k in 1:length(censuses)) {                                       # Assemble joint dist
            mins_and_maxs <- as.data.frame(rbind(sub_min[[k]], sub_max[[j]]))[,match(out[[i]][[j]][p, 1:j], var_names), drop=F]
            joint_dists[p,,k] <- smooth_dist(build_dist(sub_cen[[k]], out[[i]][[j]][p,1:j], "all"), b, mins_and_maxs, full_dist=F)
          }
        }
      }
      print(paste("Calculating generalized mutual information for all ", j, "-ensembles.", sep=""))
      if (nrow(census_combos) > 500) {                                          # If >500 image pairs, measure MI in parallel
        num_cores <- min(parallel::detectCores() - 1, cores)                    # For each combo of censuses
        print(paste("     and parallelizing over ", num_cores, " nodes.", sep=""))
        make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")        # MI for all ensembles is calculated at once
        doParallel::registerDoParallel(cl = make_cluster)
        `%dopar%` <- foreach::`%dopar%`                                         # Attempt to help R find what %dopar% means
        out[[i]][[j]][,(j + 1):(j + nrow(census_combos))] <- foreach::foreach(r = 1:nrow(census_combos), .combine=cbind) %dopar% {
          d1 <- joint_dists[,,census_combos[r,1]]
          d2 <- joint_dists[,,census_combos[r,2]]
          rowSums(d1 * log2(d1 / d2), na.rm=T)
        }                                                                       # KL is not symmetric - measure in both directions
        out[[i]][[j]][,(j + nrow(census_combos) + 1):(j + nrow(census_combos)*2)] <- foreach::foreach(r = 1:nrow(census_combos), .combine=cbind) %dopar% {
          d1 <- joint_dists[,,census_combos[r,1]]
          d2 <- joint_dists[,,census_combos[r,2]]
          rowSums(d2 * log2(d2 / d1), na.rm=T)
        }
        parallel::stopCluster(cl = make_cluster)
      } else {                                                                  # If <= 500 image pairs, measure MI serially
        for (r in 1:nrow(census_combos)) {                                      # For each combo of censuses, done in each direction
          if (dim(joint_dists)[1] > 1) {                                        # MI for all ensembles is calculated at once
            d1 <- joint_dists[,,census_combos[r,1]]
            d2 <- joint_dists[,,census_combos[r,2]]
            out[[i]][[j]][,j + r] <- rowSums(d1 * log2(d1 / d2), na.rm=T)
            out[[i]][[j]][,j + nrow(census_combos) + r] <- rowSums(d2 * log2(d2 / d1), na.rm=T)
          } else {                                                              # If there is only 1 ensemble, cannot use rowSums!
            d1 <- joint_dists[,,census_combos[r,1]]
            d2 <- joint_dists[,,census_combos[r,2]]
            out[[i]][[j]][,j + r] <- sum(d1 * log2(d1 / d2), na.rm=T)
            out[[i]][[j]][,j + nrow(census_combos) + r] <- sum(d2 * log2(d2 / d1), na.rm=T)
          }
        }
      }
    }
  }
  names(out) <- as.character(radii)                                             # Idx1 = nbhd radius, idx2 = ensemble depth
  for (i in 1:length(radii)) {                                                  # Each element is frame of bootstrap
    names(out[[i]]) <- as.character(c(1:depth))                                 # MI estimates (cols) for every ensemble (rows)
  }                                                                             # Recursively adjust for lower-order MI
  print("Averaging transMI for each pair of images.")
  for (i in 1:length(out)) {                                                    # For each radius
    for (j in 1:length(out[[i]])) {                                             # For each depth
      out[[i]][[j]][,(j+1):(nrow(census_combos)+j)] <- (out[[i]][[j]][,(j+1):(nrow(census_combos)+j)] +
                                                          out[[i]][[j]][,(nrow(census_combos)+j+1):(ncol(out[[i]][[j]]))]) / 2    # Avg MI from each direction
      out[[i]][[j]] <- out[[i]][[j]][,-((nrow(census_combos)+j+1):(ncol(out[[i]][[j]])))]
      colnames(out[[i]][[j]])[(j+1):ncol(out[[i]][[j]])] <-
        stringr::str_sub(colnames(out[[i]][[j]])[(j+1):ncol(out[[i]][[j]])], end=-3)
    }
  }
  print("Testing the modularity of transMI for significance.")
  mutual_info <- out
  remove(out)
  for (i in 1:length(mutual_info)) {
    mutual_info[[i]] <- dplyr::bind_rows(mutual_info[[i]])
    vcols <- stringr::str_count(colnames(mutual_info[[i]]), "V") == 1
    mutual_info[[i]] <- cbind(mutual_info[[i]][,vcols,drop=F], mutual_info[[i]][,!vcols])
  }
  out <- mutual_info
  for (i in 1:length(out)) {
    out[[i]] <- out[[i]][,stringr::str_count(colnames(mutual_info[[i]]), "V") == 1,drop=F]
    for (g in 1:ncol(groups)) {
      grp <- as.numeric(as.factor(groups[,g]))
      grp_nm <- colnames(groups)[g]
      num_cols <- ncol(out[[i]])
      out[[i]][,num_cols+1] <- rep(0,nrow(out[[i]]))
      null <- matrix(0, nrow=nrow(out[[i]]), ncol=bootstraps)
      for (j in 1:nrow(out[[i]])) {
        sub_mi <- as.numeric(mutual_info[[i]][j,(depth+1):ncol(mutual_info[[i]])])
        mat <- matrix(0, nrow=length(grp), ncol=length(grp))
        mat[lower.tri(mat)] <- sub_mi
        mat <- t(mat)                                                           # Rescaling KL divergence onto 0-1
        mat <- mat + t(mat)                                                     # eliminates concerns that KL is
        mat <- matrix(1, nrow=nrow(mat), ncol=ncol(mat)) - (mat / max(mat))     # inherently larger for higher-
        mat <- igraph::graph_from_adjacency_matrix(mat, weighted=T)             # order ensembles
        out[[i]][j,num_cols+1] <- igraph::modularity(mat, grp, weights=igraph::E(mat)$weight)
        for (k in 1:ncol(null)) {
          sub_mi <- sub_mi[sample(length(sub_mi))]
          mat <- matrix(0, nrow=length(grp), ncol=length(grp))
          mat[lower.tri(mat)] <- sub_mi
          mat <- t(mat)
          mat <- mat + t(mat)
          mat <- matrix(1, nrow=nrow(mat), ncol=ncol(mat)) - (mat / max(mat))
          mat <- igraph::graph_from_adjacency_matrix(mat, weighted=T)
          null[j,k] <- igraph::modularity(mat, grp, weights=igraph::E(mat)$weight)
        }
      }
      out[[i]][,num_cols+2] <- (out[[i]][,num_cols+1] - rowMeans(null)) / apply(null,1,sd)
      out[[i]][(apply(null,1,sd) == 0), num_cols+2] <- 0                        # If std dev is 0, make Z score 0
      out[[i]][,num_cols+3] <- pnorm(out[[i]][,num_cols+2], lower.tail=F)
      out[[i]][,num_cols+4] <- p.adjust(out[[i]][,num_cols+3], "BH")
      colnames(out[[i]])[(num_cols+1):(num_cols+4)] <- paste(c("TransMI_","Zscore_","Pvalue_","Padjust_"),
                                                             grp_nm, sep="")
    }
  }
  return(out)
}
