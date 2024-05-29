#' Census Images
#'
#' Draw neighborhoods on a list of images representing a single biological specimen. All images must
#' have the same dimensions. All images must be named as "O"bject or "S"calar, followed by a numeric code. If both
#' object maps and scalar sets are included, each scalar set may be independent or linked.
#' In the latter case, the object-scalar pairs of interest must be given with a link table
#' with the same name as the involved object map. For each image, seed points control
#' where neighborhoods may be centered. By default, seed points are wherever any object or scalar exists.
#' The radii and number of neighborhoods to draw at each radius should also be specified.
#'
#' @description Census neighborhoods across object maps and/or scalar sets from a single biological speciment
#' @param images        Named List: 4d array/s representing object maps and/or scalar sets
#' @param OS_pairs      Named List: link table/s that specify object-scalar pairs of interest
#' @param seed_points   Named List: vectors specifying which objects or scalars may be seed points. NULL for all is default; NA for none.
#' @param radii         Named List: vectors of X, Y, and Z pixel radii at each micron radius
#' @param sample_size   Numeric Vector: number of neighborhoods to draw at each radius
#' @param cores         Integer: number of cores to use for parallel calculations
#' @return      List: data frame of the census and nested list of the patch list
#' @export
census_image <- function(images, OS_pairs=NULL, seed_points=NULL, radii, sample_size, cores=NULL) {
  if (!is.list(images)) {
    stop("Error: Images must be provided in a named list.")
  }
  if (!is.null(OS_pairs) && !is.list(OS_pairs)) {
    stop("Error: Object-scalar pair matrices must be provided in a named list.")
  }
  if (!is.null(OS_pairs) && (!all(names(OS_pairs) %in% names(images)))) {
    stop("Error: The name of an object-scalar pair matrix does not match any of the image names.")
  }
  if (!is.null(seed_points) && !is.list(seed_points)) {
    stop("Error: The seed points must be provided in a named list.")
  }
  if (is.null(seed_points)) {
    seed_points <- vector("list", length=length(images))
    names(seed_points) <- names(images)
  }
  if (length(seed_points) != length(images)) {
    stop("Error: The number of images and the number of seed point sets do not match.")
  }
  if (!all(names(seed_points) %in% names(images))) {
    stop("Error: The name of a seed point set does not match any of the image names.")
  }
  if (length(radii) != length(sample_size)) {
    stop("Error: The number of sample sizes does not match the number of radii.")
  }
  print("Calculating expression thresholds for any scalar images.")
  image_types <- names(images)                                   # Extract image types.
  bin_thrs <- vector(mode="list", length=length(image_types))    # List of binary threshold vectors, 1 for each image
  for (i in 1:length(images)) {
    if(grepl("S",image_types[[i]])) {                            # For scalar layers, get binarization threshold for each scalar
      bin_thrs[[i]] <- rep(NA, dim(images[[i]])[4])              # This algorithm is very close to ImageJ's "IsoData" - basically:
      for (j in 1:length(bin_thrs[[i]])) {                       # Split image into bright vs. dim and find avg pixel value of each
        pix_hist <- unname(table(factor(as.vector(images[[i]][,,,j]), levels = 0:255)))
        if (sum(pix_hist > 0) == 1) {                            # If there is only one unique pixel value
          bin_thrs[[i]][j] <- which(pix_hist > 0) - 1            # That value becomes the threshold
        } else if (sum(pix_hist > 0) == 2) {                     # If there are only two unique pixel values
          bin_thrs[[i]][j] <- which(pix_hist > 0)                # The lower value + 1 becomes the threshold
        } else {
          pix_hist <- pix_hist[2:256]                            # Remove pixels with fluorescence value 0, because the number of
          thr <- min(which(pix_hist > 0)) + 1                    # these can depend on the amount of extra space in the image
          avg_lo <- sum(pix_hist[1:(thr-1)] * c(1:(thr-1))) / sum(pix_hist[1:(thr-1)])
          avg_hi <- sum(pix_hist[(thr+1):255] * c((thr+1):255)) / sum(pix_hist[(thr+1):255])
          while (round((avg_lo + avg_hi) / 2) != thr) {          # Does rounded avg of bright & dim avgs equal threshold?
            thr <- thr + 1                                       # If yes, stop. If no, increment threshold and try again.
            avg_lo <- sum(pix_hist[1:(thr-1)] * c(1:(thr-1))) / sum(pix_hist[1:(thr-1)])
            avg_hi <- sum(pix_hist[(thr+1):255] * c((thr+1):255)) / sum(pix_hist[(thr+1):255])
            if (is.nan(avg_hi)) {
              thr <- thr - 1
              break
            }
          }
          bin_thrs[[i]][j] <- thr
        }
      }
    } else {                                                     # For object layers, there is no binarization threshold
      bin_thrs[[i]] <- NA
    }
  }
  print("Choosing seedpoints for neighborhoods.")
  elig_seeds <- vector(mode="list", length=length(images))       # Maps of eligible seed points, 1 for each layer
  for (i in 1:length(images)) {
    if (grepl("O",image_types[i])) {                             # For object layers:
      if (is.null(seed_points[[i]])) {                           # NULL means all non-bkgd pixels are eligible
        elig_seeds[[i]] <- array(images[[i]] != 0, dim=dim(images[[i]])[1:3])
      } else {                                                   # Otherwise, only pix of specified object IDs are eligible
        elig_seeds[[i]] <- array(images[[i]] %in% seed_points[[i]], dim=dim(images[[i]])[1:3]) # NA means no pix are eligible
      }                                                          # Ignore linked scalars. For independent scalar layers:
    } else if (sum(substr(image_types,2,9) == substr(image_types,2,9)[i]) == 1) {
      if (is.null(seed_points[[i]])) {                           # NULL means pix pos for any scalar are eligible
        seed_points[[i]] <- 1:dim(images[[i]])[4]
      }
      if (is.na(seed_points[[i]][1])) {                          # NA means no pixels are eligible
        elig_seeds[[i]] <- array(F, dim=dim(images[[i]])[1:3])
      } else {
        elig_seeds[[i]] <- purrr::array_branch(images[[i]], margin=4)
        elig_seeds[[i]] <- mapply(function(p, v){p >= v}, p = elig_seeds[[i]], v = bin_thrs[[i]], SIMPLIFY = F)
        elig_seeds[[i]] <- elig_seeds[[i]][seed_points[[i]]]
        elig_seeds[[i]] <- Reduce("+", elig_seeds[[i]]) > 0      # Find pixels positive for at least one marker
      }
    }
  }
  elig_seeds <- purrr::compact(elig_seeds)                       # Remove NULL elig_seeds layers for linked scalars
  elig_seeds <- Reduce("+", elig_seeds) > 0                      # Pixels that are eligible in any layer are eligible overall
  elig_seeds <- which(elig_seeds, arr.ind=TRUE)                  # XYZ coordinates of eligible pixels, to sample from
  elig_seeds <- elig_seeds[sample(dim(elig_seeds)[1], max(sample_size)), ]
  slices_per_layer <- lapply(images, function(x){dim(x)[4]})     # Number of channels per layer
  image_types <- unlist(image_types)
  image_types <- paste(rep(image_types, slices_per_layer))
  images <- abind::abind(images)                                 # Merge all image layers into one 3D array
  bkgd_thrs <- unlist(bin_thrs)                                  # All binarization thresholds, substituting 1 as the
  bkgd_thrs[is.na(bkgd_thrs)] <- 1                               # threshold for object images
  glob_census <- vector("list", length(radii))
  patch_list <- vector("list", length(radii))
  for (i in 1:length(radii)) {
    print(paste("Beginning to census at radius ", names(radii)[i], " microns.", sep=""))
    num_cores <- min(parallel::detectCores() - 1, cores)
    print(paste("Preparing ", num_cores, " nodes to census in parallel.", sep=""))
    make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = make_cluster)            # For each local census, excise its bounding box
    `%dopar%` <- foreach::`%dopar%`                              # Attempt to help R find what %dopar% means
    print(paste("Censusing ", sample_size[i], " neighborhoods.", sep=""))
    nbhd_patches <- foreach::foreach (j = c(1:sample_size[i]), .packages=c('plyr','raster'), .export=c('census_patches','patch_3D')) %dopar% {
      sub_lyr <- images[max(elig_seeds[j,1] - radii[[i]][1], 1):min(elig_seeds[j,1] + radii[[i]][1], dim(images)[1]),
                        max(elig_seeds[j,2] - radii[[i]][2], 1):min(elig_seeds[j,2] + radii[[i]][2], dim(images)[2]),
                        max(elig_seeds[j,3] - radii[[i]][3], 1):min(elig_seeds[j,3] + radii[[i]][3], dim(images)[3]), , drop=F]
      incl_pix <- expand.grid(c(max(elig_seeds[j,1] - radii[[i]][1], 1):min(elig_seeds[j,1] + radii[[i]][1], dim(images)[1])),
                              c(max(elig_seeds[j,2] - radii[[i]][2], 1):min(elig_seeds[j,2] + radii[[i]][2], dim(images)[2])),
                              c(max(elig_seeds[j,3] - radii[[i]][3], 1):min(elig_seeds[j,3] + radii[[i]][3], dim(images)[3])))
      incl_pix <- t((t(incl_pix) - elig_seeds[j,]) / radii[[i]])
      incl_pix <- sqrt(rowSums((incl_pix)^2)) <= 1
      incl_pix <- array(incl_pix, dim=dim(sub_lyr)[1:3])
      census_patches(sub_lyr, incl_pix, image_types, bkgd_thrs) # Patches for just the current neighborhood
    }
    parallel::stopCluster(cl = make_cluster)
    print("Formatting the patch list and census data frame.")   # Switch nested indices to: length scale, layer, nbhd
    patch_list[[i]] <- purrr::list_transpose(nbhd_patches, template=c(1:length(nbhd_patches[[1]])), simplify=F)
    names(patch_list[[i]]) <- names(nbhd_patches[[1]])          # Restore names that are lost in list_transpose
    for (j in 1:length(patch_list[[i]])) {                      # In each layer, condense all nbhds to 1 data frame
      patch_list[[i]][[j]] <- dplyr::bind_rows(patch_list[[i]][[j]], .id = "Nbhd")
      patch_list[[i]][[j]]$Nbhd <- as.integer(patch_list[[i]][[j]]$Nbhd)
    }
    glob_census[[i]] <- summarize_patches(patch_list[[i]], OS_pairs)
    glob_census[[i]] <- cbind(glob_census[[i]], elig_seeds[1:sample_size[i],],
                              rep(as.numeric(names(radii)[i]), sample_size[i]))
    colnames(glob_census[[i]])[(ncol(glob_census[[i]])-3):ncol(glob_census[[i]])] <- c("X","Y","Z","Radius")
  }
  glob_census <- plyr::ldply(glob_census, rbind)                # Combine census frames for all local radii
  names(patch_list) <- names(radii)                             # Name the patch lists for each local radius
  return(list(glob_census, patch_list))
}

#' Census the Patches in a Neighborhood
#'
#' Find contiguous patches of pixels that belong to the same object or that are positive for the same scalar.
#' When scalars are linked to objects, the total scalar intensity on each object patch is recorded, as well.
#'
#' @description Record sizes of contiguous patches of pixels in a neighborhood, plus their object IDs and/or scalar intensities
#' @param sub_lyr      4d Numeric Array: bounding box of one neighborhood, excised from a larger image across all slices
#' @param incl_pix     3d Logical Array: mask of pixels included in a spherical neighborhood within a bounding box
#' @param image_type   Character Vector: Object map or scalar set ID of each slice of the image
#' @param bin_thr      Numeric Vector: binary threshold for each slice of the image
#' @return      Data frame: scalar intensity or object ID and pixel count for each patch
census_patches <- function(sub_lyr, incl_pix, image_type, bin_thr) {
  lyr_typ <- substr(image_type, 1, 1)                                           # Whether each layer is a "O"bject map or "S"calar set
  lyr_grp <- as.numeric(gsub("\\D", "", image_type))                            # Which grouping each layer belongs to
  out <- vector("list", length(unique(lyr_grp)))                                # Census patches separately for each grouping.
  for (i in 1:length(unique(lyr_grp))) {                                        # For each grouping:
    if (any(grepl("O", image_type[lyr_grp==i]))) {                              # If there is an object map
      obj_sub_lyr <- sub_lyr[,,,which(lyr_typ=="O" & lyr_grp==i),drop=F]        # extract it
      obj_sub_lyr <- array(obj_sub_lyr, dim=dim(incl_pix))                      # Drop 4th dim but not 3rd
      uniq_objs <- unique(obj_sub_lyr[incl_pix])                                # Find the unique objects
      uniq_objs <- uniq_objs[uniq_objs != 0]
      if (length(uniq_objs) == 0) {                                             # If there are no objects
        sub_out <- data.frame(Area = 0, V = 0)                                  # Make a dummy patch of area 0
        colnames(sub_out) <- c("Area", paste("O",i,sep=""))
        if (any(grepl("S", image_type[lyr_grp==i]))) {                          # If there are attached scalars
          sub_out <- cbind(sub_out, matrix(0, nrow = 1, ncol = sum(grepl("S", image_type[lyr_grp==i]))))
          colnames(sub_out)[3:ncol(sub_out)] <- paste("S",i,".",1:sum(grepl("S", image_type[lyr_grp==i])), sep="")
        }                                                                       # Include 0s for them, too
      } else {
        sub_out <- vector("list", length(uniq_objs))                            # Prepare to find clumps of each unique object
        for (j in c(1:length(uniq_objs))) {                                     # For each object
          clumps <- patch_3D(incl_pix & obj_sub_lyr==uniq_objs[j])              # find contiguous patches & sizes
          areas <- sapply(1:max(clumps[[1]]), function(x) {clumps[[2]][which(clumps[[1]] == x)[1]]})
          sub_out[[j]] <- data.frame(Area = areas, V = rep(uniq_objs[j], length(areas)))
          colnames(sub_out[[j]]) <- c("Area", paste("O",i,sep=""))              # Put patch areas & object in frame
          if (any(grepl("S", image_type[lyr_grp==i]))) {                        # If there are attached scalars
            sub_out[[j]] <- cbind(sub_out[[j]], matrix(0, nrow=nrow(sub_out[[j]]),# prepare extra frame columns for them
                                                       ncol=sum(grepl("S", image_type[lyr_grp==i]))))
            colnames(sub_out[[j]])[3:ncol(sub_out[[j]])] <- paste("S",i,".",1:sum(grepl("S", image_type[lyr_grp==i])), sep="")
            for (k in 1:nrow(sub_out[[j]])) {                                   # For each scalar inside each clump
              for (p in 1:sum(grepl("S", image_type[lyr_grp==i]))) {            # measure total expression
                sub_out[[j]][k,2+p] <- sum(sub_lyr[,,,which(lyr_typ=="S" & lyr_grp==i)[p]][clumps[[1]]==k])
              }
            }
          }
        }
      }
      sub_out <- dplyr::bind_rows(sub_out)
      out[[i]] <- sub_out
    } else {                                                                    # If grouping is independent scalars with no object map
      sub_out <- vector("list", sum(lyr_grp==i))
      for (j in 1:sum(lyr_grp==i)) {                                            # For each unattached scalar separately
        thr <- bin_thr[lyr_grp==i][j]                                           # Get the binarization threshold for this scalar
        brights <- NULL
        dims <- NULL
        mrk_sub_lyr <- sub_lyr[,,,which(lyr_grp==i)[j],drop=F]                  # Extract the slice for this scalar
        mrk_sub_lyr <- array(mrk_sub_lyr, dim=dim(incl_pix))                    # Drop 4th dim but not 3rd
        if (sum(incl_pix & (mrk_sub_lyr >= thr)) > 0) {                         # If there are any bright pixels
          clumps <- patch_3D(incl_pix & mrk_sub_lyr >= thr)                     # Find clumps & their areas
          areas <- sapply(1:max(clumps[[1]]), function(x) {clumps[[2]][which(clumps[[1]] == x)[1]]})
          brights <- data.frame(Area = areas, V = rep(NA, length(areas)))
          colnames(brights) <- c("Area", paste("S",i,".",j,sep=""))
          for (k in 1:nrow(brights)) {                                          # For each bright clump get total expr
            brights[k,2] <- sum(sub_lyr[,,,which(lyr_grp==i)[j]][clumps[[1]]==k])
          }
        }
        if (sum(incl_pix & (mrk_sub_lyr < thr)) > 0) {                          # If there are any dim pixels
          clumps <- patch_3D(incl_pix & mrk_sub_lyr < thr)                      # Find clumps & their areas
          areas <- sapply(1:max(clumps[[1]]), function(x) {clumps[[2]][which(clumps[[1]] == x)[1]]})
          dims <- data.frame(Area = areas, V = rep(NA, length(areas)))
          colnames(dims) <- c("Area", paste("S",i,".",j,sep=""))
          for (k in 1:nrow(dims)) {                                             # For each dim clump get total expr
            dims[k,2] <- sum(sub_lyr[,,,which(lyr_grp==i)[j]][clumps[[1]]==k])
          }
        }
        sub_out[[j]] <- dplyr::bind_rows(list(brights,dims))                    # Combine bright and dim patches into one frame
      }
      out[[i]] <- sub_out
    }
  }
  nst_lst_cts <- lapply(out, function(x) {ifelse(is.data.frame(x), 1, length(x))})
  if (any(nst_lst_cts > 1)) {
    new_out <- vector("list", length=sum(unlist(nst_lst_cts)))                  # If any of the list elements are themselves a list
    idx <- 1                                                                    # which happens for independent scalars, make each
    for (i in 1:length(out)) {                                                  # sub-entry its own element of the main list.
      if (is.data.frame(out[[i]]) || is.null(out[[i]])) {
        new_out[[idx]] <- out[[i]]
        idx <- idx + 1
      } else {
        for (j in 1:length(out[[i]])) {
          new_out[[idx]] <- out[[i]][[j]]
          idx <- idx + 1
        }
      }
    }
    out <- new_out
  }
  names(out) <- image_type[lyr_typ == "O" | (lyr_typ == "S" & !(lyr_grp %in% lyr_grp[lyr_typ=="O"]))]
  return(out)
}

#' Summarize a Patch List as a Census
#
#' From a patch list of the patches of each object and scalar inside each neighborhood of a specific radius,
#' create a census. Each row is a neighborhood and each column is a variable of interest, including
#' object abundance, avg scalar intensity at large, and avg scalar intensity attached to certain objects.
#'
#' @description Transform a patch list into a census
#' @param nbhd_patches List: Across all object maps and/or scalar sets, the contiguous patches in each neighborhood
#' @param OS_pairs     List: matrices that specify object-scalar pairs of interest
#' @return      Data frame: census of the objects and scalars in each neighborhood
summarize_patches <- function(nbhd_patches, OS_pairs) {
  num_grp <- length(nbhd_patches)                                               # Number of variable groupings included
  out <- vector("list", num_grp)
  for (i in 1:num_grp) {                                                        # For each variable grouping
    this_grp <- nbhd_patches[[i]]                                               # get just that grouping from each neighborhood
    num_nbhds <- max(this_grp$Nbhd)                                             # Record number of neighborhoods
    this_grp <- this_grp[this_grp$Area > 0, ]                                   # Get rid of rows with area = 0
    if (any(grepl("O", colnames(this_grp)))) {                                  # There will be an object column
      ocol <- which(grepl("O", colnames(this_grp)))                             # for objects & linked scalars
      ocol_nm <- colnames(this_grp)[ocol]
    } else {                                                                    # There will not be an object column
      ocol <- NULL                                                              # for independent scalars
    }
    if (!is.null(OS_pairs)) {                                                   # If link tables are given & one of
      if (names(nbhd_patches)[i] %in% names(OS_pairs)) {                        # them matches this obj set
        os_idx <- which(names(OS_pairs) == ocol_nm)
        os_pair <- as.data.frame(which(OS_pairs[[os_idx]] > 0, arr.ind=T))      # Use the provided matrix to include only
        os_pair[,1] <- rownames(OS_pairs[[os_idx]])[os_pair[,1]]                # the object-scalar pairs that are of interest
        os_pair[,2] <- colnames(OS_pairs[[os_idx]])[os_pair[,2]]
        os_pair <- unname(apply(os_pair, 1, paste, collapse="_"))
      } else {                                                                  # Otherwise, note that none of them
        os_pair <- NULL                                                         # match this set
      }
    } else {
      os_pair <- NULL
    }
    if (!is.null(ocol)) {                                                       # IF THERE ARE OBJECTS,
      colnames(this_grp)[ocol] <- "Obj"                                         # Give obj col a standardized name
      temp <- dplyr::summarise(dplyr::group_by(this_grp, Nbhd, Obj), Area = sum(Area)) # Total area per obj in each nbhd
      temp <- tidyr::pivot_wider(temp, names_from = Obj, values_from = Area)    # Make each obj its own col
      temp[is.na(temp)] <- 0                                                    # NAs are unobserved objs in each nbhd
      temp <- temp[,stringr::str_sort(colnames(temp), numeric=T)]               # Reorder columns; Write Area as % of
      temp[,1:(ncol(temp)-1)] <- 100 * (temp[,1:(ncol(temp)-1)] / rowSums(temp[,1:(ncol(temp)-1)])) # all non-bkgd area
      missnbhd <- setdiff(1:num_nbhds, temp$Nbhd)                               # Missing nbhds because no objects
      if (length(missnbhd) > 0) {                                               # If there are missing nbhds
        extra <- as.data.frame(matrix(0, nrow = length(missnbhd), ncol=ncol(temp))) # Addendum of 0s for missing nbhds
        colnames(extra) <- colnames(temp)                                       # Match the col names
        extra$Nbhd <- missnbhd                                                  # Fill in missing nbhd IDs
        temp <- rbind(temp, extra)                                              # Add to census
        temp <- temp[order(temp$Nbhd), ]                                        # Sort by nbhd
      }                                                                         # Give objs original names
      colnames(temp)[1:(ncol(temp)-1)] <- paste(ocol_nm, colnames(temp)[1:(ncol(temp)-1)], sep=".")
      if (!is.null(OS_pairs) && !is.null(os_pair)) {                            # IF THERE ARE LINKED SCALARS
        temp2 <- dplyr::summarise(dplyr::group_by(this_grp, Nbhd, Obj),         # Total scalar expr on each obj in each nbhd
                                  dplyr::across(dplyr::starts_with(c("S","A")), sum)) # plus total area of each obj
        temp2 <- dplyr::mutate(temp2, dplyr::across(dplyr::starts_with("S"), ~(. / Area))) # Avg scalar expr per pix on each obj
        temp2$Area <- NULL                                                      # Drop area - no longer necessary
        temp2 <- tidyr::pivot_wider(temp2, names_from = Obj, values_from = dplyr::starts_with("S"), # Make each obj-scl
                                    names_glue = paste(ocol_nm,".","{Obj}_{.value}",sep="")) # combo its own col
        temp2 <- temp2[,c(os_pair[os_pair %in% colnames(temp2)], "Nbhd")]       # Only present obj-scl combos of interest
        temp2[is.na(temp2)] <- 0                                                # NAs are unobserved objs
        missnbhd <- setdiff(1:num_nbhds, temp2$Nbhd)                            # Missing nbhds because no objects
        if (length(missnbhd) > 0) {                                             # If there are missing nbhds
          extra <- as.data.frame(matrix(0, nrow = length(missnbhd), ncol=ncol(temp2))) # Addendum of 0s for missing nbhds
          colnames(extra) <- colnames(temp2)                                    # Match the col names
          extra$Nbhd <- missnbhd                                                # Fill in missing nbhd IDs
          temp2 <- rbind(temp2, extra)                                          # Add to census
          temp2 <- temp2[order(temp2$Nbhd), ]                                   # Sort by nbhd
        }                                                                       # Give objs original names
        out[[i]] <- cbind(temp[,1:(ncol(temp)-1)], temp2[,1:(ncol(temp2)-1)])   # Combine obj & obj-scl, without Nbhd IDs
      } else {                                                                  # If no attached scalars, just return
        out[[i]] <- temp[,1:(ncol(temp)-1), drop=F]                             # obj without Nbhd ID
      }
    } else {                                                                    # IF THIS IS AN INDEPENDENT SCALAR
      scol <- which(grepl("S", colnames(this_grp)))                             # Col ID for the scalar measurement
      scol_nm <- colnames(this_grp)[scol]                                       # Original name of the scalar
      colnames(this_grp)[scol] <- "Scl"                                         # Give scalar col a standardized name
      temp <- dplyr::summarise(dplyr::group_by(this_grp, Nbhd),                 # Total scalar expr in each nbhd
                               dplyr::across(dplyr::starts_with(c("S","A")), sum)) # plus total area of each nbhd
      temp$Scl <- temp$Scl / temp$Area                                          # Avg scl expr per pix in each nbhd
      temp$Area <- NULL                                                         # No longer need nbhd area
      temp[is.na(temp)] <- 0                                                    # NAs shouldn't be possible anyway
      missnbhd <- setdiff(1:num_nbhds, temp$Nbhd)                               # Missing nbhds shouldn't be possible
      if (length(missnbhd) > 0) {                                               # But if there are any
        extra <- as.data.frame(matrix(0, nrow = length(missnbhd), ncol=ncol(temp))) # prepare for them anyway
        colnames(extra) <- colnames(temp)                                       # Match the col names
        extra$Nbhd <- missnbhd                                                  # Fill in missing nbhd IDs
        temp <- rbind(temp, extra)                                              # Add to census
        temp <- temp[order(temp$Nbhd), ]                                        # Sort by nbhd
      }
      colnames(temp)[2] <- scol_nm                                              # Give scalar original names
      out[[i]] <- temp[,2, drop=F]
    }
  }
  out <- dplyr::bind_cols(out)
  return(out)
}

#' Census a Table
#'
#' Draw neighborhoods on an object table to record the contents of each. The object table must give the
#' centroid coordinates in columns labeled "X", "Y", and "Z", and the numeric code for object ID in a column labeled
#' "Object". Optional further columns may give the intensity of scalars. If scalar intensities are included, the
#' object-scalar pairs of interest must be specified in a link table with the same name as the involved object map.
#' Seed points control where neighborhoods may be centered - by default, wherever any object exists.
#' The radii and number of neighborhoods should also be specified.
#'
#' @description Census neighborhoods from a table of objects
#' @param object_table   Data Frame: centroid coordinates, object type, and mean expression of scalars
#' @param OS_pairs       Named List: link table that specifies object-scalar pairs of interest
#' @param seed_points    Named List: vector specifying which objects may be seed points, or NULL for no restrictions
#' @param radii          Vector: radii of neighborhoods, in the same units as the centroid coordinates, typically microns
#' @param sample_size    Vector: number of neighborhoods to draw at each radius, no more than the number of objects
#' @param cores          Integer: number of cores to use for parallel calculations
#' @return      List: data frame of the census and nested list of the patch list
#' @export
census_table <- function(object_table, OS_pairs=NULL, seed_points=NULL, radii, sample_size, cores=NULL) {
  if (!is.data.frame(object_table) || !all(c("X","Y","Z","Object") %in% colnames(object_table))) {
    stop("Error: Object table must be a data frame with at least X, Y, Z, and Object columns.")
  }
  if (!is.null(OS_pairs) && !is.list(OS_pairs) && !is.matrix(OS_pairs)) {
    stop("Error: Object-scalar pair matrix must be provided as a named list of one link table.")
  }
  if (is.matrix(OS_pairs)) {
    OS_pairs <- list("O1" = OS_pairs)
  }
  if (is.null(seed_points)) {
    seed_points <- unique(object_table$Object)
  } else if (is.list(seed_points)) {
    seed_points <- unname(unlist(seed_points))
  }
  if (length(radii) != length(sample_size)) {
    stop("Error: The number of sample sizes does not match the number of radii.")
  }
  glob_census <- vector(mode="list", length=length(radii))                      # Set up the census
  patch_list <- vector(mode="list", length=length(radii))                       # Set up the patch list
  names(patch_list) <- radii                                                    # Name of 1st patch list idx is radius
  for (i in 1:length(radii)) {
    print(paste("Beginning to census at radius ", radii[i], " microns.", sep=""))
    seeds <- sample(which(object_table$Object %in% seed_points), sample_size[i])
    print("Finding all neighbors within the radius.")
    number_nbrs <- 10
    repeat {                                                                    # Nearest k nbrs based on guess at k
      distances <- FNN::knnx.dist(object_table[,c("X","Y","Z")], object_table[seeds,c("X","Y","Z")], k = number_nbrs)
      if (sum(distances[,number_nbrs] <= radii[i]) == 0) {                      # If kth nbr is always outside radius
        break                                                                   # this is enough info
      } else {                                                                  # Otherwise
        number_nbrs <- number_nbrs + 5                                          # increase number of nbrs and try again
      }
    }                                                                           # Then get row ID of neighbors
    neighbors <- FNN::knnx.index(object_table[,c("X","Y","Z")], object_table[seeds,c("X","Y","Z")], k = number_nbrs)
    neighbors[distances > radii[i]] <- NA                                       # Remove neighbors not within radius
    print(paste("Censusing ", sample_size[i], " neighborhoods.", sep=""))
    idxs <- as.vector(t(neighbors))                                             # Reformat row IDs to vector
    idxs <- idxs[!is.na(idxs)]
    obj_cts <- rowSums(!is.na(neighbors))                                       # Count of objects in each nbhd
    nbhd_patches <- data.frame(matrix(0, nrow=length(idxs), ncol=ncol(object_table)-2))
    colnames(nbhd_patches)[1:2] <- c("Area", "O1")
    nbhd_patches$Area <- rep(1, length(idxs))                                   # Every object has area 1
    nbhd_patches$O1 <- object_table$Object[idxs]                                # Record each object
    if (ncol(object_table) > 4) {                                               # If there are scalars too
      num_sclr <- ncol(object_table) - 4                                        # Record scalar expr for each object
      colnames(nbhd_patches)[3:(2+num_sclr)] <- paste("S1.", 1:num_sclr,sep="")
      nbhd_patches[,3:(2+num_sclr)] <- object_table[idxs,5:(4+num_sclr)]
    }
    nbhd_patches <- split(nbhd_patches, rep(1:nrow(neighbors), times=obj_cts))  # Make a list w/ 1 element per nbhd
    nbhd_patches <- lapply(nbhd_patches, function(x) {list(x)})                 # Make each list element a list itself
    print("Formatting the patch list and census data frame.")
    patch_list[[i]] <- purrr::list_transpose(nbhd_patches, template=c(1:length(nbhd_patches[[1]])), simplify=F)
    names(patch_list[[i]]) <- names(OS_pairs)                                   # Name the object set, should be "O1"
    for (j in 1:length(patch_list[[i]])) {                                      # Condense all nbhds to 1 data frame
      patch_list[[i]][[j]] <- dplyr::bind_rows(patch_list[[i]][[j]], .id = "Nbhd")
      patch_list[[i]][[j]]$Nbhd <- as.integer(patch_list[[i]][[j]]$Nbhd)
    }
    glob_census[[i]] <- summarize_patches(patch_list[[i]], OS_pairs)
    glob_census[[i]] <- cbind(glob_census[[i]], object_table[seeds,c("X","Y","Z")], rep(radii[i], sample_size[i]))
    colnames(glob_census[[i]])[(ncol(glob_census[[i]])-3):ncol(glob_census[[i]])] <- c("X","Y","Z","Radius")
  }
  glob_census <- plyr::ldply(glob_census, rbind)                                # Combine census frames for all local radii
  return(list(glob_census, patch_list))
}

#' Standardize Columns across Censuses
#'
#' For a list of censuses meant to be analyzed together, enforce that censuses contain the same columns
#' in the same order.
#'
#' @description Enforce that all censuses contain the same columns in the same order
#' @param censuses  List: Census data frames for every image in a multi-image data set
#' @return      List: Census data frames for every image in a multi-image data set, with standardized columns
#' @export
standardize_censuses <- function(censuses) {
  col_nms <- unique(unlist(lapply(censuses, colnames)))                         # Unique col names across all censuses
  col_nms <- col_nms[!(col_nms %in% c("X","Y","Z","Radius"))]                   # Remove non-variable cols temporarily
  col_nms <- stringr::str_sort(col_nms, numeric = TRUE)                         # Sort col names into an order
  col_nms <- c(col_nms, c("X","Y","Z","Radius"))                                # Re-add non-variable cols
  for (i in 1:length(censuses)) {                                               # For each census
    missing_cols <- setdiff(col_nms, colnames(censuses[[i]]))                   # Find which vars are missing
    censuses[[i]][missing_cols] <- 0                                            # Add cols of all 0s for missing vars
    censuses[[i]] <- censuses[[i]][col_nms]                                     # Put columns into standard order
  }
  return(censuses)
}
