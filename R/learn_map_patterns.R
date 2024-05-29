#' Learn a Spatial Pattern
#'
#' Draw a curve through the regions of highest data density on an ensemble's joint distribiution.
#' All neighborhoods in the census are reordered to follow the trajectory of this curve.
#' This reveals how the variables in the ensemble covary across the entire spatial extent.
#' Enrichment of features is plotted versus a reference - randomized expectations if a patch list is given,
#' or a different group of specimens if focal and reference groups are defined..
#'
#' @description Learn the spatial patterning of an ensemble by drawing a self-organized map on its joint distribution
#' @param census         Data frame: census for a single image, or list of censuses for multiple images
#' @param ensemble       Vector: variables included in the ensemble of interest
#' @param radius         Numeric: radius of neighborhoods to use
#' @param col_pal        List: named list of the color palettes encompassing all variables in the ensemble
#' @param group          Frame: data frame with grouping factors as columns and images as rows
#' @param focal          Character: group/s for which to learn the covariation pattern
#' @param reference      Character: group/s against which to compare the covariation pattern
#' @param patch_list     Nested List: patch list must be included if comparisons are made vs. random
#' @param conf_int       Numeric: confidence interval of local representation to display along latent path
#' @param obj_norm       Character: how to normalize local abundance of objects: NULL for none, "ind" for individually, "all" for collectively
#' @param scl_norm       Character: how to normalize local abundance of scalars: NULL for none, "ind" for individually, "all" for collectively
#' @param som_reps       Numeric: number of learning repetitions to run the self-organizing map
#' @param toroidal       Logical: whether the self-organizing map should be toroidal or not
#' @param smooth_window  Numeric: number of neighborhoods to smooth together
#' @param sub_sample     Logical: whether to subsample census neighborhoods for speed
#' @param plot_bkgd      Character: whether to plot on a "W"hite or "B"lack background
#' @return      Data frame: local content of each variable along the latent path, also returned graphically
#' @export
learn_pattern <- function(census, ensemble, radius, col_pal, group = NULL, focal = NULL, reference = NULL,
                          patch_list = NULL, conf_int = 0.95, obj_norm = NULL, scl_norm = "all",
                          som_reps = 50, toroidal = F, smooth_window = 100, sub_sample = F, plot_bkgd = "W") {
  if ((!is.data.frame(census) || !is.null(group)) && !is.null(patch_list)) {
    stop("Error: Provide a patch list for one census, OR provide groups for multiple censuses, but NOT both.")
  }
  if (is.data.frame(census)) {                                                  # If only one census supplied
    census <- list(census)                                                      # put it in a list anyway
  }
  for (i in 1:length(census)) {                                                 # Narrow each census to
    census[[i]] <- census[[i]][census[[i]]$Radius == radius, ]                  # just correct radius
  }
  if (length(census) > 1) {                                                     # If >1 census is given
    if (!is.null(group)) {                                                      # If groups are specified
      grp_fct <- which(apply(group, 2, function(x,y){all(y %in% x)}, y = focal) &
                       apply(group, 2, function(x,y){all(y %in% x)}, y = reference))
      if (length(grp_fct) > 1) {                                                # Infer focal & reference groups
        stop("Error: The focal and reference groups provided match multiple grouping factors.")
      } else if (length(grp_fct) == 0) {
        stop("Error: The focal and reference groups provdied do not match any grouping factors.")
      }
      foc <- dplyr::bind_rows(census[which(group[,grp_fct] %in% focal)])        # pool all focal censuses
      ref <- dplyr::bind_rows(census[which(group[,grp_fct] %in% reference)])    # and reference censuses
    } else {                                                                    # If no groups are specified
      foc <- dplyr::bind_rows(census)                                           # pool all censuses
      ref <- NULL
    }
  } else {                                                                      # If one census given
    foc <- census[[1]]                                                          # remove it from list
    if (!is.null(patch_list)) {                                                 # If a patch list is given,
      patch_list <- patch_list[[which(as.numeric(names(patch_list)) == radius)]]# only consider correct radius
      os_pairs <- colnames(foc)[which(stringr::str_count(colnames(foc), "O|S") == 2)] # extract any obj-scl pairs
      if (length(os_pairs) > 0) {
        obj_maps <- unique(gsub("\\..*","", os_pairs))                          # Object maps with attached scalars
        OS_pairs <- vector("list", length=length(obj_maps))                     # List for pairing matrices
        names(OS_pairs) <- obj_maps                                             # Name matrices by their object map
        for (i in 1:length(OS_pairs)) {                                         # For each object map represented
          sub_os_pairs <- os_pairs[which(stringr::str_count(os_pairs, obj_maps[i]) == 1)]
          sub_scalars <- sub(".*_", "", sub_os_pairs)                           # Get scalar IDs and object IDs
          sub_objs <- sub("_.*", "", sub_os_pairs)
          OS_pairs[[i]] <- matrix(0, nrow=length(unique(sub_objs)), ncol=length(unique(sub_scalars)))
          rownames(OS_pairs[[i]]) <- unique(sub_objs)                           # Create and fill in a pairing matrix
          colnames(OS_pairs[[i]]) <- unique(sub_scalars)
          for (j in 1:length(sub_os_pairs)) {
            OS_pairs[[i]][rownames(OS_pairs[[i]]) == sub_objs[j], colnames(OS_pairs[[i]]) == sub_scalars[j]] <- 1
          }
        }
      } else {
        OS_pairs <- NULL
      }
      print("Simulating a random census to use as a reference.")
      ref <- random_census(patch_list, OS_pairs)                                # Build reference random census
    } else {
      ref <- NULL
    }
  }
  foc <- foc[, ensemble, drop=F]
  if (!is.null(ref)) {
    ref <- ref[, ensemble, drop=F]
  }
  print("Constructing the latent path on the focal census(es).")
  if ((nrow(foc) > 10000) && sub_sample) {
    print(paste("     Subsampling focal neighborhoods from ", nrow(foc), " to 10000.", sep=""))
    foc <- foc[sample(nrow(foc), 10000, replace=F), ]
  }
  if (!is.null(ref)) {
    if ((nrow(ref) > 10000) && sub_sample) {
      print(paste("     Subsampling reference neighborhoods from ", nrow(ref), " to 10000.", sep=""))
      ref <- ref[sample(nrow(ref), 10000, replace=F), ]
    }
  }
  num_nodes <- min(nrow(foc), 1000)
  SOM <- kohonen::som(as.matrix(foc),                                           # Self-organizing map on focal censuses
                      grid=kohonen::somgrid(xdim=num_nodes, ydim=1, topo="rectangular",
                                            neighbourhood.fct="gaussian", toroidal=toroidal),
                      rlen=som_reps, alpha=c(0.05, 0.01))
  SOM <- as.data.frame(SOM$codes[[1]])                                          # Extract coordinates of each node
  print("Reordering neighborhoods along the latent path.")                      # Choose a starting point in the SOM curve
  if (toroidal) {                                                               # For toroidal curves:
    SOM_to_origin <- unname(sqrt(rowSums(SOM^2)))                               # Find the SOM node closest to the origin
    sto_min <- min(SOM_to_origin)                                               # Find the range of SOM nodes within 1%
    where_min <- which.min(SOM_to_origin)                                       # of distance-to-origin surrounding the closest
    where_low <- which(SOM_to_origin < sto_min + 0.01*abs(diff(range(SOM_to_origin))))
    num_before <- tail(rle(c(1, diff(where_low))[1:which(where_low==where_min)])$lengths, 1)
    num_after <- head(rle(c(1, diff(where_low))[which(where_low==where_min):length(where_low)])$lengths, 1)
    if (num_before <= num_after) {                                              # Start on shorter side of this run of nodes
      start <- where_min - num_before + 1                                       # and proceed toward the longer, where the
      if (start > 1) {                                                          # node closest to the origin is the reference
        SOM <- SOM[c(start:num_nodes, 1:(start-1)), , drop=F]                   # point somewhere in the middle
      }
    } else {
      start <- where_min + num_after - 1
      if (start < num_nodes) {
        SOM <- SOM[c(start:1, num_nodes:(start+1)), , drop=F]
      } else {
        SOM <- SOM[num_nodes:1, ,drop=F]
      }
    }
  } else {                                                                      # For non-toroidal curves:
    if (sum(SOM[1,]^2) > sum(SOM[num_nodes,]^2)) {                              # Start on the end closer to the origin
      SOM <- SOM[num_nodes:1, , drop=F]
    }
  }
  SOM_nn <- FNN::knnx.index(data=SOM, query=foc, k=1)                           # Find nearest SOM node for each nbhd, & which SOM
  foc <- foc[order(SOM_nn), ,drop=F]                                            # Reorder nbhds according to their closest SOM node
  half_window <- round(smooth_window/2)                                         # Half of smoothing window size, for rolling calcs
  out <- data.frame(X=rep(seq(0, 100, length.out = nrow(foc)), length(ensemble)),
                    Y=rep(NA, length(ensemble) * nrow(foc)), Ymin = rep(NA, length(ensemble) * nrow(foc)),
                    Ymax = rep(NA, length(ensemble) * nrow(foc)), V=rep(NA, length(ensemble) * nrow(foc)),
                    Enr = rep(NA, length(ensemble) * nrow(foc)))
  print("Calculating and smoothing rolling variable abundances.")
  for (i in 1:length(ensemble)) {                                               # For each variable, smooth its rolling conf int
    x <- c(sample(foc[,i][sort(SOM_nn) == min(SOM_nn)], half_window, replace=T), foc[,i],
           sample(foc[,i][sort(SOM_nn) == max(SOM_nn)], half_window, replace=T))
    out$Y[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <- unname(zoo::rollapply(x, 2 * half_window + 1, mean))
    out$Ymin[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <-
      unname(zoo::rollapply(x, 2 * half_window + 1, quantile, (1-conf_int)/2))
    out$Ymax[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <-
      unname(zoo::rollapply(x, 2 * half_window + 1, quantile, 1-(1-conf_int)/2))
    out$V[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <- colnames(foc)[i]
    x <- out$Ymin[((i-1) * nrow(foc) + 1):(i * nrow(foc))]
    x <- c(sample(x[1:half_window],half_window,replace=T), x, sample(x[(length(x)-half_window):length(x)],half_window,replace=T))
    out$Ymin[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <- unname(zoo::rollapply(x, 2 * half_window + 1, mean))
    x <- out$Ymax[((i-1) * nrow(foc) + 1):(i * nrow(foc))]
    x <- c(sample(x[1:half_window],half_window,replace=T), x, sample(x[(length(x)-half_window):length(x)],half_window,replace=T))
    out$Ymax[((i-1) * nrow(foc) + 1):(i * nrow(foc))] <- unname(zoo::rollapply(x, 2 * half_window + 1, mean))
  }

  if (!is.null(ref)) {                                                          # If there will be a reference
    SOM_nn_ref <- FNN::knnx.index(data=SOM, query=ref, k=1)                     # Assign each ref nbhd to a SOM node
    ref <- ref[order(SOM_nn_ref), ,drop=F]                                      # Reorder ref nbhds by nearest SOM node
    F_foc <- unname(table(factor(SOM_nn, levels = 1:num_nodes)))                # Number of focal nbhds assigned to each SOM node
    F_ref <- unname(table(factor(SOM_nn_ref, levels = 1:num_nodes)))            # Number of ref nbhds assigned to each SOM node
    SOM_foc <- SOM[sort(SOM_nn), , drop=F]                                      # SOM node corresponding to each foc nbhd
    SOM_ref <- SOM[sort(SOM_nn_ref), , drop=F]                                  # SOM node corresponding to each ref nbhd
    D_foc <- unname(sqrt(rowSums(((foc - SOM_foc)/100)^2)))                     # Dist from each foc nbhd to its SOM node
    D_ref <- unname(sqrt(rowSums(((ref - SOM_ref)/100)^2)))                     # Dist from each ref nbhd to its SOM node
    W_foc <- (sqrt(ncol(foc)) - D_foc) / sqrt(ncol(foc))                        # Weight of each foc nbhd (1 - % of max distance)
    W_ref <- (sqrt(ncol(ref)) - D_ref) / sqrt(ncol(ref))                        # Weight of each ref nbhd (1 - % of max distance)
    F_foc[F_foc > 0] <- aggregate(W_foc, by=list(sort(SOM_nn)), FUN = sum)$x    # Replace # of foc nbhds assigned to each SOM node w/ summed weights
    F_ref[F_ref > 0] <- aggregate(W_ref, by=list(sort(SOM_nn_ref)), FUN = sum)$x# Replace # of ref nbhds assigned to each SOM node w/ summed weights
    F_foc <- F_foc / nrow(foc)                                                  # Divide by total # of foc nbhds
    F_ref <- F_ref / nrow(ref)                                                  # Divide by total # of ref nbhds
    SOM_a <- log(F_foc / F_ref)                                                 # Calculate enrichment metric
    SOM_a[F_foc == 0] <- 0                                                      # If 0 foc nbhds assigned to SOM node, won't affect covariation plot
    SOM_a[F_foc > 0 & F_ref == 0] <- max(SOM_a[!is.infinite(SOM_a)])            # If 0 ref nbhds assigned to SOM node, give max observed enrichment
    x <- SOM_a[sort(SOM_nn)]                                                    # enrichment for each nbhd
    x <- c(sample(x[1:half_window],half_window,replace=T), x, sample(x[(length(x)-half_window):length(x)],half_window,replace=T))
    x <- unname(zoo::rollapply(x, 2 * half_window + 1, mean))                   # Rolling mean of enrichment
    out$Enr <- rep(x, length(ensemble))                                         # repeated for each var in ens
  } else {                                                                      # Without ref, no enrichment
    out$Enr <- rep(1, nrow(out))
  }
  out$Ymin[out$Y < out$Ymin] <- out$Y[out$Y < out$Ymin]
  out$Ymax[out$Y > out$Ymax] <- out$Y[out$Y > out$Ymax]
  scl_vars <- unique(out$V)[grepl("S",unique(out$V))]
  obj_vars <- unique(out$V)[!grepl("S",unique(out$V))]
  out$Y_norm <- out$Y
  out$Ymin_norm <- out$Ymin
  out$Ymax_norm <- out$Ymax
  if (!is.null(scl_norm) && (length(scl_vars) > 0)) {
    if (scl_norm == "all") {
      for (i in 1:length(scl_vars)) {
        out$Y_norm[out$V == scl_vars[i]] <- 100 * (out$Y[out$V == scl_vars[i]] / max(out$Ymax[out$V %in% scl_vars]))
        out$Ymin_norm[out$V == scl_vars[i]] <- 100 * (out$Ymin[out$V == scl_vars[i]] / max(out$Ymax[out$V %in% scl_vars]))
        out$Ymax_norm[out$V == scl_vars[i]] <- 100 * (out$Ymax[out$V == scl_vars[i]] / max(out$Ymax[out$V %in% scl_vars]))
      }
    } else if (scl_norm == "ind") {
      for (i in 1:length(scl_vars)) {
        out$Y_norm[out$V == scl_vars[i]] <- 100 * (out$Y[out$V == scl_vars[i]] / max(out$Ymax[out$V == scl_vars[i]]))
        out$Ymin_norm[out$V == scl_vars[i]] <- 100 * (out$Ymin[out$V == scl_vars[i]] / max(out$Ymax[out$V == scl_vars[i]]))
        out$Ymax_norm[out$V == scl_vars[i]] <- 100 * (out$Ymax[out$V == scl_vars[i]] / max(out$Ymax[out$V == scl_vars[i]]))
      }
    }
  }
  if (!is.null(obj_norm) && (length(obj_vars) > 0)) {
    if (obj_norm == "all") {
      for (i in 1:length(obj_vars)) {
        out$Y_norm[out$V == obj_vars[i]] <- 100 * (out$Y[out$V == obj_vars[i]] / max(out$Ymax[out$V %in% obj_vars]))
        out$Ymin_norm[out$V == obj_vars[i]] <- 100 * (out$Ymin[out$V == obj_vars[i]] / max(out$Ymax[out$V %in% obj_vars]))
        out$Ymax_norm[out$V == obj_vars[i]] <- 100 * (out$Ymax[out$V == obj_vars[i]] / max(out$Ymax[out$V %in% obj_vars]))
      }
    } else if (obj_norm == "ind") {
      for (i in 1:length(obj_vars)) {
        out$Y_norm[out$V == obj_vars[i]] <- 100 * (out$Y[out$V == obj_vars[i]] / max(out$Ymax[out$V == obj_vars[i]]))
        out$Ymin_norm[out$V == obj_vars[i]] <- 100 * (out$Ymin[out$V == obj_vars[i]] / max(out$Ymax[out$V == obj_vars[i]]))
        out$Ymax_norm[out$V == obj_vars[i]] <- 100 * (out$Ymax[out$V == obj_vars[i]] / max(out$Ymax[out$V == obj_vars[i]]))
      }
    }
  }

  col_pal_col <- rep(NA, length(ensemble))                                      # New col pals for col & fill
  col_pal_fil <- rep(NA, length(ensemble))
  for (i in 1:length(ensemble)) {                                               # Extract color for each var in ens
    if (grepl("_", ensemble[i])) {                                              # If var is linked obj-scl
      scl <- sub('.*_', '', ensemble[i])                                        # extract scalar
      obj <- sub('_.*', '', ensemble[i])                                        # extract object
      pal <- col_pal[[which(names(col_pal) == sub("\\..*", "", scl))]]          # Color palette for scalar
      var <- as.numeric(sub('.*\\.', '', scl))                                  # Var num for scalar
      col_pal_fil[i] <- pal[var]                                                # Scalar gives fill color
      pal <- col_pal[[which(names(col_pal) == sub("\\..*", "", obj))]]          # Color palette for object
      var <- as.numeric(sub('.*\\.', '', obj))                                  # Var num for object
      col_pal_col[i] <- pal[var]                                                # Object gives line color
    } else {                                                                    # If var is obj or scl
      pal <- col_pal[[which(names(col_pal) == sub("\\..*", "", ensemble[i]))]]  # Color palette for this var
      var <- as.numeric(sub('.*\\.', '', ensemble[i]))                          # Number of this var
      col_pal_col[i] <- pal[var]                                                # Extract this color
      col_pal_fil[i] <- pal[var]                                                # for color and fill
    }
  }
  tick_labels <- seq(0, 100, by=1)                                              # Custom tick labels to ease visual subdivisions
  tick_labels[which(tick_labels%%10 != 0)] <- ""
  len_scl <- paste(radius, " Microns", sep="")
  ttl <- paste("Spatial Pattern for ", paste(ensemble, sep="", collapse=", "), " at ", len_scl, sep="")
  if (plot_bkgd == "W") {
    if (is.null(ref)) {
      plot_A <- ggplot2::ggplot(out, ggplot2::aes(x=X, col=forcats::fct_inorder(as.factor(V)), fill=forcats::fct_inorder(as.factor(V)))) +
        ggplot2::geom_line(ggplot2::aes(y=Y_norm), linewidth=1) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Ymin_norm, ymax=Ymax_norm), alpha=0.5) +
        ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
        ggplot2::scale_discrete_manual(aesthetics=c("col"), values=col_pal_col, name="Variable") +
        ggplot2::scale_discrete_manual(aesthetics=c("fill"), values=col_pal_fil, name="Variable") +
        ggplot2::theme(text=ggplot2::element_text(size=15), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(expand=c(0,0), name="Latent Path (%)", breaks=seq(0, 100, by=1), labels=tick_labels) +
        ggplot2::scale_y_continuous(expand=c(0,0), name="Local Representation (%)")
      print(plot_A)
    } else {
      plot_A <- ggplot2::ggplot(out, ggplot2::aes(x=X, col=forcats::fct_inorder(as.factor(V)), fill=forcats::fct_inorder(as.factor(V)))) +
        ggplot2::geom_line(ggplot2::aes(y=Y_norm), linewidth=1) + ggplot2::theme_bw() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Ymin_norm, ymax=Ymax_norm), alpha=0.5) +
        ggplot2::scale_discrete_manual(aesthetics=c("col"), values=col_pal_col, name="Variable") +
        ggplot2::scale_discrete_manual(aesthetics=c("fill"), values=col_pal_fil, name="Variable") +
        ggplot2::theme(text=ggplot2::element_text(size=15), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(expand=c(0,0), name="Latent Path (%)", breaks=seq(0, 100, by=1), labels=tick_labels) +
        ggplot2::scale_y_continuous(expand=c(0,0), name="Local Representation (%)")
      plot_B <- ggplot2::ggplot(out, ggplot2::aes(x=X, y=Enr)) +
        ggplot2::geom_line(linewidth=1, color="black") + ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept=0, linetype="dashed") +
        ggplot2::theme(text=ggplot2::element_text(size=15), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank()) + ggplot2::ggtitle(ttl) +
        ggplot2::scale_x_continuous(expand=c(0,0), name=NULL, breaks=seq(0, 100, by=1), labels=NULL) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1), name="Enrichment")
      print(cowplot::plot_grid(plotlist=list(plot_B, plot_A), ncol=1, align="v", axis="lr", rel_heights=c(1, 3)))
    }
  } else if (plot_bkgd == "B") {
    if (is.null(ref)) {
      plot_A <- ggplot2::ggplot(out, ggplot2::aes(x=X, col=forcats::fct_inorder(as.factor(V)), fill=forcats::fct_inorder(as.factor(V)))) +
        ggplot2::geom_line(ggplot2::aes(y=Y_norm), linewidth=1) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Ymin_norm, ymax=Ymax_norm), alpha=0.5) +
        ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
        ggplot2::scale_discrete_manual(aesthetics=c("col"), values=col_pal_col, name="Variable") +
        ggplot2::scale_discrete_manual(aesthetics=c("fill"), values=col_pal_fil, name="Variable") +
        ggplot2::theme(text=ggplot2::element_text(size=15, color="white"), plot.background=ggplot2::element_rect(fill="black", color="black"),
                       panel.border=ggplot2::element_rect(color="white"),
                       panel.background=ggplot2::element_rect(fill="black"), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(), axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"), axis.ticks = ggplot2::element_line(colour = "white"),
                       legend.key = ggplot2::element_rect(fill = "black"), legend.background = ggplot2::element_rect(fill="black")) +
        ggplot2::scale_x_continuous(expand=c(0,0), name="Latent Path (%)", breaks=seq(0, 100, by=1), labels=tick_labels) +
        ggplot2::scale_y_continuous(expand=c(0,0), name="Local Representation (%)")
      print(plot_A)
    } else {
      plot_A <- ggplot2::ggplot(out, ggplot2::aes(x=X, col=forcats::fct_inorder(as.factor(V)), fill=forcats::fct_inorder(as.factor(V)))) +
        ggplot2::geom_line(ggplot2::aes(y=Y_norm), linewidth=1) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Ymin_norm, ymax=Ymax_norm), alpha=0.5) + ggplot2::theme_bw() +
        ggplot2::scale_discrete_manual(aesthetics=c("col"), values=col_pal_col, name="Variable") +
        ggplot2::scale_discrete_manual(aesthetics=c("fill"), values=col_pal_fil, name="Variable") +
        ggplot2::theme(text=ggplot2::element_text(size=15, color="white"), plot.background=ggplot2::element_rect(fill="black", color="black"),
                       panel.border=ggplot2::element_rect(color="white"),
                       panel.background=ggplot2::element_rect(fill="black"), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(), axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"), axis.ticks = ggplot2::element_line(colour = "white"),
                       legend.key = ggplot2::element_rect(fill = "black"), legend.background = ggplot2::element_rect(fill="black")) +
        ggplot2::scale_x_continuous(expand=c(0,0), name="Latent Path (%)", breaks=seq(0, 100, by=1), labels=tick_labels) +
        ggplot2::scale_y_continuous(expand=c(0,0), name="Local Representation (%)")
      plot_B <- ggplot2::ggplot(out, ggplot2::aes(x=X, y=Enr)) +
        ggplot2::geom_line(linewidth=1, color="white") + ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept=0, linetype="dashed", color="white") +
        ggplot2::theme(text=ggplot2::element_text(size=15, color="white"), plot.background=ggplot2::element_rect(fill="black", color="black"),
                       panel.border=ggplot2::element_rect(color="white"),
                       panel.background=ggplot2::element_rect(fill="black"), panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(), axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"), axis.ticks = ggplot2::element_line(colour = "white"),
                       legend.key = ggplot2::element_rect(fill = "black"), legend.background = ggplot2::element_rect(fill="black")) +
        ggplot2::ggtitle(ttl) +
        ggplot2::scale_x_continuous(expand=c(0,0), name=NULL, breaks=seq(0, 100, by=1), labels=NULL) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1), name="Enrichment")
      print(cowplot::plot_grid(plotlist=list(plot_B, plot_A), ncol=1, align="v", axis="lr", rel_heights=c(1, 3)))
    }
  }
  return(out)
}

#' Map a Pattern back to an Image
#'
#' Map features of a covariation plot onto the original image. The covariation plot is partitioned into
#' regions, and the census is used to determine which region each pixel in the image belongs to. New colors are
#' assigned to these regions, resulting in a new object image.
#'
#' @description Use manually-defined gates on a covariation plot to generate an image of custom-defined objects
#' @param covar_data     Data frame: covariation plot data
#' @param region_bounds  List: non-overlapping 2-entry vectors of gates on the covariation plot
#' @param img            Named List: 4d arrays of the image/s on which the covariation data is based
#' @param census         Data frame: census of the image/s from which the covariation plot was made
#' @param radius         Numeric: radius at which the covariation plot was generated
#' @param radii          Named List: for each radius in microns, named vectors of the appropriate radii in pixels in X, Y, and Z
#' @param col_pal        Vector: color palette for newly defined objects
#' @param orig_bkgd      Logical: when mapping, whether to impose the background from the original image
#' @return      List: new object image in array format and its corresponding color palette
#' @export
map_pattern <- function(covar_data, region_bounds, img, census, radius, radii, col_pal, orig_bkgd=T) {
  print("Checking boundaries of specified objects.")
  rgn_bnd_vec <- unlist(region_bounds)
  if (!(0 %in% rgn_bnd_vec)) {                                                  # If zero isn't a lower bound of a region
    rgn_bnd_vec <- c(0, rgn_bnd_vec)                                            # make it one
  }
  if (!(100 %in% rgn_bnd_vec)) {                                                # If 100 isn't an upper bound of a region
    rgn_bnd_vec <- c(rgn_bnd_vec, 100)                                          # make it one
  }
  missing_bounds <- rep(NA, length(rgn_bnd_vec)-2)                              # If there are gaps between defined regions
  for (i in 1:(length(rgn_bnd_vec)-2)) {                                        # fill them with stop-gap regions
    if (sum(rgn_bnd_vec == rgn_bnd_vec[1+i]) < 2) {
      missing_bounds[i] <- rgn_bnd_vec[i+1]
    }
  }
  missing_bounds <- missing_bounds[!is.na(missing_bounds)]
  rgn_bnd_vec <- sort(c(rgn_bnd_vec, missing_bounds))
  new_region_bounds <- unname(split(rgn_bnd_vec, rep(1:(length(rgn_bnd_vec)/2), each=2)))
  region_IDs <- rep(NA, length(new_region_bounds))                              # Distinguish originally defined regions
  region_IDs[new_region_bounds %in% region_bounds] <- c(1:length(region_bounds))# vs. catch-all label for stop-gap regions
  region_IDs[is.na(region_IDs)] <- length(region_bounds) + 1
  print("Assigning neighborhoods to specified objects.")
  var_names <- unique(covar_data$V)
  census <- census[census$Radius == radius,
                   colnames(census) %in% c(var_names,"X","Y","Z")]
  covar_data <- subset(covar_data, select=c(X,Y,V))
  covar_data <- tidyr::pivot_wider(covar_data, id_cols=X, values_from=Y, names_from=V)
  nbhd_position <- FNN::knnx.index(as.matrix(covar_data[,colnames(covar_data) %in% var_names]),
                                   as.matrix(census[,colnames(census) %in% var_names]), k=1)
  nbhd_conf <- FNN::knnx.dist(as.matrix(covar_data[,colnames(covar_data) %in% var_names]),
                              as.matrix(census[,colnames(census) %in% var_names]), k=1)
  nbhd_conf <- 1 - (nbhd_conf/(max(nbhd_conf)) + 0.0001)
  nbhd_position <- covar_data$X[nbhd_position]                                  # Find which comprehensive region each
  lower_bounds <- unlist(lapply(new_region_bounds, '[', 1))                     # neighborhood belongs to, then get the
  nbhd_region <- unlist(lapply(nbhd_position, function(x,y) {sum(x >= y)}, y=lower_bounds))
  nbhd_region <- region_IDs[nbhd_region]                                        # matching original region or the catch-all
  print("Voting to partition image into specified objects.")
  img <- abind::abind(img, along=4)                                             # Put all image layers together
  obj_img <- array(0, dim=c(dim(img)[1:3], max(region_IDs)+1))
  pix_rad <- radii[[which(as.numeric(names(radii)) == radius)]]
  for (i in c(1:nrow(census))) {                                                # For each nbhd, add votes for each region
    inc <- expand.grid(c(max(census$X[i]-pix_rad[1], 1):min(census$X[i]+pix_rad[1], dim(img)[1])),
                       c(max(census$Y[i]-pix_rad[2], 1):min(census$Y[i]+pix_rad[2], dim(img)[2])),
                       c(max(census$Z[i]-pix_rad[3], 1):min(census$Z[i]+pix_rad[3], dim(img)[3])))
    inc <- t((t(inc) - census[i,c("X","Y","Z")]) / pix_rad)
    inc <- as.integer(sqrt(rowSums((inc)^2)) <= 1)
    inc <- array(inc, dim = c(1 + min((census$X[i]+pix_rad[1]), dim(img)[1]) - max((census$X[i]-pix_rad[1]), 1),
                              1 + min((census$Y[i]+pix_rad[2]), dim(img)[2]) - max((census$Y[i]-pix_rad[2]), 1),
                              1 + min((census$Z[i]+pix_rad[3]), dim(img)[3]) - max((census$Z[i]-pix_rad[3]), 1), 1))
    obj_img[c(max(1,(census$X[i]-pix_rad[1])):min(dim(img)[1],(census$X[i]+pix_rad[1]))),
            c(max(1,(census$Y[i]-pix_rad[2])):min(dim(img)[2],(census$Y[i]+pix_rad[2]))),
            c(max(1,(census$Z[i]-pix_rad[3])):min(dim(img)[3],(census$Z[i]+pix_rad[3]))),
            nbhd_region[i]] <- obj_img[c(max(1,(census$X[i]-pix_rad[1])):min(dim(img)[1],(census$X[i]+pix_rad[1]))),
                                       c(max(1,(census$Y[i]-pix_rad[2])):min(dim(img)[2],(census$Y[i]+pix_rad[2]))),
                                       c(max(1,(census$Z[i]-pix_rad[3])):min(dim(img)[3],(census$Z[i]+pix_rad[3]))),
                                       nbhd_region[i], drop=F] + (inc*nbhd_conf[i])
  } # Each XYZ of obj_img is a pix; each layer is a region; each value is a vote count that the pixel belongs to the region
  obj_img[,,,dim(obj_img)[4]] <- 0.000001                                       # Tiny vote for bkgd in case no votes
  print("Tallying votes for final object partitioning.")
  obj_img <- apply(obj_img, c(1,2,3), which.max)                                # Which region wins vote at each pixel
  obj_img <- array(obj_img, dim=c(dim(obj_img), 1))                             # Make obj_img a 3D array with 3rd dim =1
  print("Imposing areas with no votes as background, as well as original background if requested.")
  obj_img[obj_img == (max(region_IDs)+1)] <- 0                                  # Pixels with no votes are bkgd
  if (orig_bkgd) {
    img <- array(apply(img, c(1,2,3), sum), dim=c(dim(img)[1:3],1))             # Find original content vs. bkgd
    obj_img <- obj_img * (img > 0)                                              # Reimpose original background
  }
  if (missing(col_pal)) {                                                       # If a color palette is not provided
    col_pal <- make_palette(length(region_bounds))                              # Make one for all specified zones
  }
  if (any(obj_img == (length(region_bounds) + 1))) {                            # If there any catch-all pixels in the
    col_pal <- c(col_pal, "#808080")                                            # obj map, add gray to the color palette
  }
  if (length(unique(col_pal)) < length(unique(as.vector(obj_img)))) {           # If there are more objects than colors
    for (i in 1:length(unique(col_pal))) {                                      # For each color in the palette
      which_cols <- which(col_pal == unique(col_pal)[i])                        # Which color is it
      if (length(which_cols) > 1) {                                             # If this color occurs more than once
        obj_img[obj_img %in% (which_cols - any(obj_img == 0))] <- min(which_cols) - any(obj_img == 0)
      }                                                                         # Reset all object codes of this color
    }                                                                           # to the smallest one
    col_pal <- unique(col_pal)                                                  # Reduce color palette to unique colors
  }
  print("Plotting new object image.")
  plot_palette(col_pal, any(obj_img == 0), "Object")                            # Plot the object palette
  plot_image(obj_img, "O", col_pal)                                             # Plot the objects array
  return(list(obj_img, col_pal))                                                # Return object map & color palette,
}
