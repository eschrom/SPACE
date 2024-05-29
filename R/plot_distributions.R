#' Plot a Distribution
#'
#' Plot a distribution of co/occurrence for one or two variables. Omitting a number of bins for rounding to obtain
#' a 1d density or 2d scatter plot, reflecting how the data exists in the census. Provide a number of bins for
#' rounding to obtain a probability distribution across bins, reflecting how the data is used for mutual information.
#' Omit a patch list to obtain the plot for the observed census. Provide a patch list to obtain the plot for a
#' randomized census. More than two variables and a patch list can be provided if only a randomized census is desired,
#' and no plot will be returned.
#'
#' @description Plot how often each quantitative combination of 1 or 2 variables occur in a specimen
#' @param census     Data frame: census of the specimen
#' @param ensemble   Vector: 1 or 2 variables for which to plot the distribution
#' @param radius     Numeric: radius for which to plot the distribution
#' @param bin_num    Numeric: number of rounding bins per variable; if omitted, a scatter plot is given
#' @param patch_list List: if randomizing, patch list for the corresponding image or table
#' @param plot_zoom  Logical: whether to zoom into the distribution, or keep a 0-100% range for all variables
#' @param plot_bkgd  Character: whether the plot background should be "B"lack or "W"hite
#' @return      Data frame: the joint distribution of the ensemble, and a plot
#' @export
plot_dist <- function(census, ensemble, radius, bin_num=NULL, patch_list=NULL, plot_zoom=F, plot_bkgd="W") {
  if (is.integer(ensemble)) {                                                   # If variables are specified by number,
    ensemble <- colnames(census)[ensemble]                                      # extract the names that go with them
  }                                                                             # Subset census to just correct radius
  census <- census[census$Radius == radius[1], which(stringr::str_count(colnames(census),"O|S") >= 1)]
  if (!is.null(patch_list)) {                                                   # If random nbhds desired, construct them
    plot_title <- "Randomized"                                                  # Use patches from bigger nbhds to count eligible pix
    sub_ptch_lst <- patch_list[[which(as.integer(names(patch_list)) == radius[1])]]
    os_pairs <- colnames(census)[which(stringr::str_count(colnames(census), "O|S") == 2)]
    if (length(os_pairs) > 0) {                                                 # Infer the object-scalar pairings
      obj_maps <- unique(gsub("\\..*","", os_pairs))                            # Object maps with attached scalars
      OS_pairs <- vector("list", length=length(obj_maps))                       # List for pairing matrices
      names(OS_pairs) <- obj_maps                                               # Name matrices by their obj map
      for (i in 1:length(OS_pairs)) {                                           # For each object map represented
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
    } else {
      OS_pairs <- NULL
    }
    census <- random_census(sub_ptch_lst, OS_pairs)
  } else {
    plot_title <- "Observed"
  }
  if (length(ensemble) > 2) {
    print("Distributions can only be plotted for 1 or 2 variables, but the requested census is returned.")
    return(census)
  }
  var_cols <- which(colnames(census) %in% ensemble)                             # Columns with variables in the ensemble
  col_mins <- dplyr::summarize(census, dplyr::across(all_of(var_cols), min))    # Min & max % for var
  col_maxs <- dplyr::summarize(census, dplyr::across(all_of(var_cols), max))
  if ((length(ensemble) > 1) && (which(colnames(census) == ensemble[1]) > which(colnames(census) == ensemble[2]))) {
    var_cols <- rev(var_cols)                                                   # If a 2-variable ensemble is out of order
    col_mins <- rev(col_mins)                                                   # (later variable comes first), must
    col_maxs <- rev(col_maxs)                                                   # reverse the order of these manually
  }
  if (!is.null(bin_num)) {                                                      # If rounding is desired
    adj_bin_num <- rep(NA, length(var_cols))                                    # Bin number may have to be adjusted
    for (i in 1:length(var_cols)) {                                             # For each ensemble variable, only round
      if (length(unique(census[,var_cols[i]])) > bin_num) {                     # if there are more unique values than bins
        census[,var_cols[i]] <- round_column(as.vector(census[,var_cols[i]]), as.matrix(col_mins)[,i], as.matrix(col_maxs)[,i], bin_num, T)
        adj_bin_num[i] <- bin_num                                               # In this case, keep original bin number
      } else {                                                                  # Otherwise, there's no need to round
        adj_bin_num[i] <- length(unique(census[,var_cols[i]]))                  # Record the smaller true bin number
      }
    }
  }
  if (plot_zoom) {                                                              # Axis limits depend on whether zooming
    plot_lo <- unlist(col_mins)
    plot_hi <- unlist(col_maxs)
  } else {
    plot_lo <- rep(0, length(ensemble))
    plot_hi <- rep(100, length(ensemble))
  }
  if (length(ensemble) == 1) {
    axis_label <- stringr::str_count(ensemble, "S")
    if (axis_label == 1) {
      axis_label <- "Saturation"
    } else {
      axis_label <- "Pixels"
    }
    if (!is.null(bin_num)) {
      marg_dist <- build_dist(census, ensemble)
      mins_and_maxs <- as.data.frame(rbind(col_mins, col_maxs))
      marg_dist <- smooth_dist(marg_dist, adj_bin_num, mins_and_maxs, full_dist=T)
      if (plot_bkgd == "W") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_bar(ggplot2::aes(x=marg_dist[,1], y=marg_dist[,2]), stat="identity", fill="gray67", color="black") +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15)) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble, sep="")) +
          ggplot2::scale_x_continuous(expand=c(0,0), breaks=marg_dist[,1], labels=round(marg_dist[,1]),
                                      limits=c(plot_lo - 100/bin_num, plot_hi + 100/bin_num)) +
          ggplot2::scale_y_continuous(expand=c(0,0), limits=c(0,ifelse(plot_zoom, max(marg_dist[,2]), 1))) +
          ggplot2::ylab("Frequency") + ggplot2::xlab(paste("% of ", axis_label, sep=""))
      } else if (plot_bkgd == "B") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_bar(ggplot2::aes(x=marg_dist[,1], y=marg_dist[,2]), stat="identity", fill="gray33", color="white") +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15, color="white"),
                         axis.text.y=ggplot2::element_text(color="white"), axis.text.x=ggplot2::element_text(color="white"),
                         plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black")) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble, sep="")) +
          ggplot2::scale_x_continuous(expand=c(0,0), breaks=marg_dist[,1], labels=round(marg_dist[,1]),
                                      limits=c(plot_lo - 100/bin_num, plot_hi + 100/bin_num)) +
          ggplot2::scale_y_continuous(expand=c(0,0), limits=c(0,ifelse(plot_zoom, max(marg_dist[,2]), 1))) +
          ggplot2::ylab("Frequency") + ggplot2::xlab(paste("% of ", axis_label, sep=""))
      }
    } else {
      col_id <- which(colnames(census) == ensemble)
      marg_dist <- census[,col_id]
      if (plot_bkgd == "W"){
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_density(ggplot2::aes(x=marg_dist), color="black", linewidth=2.5) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble, sep="")) +
          ggplot2::scale_x_continuous(expand=c(0.01,0.01), limits=c(plot_lo - 1, plot_hi + 1)) +
          ggplot2::scale_y_continuous(expand=c(0,0)) + ggplot2::theme_bw() +
          ggplot2::ylab("Probability") + ggplot2::xlab(paste("% of ", colnames(census)[col_id], " ", axis_label, sep="")) +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15))
      } else if (plot_bkgd == "B") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_density(ggplot2::aes(x=marg_dist), color="white", linewidth=2.5) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble, sep="")) +
          ggplot2::scale_x_continuous(expand=c(0.01,0.01), limits=c(plot_lo - 1, plot_hi + 1)) +
          ggplot2::scale_y_continuous(expand=c(0,0)) +
          ggplot2::ylab("Probability") + ggplot2::xlab(paste("% of ", colnames(census)[col_id], " ", axis_label, sep="")) +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15, color="white"), axis.text.y=ggplot2::element_text(color="white"),
                         axis.text.x=ggplot2::element_text(color="white"),
                         plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black"))
      }
    }
    print(plot_A)
    return(marg_dist)
  } else {
    axis_label <- stringr::str_count(ensemble, "S")
    for (i in 1:length(axis_label)) {
      if (axis_label[i] == 1) {
        axis_label[i] <- "Saturation"
      } else {
        axis_label[i] <- "Pixels"
      }
    }
    if (!is.null(bin_num)) {
      joint_dist <- build_dist(census, ensemble)
      mins_and_maxs <- as.data.frame(rbind(col_mins, col_maxs))
      joint_dist <- smooth_dist(joint_dist, adj_bin_num, mins_and_maxs, full_dist=T)
      joint_dist$log_freq <- log10(joint_dist$freq)                             # Calculate log probabilities for plotting
      if (plot_bkgd == "W") {
        plot_A <- ggplot2::ggplot() +                                           # White = min probability, black = max probability
          ggplot2::geom_tile(data=joint_dist, ggplot2::aes(x=joint_dist[,1], y=joint_dist[,2], fill=joint_dist[,4])) +
          ggplot2::ggtitle(paste(plot_title, " Joint Distribution of ", ensemble[1], " & ", ensemble[2], sep="")) +
          ggplot2::scale_fill_gradient(low="white", high="black") +
          ggplot2::scale_x_continuous(expand=c(0,0), limits=c(plot_lo[1] - 100/bin_num, plot_hi[1] + 100/bin_num)) +
          ggplot2::scale_y_continuous(expand=c(0,0), limits=c(plot_lo[2] - 100/bin_num, plot_hi[2] + 100/bin_num)) +
          ggplot2::theme_bw() +
          ggplot2::xlab(paste("% of ", colnames(joint_dist)[1], " ", axis_label[1], sep="")) +
          ggplot2::ylab(paste("% of ", colnames(joint_dist)[2], " ", axis_label[2], sep="")) +
          ggplot2::labs(fill="Log Prob") + ggplot2::theme(text=ggplot2::element_text(size=15))
      } else if (plot_bkgd == "B") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_tile(data=joint_dist, ggplot2::aes(x=joint_dist[,1], y=joint_dist[,2], fill=joint_dist[,4])) +
          ggplot2::ggtitle(paste(plot_title, " Joint Distribution of ", ensemble[1], " & ", ensemble[2], sep="")) +
          ggplot2::scale_fill_gradient(low="black", high="white") +
          ggplot2::scale_x_continuous(expand=c(0,0), limits=c(plot_lo[1] - 100/bin_num, plot_hi[1] + 100/bin_num)) +
          ggplot2::scale_y_continuous(expand=c(0,0), limits=c(plot_lo[2] - 100/bin_num, plot_hi[2] + 100/bin_num)) +
          ggplot2::xlab(paste("% of ", colnames(joint_dist)[1], " ", axis_label[1], sep="")) +
          ggplot2::ylab(paste("% of ", colnames(joint_dist)[2], " ", axis_label[2], sep="")) +
          ggplot2::labs(fill="Log Prob") +
          ggplot2::theme(text=ggplot2::element_text(size=15, color="white"), axis.text.y=ggplot2::element_text(color="white"),
                         axis.text.x=ggplot2::element_text(color="white"), legend.background=ggplot2::element_rect(fill="black"),
                         legend.text=ggplot2::element_text(color="white"), panel.grid.major=ggplot2::element_line(color="gray15"),
                         panel.grid.minor=ggplot2::element_line(color="gray15"),
                         plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black"))
      }
    } else {
      col_id <- which(colnames(census) %in% ensemble)
      joint_dist <- census[,col_id]
      if (plot_bkgd == "W") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x=census[,col_id[1]], y=census[,col_id[2]]),
                              color="black", alpha=0.25, size=2.5) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble[1], " & ", ensemble[2], sep="")) +
          ggplot2::ylab(paste("% of ", colnames(census)[col_id[2]], " ", axis_label[2], sep="")) +
          ggplot2::xlab(paste("% of ", colnames(census)[col_id[1]], " ", axis_label[1], sep="")) +
          ggplot2::scale_x_continuous(expand=c(0.01, 0.01), limits=c(plot_lo[1], plot_hi[1])) +
          ggplot2::scale_y_continuous(expand=c(0.01, 0.01), limits=c(plot_lo[2], plot_hi[2])) + ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15))
      } else if (plot_bkgd == "B") {
        plot_A <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x=census[,col_id[1]], y=census[,col_id[2]]),
                              color="white", alpha=0.25, size=2.5) +
          ggplot2::ggtitle(paste(plot_title, " Distribution of ", ensemble[1], " & ", ensemble[2], sep="")) +
          ggplot2::ylab(paste("% of ", colnames(census)[col_id[2]], " ", axis_label[2], sep="")) +
          ggplot2::xlab(paste("% of ", colnames(census)[col_id[1]], " ", axis_label[1], sep="")) +
          ggplot2::scale_x_continuous(expand=c(0.01, 0.01), limits=c(plot_lo[1], plot_hi[1])) +
          ggplot2::scale_y_continuous(expand=c(0.01, 0.01), limits=c(plot_lo[2], plot_hi[2])) + ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                         text=ggplot2::element_text(size=15, color="white"), axis.text.y=ggplot2::element_text(color="white"),
                         axis.text.x=ggplot2::element_text(color="white"),
                         plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black"))
      }
    }
    print(plot_A)
    if (ncol(joint_dist == 4)) {
      return(joint_dist[,-4])
    } else {
      return(joint_dist)
    }
  }
}
