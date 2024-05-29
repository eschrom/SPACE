#' Plot a Profile Table
#'
#' Plot a heatmap of object compositional profiles from a profile table. Column 1 must be labeled "Object" and give
#' numeric object IDs. Column 2 must be labeled "Count" and give the count of pixels or items belonging to each object.
#' The remaining columns give the compositional profile of scalar expression for each object.
#' Normalizing across objects is recommended: scalars are relative to their individual min vs. max levels across all objects.
#' Normalizing within objects is also possible: scalars are relative to their collective min vs. max within individual objects.
#' Normalizing uniformly is recommended, but normalizing by Z-score is also possible.
#' Bar plots show the number of pixels or items belonging to each object, and the maximum intensity for each scalar,
#' i.e. what quantitative value bright yellow indicates. For the latter plot, higher bars are better.
#' Low bars indicate that the scalar has very low and/or sparse intensity, even among the highest expressing objects.
#'
#' @description Plot a heatmap to show the compositions of each object in terms of scalars
#' @param prof_table  Data frame: profile table with object IDs, pixel/object count, and composition of each object
#' @param compare     Character: whether to normalize "A"cross objects, "W"ithin objects, or "B"oth
#' @param normalize   Character: whether to normalize "U"niformly or by "Z" score
#' @param tile_plots  Logical: whether plots should be tile or separate
#' @param plot_bkgd   Character: whether the plot should have a "W"hite or "B"lack background
#' @return            Nothing: a plot of the heatmap is generated
#' @export
plot_table <- function(prof_table, compare = "A", normalize = "U", tile_plots = F, plot_bkgd = "W") {
  prof_table_temp <- prof_table                                                 # Copy the profile table to reshape                                # Remove pixel or object counts
  prof_table_temp$Count <- NULL                                                 # Temporarily remove object or pixel count
  num_col <- dim(prof_table_temp)[2]                                            # Number of columns
  if (compare == "A") {                                                         # If comparing across object types within scalars
    if (normalize == "U") {                                                     # If using uniform normalization
      prof_table_temp[ ,2:num_col] <- scale(prof_table_temp[ ,2:num_col],
                                            center=apply(prof_table_temp[ ,2:num_col], 2, min),
                                            scale = FALSE)
      prof_table_temp[ ,2:num_col] <- scale(prof_table_temp[ ,2:num_col],
                                            center = FALSE,
                                            scale = apply(prof_table_temp[ ,2:num_col], 2, max))
      ttl <- "Uniform Normalization Across Objects"
      brk <- c(0,1)
      lbl <- c("Min", "Max")
    } else if (normalize == "Z") {                                              # If using Z-score normalization
      prof_table_temp[ ,2:num_col] <- scale(prof_table_temp[ ,2:num_col],
                                            center=T, scale=T)
      ttl <- "Z-Score Normalization Across Objects"
      brk <- seq(plyr::round_any(min(prof_table_temp[,2:num_col]), 2),
                 plyr::round_any(max(prof_table_temp[,2:num_col]), 2), 2)
      lbl <- brk
    }
  }
  if (compare == "W") {                                                         # If comparing within objects across scalars
    if (normalize == "U") {                                                     # If using uniform normalization
      prof_table_temp[ ,2:num_col] <- t(scale(t(prof_table_temp[ ,2:num_col]),
                                              center=apply(t(prof_table_temp[ ,2:num_col]), 2, min),
                                              scale = FALSE))
      prof_table_temp[ ,2:num_col] <- t(scale(t(prof_table_temp[ ,2:num_col]),
                                              center = FALSE,
                                              scale = apply(t(prof_table_temp[ ,2:num_col]), 2, max)))
      prof_table_temp[is.na(prof_table_temp)] <- 0.5
      ttl <- "Uniform Normalization Within Objects"
      brk <- c(0,1)
      lbl <- c("Min", "Max")
    } else if (normalize == "Z") {                                              # If using Z-score normalization
      prof_table_temp[ ,2:num_col] <- t(scale(t(prof_table_temp[ ,2:num_col]),
                                              center=T, scale=T))
      prof_table_temp[is.na(prof_table_temp)] <- 0
      ttl <- "Z-Score Normalization Within Objects"
      brk <- seq(plyr::round_any(min(prof_table_temp[ ,2:num_col]), 2),
                 plyr::round_any(max(prof_table_temp[ ,2:num_col]), 2), 2)
      lbl <- brk
    }
  }
  if (compare == "B") {                                                         # If comparing across the entire heatmap at once
    if (normalize == "U") {                                                     # If using uniform normalization
      prof_table_temp[ ,2:num_col] <- prof_table_temp[ ,2:num_col] - min(prof_table_temp[ ,2:num_col])
      prof_table_temp[ ,2:num_col] <- prof_table_temp[ ,2:num_col] / max(prof_table_temp[ ,2:num_col])
      ttl <- "Uniform Normalization Across & Within Objects"
      brk <- c(0,1)
      lbl <- c("Min", "Max")
    } else if (normalize == "Z") {                                              # If using Z-score normalization
      prof_table_temp[ ,2:num_col] <- (prof_table_temp[ ,2:num_col] - mean(unlist(prof_table_temp[ ,2:num_col]))) /
        sd(unlist(prof_table_temp[ ,2:num_col]))
      ttl <- "Z-Score Normalization Across & Within Objects"
      brk <- seq(plyr::round_any(min(prof_table_temp[ ,2:num_col]), 2),
                 plyr::round_any(max(prof_table_temp[ ,2:num_col]), 2), 2)
      lbl <- brk
    }
  }
  prof_table_temp <- reshape2::melt(prof_table_temp, id.vars=c("Object"))
  colnames(prof_table_temp) <- c("Object", "Component", "Amount")
  comp_max <- apply(prof_table[,3:ncol(prof_table)], 2, max)
  comp_max <- data.frame(Component = names(comp_max), max_val = unname(comp_max))
  if (tile_plots) {
    if (plot_bkgd == "W") {
      plB <- ggplot2::ggplot(data = prof_table_temp, ggplot2::aes(x = Object, y = Component, fill = Amount)) +
        ggplot2::geom_tile() +
        viridis::scale_fill_viridis(discrete=F, breaks=brk, labels=lbl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA)) +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table)))),
                                    expand=c(0,0))
    } else if (plot_bkgd == "B") {
      plB <- ggplot2::ggplot(data = prof_table_temp, ggplot2::aes(x = Object, y = Component, fill = Amount)) +
        ggplot2::geom_tile() +
        viridis::scale_fill_viridis(discrete=F, breaks=brk, labels=lbl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.title=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"),
                       axis.line = ggplot2::element_line(colour = NA),
                       legend.key=ggplot2::element_rect(fill="black"),
                       legend.background=ggplot2::element_rect(fill="black")) +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table)))),
                                    expand=c(0,0))
    }
    if (plot_bkgd == "W") {
      plC <- ggplot2::ggplot(comp_max, ggplot2::aes(x = 1:nrow(comp_max), y = max_val)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="black") + ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_continuous(expand=c(0,0)) + ggplot2::ylab("Max Amount") + ggplot2::coord_flip() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA), axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())
    } else if (plot_bkgd == "B") {
      plC <- ggplot2::ggplot(comp_max, ggplot2::aes(x = 1:nrow(comp_max), y = max_val)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="white") + ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_continuous(expand=c(0,0)) + ggplot2::ylab("Max Amount") + ggplot2::coord_flip() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.line = ggplot2::element_line(colour = NA),
                       axis.text.x=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank())
    }
    if (plot_bkgd == "W") {
      plA <- ggplot2::ggplot(prof_table, ggplot2::aes(x = Object, y = Count)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="black") + ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_continuous(expand=c(0,0)) + ggplot2::ylab("Count") + ggplot2::ggtitle(ttl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA), axis.title.x=ggplot2::element_blank())
    } else if (plot_bkgd == "B") {
      plA <- ggplot2::ggplot(prof_table, ggplot2::aes(x = Object, y = Count)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="white") + ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_continuous(expand=c(0,0)) + ggplot2::ylab("Count") + ggplot2::ggtitle(ttl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.title=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_text(color="white"),
                       axis.line = ggplot2::element_line(colour = NA),
                       axis.title.x=ggplot2::element_blank())
    }
    print((plA + patchwork::plot_spacer() + plB + plC +                         # Display all plots together
             patchwork::plot_layout(widths=c(4,1), heights=c(1,4), guides="collect",
                                    design="12
                                                  34")) &
            ggplot2::theme(plot.background = ggplot2::element_rect(fill = ifelse(plot_bkgd == "B", "black", "white"),
                                                                   color = ifelse(plot_bkgd == "B", "black", "white"))))
  } else {                                                                      # If separate, not tiled, plots are desired
    if (plot_bkgd == "W") {
      plB <- ggplot2::ggplot(data = prof_table_temp, ggplot2::aes(x = Object, y = Component, fill = Amount)) +
        ggplot2::geom_tile() +
        viridis::scale_fill_viridis(discrete=F, breaks=brk, labels=lbl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA)) + ggplot2::ggtitle(ttl) +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table)))),
                                    expand=c(0,0))
    } else if (plot_bkgd == "B") {
      plB <- ggplot2::ggplot(data = prof_table_temp, ggplot2::aes(x = Object, y = Component, fill = Amount)) +
        ggplot2::geom_tile() +
        viridis::scale_fill_viridis(discrete=F, breaks=brk, labels=lbl) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.title=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"),
                       axis.line = ggplot2::element_line(colour = NA),
                       legend.key=ggplot2::element_rect(fill="black"),
                       legend.background=ggplot2::element_rect(fill="black")) + ggplot2::ggtitle(ttl) +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table)))),
                                    expand=c(0,0))
    }
    if (plot_bkgd == "W") {
      plC <- ggplot2::ggplot(comp_max, ggplot2::aes(x = 1:nrow(comp_max), y = max_val)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="black") + ggplot2::ylab("Max Amount") +
        ggplot2::xlab("Component") +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(comp_max))), labels=colnames(prof_table)[-c(1:2)]) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA))
    } else if (plot_bkgd == "B") {
      plC <- ggplot2::ggplot(comp_max, ggplot2::aes(x = 1:nrow(comp_max), y = max_val)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="white") + ggplot2::ylab("Max Amount") +
        ggplot2::xlab("Component") +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(comp_max))), labels=colnames(prof_table)[-c(1:2)]) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.line = ggplot2::element_line(colour = NA),
                       axis.text.x=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_text(color="white"))
    }
    if (plot_bkgd == "W") {
      plA <- ggplot2::ggplot(prof_table, ggplot2::aes(x = Object, y = Count)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="black") + ggplot2::ylab("Count") +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table))))) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15), panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = NA))
    } else if (plot_bkgd == "B") {
      plA <- ggplot2::ggplot(prof_table, ggplot2::aes(x = Object, y = Count)) +
        ggplot2::geom_bar(stat="identity", color="white", fill="white") + ggplot2::ylab("Count") +
        ggplot2::scale_x_continuous(breaks=c(1:(nrow(prof_table))), labels=as.character(c(1:(nrow(prof_table))))) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       text=ggplot2::element_text(size=15, color="white"),
                       plot.background=ggplot2::element_rect(fill="black"),
                       panel.background=ggplot2::element_rect(fill="black"),
                       axis.title=ggplot2::element_text(color="white"),
                       axis.text.y=ggplot2::element_text(color="white"),
                       axis.text.x=ggplot2::element_text(color="white"),
                       axis.line = ggplot2::element_line(colour = NA))
    }
    print(plB)
    print(plA)
    print(plC)
  }
}

#' Merge Objects Together
#'
#' Merge similar objects together in a profile table, image, color palette, and/or object table that all correspond.
#'
#' @description Merge objects together in a profile table, image, color palette, and/or object table
#' @param prof_table   Data frame: profile table containing the objects to be merged
#' @param img          3d Array: image containing the objects to be merged
#' @param col_pal      Vector: Color palette for the objects including those to be merged
#' @param obj_table    Data frame: object table including the objects to be merged
#' @param obj_groups   List: vectors of numeric object IDs, grouping together the objects to be merged
#' @return      List: profile table, object image, color palette, and/or object table, but with objects merged
#' @export
merge_objects <- function(prof_table, img, col_pal, obj_table, obj_groups) {
  if (!is.list(obj_groups)) {                                                   # If only one set of objects to merge
    obj_groups <- list(obj_groups)                                              # still make it a list with 1 element
  }
  out <- vector("list", 4)
  if (!missing(prof_table)) {
    print("Merging objects in the profile table.")
    rem_prof_table <- prof_table[-(unlist(obj_groups)), ]                       # Remaining profiles w/out merged objects
    for (i in 1:length(obj_groups)) {                                           # For each set of objects to merge
      mrg_prof_table <- prof_table[obj_groups[[i]], 3:dim(prof_table)[2]]       # Original profiles for just these objects
      mrg_prof_table <- colSums(mrg_prof_table * prof_table[obj_groups[[i]], 2])# Weight orig profs by count and sum
      mrg_prof_table <- mrg_prof_table / sum(prof_table[obj_groups[[i]], 2])    # Divide by count across all objects in set
      mrg_prof_table <- c(min(obj_groups[[i]]), sum(prof_table[obj_groups[[i]], 2]), mrg_prof_table) # Append object ID & count
      rem_prof_table <- rbind(rem_prof_table, mrg_prof_table)                   # Append merged object to remaining profiles
    }
    rem_prof_table <- rem_prof_table[order(rem_prof_table[,1]), ]               # Reorder by the original object ID
    rem_prof_table[,1] <- 1:nrow(rem_prof_table)                                # Rename remaining objects to 1,2,3,...
    out[[1]] <- rem_prof_table
  }
  if (!missing(img)) {
    print("Merging objects in the image.")
    for (i in 1:length(obj_groups)) {                                           # For each set of objects to merge
      img[img %in% obj_groups[[i]]] <- min(obj_groups[[i]])                     # Reset all included objects to the smallest
    }
    rem_objs <- sort(unique(as.vector(img)))                                    # Object IDs that will not be removed
    rem_objs <- rem_objs[rem_objs != 0]
    for (i in 1:length(rem_objs)) {                                             # Reset remaining object codes to 1,2,3,...
      img[img == rem_objs[i]] <- i
    }
    out[[2]] <- img
  }
  if (!missing(col_pal)) {
    print("Merging objects in the color palette.")
    del_objs <- unlist(lapply(obj_groups, function(x) {x[!(x == min(x))]}))     # Objects to be deleted
    col_pal <- col_pal[-del_objs]
    out[[3]] <- col_pal
  }
  if (!missing(obj_table)) {
    print("Merging objects in the object table.")
    for (i in 1:length(obj_groups)) {                                           # For each set of object to merge
      obj_table$Object[obj_table$Object %in% obj_groups[[i]]] <- min(obj_groups[[i]]) # Reset all included objects to the smallest
    }
    rem_objs <- sort(unique(obj_table$Object))                                  # Object IDs that will not be removed
    for (i in 1:length(rem_objs)) {                                             # Reset remaining object codes to 1,2,3,...
      obj_table$Object[obj_table$Object == rem_objs[i]] <- i
    }
    out[[4]] <- obj_table
  }
  out <- plyr::compact(out)
  return(out)
}
