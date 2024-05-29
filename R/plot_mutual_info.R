#' Plot Mutual Information by Rank
#'
#' Compare ensembles by how non-randomly they are patterned, in the case of cisMI, or by
#' how well their patterning distinguishes specimen groups, in the case of transMI.
#'
#' @description Plot cisMI or transMI scores across ensembles, sorted from most to least significant
#' @param mi          List: cis or trans MI measurements across all radii and ensembles
#' @param radius      Numeric: radius for which to plot P values
#' @param depth       Vector: ensemble depths to include in the plot
#' @param col_pals    List: list of the palettes for each object map or scalar set represented in the ensembles
#' @param p_thr       Numeric: threshold corrected p value to consider, defaulting to 0.05
#' @param mi_thr      Numeric: threshold MI value to consider, defaulting to 0.1 for cisMI or 0 for transMI
#' @param p_adj       Logical: whether to use multiple-testing-corrected p values or raw p values
#' @param all         Vector: variables which all must be included in every ensemble of the highest depth
#' @param alo         Vector: variables of which at least one must be included in every ensemble of the highest depth
#' @param not         Vector: variables which should not be included in any ensemble
#' @param group       Character: for transMI, name of the grouping factor to plot
#' @param plot_bkgd   Character: whether the plot background should be "B"lack or "W"hite
#' @return      Data Frame: all significant ensembles and their cis/trans MI scores and a plot of this data
#' @export
plot_MI_rank <- function(mi, radius, depth=NULL, col_pals=NULL, p_thr=NULL, mi_thr=NULL, p_adj=T,
                         all=NULL, alo=NULL, not=NULL, group=NULL, plot_bkgd="W") {
  mi <- mi[[which(names(mi) == as.character(radius))]]
  v_cols <- which(stringr::str_count(colnames(mi), "V") == 1)
  if (is.null(depth)) {
    depth <- 1:length(v_cols)
  } else {
    mi <- mi[rowSums(!is.na(mi[,v_cols])) %in% depth, ]
  }
  max_depth <- max(depth)
  if (is.null(col_pals)) {
    stop("Error: color palette(s) must be provided.")
  }
  if ((length(all) + !is.null(alo)) > min(depth)) {
    stop("Error: the number of required variables exceeds ensemble depth.")
  }
  keep <- rep(T, nrow(mi))
  if (!is.null(all)) {
    keep <- keep & apply(mi[,v_cols,drop=F], 1, function(x,y) {all(y %in% x)}, y=all)
  }
  if (!is.null(alo)) {
    keep <- keep & apply(mi[,v_cols,drop=F], 1, function(x,y) {any(y %in% x)}, y=alo)
  }
  if (!is.null(not)) {
    keep <- keep & !apply(mi[,v_cols,drop=F], 1, function(x,y) {any(y %in% x)}, y=not)
  }
  mi <- mi[keep, ]
  if (is.null(p_thr)) {
    p_thr <- 0.05
  }
  if (is.null (mi_thr)) {
    if (is.null(group)) {
      mi_thr <- 0.1
    } else {
      mi_thr <- 0
    }
  }
  if (is.null(group)) {
    if (p_adj) {
      p_col <- which(stringr::str_count(colnames(mi), "Padjust") == 1)
    } else {
      p_col <- which(stringr::str_count(colnames(mi), "Pvalue") == 1)
    }
  } else {
    if (p_adj) {
      p_col <- which(colnames(mi) == paste("Padjust_", group, sep=""))
    } else {
      p_col <- which(colnames(mi) == paste("Pvalue_", group, sep=""))
    }
  }
  if (length(p_col) != 1) {
    stop("Error: exactly one grouping factor must be specified verbatim.")
  }
  if (is.null(group)) {
    m_col <- which(stringr::str_count(colnames(mi), "MI") == 1)
  } else {
    m_col <- which(colnames(mi) == paste("TransMI_", group, sep=""))
  }
  if (is.null(group)) {
    z_col <- which(stringr::str_count(colnames(mi), "Zscore") == 1)
  } else {
    z_col <- which(colnames(mi) == paste("Zscore_", group, sep=""))
  }
  mi <- mi[mi[,p_col] <= p_thr, ,drop=F]                                        # Drop p values that are too high
  if (stringr::str_count(colnames(mi)[m_col], "Cis") == 1) {                    # If this is cisMI
    mi <- mi[abs(mi[,m_col]) >= mi_thr, ,drop=F]                                # drop abs val that are too low
  } else {                                                                      # If this is transMI
    mi <- mi[mi[m_col] >= mi_thr, ,drop=F]                                      # drop vals that are too low
  }
  if (nrow(mi) == 0) {
    print("There are no significant ensembles.")
    return(NULL)
  }
  mi <- mi[order(abs(mi[,z_col]), decreasing=T), ]                              # Sort by abs Z to break P ties
  mi[,p_col] <- log10(mi[,p_col])                                               # Log transform P values
  if (any(is.infinite(mi[,p_col]))) {                                           # If any P-values are 0
    print("Warning: Exact P values for top ensembles are estimated by quadratic regression.")
    x <- c(abs(mi[!is.infinite(mi[,p_col]), z_col]))                               # Predict their log P from Z
    x2 <- c(x^2)                                                                   # Quadratic regression tends to
    y <- mi[!is.infinite(mi[,p_col]), p_col]
    mod <- lm(y ~ 1 + x + x2)
    new_x <- c(abs(mi[is.infinite(mi[,p_col]), z_col]))
    new_x2 <- c(new_x^2)
    mi[,p_col][is.infinite(mi[,p_col])] <- unname(predict(mod, newdata = data.frame(x = new_x, x2 = new_x2)))
  } else {                                                                      # If no P-values are 0, set
    mod <- NULL                                                                 # mod=NULL to remember this later
  }
  if (!is.null(col_pals) && !is.list(col_pals)) {
    col_pals <- list(col_pals)
  }
  for (i in 1:length(col_pals)) {                                               # CONSIDER REMOVING THIS in favor of
    if (col_pals[[i]][1] %in% c("#000000", "#FFFFFF")) {                        # a standardized convention of not
      col_pals[[i]] <- col_pals[[i]][-1]                                        # including bkgd in color palettes
    }
  }
  mi$size <- rep(NA, nrow(mi))                                                  # Column for point size
  ds <- unname(length(v_cols) - rowSums(is.na(mi[,v_cols,drop=F])))             # Ensemble depths
  for (i in 1:length(unique(ds))) {                                             # For each unique depth
    min_mi <- min(abs(mi[ds==unique(ds)[i], m_col]))                            # get min and max abs mi
    max_mi <- max(abs(mi[ds==unique(ds)[i], m_col]))                            # assign relative point size
    if (max_mi == min_mi) {
      mi$size[ds==unique(ds)[i]] <- 0.5
    } else {
      mi$size[ds==unique(ds)[i]] <- ((abs(mi[ds==unique(ds)[i], m_col]) - min_mi) / (max_mi - min_mi)) + 0.1
    }
  }
  if (stringr::str_count(colnames(mi)[m_col], "Cis") == 1) {                    # Title depends on cis vs. trans
    ttl <- "Ensembles that Differ from Randomized Expectations"
  } else {
    if (is.null(group)) {
      group <- substring(colnames(mi)[p_col], 9, 100)
    }
    ttl <- paste("Ensembles that Distinguish ", group, " Groups", sep="")
  }
  plot_elems <- vector("list", length(col_pals)+1)
  if (plot_bkgd == "W") {
    plot_elems[[1]] <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=c(1:nrow(mi)), y=mi[,p_col]), size=5*mi$size, color="black") +
      ggplot2::geom_hline(yintercept=log10(p_thr), linetype="dotted") +
      ggplot2::ggtitle(ttl) +
      ggplot2::scale_x_continuous(name="Ensemble", breaks=c(1:nrow(mi)), labels=NULL, expand=c(0.01,0.05,0.01,0)) +
      ggplot2::scale_y_reverse(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Corrected Log P Value") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=12), axis.title.y=ggplot2::element_text(size=15),
                     axis.text.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(),
                     plot.title=ggplot2::element_text(size=15))
    for (i in 1:length(col_pals)) {             # New layer of x axis labels for each palette
      out_cols <- mi[,v_cols,drop=F]            # Get MI columns that specify variables in each ensemble
      this_pal <- data.frame(lapply(out_cols, function(x) {grepl(paste(names(col_pals)[i],".",sep=""), x)}))
      if (sum(rowSums(this_pal)) == 0) {        # If no variables from this palette are included in any ensemble
        plot_elems[[i+1]] <- NULL               # Skip this palette and move on
        next
      }
      if (grepl("O", names(col_pals)[i])) {     # If this pal is obj, extract object variable numbers
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(names(col_pals)[i],".",sep=""),'',x)}))
        this_var <- data.frame(lapply(this_var, function(x) {gsub("_.*",'',x)}))
        this_var[as.matrix(!this_pal)] <- NA
      } else {
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(".*",names(col_pals)[i],".",sep=""),'',x)}))
        this_var[as.matrix(data.frame(lapply(this_var, function(x) {grepl("[.]",x)})))] <- NA
        this_var[as.matrix(!this_pal)] <- NA
      }
      this_lnk <- as.matrix(data.frame(lapply(out_cols, function(x) {grepl("_",x)}))) * 1
      this_lnk <- t(apply(this_lnk, 1, cumsum)) * this_lnk
      out_syms <- data.frame(X = rep(1:nrow(mi), times=rowSums(this_pal)),
                             Y = unlist(apply(this_pal, 1, which)),
                             C = as.numeric(t(as.matrix(this_var)))[as.vector(t(as.matrix(this_pal)))])
      xlab <- ifelse(grepl("O", names(col_pals)[i], fixed=T), "Object Map", "Scalar Set")
      plot_elems[[i+1]] <- ggplot2::ggplot(out_syms) +
        ggplot2::geom_point(ggplot2::aes(x=X, y=Y, color=as.factor(C)), shape=15, size=5, show.legend=F) +
        ggplot2::geom_text(ggplot2::aes(x=X, y=Y, label=as.character(C)), size=3, color="white") +
        ggplot2::scale_color_manual(values = col_pals[[i]][sort(unique(out_syms$C))]) + ggplot2::theme_bw() +
        ggplot2::xlab(paste(xlab, names(col_pals)[i])) +
        ggplot2::scale_x_continuous(expand=c(0.01,0.05,0.01,0), limits = c(1,nrow(mi))) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1), limits = c(1,ncol(out_cols))) +
        ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                       panel.background=ggplot2::element_blank(), axis.text.y=ggplot2::element_text(size=13, color="white"),
                       axis.title.y=ggplot2::element_text(size=15, color="white"), axis.ticks.y=ggplot2::element_line(color="white"),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),panel.border=ggplot2::element_blank())
      remove(out_cols, this_pal, this_var, xlab)
    }
    plot_elems <- purrr::compact(plot_elems)
    print(cowplot::plot_grid(plotlist=plot_elems, ncol=1, align="v",
                             rel_heights=c(4, rep(1,length(plot_elems)-1))))
    remove(plot_elems)
  } else if (plot_bkgd == "B") {
    plot_elems[[1]] <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=c(1:nrow(mi)), y=mi[,p_col]), size=5*mi$size, color="white") +
      ggplot2::ggtitle(ttl) +
      ggplot2::geom_hline(yintercept=log10(p_thr), color="white", linetype="dotted") +
      ggplot2::scale_x_continuous(name="Ensemble", breaks=c(1:nrow(mi)), labels=NULL, expand = c(0.01,0.05,0.01,0)) +
      ggplot2::scale_y_reverse(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Corrected Log P Value") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     plot.background=ggplot2::element_rect(fill="black",color="black"), panel.background=ggplot2::element_rect(fill="black",color="black"),
                     axis.text.y=ggplot2::element_text(size=12, color="white"),axis.title.y=ggplot2::element_text(size=15, color="white"),
                     axis.text.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_line(color="white"),
                     panel.border=ggplot2::element_rect(color="white"), axis.ticks.y=ggplot2::element_line(color="white"),
                     plot.title=ggplot2::element_text(size=15, color="white"))
    for (i in 1:length(col_pals)) {             # New layer of x axis labels for each palette
      out_cols <- mi[,v_cols,drop=F]            # Which variables are in the current palette
      this_pal <- data.frame(lapply(out_cols, function(x) {grepl(paste(names(col_pals)[i],".",sep=""), x)}))
      if (sum(rowSums(this_pal)) == 0) {        # If no variables from this palette are included in any ensemble
        plot_elems[[i+1]] <- NULL               # Skip this palette and move on
        next
      }
      if (grepl("O", names(col_pals)[i])) {     # If this pal is obj, extract object variable numbers
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(names(col_pals)[i],".",sep=""),'',x)}))
        this_var <- data.frame(lapply(this_var, function(x) {gsub("_.*",'',x)}))
        this_var[as.matrix(!this_pal)] <- NA
      } else {
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(".*",names(col_pals)[i],".",sep=""),'',x)}))
        this_var[as.matrix(data.frame(lapply(this_var, function(x) {grepl("[.]",x)})))] <- NA
        this_var[as.matrix(!this_pal)] <- NA
      }
      this_lnk <- as.matrix(data.frame(lapply(out_cols, function(x) {grepl("_",x)}))) * 1
      this_lnk <- t(apply(this_lnk, 1, cumsum)) * this_lnk
      out_syms <- data.frame(X = rep(1:nrow(mi), times=rowSums(this_pal)),
                             Y = unlist(apply(this_pal, 1, which)),
                             C = as.numeric(t(as.matrix(this_var)))[as.vector(t(as.matrix(this_pal)))])
      xlab <- ifelse(grepl("O", names(col_pals)[i], fixed=T), "Object Map", "Scalar Set")
      plot_elems[[i+1]] <- ggplot2::ggplot(out_syms) +
        ggplot2::geom_point(ggplot2::aes(x=X, y=Y, color=as.factor(C)), shape=15, size=5, show.legend=F) +
        ggplot2::geom_text(ggplot2::aes(x=X, y=Y, label=as.character(C)), size=3, color="white") +
        ggplot2::scale_color_manual(values = col_pals[[i]][sort(unique(out_syms$C))]) + ggplot2::theme_bw() +
        ggplot2::xlab(paste(xlab, names(col_pals)[i])) +
        ggplot2::scale_x_continuous(expand=c(0.01,0.05,0.01,0), limits = c(1,nrow(mi))) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1), limits = c(1,ncol(out_cols))) +
        ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_text(size=13, color="black"),
                       plot.background=ggplot2::element_rect(fill="black",color="black"), panel.background=ggplot2::element_rect(fill="black",color="black"),
                       axis.title.y=ggplot2::element_text(size=15, color="black"), axis.ticks.y=ggplot2::element_line(color="black"),
                       axis.text.x=ggplot2::element_blank(),axis.title.x=ggplot2::element_text(color="white"),
                       axis.ticks.x=ggplot2::element_blank(),panel.border=ggplot2::element_blank())
      remove(out_cols, this_pal, this_var, xlab)
    }
    plot_elems <- purrr::compact(plot_elems)
    print(cowplot::plot_grid(plotlist=plot_elems, ncol=1, align="v",
                             rel_heights=c(4, rep(1,length(plot_elems)-1))))
    remove(plot_elems)
  }
  uniq_vars <- unique(unlist(mi[,v_cols,drop=F]))
  uniq_vars <- uniq_vars[!is.na(uniq_vars)]
  p_agg <- data.frame(V = uniq_vars, Zscore = rep(0, length(uniq_vars)), size = rep(0, length(uniq_vars)))
  colnames(p_agg)[2] <- colnames(mi)[stringr::str_count(colnames(mi), "Zscore") == 1]
  for (i in 1:nrow(p_agg)) {
    sub_mi <- mi[rowSums(mi[,v_cols,drop=F] == p_agg$V[i], na.rm=T) == 1, ]
    p_agg[i,2] <- sum(abs(sub_mi[,z_col])) / nrow(mi)
    p_agg[i,3] <- sum(sub_mi$size) / nrow(mi)
  }
  p_agg <- p_agg[order(p_agg[,2], decreasing=T), ]
  p_agg$size <- ((p_agg$size - min(p_agg$size)) / (max(p_agg$size) - min(p_agg$size))) + 0.1
  if (stringr::str_count(colnames(mi)[m_col], "Cis") == 1) {                    # Title depends on cis vs. trans
    ttl <- "Variables that Differ from Randomized Expectations"
  } else {
    if (is.null(group)) {
      group <- substring(colnames(mi)[p_col], 9, 100)
    }
    ttl <- paste("Variables that Distinguish ", group, " Groups", sep="")
  }
  plot_elems <- vector("list", length(col_pals)+1)
  if (plot_bkgd == "W") {
    plot_elems[[1]] <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=c(1:nrow(p_agg)), y=p_agg[,2]), size=5*p_agg$size, color="black") +
      ggplot2::ggtitle(ttl) +
      ggplot2::scale_x_continuous(name="Variable", breaks=c(1:nrow(mi)), labels=NULL, expand=c(0.01,0.05,0.01,0)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Average Z Score") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=12), axis.title.y=ggplot2::element_text(size=15),
                     axis.text.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(),
                     plot.title=ggplot2::element_text(size=15))
    for (i in 1:length(col_pals)) {             # New layer of x axis labels for each palette
      out_cols <- p_agg[,1,drop=F]              # Which variables are in the current palette
      this_pal <- data.frame(lapply(out_cols, function(x) {grepl(paste(names(col_pals)[i],".",sep=""), x)}))
      if (sum(rowSums(this_pal)) == 0) {        # If no variables from this palette are included in any ensemble
        plot_elems[[i+1]] <- NULL               # Skip this palette and move on
        next
      }
      if (grepl("O", names(col_pals)[i])) {     # If this pal is obj, extract object variable numbers
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(names(col_pals)[i],".",sep=""),'',x)}))
        this_var <- data.frame(lapply(this_var, function(x) {gsub("_.*",'',x)}))
        this_var[as.matrix(!this_pal)] <- NA
      } else {
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(".*",names(col_pals)[i],".",sep=""),'',x)}))
        this_var[as.matrix(data.frame(lapply(this_var, function(x) {grepl("[.]",x)})))] <- NA
        this_var[as.matrix(!this_pal)] <- NA
      }
      this_lnk <- as.matrix(data.frame(lapply(out_cols, function(x) {grepl("_",x)}))) * 1
      out_syms <- data.frame(X = rep(1:nrow(p_agg), times=rowSums(this_pal)),
                             Y = unlist(apply(this_pal, 1, which)),
                             C = as.numeric(t(as.matrix(this_var)))[as.vector(t(as.matrix(this_pal)))])
      xlab <- ifelse(grepl("O", names(col_pals)[i], fixed=T), "Object Map", "Scalar Set")
      plot_elems[[i+1]] <- ggplot2::ggplot(out_syms) +
        ggplot2::geom_point(ggplot2::aes(x=X, y=Y, color=as.factor(C)), shape=15, size=5, show.legend=F) +
        ggplot2::geom_text(ggplot2::aes(x=X, y=Y, label=as.character(C)), size=3, color="white") +
        ggplot2::scale_color_manual(values = col_pals[[i]][sort(unique(out_syms$C))]) + ggplot2::theme_bw() +
        ggplot2::xlab(paste(xlab, names(col_pals)[i])) +
        ggplot2::scale_x_continuous(expand=c(0.01,0.05,0.01,0), limits = c(1,nrow(p_agg))) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1)) +
        ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                       panel.background=ggplot2::element_blank(), axis.text.y=ggplot2::element_text(size=13, color="white"),
                       axis.title.y=ggplot2::element_text(size=15, color="white"), axis.ticks.y=ggplot2::element_line(color="white"),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),panel.border=ggplot2::element_blank())
      remove(out_cols, this_pal, this_var, xlab)
    }
    plot_elems <- purrr::compact(plot_elems)
    print(cowplot::plot_grid(plotlist=plot_elems, ncol=1, align="v",
                             rel_heights=c(4, rep(1,length(plot_elems)-1))))
    remove(plot_elems)
  } else if (plot_bkgd == "B") {
    plot_elems[[1]] <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=c(1:nrow(p_agg)), y=p_agg[,2]), size=5*p_agg$size, color="white") +
      ggplot2::ggtitle(ttl) +
      ggplot2::scale_x_continuous(name="Variable", breaks=c(1:nrow(mi)), labels=NULL, expand = c(0.01,0.05,0.01,0)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Average Z Score") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     plot.background=ggplot2::element_rect(fill="black",color="black"), panel.background=ggplot2::element_rect(fill="black",color="black"),
                     axis.text.y=ggplot2::element_text(size=12, color="white"),axis.title.y=ggplot2::element_text(size=15, color="white"),
                     axis.text.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_line(color="white"),
                     panel.border=ggplot2::element_rect(color="white"), axis.ticks.y=ggplot2::element_line(color="white"),
                     plot.title=ggplot2::element_text(size=15, color="white"))
    for (i in 1:length(col_pals)) {             # New layer of x axis labels for each palette
      out_cols <- p_agg[,1,drop=F]              # Which variables are in the current palette
      this_pal <- data.frame(lapply(out_cols, function(x) {grepl(paste(names(col_pals)[i],".",sep=""), x)}))
      if (sum(rowSums(this_pal)) == 0) {        # If no variables from this palette are included in any ensemble
        plot_elems[[i+1]] <- NULL               # Skip this palette and move on
        next
      }
      if (grepl("O", names(col_pals)[i])) {     # If this pal is obj, extract object variable numbers
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(names(col_pals)[i],".",sep=""),'',x)}))
        this_var <- data.frame(lapply(this_var, function(x) {gsub("_.*",'',x)}))
        this_var[as.matrix(!this_pal)] <- NA
      } else {
        this_var <- data.frame(lapply(out_cols, function(x) {gsub(paste(".*",names(col_pals)[i],".",sep=""),'',x)}))
        this_var[as.matrix(data.frame(lapply(this_var, function(x) {grepl("[.]",x)})))] <- NA
        this_var[as.matrix(!this_pal)] <- NA
      }
      this_lnk <- as.matrix(data.frame(lapply(out_cols, function(x) {grepl("_",x)}))) * 1
      out_syms <- data.frame(X = rep(1:nrow(p_agg), times=rowSums(this_pal)),
                             Y = unlist(apply(this_pal, 1, which)),
                             C = as.numeric(t(as.matrix(this_var)))[as.vector(t(as.matrix(this_pal)))])
      xlab <- ifelse(grepl("O", names(col_pals)[i], fixed=T), "Object Map", "Scalar Set")
      plot_elems[[i+1]] <- ggplot2::ggplot(out_syms) +
        ggplot2::geom_point(ggplot2::aes(x=X, y=Y, color=as.factor(C)), shape=15, size=5, show.legend=F) +
        ggplot2::geom_text(ggplot2::aes(x=X, y=Y, label=as.character(C)), size=3, color="white") +
        ggplot2::scale_color_manual(values = col_pals[[i]][sort(unique(out_syms$C))]) + ggplot2::theme_bw() +
        ggplot2::xlab(paste(xlab, names(col_pals)[i])) +
        ggplot2::scale_x_continuous(expand=c(0.01,0.05,0.01,0), limits = c(1,nrow(p_agg))) +
        ggplot2::scale_y_continuous(expand=c(0.1,0.1)) +
        ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_text(size=13, color="black"),
                       plot.background=ggplot2::element_rect(fill="black",color="black"), panel.background=ggplot2::element_rect(fill="black",color="black"),
                       axis.title.y=ggplot2::element_text(size=15, color="black"), axis.ticks.y=ggplot2::element_line(color="black"),
                       axis.text.x=ggplot2::element_blank(),axis.title.x=ggplot2::element_text(color="white"),
                       axis.ticks.x=ggplot2::element_blank(),panel.border=ggplot2::element_blank())
      remove(out_cols, this_pal, this_var, xlab)
    }
    plot_elems <- purrr::compact(plot_elems)
    print(cowplot::plot_grid(plotlist=plot_elems, ncol=1, align="v",
                             rel_heights=c(4, rep(1,length(plot_elems)-1))))
    remove(plot_elems)
  }
  return(list(mi, p_agg))
}

#' Plot Mutual Information by Radius
#'
#' Compare across length scales how non-randomly one ensemble is patterned, in the case of cisMI,
#' or how well its patterning distinguishes specimen groups, in the case of transMI.
#'
#' @description Plot cisMI or transMI scores across radii, for a single ensemble
#' @param mi          List: cis or trans MI measurements across all radii and ensembles
#' @param ensemble    Numeric: vector of the variables included in the ensemble of interest
#' @param p_thr       Numeric: threshold corrected p-value, defaulting to 0.05
#' @param p_adj       Logical: whether to use multiple-testing-corrected p values or raw p values
#' @param group       Character: for transMI, name of the grouping factor to consider if there are more than 1
#' @param plot_bkgd   Character: whether the plot background should be "B"lack or "W"hite
#' @return      Data Frame: the focal ensemble's cis/trans MI scores and a plot of this data
#' @export
plot_MI_radius <- function(mi, ensemble, p_thr=NULL, p_adj=F, group=NULL, plot_bkgd="W") {
  out <- data.frame(radius = as.numeric(names(mi)), MI = rep(NA, length(mi)), Padjust = rep(NA, length(mi)),
                    size = rep(NA,length(mi)))
  if (is.null(p_thr)) {
    p_thr <- 0.05
  }
  for (i in 1:length(mi)) {
    v_cols <- which(stringr::str_count(colnames(mi[[i]]), "V") == 1)
    if (is.null(group)) {
      m_col <- which(stringr::str_count(colnames(mi[[i]]), "MI") == 1)
    } else {
      m_col <- which(colnames(mi[[i]]) == paste("TransMI_", group, sep=""))
    }
    if (is.null(group)) {
      if (p_adj) {
        p_col <- which(stringr::str_count(colnames(mi[[i]]), "Padjust") == 1)
      } else {
        p_col <- which(stringr::str_count(colnames(mi[[i]]), "Pvalue") == 1)
      }
    } else {
      if (p_adj) {
        p_col <- which(colnames(mi[[i]]) == paste("Padjust_", group, sep=""))
      } else {
        p_col <- which(colnames(mi[[i]]) == paste("Pvalue_", group, sep=""))
      }
    }
    if (length(p_col) != 1) {
      stop("Error: exactly one grouping factor must be specified verbatim.")
    }
    if (i == 1) {
      colnames(out)[2] <- colnames(mi[[i]])[m_col]
    }
    row <- apply(mi[[i]][,v_cols,drop=F], 1, function(x,y) {sum(x %in% y) == length(ensemble)}, y=ensemble)
    if (sum(row) > 1) {
      row <- row & (rowSums(is.na(mi[[i]][,v_cols,drop=F])) == (length(v_cols) - length(ensemble)))
    }
    row <- which(row)
    out[i,2] <- mi[[i]][row,m_col]
    mi[[i]][,p_col] <- log10(mi[[i]][,p_col])
    if (is.infinite(mi[[i]][row,p_col])) {
      if (is.null(group)) {
        z_col <- which(stringr::str_count(colnames(mi[[i]]), "Zscore") == 1)
      } else {
        z_col <- which(colnames(mi[[i]]) == paste("Zscore_", group, sep=""))
      }
      print("Warning: Exact P values are estimated by quadratic regression.")
      x <- c(abs(mi[[i]][!is.infinite(mi[[i]][,p_col]), z_col]))                   # Predict their log P from Z
      x2 <- c(x^2)                                                                 # via quadratic regression
      y <- mi[[i]][!is.infinite(mi[[i]][,p_col]), p_col]
      mod <- lm(y ~ 1 + x + x2)
      new_x <- abs(mi[[i]][row, z_col])
      new_x2 <- new_x^2
      out[i,3] <- unname(predict(mod, newdata = data.frame(x = new_x, x2 = new_x2)))
    } else {
      out[i,3] <- mi[[i]][row,p_col]
    }
  }
  if (nrow(out) > 1) {
    out$size <- ((abs(out[,2]) - min(abs(out[,2]))) / (max(abs(out[,2])) - min(abs(out[,2])))) + 0.5
  } else {
    out$size <- 1
  }
  if (colnames(out)[2] == "CisMI") {                                            # Title depends on cis vs. trans
    ttl <- "Difference from Randomized Expectations across Length Scales"
  } else {
    if (is.null(group)) {
      group <- substring(colnames(mi)[p_col], 9, 100)
    }
    ttl <- paste("Distinction among ", group, " Groups across Length Scales", sep="")
  }
  if (plot_bkgd == "W") {
    plA <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept=log10(p_thr), linetype="dotted") +
      ggplot2::geom_line(ggplot2::aes(x=out$radius, y=out[,3]), size=0.25, color="black") +
      ggplot2::geom_point(ggplot2::aes(x=out$radius, y=out[,3]), size=5*out$size, color="black") +
      ggplot2::ggtitle(ttl) +
      ggplot2::xlab("Length Scale (um)") + ggplot2::scale_y_reverse(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Corrected Log P Value") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=12), axis.title.y=ggplot2::element_text(size=15),
                     axis.text.x=ggplot2::element_text(size=12), axis.title.x=ggplot2::element_text(size=15),
                     plot.title=ggplot2::element_text(size=15))
  } else if (plot_bkgd == "B") {
    plA <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept=log10(p_thr), linetype="dotted", color="white") +
      ggplot2::geom_line(ggplot2::aes(x=out$radius, y=out[,3]), size=0.25, color="white") +
      ggplot2::geom_point(ggplot2::aes(x=out$radius, y=out[,3]), size=5*out$size, color="white") +
      ggplot2::ggtitle(ttl) +
      ggplot2::xlab("Length Scale (um)") + ggplot2::scale_y_reverse(expand = c(0.05, 0.05)) +
      ggplot2::theme_bw() + ggplot2::ylab("Corrected Log P Value") +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                     plot.background=ggplot2::element_rect(fill="black",color="black"), panel.background=ggplot2::element_rect(fill="black",color="black"),
                     axis.text.y=ggplot2::element_text(size=12, color="white"),axis.title.y=ggplot2::element_text(size=15, color="white"),
                     axis.text.x=ggplot2::element_text(size=12, color="white"),axis.title.x=ggplot2::element_text(size=15, color="white"), axis.ticks.x=ggplot2::element_line(color="white"),
                     panel.border=ggplot2::element_rect(color="white"), axis.ticks.y=ggplot2::element_line(color="white"),
                     plot.title=ggplot2::element_text(size=15, color="white"))
  }
  print(plA)
  return(out)
}
