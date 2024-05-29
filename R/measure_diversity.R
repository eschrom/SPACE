#' Measure Composition of a Full Image
#'
#' Quantify and plot the fraction of an image that belongs to each object or that is positive for each scalar.
#' The number of different object or scalars and how evenly represented they are throughout the image
#' determines the image's alpha diversity. Alpha diversity is not available for object-specific scalars.
#'
#' @description Quantify the object or scalar composition of an entire image
#' @param img       Named List: 4d array of the image whose composition is to be measured
#' @param img_type  Character: whether the image represents "O"bjects or "S"calars
#' @param col_pal   Named List: color palette for the image
#' @param plot_bkgd Character: whether the plot background should be "W"hite or "B"lack
#' @return      Data Frame: Percentage of image belonging to each variable and a plot of this data
#' @export
alpha_diversity <- function(img, img_type, col_pal, plot_bkgd="W") {
  if (is.list(img) && (length(img) > 1)) {
    stop("Error: Alpha diversity can be measured for only one image at a time.")
  } else if (is.list(img)) {
    img <- img[[1]]
  }
  if (is.list(col_pal) && (length(col_pal) > 1)) {
    stop("Error: Alpha diversity can be measured for only one image with one color palette at a time.")
  } else if (is.list(col_pal)) {
    col_pal <- col_pal[[1]]
  }
  out <- rep(NA, length(col_pal))
  if (img_type == "O") {
    for (i in 1:length(col_pal)) {
      out[i] <- sum(img == i)
    }
  } else if (img_type == "S") {
    for (i in 1:length(col_pal)) {
      out[i] <- Reduce("+",img[,,,i])
    }
  }
  out <- out/sum(out)
  shan_ent <- -sum(out[out>0]*log2(out[out>0]))
  out <- data.frame(V = 1:length(col_pal), P = 100*out)
  if (plot_bkgd == "W") {
    pl_A <- ggplot2::ggplot(out, ggplot2::aes(x=V, y=P)) + ggplot2::geom_bar(stat="identity", fill=col_pal) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                     text=ggplot2::element_text(size=15)) +
      ggplot2::scale_x_continuous(breaks=out$V, labels=as.character(out$V)) +
      ggplot2::ggtitle(paste("Image Composition of ", ifelse(img_type=="O","Objects ", "Scalars "),
                             "\nAlpha Diversity: ", round(shan_ent, 2), " bits", sep="")) +
      ggplot2::ylab("Percentage (%)") + ggplot2::xlab(ifelse(img_type=="O","Objects", "Scalars"))
  } else if (plot_bkgd == "B") {
    pl_A <- ggplot2::ggplot(out, ggplot2::aes(x=V, y=P)) + ggplot2::geom_bar(stat="identity", fill=col_pal) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                     text=ggplot2::element_text(size=15, color="white"),
                     axis.text.y=ggplot2::element_text(color="white"), axis.text.x=ggplot2::element_text(color="white"),
                     plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black")) +
      ggplot2::scale_x_continuous(breaks=out$V, labels=as.character(out$V)) +
      ggplot2::ggtitle(paste("Image Composition of ", ifelse(img_type=="O","Objects", "Scalars"),
                             "\nAlpha Diversity: ", round(shan_ent, 2), " bits", sep="")) +
      ggplot2::ylab("Percentage (%)") + ggplot2::xlab(ifelse(img_type=="O","Objects", "Scalars"))
  }
  plot(pl_A)
  return(list(out, shan_ent))
}

#' Measure Compositions of Objects Separately
#'
#' Quantify and plot the fraction of each parent object that belongs to each sub-object or that is positive for each
#' scalar. The difference in composition across the parent objects determines the image's beta diversity.
#' The number of different objects or scalars and how evenly represented they are throughout each parent object
#' determines the parent object's alpha diversity, which is also reported.
#'
#' @description Quantify the object or scalar composition of each parent object in an image separately
#' @param img          Named List: 4d arrays of the images of the parent objects and the constituent objects or scalars
#' @param img_type     Character: whether the constituent variables are "O"bjects or "S"calars
#' @param col_pal      Named List: color palettes for the object images and the constituent objects or scalars
#' @param plot_bkgd    Character: whether the plot background should be "W"hite or "B"lack
#' @return      Data frame: object or scalar composition of each parent object separately and a plot of this data
#' @export
beta_diversity <- function(img, img_type, col_pal, plot_bkgd="W") {
  if (!is.list(img) || (length(img) != 2)) {
    stop("Error: Exactly two images must be provided in a list: the parent objects and the constituent variables.")
  }
  if (!is.list(col_pal) || (length(col_pal) != 2)) {
    stop("Error: Exactly two color palettes must be provided in a list: the parent objects and the constituent variables.")
  }
  obj_img <- img[[1]]
  img <- img[[2]]
  obj_col_pal <- col_pal[[1]]
  col_pal <- col_pal[[2]]
  objects <- 1:length(obj_col_pal)
  if (img_type == "O") {
    vars <- 1:length(col_pal)
  } else if (img_type == "S") {
    vars <- 1:dim(img)[4]
  }
  out <- as.data.frame(matrix(NA, nrow = length(objects), ncol = (length(vars) + 1)))
  colnames(out) <- c("O", paste("V", vars, sep=""))
  out$O <- objects
  if (img_type == "O") {
    for (i in 1:nrow(out)) {
      for (j in 1:length(vars)) {
        if (sum(obj_img == i) == 0) {
          out[i, j+1] <- 0
        } else {
          out[i, j+1] <- sum(img[obj_img == i] == j)
        }
      }
    }
  } else if (img_type == "S") {
    for (i in 1:nrow(out)) {
      for (j in 1:length(vars)) {
        if (sum(obj_img == i) == 0) {
          out[i, j+1] <- 0
        } else {
          out[i, j+1] <- sum(img[,,,j,drop=F][obj_img == i])
        }
      }
    }
  }
  pos_rows <- rowSums(out[,2:ncol(out)]) > 0
  out[pos_rows, 2:ncol(out)] <- out[pos_rows, 2:ncol(out)] / rowSums(out[pos_rows, 2:ncol(out)])
  out$E <- rep(0, nrow(out))
  out$E[pos_rows] <- apply(out[pos_rows,2:ncol(out)], 1, function(x) {-sum(x[x!=0] * log2(x[x!=0]))})
  out_beta <- out[pos_rows, ]
  region_wts <- rep(NA, nrow(out_beta))                               # Beta diversity across regions is
  for (i in 1:length(region_wts)) {                                   # avg KL div of each regional
    region_wts[i] <- sum(obj_img == out_beta$O[i])                    # comp from total avg comp,
  }                                                                   # where both avgs are weighted by
  region_wts <- region_wts / sum(region_wts)                          # relative prop of tissue belonging
  wt_avg_dist <- colSums(out_beta[,2:(ncol(out_beta)-1)] * region_wts)# to each region
  wt_avg_dist <- wt_avg_dist / sum(wt_avg_dist)
  tiny_prob <- 0.000001 / length(wt_avg_dist)
  if (any(wt_avg_dist == 0)) {
    wt_avg_dist[wt_avg_dist == 0] <- tiny_prob
    wt_avg_dist <- wt_avg_dist / sum(wt_avg_dist)
  }
  KL_divs <- rep(NA, nrow(out_beta))
  for (i in 1:nrow(out_beta)) {
    region_dist <- out_beta[i,2:(ncol(out_beta)-1)]
    if (any(region_dist == 0)) {
      region_dist[region_dist == 0] <- tiny_prob
      region_dist <- region_dist / sum(region_dist)
    }
    KL_divs[i] <- sum(region_dist * log2(region_dist / wt_avg_dist))
  }
  beta_div <- sum(KL_divs * region_wts)
  out[,2:(ncol(out)-1)] <- 100*out[,2:(ncol(out)-1)]
  out_reformat <- tidyr::pivot_longer(out, 2:(ncol(out)-1), names_to="V", values_to="P")
  if (plot_bkgd == "W") {
    pl_A <- ggplot2::ggplot(out_reformat, ggplot2::aes(fill=factor(V, levels=unique(V)), y=P, x=O)) +
      ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                     text=ggplot2::element_text(size=15), axis.text.x=ggplot2::element_text(color=obj_col_pal)) +
      ggplot2::scale_fill_manual(values=col_pal, name=ifelse(img_type=="O", "Object", "Scalar"),
                                 labels=substring(out_reformat$V,2)) +
      ggplot2::scale_x_continuous(breaks=unique(out_reformat$O), labels=as.character(unique(out_reformat$O)),
                                  name="Parent Object") +
      ggplot2::scale_y_continuous(breaks=seq(0,100,10), limits=c(0,105), labels=seq(0,100,10), name="Percentage (%)") +
      ggplot2::ggtitle(paste("Separate Compositions of Parent Objects",
                             "\nBeta Diversity: ", round(beta_div, 2), " bits, Alpha Diversity Below", sep="")) +
      ggplot2::annotate(geom="text", x=c(out$O), y=105, label=round(out$E,2),
                        size=5, color=obj_col_pal)
  } else if (plot_bkgd == "B") {
    pl_A <- ggplot2::ggplot(out_reformat, ggplot2::aes(fill=factor(V, levels=unique(V)), y=P, x=O)) +
      ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                     text=ggplot2::element_text(size=15, color="white"), legend.background=ggplot2::element_rect(fill="black"),
                     axis.text.y=ggplot2::element_text(color="white"), axis.text.x=ggplot2::element_text(color=obj_col_pal),
                     plot.background=ggplot2::element_rect(fill="black"), panel.background=ggplot2::element_rect(fill="black"),
                     legend.key=ggplot2::element_rect(color="black",fill="black")) +
      ggplot2::scale_fill_manual(values=col_pal, name=ifelse(img_type=="O", "Object", "Scalar"),
                                 labels=substring(out_reformat$V,2)) +
      ggplot2::scale_x_continuous(breaks=unique(out_reformat$O), labels=as.character(unique(out_reformat$O)),
                                  name="Parent Object") +
      ggplot2::scale_y_continuous(breaks=seq(0,100,10), limits=c(0,105), labels=seq(0,100,10), name="Percentage (%)") +
      ggplot2::ggtitle(paste("Separate Compositions of Parent Objects",
                             "\nBeta Diversity: ", round(beta_div, 2), " bits, Alpha Diversity Below", sep="")) +
      ggplot2::annotate(geom="text", x=c(out$O), y=105, label=round(out$E,2),
                        size=5, color=obj_col_pal)
  }
  print(pl_A)
  return(list(out, beta_div))
}
