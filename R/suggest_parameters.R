#' Suggest Radii of Neighborhoods to Census
#'
#' Calculate the neighborhood radii in pixels based on desired radii in microns and the image resolution
#' in X, Y, and Z.
#'
#' @description Suggest radii in pixels of neighborhoods to collect for the census, based on desired radii in microns
#' @param target       Vector: desired radii in microns
#' @param pix_res      Named Vector: resolution of the image in microns/pixel in the X, Y, and Z dimensions
#' @return      Vector: pixel radii to use for censusing
#' @export
suggest_radii <- function(target, pix_res) {
  plot_range <- floor(max(min(target) - 10, 0)):ceiling(max(target) + 10)
  out <- vector("list", length(target))
  names(out) <- target
  for (i in 1:length(target)) {
    out[[i]] <- c(round(target[i]/pix_res[1]), round(target[i]/pix_res[2]), round(target[i]/pix_res[3]))
    out[[i]][out[[i]] < 1] <- 1
  }
  for (i in 1:3) {
    if (i == 1) {
      ttl <- "X"
    } else if (i ==2) {
      ttl <- "Y"
    } else {
      ttl <- "Z"
    }
    pl <- ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x=target, xend=target, y=0, yend=sapply(out,"[[",i)),
                            color="red", linewidth=1.5) +
      ggplot2::geom_segment(ggplot2::aes(x=0, xend=target, y=sapply(out,"[[",i), yend=sapply(out,"[[",i)),
                            color="red", linewidth=1.5) +
      ggplot2::geom_line(ggplot2::aes(x=plot_range, y=plot_range/pix_res[i]), col="black", linewidth=3) +
      ggplot2::xlab("Neighborhood Radius (microns)") + ggplot2::ylab ("Neighborhood Radius (pixels)") +
      ggplot2::theme_bw() + ggplot2::ggtitle(paste("Micron to Pixel Conversion in ", ttl, sep="")) +
      ggplot2::theme(axis.title=ggplot2::element_text(size=20), axis.text.x=ggplot2::element_text(size=15),
                     axis.text.y=ggplot2::element_text(size=15), plot.title=ggplot2::element_text(size=20)) +
      ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))
    print(pl)
  }
  return(out)
}

#' Suggest Number of Neighborhoods to Census
#'
#' Calculate the number of neighborhoods to census, based on neighborhood radii and the desired coverage.
#' Coverage is the number of neighborhoods each non-background pixel is included in, on average.
#' Coverage should not exceed 5x, to prevent pseudo-replication.
#'
#' @description Suggest the number of neighborhoods of each radius to census, based on a desired coverage
#' @param coverage      Vector: desired image coverage/s
#' @param radii         List: named vectors of neighborhood radii for the X, Y, and Z dimensions in pixels
#' @param images        List: 3d arrays of all images to be analyzed
#' @return      Vector: Number of neighborhoods to census at each radius
#' @export
suggest_number <- function(coverage, radii, images) {
  if (coverage > 5) {
    print("WARNING: Coverage exceeding 5x will yield invalid inference due to pseudo-replication.")
  }
  if (is.list(images)) {                                                        # Remove linked scalar images
    images <- images[grepl("O",names(images)) | (purrr::map_dbl(substr(names(images),2,9), ~sum(.x==substr(names(images),2,9))) == 1)]
    for (i in 1:length(images)) {
      images[[i]] <- apply(images[[i]], c(1,2,3), sum)
    }
    images <- simplify2array(images)
  }
  images <- apply(images, c(1,2,3), sum)
  num_pix <- sum(images > 0)
  nbhd_vol <- calc_vols(radii, dim(images))
  plot_range <- c(max(min(coverage) - 1, 0.001), coverage, max(coverage) + 1)
  out <- data.frame(X = rep(plot_range, each=length(radii)),
                    L = rep(as.numeric(names(radii)), length(plot_range)),
                    N = rep(NA, length(plot_range)*length(radii)))
  for (i in 1:length(nbhd_vol)) {
    out$N[out$L == names(radii)[i]] <- ceiling(plot_range*num_pix/nbhd_vol[i])
  }
  pl <- ggplot2::ggplot() +
    ggplot2::geom_line(data=out, ggplot2::aes(x=X, y=N, color=factor(L)), linewidth=3) +
    ggplot2::xlab("Coverage") + ggplot2::ylab("Number of Neighborhoods") + ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::theme(axis.title=ggplot2::element_text(size=20), axis.text.x=ggplot2::element_text(size=15),
                   axis.text.y=ggplot2::element_text(size=15)) +
    ggplot2::scale_color_grey(name="Neighborhood Radius (microns)")
  print(pl)
  out <- out$N[out$X == coverage]
  return(out)
}
