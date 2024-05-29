#' Plot an Image
#'
#' Visualize either the entirety of a 2d image or a slice of a 3d image. Subsets of the image's scalars or objects
#' may be specified to reduce what is plotted. The plotted image is also saved as a .tif file.
#'
#' @description Plot an image from its array representation
#' @param img       Array 4d: the image to plot
#' @param img_type  Character: whether the image represents "O"bjects or "S"calars
#' @param col_pal   Vector: color palette for the objects or scalars
#' @param slice     Named Vector: "X", "Y", or "Z" coordinate at which to slice
#' @param objects   Vector: for object images, numeric codes of the objects to visualize, or all of them with NULL
#' @param scalars   Vector: for scalar images, numeric codes of the scalars to visualize, or all of them with NULL
#' @param enh_cnt   Numeric: quantile of pixels to saturate for each scalar
#' @return      Array 4d: array representation of the image that was plotted.
#' @export
plot_image <- function(img, img_type, col_pal, slice = c("Z" = 1), objects=NULL, scalars=NULL, enh_cnt=0) {
  if (img_type == "O") {                                                        # If it's an object image
    if (names(slice) == "X") {                                                  # Excise the correct slice
      img <- img[slice,,,1]
    } else if (names(slice) == "Y") {
      img <- img[,slice,,1]
    } else if (names(slice) == "Z") {
      img <- img[,,slice,1]
    } else {
      stop("Slice must be named X, Y, or Z")
    }
    if (is.null(objects)) {                                                     # If all objs are desired (default)
      objects <- as.numeric(names(table(img)))                                  # get numeric codes
      objects <- objects[objects != 0]                                          # Don't include background
    }
    any_bkgd <- any(img == 0)                                                   # Note whether there is any bkgd
    objects <- sort(objects)                                                    # Put desired objects in order
    col_pal <- col_pal[objects]                                                 # Reduce col pal to just desired objs
    img[!(img %in% c(0,objects))] <- max(objects) + 1                           # Any undesired obj set to new obj
    if (sum(img == (max(objects) + 1)) > 0) {                                   # If any undesired objs,
      objects <- c(objects, max(objects)+1)                                     # include it among the objects and
      col_pal <- c(col_pal, "#808080")                                          # give it gray in the palette
    }
    for (i in c(1:length(objects))) {                                           # Reduce obj codes to 1, 2, 3, ...
      img[img == objects[i]] <- i                                               # where 0 is bkgd
    }
    if (any_bkgd) {                                                             # If bkgd, include black in pal
      col_pal <- c("#000000", col_pal)
    }
    img_plot <- array(0, dim=c(dim(img)[1], dim(img)[2], 3))                    # 3D array for RGB at each pix
    for (j in 1:length(col_pal)) {                                              # For each object to be plotted
      for (i in 1:3) {                                                          # Specify R, G, & B separately
        img_plot[,,i][img==(j-1)] <- col2rgb(col_pal[j])[i]
      }                                                                         # for each scalar's contribution
    }
    par(mar=c(0,0,0,0))                                                         # No margins = bigger picture
    raster::plotRGB(raster::brick(img_plot, xmn=0, xmx=dim(img)[2], ymn=0, ymx=dim(img)[1]),
                    maxpixels=dim(img)[1]*dim(img)[2])
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    raster::writeRaster(raster::brick(img_plot, xmn=0, xmx=dim(img)[2], ymn=0, ymx=dim(img)[1]),
                        filename="Plot.tif", overwrite=T, datatype="INT1U",
                        options="TFW=YES", format="GTiff")
    return(img_plot)
  } else if (img_type == "S") {                                                 # If it's a scalar image
    if (is.null(scalars)) {
      scalars <- c(1:dim(img)[4])
    }
    if (names(slice) == "X") {                                                  # Excise the correct slice
      img <- abind::adrop(img[slice,,,,drop=FALSE], drop=1)
    } else if (names(slice) == "Y") {
      img <- abind::adrop(img[,slice,,,drop=FALSE], drop=2)
    } else if (names(slice) == "Z") {
      img <- abind::adrop(img[,,slice,,drop=FALSE], drop=3)
    } else {
      stop("Slice must be named X, Y, or Z")
    }
    img_plot <- array(0, dim=c(dim(img)[1], dim(img)[2], 3))                    # 3D array for RGB at each pix
    for (j in 1:length(scalars)) {                                              # For each scalar to be plotted
      for (i in 1:3) {                                                          # Increase R, G, & B separately
        orig_max <- max(img[,,scalars[j]])
        img_plot[,,i] <- img_plot[,,i] + (img[,,scalars[j]]/quantile(img[,,scalars[j]], 1-enh_cnt)) * (col2rgb(col_pal[scalars[j]])[i])
        img_plot[,,i][img_plot[,,i] > orig_max] <- orig_max
      }                                                                         # for each scalar's contribution
    }
    par(mar=c(0,0,0,0))                                                         # No margins = bigger picture
    raster::plotRGB(raster::brick(img_plot, xmn=0, xmx=dim(img)[2], ymn=0, ymx=dim(img)[1]),
                    maxpixels=dim(img)[1]*dim(img)[2])
    par(mar = c(5.1, 4.1, 4.1, 2.1))                                            # Reset to default margins
    raster::writeRaster(raster::brick(img_plot, xmn=0, xmx=dim(img)[2], ymn=0, ymx=dim(img)[1]),
                        filename="Plot.tif", overwrite=T, datatype="INT1U",
                        options="TFW=YES", format="GTiff")
    return(img_plot)
  }                                                                             # Plot the resulting image
}

#' Make a Color Palette
#'
#' Generate a set of color with maximal contrast. Beyond 12 colors, visual distinction can still be difficult.
#'
#' @description Randomly select contrasting colors to form a palette of a specified size
#' @param num_cols  Numeric: number of colors
#' @return      Vector: hex codes for each color in the palette
#' @export
make_palette <- function(num_cols) {
  cols <- matrix(rep(NA, num_cols*3), nrow=num_cols)                     # Matrix of hsv representations of colors
  hues <- seq(0, (1-1/num_cols), (1/num_cols)) + runif(1,0,(1/num_cols)) # Hue of each color, maximally & evenly spaced
  cols[,1] <- sample(hues)                                               # Randomly reorder the hues
  cols[,2] <- rbeta(num_cols,1.5,1)                                      # Sample saturations from beta dist to skew toward 1
  cols[,3] <- rbeta(num_cols,4.5,1)                                      # Sample values from beta dist to skew toward 1
  col_pal <- rep(NA, num_cols)                                           # Initialize hex code color palette
  for (i in 1:num_cols) {
    col_pal[i] <- hsv(h=cols[i,1], s=cols[i,2], v=cols[i,3])             # Translate each row into a single hex color
  }
  return(col_pal)
}

#' Plot a Color Palette
#'
#' Plot a color palette, optionally including short names for what each  variable represents biologically.
#'
#' @description Plot all colors in a palette with labels
#' @param col_pal    Vector: color palette to plot
#' @param axis_label Character: label for the entire palette of colors
#' @param col_labels Vector: labels for individual colors
#' @param vertical   Logical: whether to plot the colors vertically or horizontally
#' @param plot_bkgd  Character: whether the plot background should be "W"hite or "B"lack
#' @return      Nothing: only a plot of the color palette is created
#' @export
plot_palette <- function(col_pal, axis_label, col_labels, vertical=F, plot_bkgd="W") {
  if (missing(axis_label)) {
    axis_label <- ""
  }
  if (missing(col_labels)) {
    col_labels <- as.character(c(1:(length(col_pal))))
  }
  if (vertical) {
    col_pal <- rev(col_pal)
    col_labels <- rev(col_labels)
  }
  par(bg = ifelse(plot_bkgd == "B", "black", "white"))
  if (vertical) {
    par(mar = c(5.1, 4.1 + 0.5*(max(nchar(col_labels))-4), 4.1, 2.1))           # Customize margins to accommodate names
    barplot(rep(1,length(col_pal)), col=col_pal, axes=F, names=col_labels, horiz=T, las=T,
            border=ifelse(plot_bkgd == "B", "white", "black"), col.axis=ifelse(plot_bkgd == "B", "white", "black"),
            col.lab=ifelse(plot_bkgd == "B", "white", "black"))
    mtext(side = 2, text = axis_label, line = (max(nchar(col_labels)))/2, col=ifelse(plot_bkgd == "B", "white", "black"))
  } else {
    barplot(rep(1,length(col_pal)), col=col_pal, axes=F, xlab=axis_label, names=col_labels,
            border=ifelse(plot_bkgd == "B", "white", "black"), col.axis=ifelse(plot_bkgd == "B", "white", "black"),
            col.lab=ifelse(plot_bkgd == "B", "white", "black"))
  }
  par(mar = c(5.1, 4.1, 4.1, 2.1), bg = "white")                                # Reset to default margins
}
