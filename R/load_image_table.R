#' Load an Image
#'
#' Load an image into array format. Object images must be RGB color .tif. Scalar images must be 8-bit grayscale .tif
#' possibly including multiple channels. 3d or 2d images are permitted.
#'
#' @description Load an object or scalar image and transform it into a 4d array
#' @param in_file    Character: path to the input file
#' @param img_type   Character: whether the image is a "O"bject or "S"calar image
#' @param bkgd_col   Character: for object images, the background color's hex code/nickname, or NULL if none
#' @param num_chs    Numeric: for scalar images, the number of channels in the image
#' @param keep_chs   Vector: for scalar images, the channels to keep upon loading, or NULL if all
#' @return      List: 4d array representing the image: x,y,z,c. And for object maps, the extracted color palette
#' @export
load_image <- function(in_file, img_type, bkgd_col, num_chs, keep_chs=NULL) {
  img <- tiff::readTIFF(in_file, all = T, as.is = F)                            # Read in image
  xy_dim <- c(dim(img[[1]])[1], dim(img[[1]])[2])                               # XY dimensions of image
  zc_dim <- length(img)                                                         # Total number of Zs and Cs
  img <- round(simplify2array(img)*255)                                         # Collapse to array: x,y,rgb,z
  if (img_type == "O") {                                                        # If the image is an object map:
    if (zc_dim == 1) {                                                          # If only one z, still make a
      img <- array(img, dim = c(dim(img),1))                                    # 4th dim for the array
    }                                                                           # Collapse rgb cols to hex codes
    img <- apply(img, c(1,2,4), function(x) {rgb(x[1],x[2],x[3],maxColorValue=255)})
    col_pal <- unique(as.vector(img))                                           # Unique colors
    if (!is.null(bkgd_col)) {                                                   # If there is background
      if (substr(bkgd_col,1,1) != "#") {                                        # If bkgd isn't given as hex
        bkgd_col <- col2rgb(bkgd_col)                                           # Convert to rgb then to hex
        bkgd_col <- rgb(bkgd_col[1], bkgd_col[2], bkgd_col[3], maxColorValue=255)
      }
      col_pal <- c(bkgd_col, col_pal[col_pal != bkgd_col])                      # Reorder to put bkgd 1st
    }
    img <- array(match(img, col_pal), dim=dim(img))                             # Match each pix to its color code
    if (!is.null(bkgd_col)) {                                                   # If there is background,
      img <- img - as.integer(1)                                                # Subtract 1 from every code
      col_pal <- col_pal[-1]                                                    # Remove background from palette
    }
    img <- array(as.integer(img), dim=c(dim(img),1))                            # Make 4th dim, but just 1 channel
    return(list(img, col_pal))
  } else if (img_type == "S") {                                                 # If the image is a scalar set:
    if (zc_dim == num_chs) {                                                    # If there is only one z slice
      img <- array(as.integer(img), dim=c(dim(img)[1:2], 1, dim(img)[3]))       # Insert a single z slice
    } else {                                                                    # If there are multiple z slices
      num_zs <- zc_dim / num_chs                                                # Calculate how many
      img_new <- array(0, dim=c(dim(img)[1:2], num_zs, num_chs))                # Current order: all chs at z1,
      for (i in 1:num_chs) {                                                    # all chs at z2, and so on
        img_new[,,,i] <- img[,,(seq(1,zc_dim-num_chs+1,by=num_chs)+(i-1))]      # Separate zs and chs
      }
      img <- img_new                                                            # Overwrite old image
    }
    if (!is.null(keep_chs)) {                                                   # If a subset of channels is given
      img <- img[,,,keep_chs]                                                   # keep only these
    }
    img <- array(as.integer(img), dim=dim(img))
    return(img)
  }
}

#' Load a Table
#'
#' Load a table in .csv format. The table can be a profile, object, or link table.
#' For profile tables, column 1 must be labeled "Object" and give numeric object IDs.
#' Column 2 must be labeled "Count" and give the count of pixels or items belonging to each object.
#' Optionally, columns 3-5 can labeled "X", "Y", and "Z" and give the coordinates of one pixel belonging
#' to each object in the image, to speed the matching of objects between the table and the corresponding image.
#' The remaining columns give the compositional profile of each object, in terms of scalars or sub-objects.
#' For object tables, columns 1-3 must be labeled "X", "Y", and "Z" and give the centroid coordinates for each object.
#' Column 4 must be labeled "Object" and give the numeric code for each object ID. Optional further columns
#' may give scalar intensity values. For link tables, the rows must be labeled with the alpha-numeric codes of the
#' objects, e.g. "O1.1", and the columns must be labeled with the alpha-numeric codes of the scalars, e.g. "S1.1".
#' The object set and scalar set must share the same numeric code. Each entry is 1 or 0, indicating whether the
#' indicated combination should be analyzed or ignored.
#'
#' @description Load a profile, object, or link table
#' @param in_file    Character: path to the table file
#' @param table_type Character: whether the input is a "P"rofile, "O"bject, or "L"ink table
#' @param img        Array 4d: for profile tables, the object image to which the table corresponds
#' @param col_pal    Vector: for profile tables, the color palette to which the table corresponds
#' @return      Data frame: the profile table where objects are reordered to match the corresponding object image, or the object or link table exactly as in the source file
#' @export
load_table <- function(in_file, table_type, img, col_pal) {
  if (!(table_type %in% c("P","O","L"))) {
    stop("Error: table_type must be P, O, or L.")
  }
  if (table_type == "L") {                                                      # For link tables,
    out <- as.matrix(read.csv(in_file, header=T, row.names=1))                  # assume row & col names
    return(out)
  }
  if (table_type == "O") {                                                      # For object tables,
    out <- read.csv(in_file, header=T)                                          # assume only col names
    if (!("Z" %in% colnames(out))) {                                            # If there is no Z column,
      out$Z <- rep(1, nrow(out))                                                # add it
      out <- out[,c("X","Y","Z",colnames(out)[!(colnames(out) %in% c("X","Y","Z"))])]
    }
    return(out)
  }                                                                             # Anything else is a profile table
  out <- read.csv(in_file, header=T)                                            # Assume only col names
  out <- out[order(out$Object), ]                                               # Put objs in ascending order.                                                                            # Can trust the ascending order
  if (missing(img) && missing(col_pal)) {                                       # if no image or col pal is given
    warning("Without providing an image and color palette, the objects might not match the order of the color palette.")
    return(out)                                                                 # Otherwise, be careful to match
  }                                                                             # objs in table to order of col pal
  if ("X" %in% colnames(out) && "Y" %in% colnames(out) && "Z" %in% colnames(out)) { # If coors of 1pix/object given
    any_bkgd <- 1*any(img == 0)                                                 # Whether there is any background
    objID_reorder <- rep(NA, length(col_pal) - any_bkgd)                        # Account for all non-bkgd objects
    for (i in 1:nrow(out)) {                                                    # Object reordering from what's in csv to img
      object <- img[out$Row[i], out$Column[i], ]                                # index = img objID, value = .csv objID
      if (object > 0) {                                                         # As long as the image object is not bkgd
        objID_reorder[object] <- i                                              # Record at the img object index what the
      }                                                                         # csv objID is
    }
    out <- out[objID_reorder, ]                                                 # Reorder rows to match img objID
    out$Object <- c(1:(length(col_pal) - any_bkgd))                             # Rename to img objIDs
    out$X <- NULL                                                               # Remove representative row & col pixels
    out$Y <- NULL
    out$Z <- NULL
    return(out)
  } else {
    if (length(unique(out[,2])) < length(out[,2])) {
      stop("Cannot match objects from image to table based on counts alone: specify a pixel coordinate for each object.")
    }
    if (sum(out$Count) >= sum(img > 0)) {                                       # In pix mode, match pixel counts
      any_bkgd <- 1*any(img == 0)                                               # Whether there is any background
      pix_cnt_img <- rep(NA, length(col_pal))                                   # Pixel count for each non-bkgd,
      for (i in c(1:length(pix_cnt_img))) {                                     # counted directly from the image
        pix_cnt_img[i] <- sum(img == i)
      }
      objID_reorder <- rep(NA, length(pix_cnt_img))                             # Object reordering from what's in csv to img
      for (i in c(1:length(pix_cnt_img))) {                                     # index = img objID, value = .csv objID
        objID_reorder[i] <- which(out$Count == pix_cnt_img[i])
      }
      out <- out[objID_reorder, ]                                               # Reorder rows to match img objID
      out$Object <- c(1:length(col_pal))                                        # Rename to img objIDs; 0 is background
      return(out)
    } else {                                                                    # In obj mode, match object counts
      any_bkgd <- 1*any(img == 0)                                               # Whether there is any background
      obj_cnt_img <- rep(NA, length(col_pal))
      for (i in c(1:length(obj_cnt_img))) {                                     # For each objID in image, count patches
        patch_img <- raster::raster(matrix(img == i, nrow = dim(img)[1], ncol = dim(img)[2]))
        obj_cnt_img[i] <- max(raster::freq(raster::clump(patch_img))[,1], na.rm = T)
      }
      objID_reorder <- rep(NA, length(obj_cnt_img))                             # Object reordering from what's in csv to img
      for (i in c(1:length(obj_cnt_img))) {                                     # index = img objID, value = .csv objID
        objID_reorder[i] <- which(out$Count == obj_cnt_img[i])
      }
      out <- out[objID_reorder, ]                                               # Reorder rows to match img objID
      out$Object <- c(1:length(col_pal))                                        # Rename to img objIDs
      return(out)
    }
  }
}
