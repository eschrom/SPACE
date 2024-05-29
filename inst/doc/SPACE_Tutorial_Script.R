library(SPACE)
setwd("C:/Users/schromec/Desktop/SPACE/SPACE_Tutorial")
load(".Rdata")
##### Import and Visualize Input Data ##########################################
# Fig. 1. Object image from cell segmentation and classification
IMG_obj1 <- load_image(in_file = "HuLN_Object1.tif", img_type = "O", bkgd_col = "black")
PAL_obj1 <- IMG_obj1[[2]]
IMG_obj1 <- IMG_obj1[[1]]
plot_image(img = IMG_obj1, img_type = "O", col_pal = PAL_obj1)
# Fig. 2. Object image from pixel classification
IMG_obj2 <- load_image(in_file = "HuLN_Object2.tif", img_type = "O", bkgd_col = "black")
PAL_obj2 <- IMG_obj2[[2]]
IMG_obj2 <- IMG_obj2[[1]]
plot_image(img = IMG_obj2, img_type = "O", col_pal = PAL_obj2)
# Fig 3. Color palette for the pixel-based object image
plot_palette(col_pal = PAL_obj2, axis_label = "Object Types", plot_bkgd = "W")
# No Fig. Load profile table for pixel-based objects
TBL_obj2_prof <- load_table(in_file = "HuLN_Object2_ProfTable.csv", table_type = "P", img = IMG_obj2, 
                            col_pal = PAL_obj2)
head(TBL_obj2_prof)
# Fig. 4. Profiles of biomolecule expression for the pixel-based objects
plot_table(prof_table = TBL_obj2_prof, compare = "A", normalize = "U", tile_plots = T, plot_bkgd = "W")
# Fig. 5. Color palette for the pixel-based object image with descriptive names
NMS_obj2 <- c("Matrix","Blood Vessel","Mesenchyma","Smooth Muscle","FDC Network","Lymphatic Vessel")
plot_palette(col_pal = PAL_obj2, axis_label = "Object Type", col_labels = NMS_obj2, vertical = T, plot_bkgd = "W")
# Fig. 6. Color palette for the segmentation-based object image
plot_palette(col_pal = PAL_obj1, axis_label = "Object Type", plot_bkgd = "W")
# Fig. 7. Profiles of biomolecule expression for the segmentation-based objects
TBL_obj1_prof <- load_table(in_file = "HuLN_Object1_ProfTable.csv", table_type = "P", img = IMG_obj1, 
                            col_pal = PAL_obj1)
plot_table(prof_table = TBL_obj1_prof, compare = "A", normalize = "U", tile_plots = T, plot_bkgd = "W")
# No Fig. Merging segmented objects that are too similar.
MERGE <- merge_objects(prof_table = TBL_obj1_prof, img = IMG_obj1, col_pal = PAL_obj1, 
                       obj_groups = list(c(7,8), c(9,14,15,16,17,18)))
TBL_obj1_prof <- MERGE[[1]]
IMG_obj1 <- MERGE[[2]]
PAL_obj1 <- MERGE[[3]]
remove(MERGE)
# Fig. 8. Segmentation-based object image after merging
plot_image(img = IMG_obj1, img_type = "O", col_pal = PAL_obj1)
# Fig. 9. Profiles of biomolecule expression for the segmentation-based objects after merging
plot_table(prof_table = TBL_obj1_prof, tile_plots = T)
# Fig. 10. Color palette for the segmentation-based object image with descriptive names
NMS_obj1 <- c("T-reg","CD8+ T","B-Eng CD4_ T","CD4+ T","APC-Eng CD4+ T","GC B","DC","Foll B",
              "M1 Mac","Coll Mac","M2 Mac","Plasma","MAIT","Reg Mac")
plot_palette(col_pal = PAL_obj1, axis_label = "Object Type", col_labels = NMS_obj1, vertical=T, plot_bkgd = "W")
# No Fig. Create new palette with more distinguishable colors
make_palette(num_cols = 14)
load("PAL_obj1_alt.Rdata")
# Fig. 11. New color palette for the segmentation-based object image
plot_palette(col_pal = PAL_obj1_alt, axis_label = "Object Type", col_labels = NMS_obj1, vertical = T, 
             plot_bkgd = "W")
# Fig. 12. Segmentation-based object image with new color palette
plot_image(img = IMG_obj1, img_type = "O", col_pal = PAL_obj1_alt)
# No Fig. Import table of centroids and object types for segmentation-based objects
TBL_obj1_obj <- load_table("HuLN_Object1_ObjTable_Types.csv", table_type = "O")
head(TBL_obj1_obj)
# No Fig. Merge objects in a table
MERGE <- merge_objects(obj_table = TBL_obj1_obj, obj_groups = list(c(8,15), c(3,5,6,7,9,11)))
TBL_obj1_obj <- MERGE[[1]]
remove(MERGE)
# No Fig. Import table of centroids and expression for segmentation-based objects
TBL_obj1_xpr <- load_table("HuLN_Object1_ObjTable_Exprs.csv", table_type = "O")
head(TBL_obj1_xpr)
# No Fig. Learn objects from expression of segmentation-based objects
TBL_obj1_alt <- learn_objects(input_data = TBL_obj1_xpr, input_type = "T", k_value = 20) # Produces error!!
head(TBL_obj1_alt[[1]])
# Fig. 13. Scalar image of CXCL13 (blue), CXCL12 (red), and Fibronectin (yellow), to be used independently
IMG_scl3 <- load_image(in_file = "HuLN_Scalar3.tif", img_type = "S", num_chs = 3)
PAL_scl3 <- c("blue", "red", "yellow")
plot_image(img = IMG_scl3, img_type = "S", col_pal = PAL_scl3, slice = c("Z" = 1))
# Fig. 14. Scalar image of Ki67 (green) and PD-1 (magenta), to be linked to objects
IMG_scl1 <- load_image(in_file = "HuLN_Scalar1.tif", img_type = "S", num_chs = 2)
PAL_scl1 <- c("green", "magenta")
plot_image(img = IMG_scl1, img_type = "S", col_pal = PAL_scl1, slice = c("Z" = 1))
# NO Fig. Import link table for scalar set 1 and segmentation-based objects
TBL_link <- load_table("HuLN_LinkTable.csv", table_type = "L")
TBL_link
##### Census Input Data ########################################################
# Fig. 15. Translation of desired radii in microns to radii in pixels
PAR_radii <- suggest_radii(target = seq(10, 50, by=10), pix_res = c(X=0.284, Y=0.284, Z=1))
# Fig. 16. Translation of desired coverage into required sample size
PAR_number <- suggest_number(coverage = 5, radii = PAR_radii, 
                             images = list(S1 = IMG_scl1, S3 = IMG_scl3, O1 = IMG_obj1, O2 = IMG_obj2))
# No Fig. Census the image data.
CEN_images <- census_image(images = list(O1 = IMG_obj1, O2 = IMG_obj2, S1 = IMG_scl1, S3 = IMG_scl3),
                           OS_pairs = list(O1 = TBL_link), radii = PAR_radii, sample_size = PAR_number)
PLS_images <- CEN_images[[2]]
CEN_images <- CEN_images[[1]]
#save(CEN_images, file="CEN_images.Rdata")
#save(PLS_images, file="PLS_images.Rdata")
# No Fig. Census the table data
CEN_table <- census_table(object_table = TBL_obj1_obj, radii = seq(10,50,by=10), sample_size = PAR_number)
PLS_table <- CEN_table[[2]]
CEN_table <- CEN_table[[1]]
#save(CEN_table, file="CEN_table.Rdata")
#save(PLS_table, file="PLS_table.Rdata")
##### Visualize a Census Directly ##############################################
# Fig. 17. Distribution of GC B cell (O1.6) abundance among 10-micron neighborhoods
plot_dist(census = CEN_images, ensemble = "O1.6", radius = 10, bin_num = 11)
# Fig. 18. Joint distribution of GC B cells (O1.6) and FDC network (O2.5) across observed 10-micron neighborhoods
plot_dist(census = CEN_images, ensemble = c("O1.6","O2.5"), radius = 10, bin_num = 11)
# Fig. 19. Scatter plot of GC B cell (O1.6) and FDC network (O2.5) across observed 10-micron neighborhoods
plot_dist(census = CEN_images, ensemble = c("O1.6","O2.5"), radius = 10)
# Fig. 20. Joint distribution of GC B cells (O1.6) and FDC network (O2.5) across randomized 10-micron neighborhoods
plot_dist(census = CEN_images, ensemble = c("O1.6","O2.5"), radius = 10, bin_num = 11, patch_list = PLS_images)
##### Discover Spatial Patterns within Single Biological Samples ###############
MIS_images_10um <- measure_cisMI(census = CEN_images, patch_list = PLS_images, depth = 3, radii = 10, 
                                 bootstraps = 100)
#save(MIS_images_10um, file="MIS_images_10um.Rdata")
MIS_images_all <- measure_cisMI(census = CEN_images, patch_list = PLS_images, depth = 3, radii = NULL, 
                                bootstraps = 100)
#save(MIS_images_all, file="MIS_images_all.Rdata")
# Fig. 21. Significant ensembles at a length scale of 10um, ranked by P value.
# Fig. 22. Variables involved in significant ensembles at a length scale of 10um, ranked by average Z score.
MIS_plot_20um <- plot_MI_rank(mi = MIS_images_all, radius = 20, 
                              col_pals = list(O1=PAL_obj1_alt, O2=PAL_obj2, S1=PAL_scl1, S3=PAL_scl3))
# Fig. 23. Significance and magnitude of cisMI for O1.8 and O2.5 across length scales.
MIS_plot_ens <- plot_MI_radius(mi = MIS_images_all, ensemble = c("O1.8","O2.5"))
##### Describe Spatial Patterns within Single Biological Samples ###############
# Fig. 24. Covariation plot for O1.8 and O2.5 at a length scale of 20um.
CVP_ens <- learn_pattern(census = CEN_images, ensemble = c("O1.8", "O2.5"), radius = 20,
                         col_pal = list(O1 = PAL_obj1_alt, O2 = PAL_obj2), patch_list = PLS_images) 
# Fig. 25. New object image of the broad zonation created by O1.8 and O2.5.
MAP_img <- map_pattern(covar_data = CVP_ens, region_bounds = list(c(0,15),c(15,55),c(55,67),c(67,100)),
                       img = list(O1 = IMG_obj1, O2 = IMG_obj2), census = CEN_images, radius = 20, radii = PAR_radii)
MAP_pal <- MAP_img[[2]]
MAP_img <- MAP_img[[1]]
#save(MAP_img, file="MAP_img.Rdata")
#save(MAP_pal, file="MAP_pal.Rdata")
##### Discover Spatial Patterns that Distinguish Groups of Biological Samples ##
load("TB_CEN.Rdata")
load("TB_PLS.Rdata")
load("TB_CPM.Rdata")
load("TB_MPM.Rdata")
load("TB_TIF.Rdata")
load("TB_CPM_NMS.Rdata")
load("TB_CPM_PAL.Rdata")
load("TB_MPM_NMS.Rdata")
load("TB_MPM_PAL.Rdata")
load("TB_TIF_NMS.Rdata")
load("TB_TIF_PAL.Rdata")
load("TB_PAR_RAD.Rdata")
load("TB_PAR_NUM.Rdata")
load("TB_TMI.Rdata")
TB_GRP <- data.frame(Status = c(rep("Resection",6), rep("Postmortem",6), rep("Biopsy",18)))
TB_TMI <- measure_transMI(censuses = TB_CEN, groups = TB_GRP, depth = 3, radii = 10)
#save(TB_TMI, file = "TB_TMI.Rdata")
# Fig. 26. Ensembles that distinguish among sample groups at a length scale of 10um, ranked by P value
# Fig. 27. Variables involved in distinguishing ensembles at a length scale of 10um, ranked by average Z score
TB_TMI_plot <- plot_MI_rank(mi = TB_TMI, radius = 10, col_pals = list(O1=TB_CPM_PAL, S1=TB_TIF_PAL), p_thr = 1e-7,
                            group = "Status")
##### Describe Spatial Patterns that Distinguish Groups of Biological Samples ##
# Fig. 28. Covariation plot for O1.7_S1.3, O1.10, and O1.13_S1.9 at length scale 10um, comparing biopsy versus 
# resection and post-mortem samples.
TB_CVP <- learn_pattern(census = TB_CEN, ensemble = c("O1.7_S1.3", "O1.10", "O1.13_S1.9"), radius = 10, 
                        col_pal = list(O1 = TB_CPM_PAL, S1 = TB_TIF_PAL), group = TB_GRP, 
                        focal = c("Resection", "Postmortem"), reference = "Biopsy", smooth_window=500)
##### Measure Diversity of Spatial Elements ####################################
# Fig. 29. Alpha diversity of segmented cell types in the 27th TB granuloma.
AD <- alpha_diversity(TB_CPM[[27]], "O", list(O1 = TB_CPM_PAL))
# Fig. 30. Beta diversity of microenvironments with respect to the segmented cell types in the 27th TB granuloma.
BD <- beta_diversity(list(O2 = TB_MPM[[27]], O1 = TB_CPM[[27]]), "O", list(O2 = TB_MPM_PAL, O2 = TB_CPM_PAL))