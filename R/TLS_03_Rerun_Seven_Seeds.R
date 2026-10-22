# Repair the seven visible stem seeds in the plot167 central pilot.
# In RStudio: open this file and click Source, or run source(file.choose()).
# Uses tls_geometry if it is in the current session. Otherwise uses the saved
# pilot_normalised.rds checkpoint. It never reads the full plot LAZ.
#
# The seven centers were fitted to TLS_DBH_Slice_ID_Check.csv (2,834 points at
# 1.20-1.40 m HAG). XY clusters were reviewed against candidate_stem_map.png;
# RANSAC circle fits used a 0.02 m radial tolerance followed by robust refinement.
# IDs 1-4 retain their identities; 5 = northwest, 6 = northeast, 7 = west.
# Z = 1.30 m ABOVE GROUND, matching the already-normalized processing cloud.
# These are seed estimates for this pilot. They are not a full-plot inventory.
# The circle radii below support ID checks only; they are not validated DBH.
#
# References checked 2026-09-24:
# https://github.com/JulFrey/CspStandSegmentation
# https://r-lidar.r-universe.dev/lidR/doc/manual.html
# The CSV fits were computed and checked. This R script was not executed here.

for (pkg in c("lidR", "CspStandSegmentation", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is missing: ", pkg)
  }
}
library(lidR)

# ---- Corrected seeds ---------------------------------------------------------
# Coordinates are embedded so no separate seed-file selection is necessary.
tls_seeds_v2 <- read.csv(text = "X,Y,Z,TreeID
420409.339,6481818.000,1.300,1
420410.381,6481816.073,1.300,2
420411.695,6481813.066,1.300,3
420406.740,6481811.905,1.300,4
420403.218,6481819.670,1.300,5
420410.515,6481819.733,1.300,6
420404.018,6481814.116,1.300,7")
tls_seeds_v2$TreeID <- as.integer(tls_seeds_v2$TreeID)
tls_reference_radii <- c(0.1261, 0.1054, 0.1352, 0.1584,
                         0.2225, 0.1480, 0.1827)

# ---- Resume the pilot --------------------------------------------------------
tls_parent_dir_v2 <- if (exists("tls_output_dir", inherits = TRUE)) {
  tls_output_dir
} else {
  Sys.getenv("TLS_PILOT_DIR", unset = "data/TLS/TLS_pilot")
}

if (exists("tls_geometry", inherits = TRUE) && inherits(tls_geometry, "LAS")) {
  # Deep copy because some R/data.table operations can modify objects by reference.
  tls_geometry_v2 <- data.table::copy(tls_geometry)
  message("Reusing the pilot geometry from this R session.")
} else {
  tls_resume_path_v2 <- file.path(tls_parent_dir_v2, "pilot_normalised.rds")
  if (!file.exists(tls_resume_path_v2)) {
    message("Select pilot_normalised.rds from your saved pilot folder.")
    tls_resume_path_v2 <- file.choose()
    tls_parent_dir_v2 <- dirname(tls_resume_path_v2)
  }
  tls_resume_cloud_v2 <- readRDS(tls_resume_path_v2)
  if (!inherits(tls_resume_cloud_v2, "LAS") ||
      !all(c("Classification", "Zref") %in% names(tls_resume_cloud_v2@data))) {
    stop("Select the saved, ground-classified and normalized pilot LAS object.")
  }
  tls_geometry_v2 <- lidR::filter_poi(
    tls_resume_cloud_v2, Classification != 2L & Z > 0.20
  )
  rm(tls_resume_cloud_v2)
  if (lidR::is.empty(tls_geometry_v2)) stop("The pilot vegetation cloud is empty.")
  tls_geometry_v2 <- CspStandSegmentation::add_geometry(
    tls_geometry_v2, k = 10L, n_cores = 1
  )
}

if (lidR::is.empty(tls_geometry_v2)) stop("The pilot geometry is empty.")
if (!"Zref" %in% names(tls_geometry_v2@data)) {
  stop("This workflow expects the normalized pilot with its Zref column.")
}
# Catch accidental use on a different plot or the full cloud.
if (min(tls_geometry_v2$X) < 420402.55 ||
    max(tls_geometry_v2$X) > 420412.65 ||
    min(tls_geometry_v2$Y) < 6481809.95 ||
    max(tls_geometry_v2$Y) > 6481820.05) {
  stop("The cloud does not match the 10 x 10 m pilot used to fit these seeds.")
}
if (!dir.exists(tls_parent_dir_v2)) {
  stop("Pilot output folder not found: ", tls_parent_dir_v2)
}
for (nm in intersect(c("TreeID", "treeID", "x_vox", "y_vox", "z_vox"),
                     names(tls_geometry_v2@data))) {
  tls_geometry_v2 <- lidR::remove_lasattribute(tls_geometry_v2, nm)
}

# ---- Rerun the same segmentation with seven seeds -----------------------------
tls_output_dir_v2 <- tempfile("seven_seeds_", tmpdir = tls_parent_dir_v2)
if (!dir.create(tls_output_dir_v2)) stop("Cannot create the rerun output folder.")
write.csv(tls_seeds_v2,
          file.path(tls_output_dir_v2, "corrected_stem_positions.csv"),
          row.names = FALSE)
cat("Seed positions for this rerun:\n")
print(tls_seeds_v2, digits = 10, row.names = FALSE)

tls_trees_v2 <- CspStandSegmentation::csp_cost_segmentation(
  tls_geometry_v2,
  map = tls_seeds_v2,
  Voxel_size = 0.30,
  V_w = 0.50,
  N_cores = 1
)
if (!inherits(tls_trees_v2, "LAS") ||
    !"TreeID" %in% names(tls_trees_v2@data)) {
  stop("Segmentation did not return a LAS object with TreeID.")
}
saveRDS(tls_trees_v2, file.path(tls_output_dir_v2, "pilot_tree_ids_v2.rds"),
        compress = FALSE)
saveRDS(tls_seeds_v2, file.path(tls_output_dir_v2, "corrected_stem_positions.rds"))

# ---- Check actual point IDs around each fitted ring ---------------------------
tls_slice_v2 <- lidR::filter_poi(tls_trees_v2, Z >= 1.20 & Z <= 1.40)
if (lidR::is.empty(tls_slice_v2)) stop("No points remain in the 1.20-1.40 m band.")
tls_slice_df_v2 <- data.frame(
  X = tls_slice_v2$X, Y = tls_slice_v2$Y,
  Height_m = tls_slice_v2$Z, TreeID = tls_slice_v2$TreeID
)
write.csv(tls_slice_df_v2,
          file.path(tls_output_dir_v2, "TLS_DBH_Slice_ID_Check_v2.csv"),
          row.names = FALSE)

tls_id_check_v2 <- do.call(rbind, lapply(seq_len(nrow(tls_seeds_v2)), function(i) {
  radial_distance <- sqrt((tls_slice_v2$X - tls_seeds_v2$X[i])^2 +
                          (tls_slice_v2$Y - tls_seeds_v2$Y[i])^2)
  # Use a 3 cm annulus around the reviewed ring, not every branch in the cluster.
  observed <- tls_slice_v2$TreeID[
    abs(radial_distance - tls_reference_radii[i]) <= 0.03
  ]
  counts <- table(observed, useNA = "ifany")
  n <- length(observed)
  expected_fraction <- if (n > 0L) {
    sum(observed == tls_seeds_v2$TreeID[i], na.rm = TRUE) / n
  } else NA_real_
  dominant <- if (length(counts)) {
    suppressWarnings(as.integer(names(counts)[which.max(counts)]))
  } else NA_integer_
  data.frame(
    ExpectedTreeID = tls_seeds_v2$TreeID[i], RingPoints = n,
    DominantTreeID = dominant, ExpectedFraction = expected_fraction,
    ObservedIDs = paste(names(counts), collapse = ";"),
    SliceCheck = if (n >= 20L && is.finite(expected_fraction) &&
                    expected_fraction >= 0.95) "consistent" else "inspect"
  )
}))
write.csv(tls_id_check_v2, file.path(tls_output_dir_v2, "stem_id_check_v2.csv"),
          row.names = FALSE)
cat("\nID check in the breast-height slice (not a crown or DBH validation):\n")
print(tls_id_check_v2, row.names = FALSE)

# ---- Plot with a separate color for each numeric ID ---------------------------
tls_palette_v2 <- c("0" = "grey75",
                    setNames(grDevices::hcl.colors(7, "Dark 3"), as.character(1:7)))
tls_point_col_v2 <- unname(tls_palette_v2[as.character(tls_slice_v2$TreeID)])
tls_point_col_v2[is.na(tls_point_col_v2)] <- "grey75"
tls_draw_slice_v2 <- function() {
  plot(tls_slice_v2$X, tls_slice_v2$Y, asp = 1,
       col = tls_point_col_v2, pch = 16, cex = 0.35,
       xlab = "X (m)", ylab = "Y (m)",
       main = "Rerun: point IDs at 1.20-1.40 m above ground")
  points(tls_seeds_v2$X, tls_seeds_v2$Y, pch = 3, cex = 1.3, lwd = 2)
  text(tls_seeds_v2$X, tls_seeds_v2$Y, labels = tls_seeds_v2$TreeID, pos = 3)
  legend("bottomleft", legend = paste("ID", 0:7),
         col = unname(tls_palette_v2), pch = 16, bty = "n", cex = 0.8)
}
grDevices::png(file.path(tls_output_dir_v2, "stem_ids_v2.png"),
               width = 1600, height = 1400, res = 180)
tryCatch(tls_draw_slice_v2(), finally = grDevices::dev.off())
if (interactive()) tls_draw_slice_v2()

capture.output(sessionInfo(), file = file.path(tls_output_dir_v2, "sessionInfo.txt"))
cat("\nSaved rerun outputs to:", tls_output_dir_v2, "\n")
cat("Expected: seven rings with dominant IDs 1-7 matching their seed labels.\n")
cat("Inspect any row marked 'inspect' before trunk extraction or DBH fitting.\n")
cat("Share stem_ids_v2.png and stem_id_check_v2.csv for the next step.\n")

if (interactive()) {
  tryCatch({
    set.seed(260924)
    tls_view_v2 <- lidR::decimate_points(
      tls_trees_v2, lidR::random_per_voxel(res = 0.05, n = 1)
    )
    plot(tls_view_v2, color = "TreeID", pal = unname(tls_palette_v2),
         breaks = seq(-0.5, 7.5, by = 1), NAcol = "grey75", size = 2)
  }, error = function(e) message("3D display unavailable; saved results are intact: ",
                                 conditionMessage(e)))
}
