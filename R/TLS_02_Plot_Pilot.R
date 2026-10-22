# TLS plot pilot and candidate stem detection
# Prepared for a LAZ containing 290,336,075 points over approximately 30 x 30 m.
# Assumption: X, Y and Z use metres. Verify this in the scan/export settings.
# This is a diagnostic pilot, not a final trunk extraction or DBH workflow.
# The input LAZ is never overwritten. No EPSG code is assigned.
#
# HOW TO RUN
# 1. Open this file in RStudio and click Source. Select the plot LAZ when asked.
#    This loads a sampled central crop, prints diagnostics, and saves a checkpoint.
# 2. In the Console, run:
#      tls_ground <- tls_classify_pilot_ground(tls_pilot)
#    Inspect the displayed ground points. They should follow the forest floor.
# 3. Once the ground classification looks sensible, AND if input Z is not already
#    height above ground, run:
#      tls_stems <- tls_detect_pilot_stems(tls_ground, tls_output_dir)
#    This normalizes height, detects candidate stems, and saves an overlay PNG.
#    If your LAZ is already normalized, stop after step 2 and report that fact;
#    do not normalize it again without reviewing its height reference.
#
# A small crop can truncate crowns and contain branches from trees outside it.
# Its tree count, heights and candidate positions are not a plot inventory.
# Random sampling is for this diagnostic pilot only. Final stem fitting should
# revisit the original stem neighbourhoods at higher resolution.
#
# Documentation checked 2026-09-24:
# https://r-lidar.r-universe.dev/lidR/doc/manual.html
# https://github.com/JulFrey/CspStandSegmentation

if (!requireNamespace("lidR", quietly = TRUE)) {
  stop("lidR is not installed in this R library.")
}
library(lidR)

# ---- Configuration -----------------------------------------------------------
# To avoid the file chooser, replace NA_character_ with "data/TLS/your_plot.laz".
tls_input <- NA_character_
tls_roi <- c(xmin = 420402.6, ymin = 6481810,
             xmax = 420412.6, ymax = 6481820)
tls_sample_fraction <- 0.02  # retain about 2% within the crop during reading
tls_voxel_m <- 0.02         # then retain one observed point per 2 cm voxel

# ---- Helper functions --------------------------------------------------------
tls_plot_3d <- function(cloud, colour = "Z") {
  # Limit the graphical workload independently of the analysis cloud.
  set.seed(260924)
  view_cloud <- lidR::decimate_points(
    cloud, lidR::random_per_voxel(res = 0.05, n = 1)
  )
  tryCatch(
    plot(view_cloud, color = colour, size = 2),
    error = function(e) message("3D display unavailable: ", conditionMessage(e))
  )
  invisible(NULL)
}

tls_classify_pilot_ground <- function(cloud, cloth_m = 0.20,
                                      ground_tolerance_m = 0.10,
                                      rigidness = 2L) {
  if (!requireNamespace("RCSF", quietly = TRUE)) {
    stop("RCSF is required. Install it with install.packages('RCSF').")
  }
  # These are trial settings. Inspect the result before height normalization.
  classified <- lidR::classify_ground(
    cloud,
    lidR::csf(cloth_resolution = cloth_m,
              class_threshold = ground_tolerance_m,
              rigidness = rigidness),
    last_returns = FALSE
  )
  ground <- lidR::filter_poi(classified, Classification == 2L)
  if (lidR::is.empty(ground)) {
    stop("No ground was detected. Review the cloud and ground settings.")
  }
  cat("Ground points:", lidR::npoints(ground), "\n")
  cat("Inspect these points: they should follow the forest floor.\n")
  tls_plot_3d(ground)
  classified
}

tls_detect_pilot_stems <- function(classified, output_dir,
                                   zmin_m = 0.5, zmax_m = 2,
                                   raster_m = 0.10,
                                   density_quantile = 0.975,
                                   merge_radius_m = 0.20) {
  if (!requireNamespace("CspStandSegmentation", quietly = TRUE)) {
    stop("CspStandSegmentation is not installed. Complete its installation first.")
  }
  if (!any(classified$Classification == 2L, na.rm = TRUE)) {
    stop("Run and inspect tls_classify_pilot_ground() first.")
  }
  # Rescale this processing copy to millimetre coordinate storage. The original
  # LAZ and the saved pilot_input.rds retain their original coordinates.
  work <- lidR::las_rescale(classified, 0.001, 0.001, 0.001)
  normalised <- lidR::normalize_height(work, lidR::tin())
  # lidR retains pre-normalization Z in Zref. CspStandSegmentation recognizes it.
  if (!"Zref" %in% names(normalised@data)) {
    stop("Normalization did not retain Zref; inspect the installed lidR version.")
  }
  band <- lidR::filter_poi(normalised, Z > zmin_m & Z < zmax_m)
  if (lidR::is.empty(band)) {
    stop("No points in the stem detection band. Review the ground model.")
  }
  seeds <- CspStandSegmentation::find_base_coordinates_raster(
    normalised, res = raster_m, zmin = zmin_m, zmax = zmax_m,
    q = density_quantile, eps = merge_radius_m
  )
  seeds <- as.data.frame(seeds)
  if (nrow(seeds) == 0L) {
    stop("No candidate stems found. Review the height band and detection settings.")
  }
  required <- c("X", "Y", "Z", "TreeID")
  if (!all(required %in% names(seeds))) {
    stop("Unexpected stem-map columns: ", paste(names(seeds), collapse = ", "))
  }
  if (any(!is.finite(as.matrix(seeds[, required])))) {
    stop("Non-finite values in the stem map. Inspect ground support and detection.")
  }
  # Z is the algorithm's seed height, not DBH or tree height.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(seeds, file.path(output_dir, "candidate_stem_positions.csv"),
            row.names = FALSE)
  saveRDS(normalised, file.path(output_dir, "pilot_normalised.rds"),
          compress = FALSE)
  saveRDS(seeds, file.path(output_dir, "candidate_stem_positions.rds"))

  set.seed(260924)
  ids <- sample.int(lidR::npoints(band),
                    min(100000L, lidR::npoints(band)))
  draw_overlay <- function() {
    plot(band$X[ids], band$Y[ids], asp = 1, pch = 16, cex = 0.20,
         col = "grey60", xlab = "X (m)", ylab = "Y (m)",
         main = "Candidate stems in the 0.5-2 m height band")
    points(seeds$X, seeds$Y, pch = 3, col = "red3", cex = 1.4, lwd = 2)
    text(seeds$X, seeds$Y, labels = seeds$TreeID, pos = 3,
         cex = 0.8, col = "red3")
  }
  png_path <- file.path(output_dir, "candidate_stem_map.png")
  grDevices::png(png_path, width = 1600, height = 1400, res = 180)
  tryCatch(draw_overlay(), finally = grDevices::dev.off())
  draw_overlay()
  cat("Candidate stems:", nrow(seeds), "\n")
  print(seeds)
  cat("Overlay saved to:", png_path, "\n")
  cat("Review omissions, duplicate seeds and false detections before segmentation.\n")
  invisible(list(cloud = normalised, stem_positions = seeds))
}

# ---- Read and inspect a pilot -------------------------------------------------
if (is.na(tls_input) || !nzchar(tls_input)) tls_input <- file.choose()
tls_input <- normalizePath(tls_input, winslash = "/", mustWork = TRUE)
if (!grepl("\\.la[sz]$", tls_input, ignore.case = TRUE)) {
  stop("Select a LAS or LAZ point cloud.")
}
stopifnot(tls_sample_fraction > 0, tls_sample_fraction <= 1, tls_voxel_m > 0,
          tls_roi["xmax"] > tls_roi["xmin"],
          tls_roi["ymax"] > tls_roi["ymin"])
tls_output_dir <- tempfile(pattern = "TLS_pilot_", tmpdir = dirname(tls_input))
if (!dir.create(tls_output_dir, recursive = TRUE)) {
  stop("Cannot create the output folder beside the LAZ. Choose a writable folder.")
}
tls_filter <- sprintf(
  "-inside %.3f %.3f %.3f %.3f -keep_random_fraction %.5f",
  tls_roi["xmin"], tls_roi["ymin"], tls_roi["xmax"], tls_roi["ymax"],
  tls_sample_fraction
)
cat("Reading a sampled crop. A large LAZ may still take time to scan.\n")
tls_pilot <- lidR::readLAS(tls_input, select = "xyzc", filter = tls_filter)
if (lidR::is.empty(tls_pilot)) {
  stop("The crop is empty. Check the selected file and its coordinate bounds.")
}
tls_count_before <- lidR::npoints(tls_pilot)
set.seed(260924)
tls_pilot <- lidR::decimate_points(
  tls_pilot, lidR::random_per_voxel(res = tls_voxel_m, n = 1)
)
saveRDS(tls_pilot, file.path(tls_output_dir, "pilot_input.rds"), compress = FALSE)
tls_diagnostics <- capture.output({
  cat("Input:", tls_input, "\n")
  cat("Read filter:", tls_filter, "\n")
  cat("Sampled points before voxel thinning:", tls_count_before, "\n")
  cat("Pilot points after voxel thinning:", lidR::npoints(tls_pilot), "\n")
  cat("Pilot object size:", format(object.size(tls_pilot), units = "auto"), "\n")
  cat("Input Z quantiles (height reference not yet established):\n")
  print(quantile(tls_pilot$Z, c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)))
  cat("Existing classification counts:\n")
  print(table(tls_pilot$Classification, useNA = "ifany"))
  cat("CRS reported by sf:\n")
  print(sf::st_crs(tls_pilot))
  cat("lidR version:", as.character(packageVersion("lidR")), "\n")
  cat("CspStandSegmentation available:",
      requireNamespace("CspStandSegmentation", quietly = TRUE), "\n")
})
writeLines(tls_diagnostics, file.path(tls_output_dir, "pilot_diagnostics.txt"))
cat(paste(tls_diagnostics, collapse = "\n"), "\n")
capture.output(sessionInfo(), file = file.path(tls_output_dir, "sessionInfo.txt"))
tls_plot_3d(tls_pilot)
cat("\nPilot outputs:", tls_output_dir, "\n")
cat("Next: tls_ground <- tls_classify_pilot_ground(tls_pilot)\n")
cat("Inspect ground before running tls_detect_pilot_stems().\n")
# To resume in a later R session, source the function definitions above and use:
# tls_pilot <- readRDS("data/TLS/TLS_pilot_xxxx/pilot_input.rds")
# LASlib's random read filter is not seeded by set.seed(). Saving pilot_input.rds
# makes the exact pilot reusable; do not assume a fresh import is identical.
