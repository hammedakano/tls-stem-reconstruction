# Inspect the rejected cross sections from step 4 without changing the extraction.
# Open in RStudio and click Source, or use source(file.choose()).
# Uses the saved trunk_slice_checks.csv and trace_settings.rds together.
# Uses tls_trees_v2 if available; otherwise select pilot_tree_ids_v2.rds.
# Export includes each rejected slice and its preceding accepted reference slice.
# X/Y are projected coordinates in metres; HAG_m is height above ground.
# InsideOriginalSearch identifies the step 4 candidate pool on rejected rows.
# Reference rows use the same context aperture and have NA search flags/counts.
# Additional nearby points provide context. They are not newly accepted stems.
# A row with NA coordinates explicitly records an empty context window.
# R execution on the user's pilot is required; no fitting is performed here.

if (!requireNamespace("lidR", quietly=TRUE)) stop("lidR is required.")

tls05_checks_file <- if (exists("tls_trunk_output_dir",inherits=TRUE)) {
  file.path(tls_trunk_output_dir,"trunk_slice_checks.csv")
} else ""
if (!file.exists(tls05_checks_file)) {
  message("Select trunk_slice_checks.csv from the step 4 output folder.")
  tls05_checks_file <- file.choose()
}
tls05_parent <- dirname(tls05_checks_file)
tls05_settings_file <- file.path(tls05_parent,"trace_settings.rds")
if (!file.exists(tls05_settings_file)) {
  stop("trace_settings.rds must be beside trunk_slice_checks.csv.")
}
tls05_models <- read.csv(tls05_checks_file)
tls05_p <- readRDS(tls05_settings_file)
if (!all(c("TreeID","Bin","Height_m","Accepted","X","Y","Radius_m") %in%
         names(tls05_models))) stop("Unexpected slice-check columns.")
if (!is.logical(tls05_models$Accepted)) stop("Accepted must contain TRUE/FALSE.")
tls05_failed <- tls05_models[!is.na(tls05_models$Accepted) & !tls05_models$Accepted,,drop=FALSE]
if (!nrow(tls05_failed)) stop("There are no rejected slices to inspect.")

tls05_las <- if (exists("tls_trees_v2",inherits=TRUE) && inherits(tls_trees_v2,"LAS")) {
  tls_trees_v2
} else {
  message("Select pilot_tree_ids_v2.rds from the step 3 output folder.")
  readRDS(file.choose())
}
if (!inherits(tls05_las,"LAS") ||
    !all(c("TreeID","Zref") %in% names(tls05_las@data))) {
  stop("Use the normalized, seven-ID pilot_tree_ids_v2.rds.")
}
tls05_d <- data.frame(X=tls05_las$X,Y=tls05_las$Y,
                       HAG_m=tls05_las$Z,TreeID=tls05_las$TreeID)
if (any(!is.finite(as.matrix(tls05_d[,c("X","Y","HAG_m")])))) {
  stop("Non-finite point coordinates.")
}
if (min(tls05_d$X)<420402.55 || max(tls05_d$X)>420412.65 ||
    min(tls05_d$Y)<6481809.95 || max(tls05_d$Y)>6481820.05) {
  stop("These checks apply to the central 10 x 10 m pilot.")
}
tls05_d$Bin <- floor((tls05_d$HAG_m-(tls05_p$anchor_h-tls05_p$step/2))/tls05_p$step+1e-8)

tls05_points <- list()
tls05_index <- list()
for (i in seq_len(nrow(tls05_failed))) {
  rejected <- tls05_failed[i,,drop=FALSE]
  # The tracing algorithm starts at the anchor and walks outwards in height.
  reference <- tls05_models[which(tls05_models$Accepted &
    tls05_models$TreeID == rejected$TreeID),,drop=FALSE]
  if (rejected$Height_m > tls05_p$anchor_h) {
    reference <- reference[reference$Height_m < rejected$Height_m,,drop=FALSE]
    reference <- reference[order(reference$Height_m,decreasing=TRUE),,drop=FALSE]
  } else if (rejected$Height_m < tls05_p$anchor_h) {
    reference <- reference[reference$Height_m > rejected$Height_m,,drop=FALSE]
    reference <- reference[order(reference$Height_m),,drop=FALSE]
  } else stop("An anchor fit failed; inspect that tree separately.")
  if (!nrow(reference)) stop("No preceding accepted fit for tree ",rejected$TreeID)
  ref <- reference[1,,drop=FALSE]
  nstep <- max(1,abs(rejected$Height_m-ref$Height_m)/tls05_p$step)
  shift <- tls05_p$center_step*nstep
  rlo <- max(tls05_p$min_radius,ref$Radius_m-tls05_p$radius_step*nstep)
  rhi <- min(tls05_p$max_radius,ref$Radius_m+tls05_p$radius_step*nstep)
  search_r <- rhi+shift+tls05_p$tolerance
  context_r <- max(0.60,search_r+0.20)
  case_id <- sprintf("T%02d_B%03d",rejected$TreeID,rejected$Bin)
  for (role in c("rejected","reference")) {
    target <- if (role=="rejected") rejected else ref
    rows <- which(tls05_d$TreeID==target$TreeID & tls05_d$Bin==target$Bin &
      (tls05_d$X-ref$X)^2+(tls05_d$Y-ref$Y)^2 <= context_r^2)
    pts <- tls05_d[rows,c("X","Y","HAG_m"),drop=FALSE]
    original_search <- if (nrow(pts)) {
      (pts$X-ref$X)^2+(pts$Y-ref$Y)^2 <= search_r^2
    } else logical()
    n_context <- nrow(pts)
    n_search <- if (role=="rejected") sum(original_search) else NA_integer_
    if (!nrow(pts)) {
      pts <- data.frame(X=NA_real_,Y=NA_real_,HAG_m=NA_real_)
      original_search <- NA
    }
    if (role=="reference") original_search <- rep(NA,nrow(pts))
    # Limits below describe the rejected fit, including on reference-context rows.
    metadata <- data.frame(CaseID=case_id,TreeID=rejected$TreeID,Role=role,
      RejectedHeight_m=rejected$Height_m,SliceHeight_m=target$Height_m,
      Bin=target$Bin,SliceWidth_m=tls05_p$step,
      RefHeight_m=ref$Height_m,RefX=ref$X,RefY=ref$Y,RefRadius_m=ref$Radius_m,
      CenterShiftLimit_m=shift,RadiusMin_m=rlo,RadiusMax_m=rhi,
      FitTolerance_m=tls05_p$tolerance,MinInliers=tls05_p$min_points,
      MinInlierFraction=tls05_p$min_fraction,MinSectorsOf12=tls05_p$min_sectors,
      MaxEmptyAngle_deg=tls05_p$max_empty_angle,MaxRMSE_m=tls05_p$max_rmse,
      SearchRadius_m=search_r,ContextRadius_m=context_r,
      ContextPoints=n_context,SearchPoints=n_search)
    block <- cbind(metadata[rep(1L,nrow(pts)),,drop=FALSE],
                   InsideOriginalSearch=original_search,pts)
    tls05_points[[length(tls05_points)+1L]] <- block
    tls05_index[[length(tls05_index)+1L]] <- metadata
  }
}
tls05_out <- tempfile("rejected_slices_",tmpdir=tls05_parent)
if (!dir.create(tls05_out)) stop("Cannot create diagnostic folder.")
tls05_export <- do.call(rbind,tls05_points)
tls05_inventory <- do.call(rbind,tls05_index)
write.csv(tls05_export,file.path(tls05_out,"rejected_slice_points.csv"),row.names=FALSE)
write.csv(tls05_inventory,file.path(tls05_out,"rejected_slice_index.csv"),row.names=FALSE)
print(tls05_inventory[tls05_inventory$Role=="rejected",
                     c("CaseID","TreeID","RejectedHeight_m","SearchPoints","ContextPoints")],
      row.names=FALSE)
cat("\nRejected slices exported:",nrow(tls05_failed),"\n")
cat("Saved to:",tls05_out,"\n")
cat("Upload rejected_slice_points.csv for inspection.\n")
