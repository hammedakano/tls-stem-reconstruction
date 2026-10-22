# Step 6: continuity-guided circle fits with explicitly flagged coverage trials.
# Run with source(file.choose()), selecting THIS script.
# Uses tls_trees_v2 in memory or asks for pilot_tree_ids_v2.rds.
# Writes a NEW folder. Previous results and input cloud are preserved.
#
# Evidence: 21 uploaded rejected slices, 2026-09-24. Several trunk rings
# failed angular coverage despite small residuals. A previous-circle start
# recovered strict fits that the original unconstrained proposals missed.
#
# TrunkFitQuality: 1 = passes original checks; 2 = provisional coverage trial.
# Coverage trials keep >=20 inliers and >=35% ring fraction, allow >=7/12
# sectors and <=150 degree gap, but tighten RMSE to <=8 mm, centre movement
# to <=25 mm and radius change to <=20 mm per 0.2 m height step.
# These are pilot sensitivity settings, NOT validated measurement thresholds.
# Coverage fallback is disabled at/below 2 m HAG and after three consecutive
# coverage-only slices. A strict slice resets this consecutive counter.
# PathUsesCoverage flags downstream slices whose trace used a coverage trial.
# Inspect provisional segments before using them in an inventory.
# Set allow_coverage=FALSE below for the improved fitter with original gates.
# R is unavailable in the assistant workspace; run this in your R session.
# HOW TO RUN: open this file in RStudio and click Source.
# Uses tls_trees_v2 from step 3. If absent, select pilot_tree_ids_v2.rds.
#
# Method: start at 1.3 m HAG, fit circles in adjacent 0.2 m horizontal slices,
# and track their centers/radii upward and downward. Two consecutive failed
# slices stop a direction. Failed slices remain gaps; no points are invented.
# This is a conservative, custom pilot workflow, not a validated semantic
# classifier. Junction points may remain; leaning, irregular or occluded stems
# may stop early. Reported trace extent is NOT total tree or trunk height.
# Circle diameters are provisional horizontal diameters on the thinned pilot,
# not final DBH. Final measurement should revisit denser original neighborhoods
# and address ground reference, lean, buttresses and the measurement protocol.
#
# Processing RDS: Z = height above ground; Zref = pre-normalization elevation.
# Exported LAZ: Z = Zref; HeightAboveGround and TreeID are extra attributes.
# LAZ coordinate storage is 0.001 m. No CRS is guessed.
#
# Checked API references, 2026-09-24:
# https://julfrey.r-universe.dev/CspStandSegmentation/doc/manual.html
# https://rdrr.io/cran/CspStandSegmentation/man/ransac_circle_fit.html
# https://r-lidar.r-universe.dev/lidR/doc/manual.html
# R execution and extraction on the full pilot must be checked in your session.

for (pkg in c("lidR", "CspStandSegmentation", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
library(lidR)
if (!"ransac_circle_fit" %in% getNamespaceExports("CspStandSegmentation")) {
  stop("The installed CspStandSegmentation does not export ransac_circle_fit().")
}

# Trial settings for this sampled, 2 cm voxel pilot. Units are metres.
tls06_settings <- list(
  step = 0.20, anchor_h = 1.30, lowest_center = 0.30,
  tolerance = 0.02, keep_tolerance = 0.025,
  min_points = 20L, min_fraction = 0.35, min_sectors = 8L,
  max_empty_angle = 120, max_rmse = 0.012,
  min_radius = 0.04, max_radius = 0.60,
  center_step = 0.10, radius_step = 0.04,
  max_failures = 2L, attempts = 3L, iterations = 300L,
  allow_coverage = TRUE, review_min_h = 2,
  review_min_sectors = 7L, review_max_gap = 150,
  review_max_rmse = 0.008, review_center_step = 0.025,
  review_radius_step = 0.020, max_coverage_run = 3L
)
tls06_seeds <- data.frame(
  TreeID = 1:7,
  X = c(420409.339,420410.381,420411.695,420406.740,
        420403.218,420410.515,420404.018),
  Y = c(6481818.000,6481816.073,6481813.066,6481811.905,
        6481819.670,6481819.733,6481814.116),
  Radius = c(0.1261,0.1054,0.1352,0.1584,0.2225,0.1480,0.1827)
)

tls06_fit_strict <- function(d, bin, ref, p) {
  h <- p$anchor_h + bin * p$step
  nstep <- max(1, abs(h - ref$h) / p$step)
  shift_limit <- p$center_step * nstep
  rlo <- max(p$min_radius, ref$r - p$radius_step * nstep)
  rhi <- min(p$max_radius, ref$r + p$radius_step * nstep)
  sel <- which(d$Bin == bin &
    (d$x-ref$x)^2 + (d$y-ref$y)^2 <= (rhi+shift_limit+p$tolerance)^2)
  failure <- function(reason) list(ok = FALSE, reason = reason, bin = bin, h = h)
  if (length(sel) < p$min_points) return(failure("too_few_local_points"))
  xy <- cbind(d$x[sel]-ref$x, d$y[sel]-ref$y)
  best <- NULL
  radial <- function(cc) sqrt((xy[,1]-cc[1])^2+(xy[,2]-cc[2])^2)-cc[3]
  for (attempt in seq_len(p$attempts)) {
    # Current lidR exposes fit_circle(); retain a fallback for older builds.
    cc <- tryCatch({
      if ("fit_circle" %in% getNamespaceExports("lidR")) {
        fit <- lidR::fit_circle(cbind(xy,d$HAG[sel]),
          num_iterations=p$iterations,inlier_threshold=p$tolerance)
        as.numeric(c(fit$center_x,fit$center_y,fit$radius))
      } else {
        fit <- CspStandSegmentation::ransac_circle_fit(
          xy,n_iterations=p$iterations,distance_threshold=p$tolerance,
          min_inliers=p$min_points)
        as.numeric(fit$circle)
      }
    },error=function(e) numeric())
    if (length(cc) != 3L || any(!is.finite(cc))) next
    if (sqrt(sum(cc[1:2]^2)) > shift_limit || cc[3] < rlo || cc[3] > rhi) next
    # Refine on ring inliers only. Work near (0,0) for numeric stability.
    for (j in 1:2) {
      inside <- abs(radial(cc)) <= p$tolerance
      if (sum(inside) < p$min_points) break
      refined <- tryCatch(stats::optim(cc, function(v) {
        sum((sqrt((xy[inside,1]-v[1])^2+(xy[inside,2]-v[2])^2)-v[3])^2)
      }, method = "L-BFGS-B",
      lower = c(-shift_limit,-shift_limit,rlo),
      upper = c(shift_limit,shift_limit,rhi)), error = function(e) NULL)
      if (is.null(refined) || refined$convergence != 0L) break
      cc <- refined$par
    }
    if (sqrt(sum(cc[1:2]^2)) > shift_limit) next
    residual <- radial(cc)
    inside <- abs(residual) <= p$tolerance
    n <- sum(inside)
    if (n < p$min_points || n / nrow(xy) < p$min_fraction) next
    angles <- sort(atan2(xy[inside,2]-cc[2],xy[inside,1]-cc[1]) %% (2*pi))
    sectors <- length(unique(floor(angles / (2*pi) * 12)))
    empty_angle <- max(diff(c(angles,angles[1]+2*pi))) * 180/pi
    rmse <- sqrt(mean(residual[inside]^2))
    if (sectors < p$min_sectors || empty_angle > p$max_empty_angle ||
        rmse > p$max_rmse) next
    candidate <- list(ok = TRUE, reason = "accepted", bin = bin, h = h,
      x = cc[1]+ref$x, y = cc[2]+ref$y, r = cc[3],
      n = n, fraction = n/nrow(xy), sectors = sectors,
      empty_angle = empty_angle, rmse = rmse,
      keep = d$Row[sel[abs(residual) <= p$keep_tolerance]])
    if (is.null(best) || n > best$n || (n == best$n && rmse < best$rmse)) {
      best <- candidate
    }
  }
  if (is.null(best)) failure("circle_support_or_continuity_failed") else best
}


# First try the preceding circle as an initial model. Refinement uses only
# points within the original 2 cm annulus and does not manufacture returns.
tls06_fit_slice <- function(d, bin, ref, p) {
  h <- p$anchor_h + bin*p$step
  nstep <- max(1, abs(h-ref$h)/p$step)
  shift <- p$center_step*nstep
  rlo <- max(p$min_radius,ref$r-p$radius_step*nstep)
  rhi <- min(p$max_radius,ref$r+p$radius_step*nstep)
  sel <- which(d$Bin == bin &
    (d$x-ref$x)^2+(d$y-ref$y)^2 <= (rhi+shift+p$tolerance)^2)
  xy <- cbind(d$x[sel]-ref$x,d$y[sel]-ref$y)
  cc <- c(0,0,ref$r)
  radial <- function(v) sqrt((xy[,1]-v[1])^2+(xy[,2]-v[2])^2)-v[3]
  fit_error <- ""
  if(length(sel)>=6L) for(j in 1:3) {
    inside <- abs(radial(cc)) <= p$tolerance
    if(sum(inside)<6L) break
    opt <- tryCatch(stats::optim(cc,function(v) {
      sum((sqrt((xy[inside,1]-v[1])^2+(xy[inside,2]-v[2])^2)-v[3])^2)
    },method="L-BFGS-B",lower=c(-shift,-shift,rlo),upper=c(shift,shift,rhi),
    control=list(factr=1e5,pgtol=1e-9,ndeps=rep(1e-5,3))),
    error=function(e) {fit_error <<- conditionMessage(e); NULL})
    if(is.null(opt)) break
    if(opt$convergence != 0L) {fit_error <- paste("optim code",opt$convergence); break}
    cc <- opt$par
  }
  residual <- radial(cc)
  inside <- abs(residual)<=p$tolerance
  n <- sum(inside)
  fraction <- if(length(sel)) n/length(sel) else 0
  if(n) {
    angles <- sort(atan2(xy[inside,2]-cc[2],xy[inside,1]-cc[1]) %% (2*pi))
    sectors <- length(unique(floor(angles/(2*pi)*12)))
    gap <- max(diff(c(angles,angles[1]+2*pi)))*180/pi
    rmse <- sqrt(mean(residual[inside]^2))
  } else {sectors <- 0L; gap <- 360; rmse <- NA_real_}
  move <- sqrt(sum(cc[1:2]^2))
  common <- character()
  if(length(sel)<p$min_points) common <- c(common,"search_points_below_min")
  if(n<p$min_points) common <- c(common,"ring_points_below_min")
  if(fraction<p$min_fraction) common <- c(common,"ring_fraction_below_min")
  if(move>shift || cc[3]<rlo || cc[3]>rhi) common <- c(common,"continuity_limit")
  if(!is.finite(rmse) || rmse>p$max_rmse) common <- c(common,"radial_rmse_limit")
  if(nzchar(fit_error)) common <- c(common,"prior_fit_optim_error")
  strict_fail <- common
  if(sectors<p$min_sectors) strict_fail <- c(strict_fail,"occupied_sectors_below_min")
  if(gap>p$max_empty_angle) strict_fail <- c(strict_fail,"angular_gap_limit")
  out <- list(ok=FALSE,reason=paste(strict_fail,collapse=";"),bin=bin,h=h,
    x=cc[1]+ref$x,y=cc[2]+ref$y,r=cc[3],n=n,fraction=fraction,
    sectors=sectors,empty_angle=gap,rmse=rmse,search_n=length(sel),
    qc=0L,method="previous_circle",fit_error=fit_error,
    strict_fail=paste(strict_fail,collapse=";"),coverage_run=0L,
    path_review=isTRUE(ref$path_review),
    keep=d$Row[sel[abs(residual)<=p$keep_tolerance]])
  if(!length(strict_fail)) {
    out$ok <- TRUE; out$reason <- "accepted_strict"; out$qc <- 1L
    return(out)
  }
  # Retain the original RANSAC search as another opportunity for a strict fit.
  strict <- tls06_fit_strict(d,bin,ref,p)
  if(strict$ok) {
    strict$qc <- 1L; strict$reason <- "accepted_strict"
    strict$method <- "ransac"; strict$search_n <- length(sel)
    strict$fit_error <- ""; strict$strict_fail <- ""; strict$coverage_run <- 0L
    strict$path_review <- isTRUE(ref$path_review)
    return(strict)
  }
  run <- if(is.null(ref$coverage_run)) 0L else ref$coverage_run
  coverage_ok <- isTRUE(p$allow_coverage) && h>p$review_min_h &&
    run<p$max_coverage_run && !length(common) &&
    sectors>=p$review_min_sectors && gap<=p$review_max_gap &&
    rmse<=p$review_max_rmse && move<=p$review_center_step*nstep &&
    abs(cc[3]-ref$r)<=p$review_radius_step*nstep
  if(coverage_ok) {
    out$ok <- TRUE; out$qc <- 2L; out$reason <- "accepted_coverage_review"
    out$coverage_run <- run+1L; out$path_review <- TRUE
  } else {
    # These are the prior-guided model's measured failures. The old RANSAC
    # routine also failed to return a strict fit; its trial models were not logged.
    if(run>=p$max_coverage_run) out$reason <- paste(out$reason,"coverage_run_limit",sep=";")
    out$keep <- integer()
  }
  out
}

tls06_trace_one <- function(d, seed, p) {
  d$x <- d$X - seed$X
  d$y <- d$Y - seed$Y
  # Half-open, non-overlapping bins, anchored on [1.2,1.4).
  d$Bin <- floor((d$HAG-(p$anchor_h-p$step/2))/p$step + 1e-8)
  ref <- list(x=0,y=0,r=seed$Radius,h=p$anchor_h,coverage_run=0L,path_review=FALSE)
  anchor <- tls06_fit_slice(d, 0L, ref, p)
  fits <- list(anchor)
  stops <- c(down="anchor_failed",up="anchor_failed")
  if (anchor$ok) {
    kmin <- ceiling((p$lowest_center-p$anchor_h)/p$step - 1e-8)
    kmax <- max(d$Bin)
    walks <- list(down=if(kmin <= -1) seq.int(-1L,kmin) else integer(),
                  up=if(kmax >= 1) seq_len(kmax) else integer())
    for (direction in names(walks)) {
      ref <- anchor
      failed <- 0L
      stops[direction] <- "reached_search_limit"
      for (bin in walks[[direction]]) {
        fit <- tls06_fit_slice(d, bin, ref, p)
        fits[[length(fits)+1L]] <- fit
        if (fit$ok) {
          ref <- fit
          failed <- 0L
        } else {
          failed <- failed+1L
          if (failed >= p$max_failures) {
            stops[direction] <- paste0("stopped_after_",failed,"_failed_slices")
            break
          }
        }
      }
    }
  }
  models <- do.call(rbind,lapply(fits,function(f) data.frame(
    TreeID=seed$TreeID, Bin=f$bin, Height_m=f$h, Accepted=f$ok, Reason=f$reason,
    FitQuality=f$qc, FitMethod=f$method, PathUsesCoverage=f$path_review,
    StrictChecksFailed=f$strict_fail, FitError=f$fit_error, SearchPoints=f$search_n,
    # Rejected fitted coordinates are diagnostics, never accepted trunk models.
    X=f$x+seed$X,Y=f$y+seed$Y,Radius_m=f$r,RingPoints=f$n,
    RingFraction=f$fraction,RadialRMSE_m=f$rmse,
    OccupiedSectorsOf12=f$sectors,MaxEmptyAngle_deg=f$empty_angle
  )))
  kept <- unique(unlist(lapply(fits,function(f) if(f$ok) f$keep else integer()),
                        use.names=FALSE))
  heights <- d$HAG[match(kept,d$Row)]
  summary <- data.frame(TreeID=seed$TreeID, CandidatePoints=length(kept),
    AcceptedSlices=sum(models$Accepted), FailedSlices=sum(!models$Accepted),
    StrictSlices=sum(models$FitQuality==1L),
    CoverageReviewSlices=sum(models$FitQuality==2L),
    PathUsesCoverageSlices=sum(models$Accepted & models$PathUsesCoverage),
    TracedMinHAG_m=if(length(kept)) min(heights) else NA_real_,
    TracedMaxHAG_m=if(length(kept)) max(heights) else NA_real_,
    PilotDiameterAt1p3_cm=if(anchor$ok) 200*anchor$r else NA_real_,
    LowerStop=unname(stops["down"]), UpperStop=unname(stops["up"]))
  review_kept <- unique(unlist(lapply(fits,function(f)
    if(f$ok && f$qc==2L) f$keep else integer()),use.names=FALSE))
  dependent_kept <- unique(unlist(lapply(fits,function(f)
    if(f$ok && f$path_review) f$keep else integer()),use.names=FALSE))
  list(models=models,keep=kept,review_keep=review_kept,
       dependent_keep=dependent_kept,summary=summary)
}

tls06_draw_profiles <- function(d, keep, seeds, palette, review=NULL) {
  oldpar <- par(mfrow=c(2,7),mar=c(3.2,3.1,2.2,0.5),oma=c(1,1,3,0))
  on.exit(par(oldpar))
  set.seed(260924)
  for (axis in c("X","Y")) for (i in seq_len(nrow(seeds))) {
    id <- seeds$TreeID[i]
    rows <- which(d$TreeID == id)
    retained <- rows[keep[rows]]
    background <- if(length(rows)>8000L) sample(rows,8000L) else rows
    provisional <- if(is.null(review)) integer() else rows[review[rows]]
    offset <- seeds[[axis]][i]
    xr <- if(length(retained)) range(d[[axis]][retained]-offset)+c(-.4,.4) else c(-1,1)
    zr <- c(0,max(d$HAG[rows]))
    plot(d[[axis]][background]-offset,d$HAG[background],xlim=xr,ylim=zr,
      pch=16,cex=.15,col="grey80",xlab=paste(axis,"offset (m)"),
      ylab="Height above ground (m)",main=paste("Tree",id,axis,"view"))
    if(length(retained)) points(d[[axis]][retained]-offset,d$HAG[retained],
                               pch=16,cex=.22,col=palette[as.character(id)])
    if(length(provisional)) points(d[[axis]][provisional]-offset,d$HAG[provisional],
                                  pch=16,cex=.25,col="black")
  }
  mtext("Candidate trunks: grey context, tree colours strict fits, black provisional coverage",
        outer=TRUE,cex=1.05,line=1)
}

tls06_run <- function(las, seeds, parent_dir, p=tls06_settings, show_3d=interactive()) {
  if (!inherits(las,"LAS") ||
      !all(c("TreeID","Zref") %in% names(las@data))) {
    stop("Use the normalized pilot_tree_ids_v2.rds from step 3.")
  }
  if (min(las$X)<420402.55 || max(las$X)>420412.65 ||
      min(las$Y)<6481809.95 || max(las$Y)>6481820.05) {
    stop("These seeds apply only to the plot167 central 10 x 10 m pilot.")
  }
  if (!dir.exists(parent_dir)) stop("Output parent folder does not exist.")
  needed <- c("data","n_iterations","distance_threshold","min_inliers")
  if (!all(needed %in% names(formals(CspStandSegmentation::ransac_circle_fit)))) {
    stop("Unexpected ransac_circle_fit() arguments; report packageVersion('CspStandSegmentation').")
  }
  d <- data.frame(Row=seq_len(lidR::npoints(las)),X=las$X,Y=las$Y,
                  HAG=las$Z,Zref=las$Zref,TreeID=las$TreeID)
  if (any(!is.finite(as.matrix(d[,c("X","Y","HAG","Zref")])))) {
    stop("Non-finite coordinates in the processing cloud.")
  }
  if (!all(seeds$TreeID %in% d$TreeID)) stop("One or more of the seven tree IDs is missing.")
  out <- tempfile("trunk_trace_v3_",tmpdir=parent_dir)
  if(!dir.create(out)) stop("Cannot create output folder.")
  message("Saving this trace to: ",out)
  saveRDS(p,file.path(out,"trace_settings.rds"))
  write.csv(seeds,file.path(out,"trace_seeds.csv"),row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  results <- vector("list",nrow(seeds))
  keep <- rep(FALSE,nrow(d))
  review <- dependent <- rep(FALSE,nrow(d))
  set.seed(260924)
  for(i in seq_len(nrow(seeds))) {
    id <- seeds$TreeID[i]
    message("Tracing tree ",id," of ",nrow(seeds)," ...")
    tree <- d[which(d$TreeID == id),,drop=FALSE]
    results[[i]] <- tls06_trace_one(tree,seeds[i,,drop=FALSE],p)
    keep[results[[i]]$keep] <- TRUE
    review[results[[i]]$review_keep] <- TRUE
    dependent[results[[i]]$dependent_keep] <- TRUE
    # Save each completed trace so progress survives an interrupted later tree.
    saveRDS(results[[i]],file.path(out,sprintf("trace_tree_%02d.rds",id)))
  }
  models <- do.call(rbind,lapply(results,`[[`,"models"))
  models <- models[order(models$TreeID,models$Height_m),]
  summary <- do.call(rbind,lapply(results,`[[`,"summary"))
  write.csv(models,file.path(out,"trunk_slice_checks.csv"),row.names=FALSE)
  write.csv(summary,file.path(out,"trunk_trace_summary.csv"),row.names=FALSE)
  palette <- setNames(grDevices::hcl.colors(7,"Dark 3"),as.character(1:7))
  grDevices::png(file.path(out,"trunk_profiles.png"),width=3200,height=1800,res=160)
  tryCatch(tls06_draw_profiles(d,keep,seeds,palette,review),finally=grDevices::dev.off())
  if(interactive()) tls06_draw_profiles(d,keep,seeds,palette,review)

  trunks <- NULL
  export_status <- list()
  if(any(keep)) {
    # Filtering preserves LAS metadata and point attributes.
    tagged <- data.table::copy(las)
    quality <- integer(nrow(d)); quality[keep] <- 1L; quality[review] <- 2L
    tagged <- lidR::add_lasattribute(tagged,quality,"TrunkFitQuality",
                                   "1 strict; 2 coverage review")
    tagged <- lidR::add_lasattribute(tagged,as.integer(dependent),"PathUsesCoverage",
                                   "Trace used coverage trial")
    trunks <- lidR::filter_poi(tagged, keep)
    saveRDS(trunks,file.path(out,"candidate_trunks_HAG.rds"),compress=FALSE)
    # Restore the input elevation on a copy. Starting from the existing LAS
    # avoids the default coordinate quantization of a new LAS constructor.
    selected <- d[keep,,drop=FALSE]
    export_prep_error <- NULL
    exchange <- tryCatch({
      x <- lidR::unnormalize_height(data.table::copy(trunks))
      x <- lidR::las_rescale(x,0.001,0.001,0.001)
      x <- lidR::add_lasattribute(x,selected$HAG,
                                 "HeightAboveGround","Height above ground (m)")
      CspStandSegmentation::add_las_attributes(x)
    },error=function(e) {
      export_prep_error <<- conditionMessage(e)
      NULL
    })
    export_one <- function(cloud,filename) {
      tryCatch({
        lidR::writeLAS(cloud,file.path(out,filename))
        data.frame(File=filename,Status="written",Error="")
      },error=function(e) data.frame(File=filename,Status="failed",Error=conditionMessage(e)))
    }
    if(is.null(exchange)) {
      export_status[[1]] <- data.frame(File="LAZ preparation",Status="failed",
                                      Error=export_prep_error)
    } else {
      export_status[[1]] <- export_one(exchange,"candidate_trunks_originalZ.laz")
      for(id in sort(unique(selected$TreeID))) {
        one <- lidR::filter_poi(exchange, TreeID == id)
        export_status[[length(export_status)+1L]] <- export_one(
          one,sprintf("stem_candidate_%02d_originalZ.laz",id))
      }
    }
    write.csv(do.call(rbind,export_status),file.path(out,"LAZ_export_status.csv"),row.names=FALSE)
    if(show_3d) tryCatch({
      # Explicit bins ensure that each numeric TreeID has its own color.
      plot(trunks,color="TreeID",pal=unname(palette),breaks=seq(.5,7.5,1),size=2)
    },error=function(e) message("3D display unavailable: ",conditionMessage(e)))
  }
  cat("\nPilot trunk trace summary:\n")
  print(summary,row.names=FALSE,digits=4)
  if(length(export_status)) print(do.call(rbind,export_status),row.names=FALSE)
  cat("\nOutputs:",out,"\n")
  cat("TracedMaxHAG_m is retained extent, not total tree/trunk height.\n")
  cat("Black profile points and TrunkFitQuality=2 require coverage review.\n")
  cat("PathUsesCoverage=1 also flags downstream dependence on coverage trials.\n")
  cat("Share trunk_profiles.png, trunk_trace_summary.csv and trunk_slice_checks.csv.\n")
  invisible(list(trunks=trunks,models=models,summary=summary,output_dir=out,
                 export_status=export_status))
}

# ---- Run on the corrected pilot ----------------------------------------------
if(exists("tls_trees_v2",inherits=TRUE) && inherits(tls_trees_v2,"LAS")) {
  tls06_parent <- if(exists("tls_output_dir_v2",inherits=TRUE)) tls_output_dir_v2 else NULL
  if(is.null(tls06_parent) || !dir.exists(tls06_parent)) {
    message("Select pilot_tree_ids_v2.rds to locate its output folder.")
    tls06_checkpoint <- file.choose()
    tls06_parent <- dirname(tls06_checkpoint)
  }
  tls06 <- tls06_run(tls_trees_v2,tls06_seeds,tls06_parent)
} else {
  message("Select pilot_tree_ids_v2.rds from the step 3 output folder.")
  tls06_checkpoint <- file.choose()
  tls06 <- tls06_run(readRDS(tls06_checkpoint),tls06_seeds,dirname(tls06_checkpoint))
}
tls_trunks_v3 <- tls06$trunks
tls_trunk_summary_v3 <- tls06$summary
tls_trunk_models_v3 <- tls06$models
tls_trunk_output_dir_v3 <- tls06$output_dir
