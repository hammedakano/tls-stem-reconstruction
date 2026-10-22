# Step 4: trace candidate trunk surfaces in the seven-tree plot167 pilot.
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
tls04_settings <- list(
  step = 0.20, anchor_h = 1.30, lowest_center = 0.30,
  tolerance = 0.02, keep_tolerance = 0.025,
  min_points = 20L, min_fraction = 0.35, min_sectors = 8L,
  max_empty_angle = 120, max_rmse = 0.012,
  min_radius = 0.04, max_radius = 0.60,
  center_step = 0.10, radius_step = 0.04,
  max_failures = 2L, attempts = 3L, iterations = 300L
)
tls04_seeds <- data.frame(
  TreeID = 1:7,
  X = c(420409.339,420410.381,420411.695,420406.740,
        420403.218,420410.515,420404.018),
  Y = c(6481818.000,6481816.073,6481813.066,6481811.905,
        6481819.670,6481819.733,6481814.116),
  Radius = c(0.1261,0.1054,0.1352,0.1584,0.2225,0.1480,0.1827)
)

tls04_fit_slice <- function(d, bin, ref, p) {
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

tls04_trace_one <- function(d, seed, p) {
  d$x <- d$X - seed$X
  d$y <- d$Y - seed$Y
  # Half-open, non-overlapping bins, anchored on [1.2,1.4).
  d$Bin <- floor((d$HAG-(p$anchor_h-p$step/2))/p$step + 1e-8)
  ref <- list(x=0,y=0,r=seed$Radius,h=p$anchor_h)
  anchor <- tls04_fit_slice(d, 0L, ref, p)
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
        fit <- tls04_fit_slice(d, bin, ref, p)
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
    X=if(f$ok) f$x+seed$X else NA_real_,
    Y=if(f$ok) f$y+seed$Y else NA_real_,
    Radius_m=if(f$ok) f$r else NA_real_,
    RingPoints=if(f$ok) f$n else 0L,
    RingFraction=if(f$ok) f$fraction else NA_real_,
    RadialRMSE_m=if(f$ok) f$rmse else NA_real_,
    OccupiedSectorsOf12=if(f$ok) f$sectors else NA_integer_,
    MaxEmptyAngle_deg=if(f$ok) f$empty_angle else NA_real_
  )))
  kept <- unique(unlist(lapply(fits,function(f) if(f$ok) f$keep else integer()),
                        use.names=FALSE))
  heights <- d$HAG[match(kept,d$Row)]
  summary <- data.frame(TreeID=seed$TreeID, CandidatePoints=length(kept),
    AcceptedSlices=sum(models$Accepted), FailedSlices=sum(!models$Accepted),
    TracedMinHAG_m=if(length(kept)) min(heights) else NA_real_,
    TracedMaxHAG_m=if(length(kept)) max(heights) else NA_real_,
    PilotDiameterAt1p3_cm=if(anchor$ok) 200*anchor$r else NA_real_,
    LowerStop=unname(stops["down"]), UpperStop=unname(stops["up"]))
  list(models=models,keep=kept,summary=summary)
}

tls04_draw_profiles <- function(d, keep, seeds, palette) {
  oldpar <- par(mfrow=c(2,7),mar=c(3.2,3.1,2.2,0.5),oma=c(1,1,3,0))
  on.exit(par(oldpar))
  set.seed(260924)
  for (axis in c("X","Y")) for (i in seq_len(nrow(seeds))) {
    id <- seeds$TreeID[i]
    rows <- which(d$TreeID == id)
    retained <- rows[keep[rows]]
    background <- if(length(rows)>8000L) sample(rows,8000L) else rows
    offset <- seeds[[axis]][i]
    xr <- if(length(retained)) range(d[[axis]][retained]-offset)+c(-.4,.4) else c(-1,1)
    zr <- c(0,max(d$HAG[rows]))
    plot(d[[axis]][background]-offset,d$HAG[background],xlim=xr,ylim=zr,
      pch=16,cex=.15,col="grey80",xlab=paste(axis,"offset (m)"),
      ylab="Height above ground (m)",main=paste("Tree",id,axis,"view"))
    if(length(retained)) points(d[[axis]][retained]-offset,d$HAG[retained],
                               pch=16,cex=.22,col=palette[as.character(id)])
  }
  mtext("Candidate trunk surfaces in local neighborhoods: gray context, color retained",
        outer=TRUE,cex=1.05,line=1)
}

tls04_run <- function(las, seeds, parent_dir, p=tls04_settings, show_3d=interactive()) {
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
  out <- tempfile("trunk_trace_",tmpdir=parent_dir)
  if(!dir.create(out)) stop("Cannot create output folder.")
  message("Saving this trace to: ",out)
  saveRDS(p,file.path(out,"trace_settings.rds"))
  write.csv(seeds,file.path(out,"trace_seeds.csv"),row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  results <- vector("list",nrow(seeds))
  keep <- rep(FALSE,nrow(d))
  set.seed(260924)
  for(i in seq_len(nrow(seeds))) {
    id <- seeds$TreeID[i]
    message("Tracing tree ",id," of ",nrow(seeds)," ...")
    tree <- d[which(d$TreeID == id),,drop=FALSE]
    results[[i]] <- tls04_trace_one(tree,seeds[i,,drop=FALSE],p)
    keep[results[[i]]$keep] <- TRUE
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
  tryCatch(tls04_draw_profiles(d,keep,seeds,palette),finally=grDevices::dev.off())
  if(interactive()) tls04_draw_profiles(d,keep,seeds,palette)

  trunks <- NULL
  export_status <- list()
  if(any(keep)) {
    # Filtering preserves LAS metadata and point attributes.
    trunks <- lidR::filter_poi(las, keep)
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
  cat("Share trunk_profiles.png and trunk_trace_summary.csv for review.\n")
  invisible(list(trunks=trunks,models=models,summary=summary,output_dir=out,
                 export_status=export_status))
}

# ---- Run on the corrected pilot ----------------------------------------------
if(exists("tls_trees_v2",inherits=TRUE) && inherits(tls_trees_v2,"LAS")) {
  tls04_parent <- if(exists("tls_output_dir_v2",inherits=TRUE)) tls_output_dir_v2 else NULL
  if(is.null(tls04_parent) || !dir.exists(tls04_parent)) {
    message("Select pilot_tree_ids_v2.rds to locate its output folder.")
    tls04_checkpoint <- file.choose()
    tls04_parent <- dirname(tls04_checkpoint)
  }
  tls04 <- tls04_run(tls_trees_v2,tls04_seeds,tls04_parent)
} else {
  message("Select pilot_tree_ids_v2.rds from the step 3 output folder.")
  tls04_checkpoint <- file.choose()
  tls04 <- tls04_run(readRDS(tls04_checkpoint),tls04_seeds,dirname(tls04_checkpoint))
}
tls_trunks <- tls04$trunks
tls_trunk_summary <- tls04$summary
tls_trunk_models <- tls04$models
tls_trunk_output_dir <- tls04$output_dir
