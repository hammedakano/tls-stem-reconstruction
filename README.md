# tls-stem-reconstruction

R workflow for segmenting tree stems from terrestrial laser scanning (TLS) point clouds and reconstructing each stem slice by slice. The output is a diameter profile, DBH at 1.3 m and a quality log for every slice.

**Skills shown:** LiDAR processing, automated tree detection and segmentation, stem reconstruction, scalable processing of large point clouds, honest quality reporting.

## Workflow

| Script | Purpose |
|---|---|
| `R/TLS_02_Plot_Pilot.R` | Reads a cropped, sampled and voxel-thinned subset of a large LAZ file, classifies ground with cloth simulation (RCSF), normalises height and detects candidate stems. Saves checkpoints so later steps never reread the full cloud. |
| `R/TLS_03_Rerun_Seven_Seeds.R` | Re-segments trees from reviewed stem seeds (CspStandSegmentation) and checks point IDs around each ring at 1.2 to 1.4 m. |
| `R/TLS_04_Trace_Candidate_Trunks.R` | Fits RANSAC circles in 0.2 m slices from 1.3 m upward and downward and tracks centre and radius continuity. Two consecutive failed slices stop a trace. Gaps are left as gaps. |
| `R/TLS_05_Inspect_Rejected_Slices.R` | Exports each rejected slice with its preceding accepted slice for diagnosis. |
| `R/TLS_06_Trace_With_Coverage_QC.R` | Continuity-guided refitting with explicitly flagged coverage trials (tighter RMSE, centre and radius limits). Every slice is labelled strict, coverage trial or rejected. |
| `python/make_tls_figures.py` | Rebuilds the pilot figures and summary table from the slice logs. |

Run the R scripts in order in RStudio (open, then Source). Each script reuses objects from the previous step if they are in memory, or asks for the saved `.rds` checkpoint. Outputs go to a new folder each time and inputs are never overwritten.

## Scale

The pipeline was written for a plot of about 290 million points (about 2 GB LAZ). It runs on a laptop because it crops, samples and thins before analysis and saves checkpoints at each stage. Final fitting can revisit full-density neighbourhoods around each stem.

## Pilot result

Seven trees in the central 10 x 10 m of one plot.

- 540 of 563 tested slices retained (95.9%)
- traced stem extent 13.4 to 18.2 m
- fitted diameter at 1.3 m 21.1 to 44.4 cm
- median circle-fit residual 5.3 mm

These are internal fit diagnostics, not independent validation against tape DBH. Per-tree values are in `docs/tls_pilot_operational_summary.csv`.

## Installation

```r
source("install_packages.R")   # lidR, RCSF, data.table, CspStandSegmentation
```
For the figures, install the Python requirements with `pip install -r python/requirements.txt` and place the slice log tables in `python/inputs/`.

## Data

Point clouds are not included. `.gitignore` excludes `data/`, `outputs/` and all point cloud and checkpoint formats.

## Related repositories

- [tls-stem-reconstruction](../tls-stem-reconstruction): TLS stem segmentation and diameter profiles
- [uav-sfm-tree-metrics](../uav-sfm-tree-metrics): drone point cloud tree detection and crown metrics
- [forest-model-transfer-uncertainty](../forest-model-transfer-uncertainty): site-transfer validation and prediction intervals

## Author

Hammed Akano, PhD candidate, Deakin University. Contact hammedaakano@gmail.com.

## Licence

MIT. See `LICENSE`.
