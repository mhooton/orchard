# Orchard

![Alt text](orchard.png)

Automated photometry pipeline for the SPECULOOS telescope network. Processes raw astronomical images through calibration, astrometry, and differential photometry to produce publication-quality light curves for exoplanet transit detection and characterisation.

## Overview

Orchard is a modernised, dockerised photometry pipeline originally developed for the SPECULOOS project. The pipeline handles the complete data reduction workflow from raw CCD images to calibrated light curves, supporting multiple telescope/instrument configurations through a centralised JSON configuration system.

**Key Features:**
- Multi-instrument support (ANDOR, SPIRIT NIR)
- Automated nightly processing via cron jobs
- Local Gaia DR3 catalogue for fast astrometric calibration
- Proper motion corrections for multi-year baselines
- PWV atmospheric correction for cool stars at Paranal
- Comprehensive quality monitoring and reporting
- Dockerised deployment for reproducibility

**Supported Telescopes:**
- SPECULOOS Southern Observatory: Io, Europa, Ganymede, Callisto (SSO)
- SPECULOOS Northern Observatory: Artemis
- SAINT-EX

## Architecture

### Technology Stack

- **Core Language:** Python 3.9+
- **Key Dependencies:**
  - `casutools` - CASU pipeline tools (imcore, imstack, wcsfit)
  - `astropy` - FITS handling and coordinate transformations
  - `photutils` - Source detection and photometry
  - `twirl` - Fast pattern matching for plate solving
- **Deployment:** Docker with docker-compose
- **Catalogue:** Local Gaia DR3 + 2MASS cross-match database (~50GB)

### Pipeline Philosophy

Orchard follows a staged processing model where each stage can be run independently or skipped via command-line flags. This design enables:
- Efficient debugging of specific stages
- Reprocessing with different parameters
- Recovery from partial failures
- Flexible deployment for different use cases

The pipeline emphasises **robustness over optimisation** - it's designed to handle real-world observing conditions including clouds, tracking errors, focus issues, and occasional corrupted frames.

## Data Directory Structure

The pipeline expects data organised as follows:

```
/data/SPECULOOSPipeline/
├── Observations/
│   └── {TELESCOPE}/              # e.g., Io, Europa, Artemis
│       └── images/
│           └── {YYYYMMDD}/       # Observation date
│               ├── IMAGE*.fits   # Raw science frames
│               ├── IMAGE*.fts    # Alternative extension
│               └── *.log         # Acquisition logs
│
└── PipelineOutput/
    └── {VERSION}/                # e.g., v2, v3
        └── {TELESCOPE}/
            ├── output/
            │   ├── {YYYYMMDD}/
            │   │   ├── reduction/              # Master calibration frames
            │   │   │   ├── {run}_MasterBias.fits
            │   │   │   ├── {run}_MasterDark*.fits
            │   │   │   ├── {run}_MasterFlat_*.fits
            │   │   │   ├── {run}_BadPixelMap.fits
            │   │   │   ├── readoutnoise.dat
            │   │   │   ├── darkcurrent.dat
            │   │   │   └── *.list              # Image classification lists
            │   │   │
            │   │   └── {TARGET}/               # Per-target processing
            │   │       └── {run}/
            │   │           ├── proc*.fits      # Reduced science frames
            │   │           ├── *_processed.dat # Frame lists by filter
            │   │           ├── *.phot          # Photometry catalogues
            │   │           └── lightcurves/    # Differential photometry
            │   │
            │   └── StackImages/                # Reference catalogues
            │       ├── {GAIA_ID}_outstack_*.fits
            │       └── {GAIA_ID}_stack_catalogue_*.fits
            │
            ├── reports/                        # Quality control reports
            │   └── QC/
            │
            └── logs/                           # Pipeline execution logs
                └── {YYYYMMDD}_*.log
```

### Key Directories

- **`Observations/{TELESCOPE}/images/{DATE}/`**: Raw data from telescopes (read-only during processing)
- **`PipelineOutput/{VERSION}/{TELESCOPE}/output/`**: All pipeline products
- **`StackImages/`**: Deep reference images and source catalogues (reused across nights)
- **`reports/QC/`**: Time-series quality metrics (readout noise, zero point, dark current, etc.)

## Instrument Configuration

Telescope and instrument parameters are defined in `src/calibration/instrument_config.json`. This centralised configuration drives all instrument-specific processing including:

- Gain and saturation levels
- CCD geometry (overscan, trim regions)
- Dark frame handling (exposure-matched vs combined)
- Bad pixel detection thresholds
- Header keyword mappings

### Example Configuration

```json
{
  "telescopes": {
    "SPECULOOS-CALLISTO": {
      "andor": {
        "gain": 1.0029,
        "saturation_threshold": 62000,
        "naxis2": 2088,
        "dark_matching_exptime": false,
        "overscan": {
          "top_row": -15,
          "bottom_row": null,
          "left_col": 4,
          "right_col": null
        },
        "trim": {
          "top_row": 22,
          "bottom_row": 2066,
          "left_col": 2,
          "right_col": 2048
        }
      },
      "spirit": {
        "gain": 5.09,
        "saturation_threshold": 11000,
        "naxis2": 1280,
        "dark_matching_exptime": true,
        "bad_pixel_correction": {
          "preferred_flat_filter": "zYJ",
          "hot_sigma_threshold": 5,
          "cold_sigma_threshold": 5,
          "flat_threshold": 0.1
        }
      }
    }
  }
}
```

### Instrument-Specific Features

**ANDOR (Optical):**
- Combined dark frame (scaled by exposure time)
- Overscan correction with sigma clipping
- 2044×2048 pixel trimmed output
- No bad pixel correction

**SPIRIT (NIR):**
- Exposure-matched dark frames
- Hot/cold/flatbad pixel detection and interpolation
- 1280×1280 pixel output
- Essential for source detection in noisy NIR images

## Usage

### Docker Setup

The pipeline runs inside a Docker container for reproducibility and portability.

**Starting the Container:**

```bash
docker-compose -f /appct/data/speculoos/orchard/docker/docker-compose.server.yml up -d
```

This starts the `orchard-server` container with:
- Source code mounted from `/appct/data/speculoos/orchard/src`
- Data access to `/export/data/SPECULOOSPipeline`
- Gaia database mounted at `/gaia_database`
- Environment variables: `N_CORES=20`, `GAIADATABASEPATH`, etc.

**Entering the Container:**

```bash
docker exec -it orchard-server bash
```

**Example with Pipeline Execution:**

```bash
docker exec orchard-server bash -c "./main/ZLP_pipeline.sh --force-platesolve 1 /data/SPECULOOSPipeline 20260128 8 2 Io"
```

**Note:** All subsequent commands assume you're running inside the container. The `docker exec orchard-server bash -c` prefix is omitted for clarity.

---

### Basic Usage

The main pipeline script is `ZLP_pipeline.sh` located in `src/main/`.

**Syntax:**

```bash
./main/ZLP_pipeline.sh [OPTIONS] <runname> <root-directory> <dates> <c_thresh> <s_thresh> <telescope> [input_targets]
```

**Mandatory Positional Arguments:**

| Argument | Description | Example |
|----------|-------------|---------|
| `runname` | Identifier for this pipeline run | `1` |
| `root-directory` | Base directory containing Observations and PipelineOutput | `/data/SPECULOOSPipeline` |
| `dates` | Space-separated observation dates (YYYYMMDD format) | `20260128` or `20260128 20260129 20260130` |
| `c_thresh` | Detection threshold for catalogue creation (σ) | `2` (default for deep stacks) |
| `s_thresh` | Detection threshold for plate solving (σ) | `20` (higher for reliable WCS) |
| `telescope` | Telescope name | `Io`, `Europa`, `Ganymede`, `Callisto`, `Artemis` |
| `input_targets` | Optional: specific targets to process | Target names (space-separated) |

**Optional Parameters:**

| Flag | Description | Default |
|------|-------------|---------|
| `--cores N` | Number of parallel processes | `N_CORES` env var (20), or 1 |
| `--force-platesolve` | Re-solve images that already have WCS | Disabled |
| `--no_T1` through `--no_T12` | Skip specific pipeline stage | All enabled |
| `--only_T1` through `--only_T12` | Run only specific stage | All enabled |

---

### Common Workflows

#### 1. Process Single Night (Full Pipeline)

Process all targets observed on 2026-01-28 with Io telescope:

```bash
./main/ZLP_pipeline.sh 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

**What Happens:**
- T1: Classifies all images from `Observations/Io/images/20260128/`
- T2-T5: Creates master calibration frames
- T6: Reduces all science images
- T7: Plate solves reduced images
- T8: Creates/updates stack catalogues for each target
- T9: Performs aperture photometry (apertures 3-8)
- T10: Generates differential photometry light curves
- T11: Creates PDF report
- T12: Migrates successful results to production

**Output:** Complete processing for all targets, ready for science analysis.

---

#### 2. Process Date Range

Process three consecutive nights:

```bash
./main/ZLP_pipeline.sh 1 /data/SPECULOOSPipeline "20260128 20260129 20260130" 2 20 Europa
```

The pipeline loops over each date sequentially.

---

#### 3. Process Specific Targets

Process only two targets from a night's observations:

```bash
./main/ZLP_pipeline.sh 1 /data/SPECULOOSPipeline 20260128 2 20 Ganymede "TRAPPIST-1 TOI-1234"
```

**Note:** Target names with spaces must use double hyphens in file lists (handled automatically by `createlists.py`).

---

#### 4. Reprocess with Forced Plate Solving

Re-run plate solving even for images with existing WCS solutions:

```bash
./main/ZLP_pipeline.sh --force-platesolve 1 /data/SPECULOOSPipeline 20260128 2 20 Callisto
```

**Use Case:** After Gaia database updates, plate scale changes, or debugging WCS issues.

---

#### 5. Skip Calibration (Reprocess Photometry Only)

Reprocess photometry using existing calibration frames:

```bash
./main/ZLP_pipeline.sh --no_T1 --no_T2 --no_T3 --no_T4 --no_T5 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

Runs T6 onward, using existing master bias/dark/flat from previous run.

---

#### 6. Run Single Stage (Debugging)

Create stack catalogues only (useful for testing catalogue matching):

```bash
./main/ZLP_pipeline.sh --only_T8 1 /data/SPECULOOSPipeline 20260128 2 20 Artemis
```

**Available Stages:**
- `--only_T1`: Image classification
- `--only_T6`: Science reduction
- `--only_T7`: Plate solving
- `--only_T8`: Stack creation
- `--only_T9`: Aperture photometry
- `--only_T10`: Differential photometry

---

#### 7. Serial Processing (Debugging)

Run with single core for clearer error messages:

```bash
./main/ZLP_pipeline.sh --cores 1 1 /data/SPECULOOSPipeline 20260128 2 20 Europa
```

**Recommendation:** Always debug with `--cores 1` to avoid interleaved log output from parallel processes.

---

#### 8. Multi-Core Production Run

Process with maximum parallelisation:

```bash
./main/ZLP_pipeline.sh --cores 20 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

Significant speedup for:
- T6: Science reduction (parallel processing of frames)
- T7: Plate solving (parallel WCS fitting)
- T8: WCS solving during stack creation
- T9: Aperture photometry (parallel processing)

**Typical Performance:**
- Single core: ~30-60 minutes per night
- 20 cores: ~5-15 minutes per night (depending on target count)

---

### Advanced Usage

#### Custom Aperture Photometry

The pipeline processes apertures 3-8 by default (in `SPlightcurve.py`). To process a single aperture:

```bash
cd /opt/orchard/src
python diff_photometry/SPlightcurve.py \
  --date 20260128 \
  --targ_gaia 1234567890123456789 \
  --ap 6 \
  --filt I+z \
  --outfits /data/SPECULOOSPipeline/PipelineOutput/v3/Io/output/20260128/1234567890123456789/1234567890123456789_Iz_20260128_output.fits \
  --goutfits /data/SPECULOOSPipeline/PipelineOutput/v3/Io/output/1234567890123456789_Iz_output.fits \
  --outdir /data/SPECULOOSPipeline/PipelineOutput/v3/Io/output/20260128/1234567890123456789/ \
  --lcdir /data/SPECULOOSPipeline/PipelineOutput/v3/Io/output/20260128/1234567890123456789/lightcurves/ \
  --version v3
```

**Aperture Radii:**
- Aperture 3: 1.0 × rcore (4 pixels = 1.4" for SPIRIT)
- Aperture 4: √2 × rcore (5.66 pixels = 2.0")
- Aperture 5: 2.0 × rcore (8 pixels = 2.8")
- Aperture 6: 2√2 × rcore (11.3 pixels = 4.0") ← **Typical optimal**
- Aperture 7: 4.0 × rcore (16 pixels = 5.6")
- Aperture 8: 5.0 × rcore (20 pixels = 7.0")

**Selection Criterion:** Choose aperture that minimises point-to-point scatter while avoiding significant sky background contribution.

---

#### PWV Correction (Paranal Only)

Apply water vapour correction to cool star observations:

```bash
cd /opt/orchard/src
python water_vapour/water_vapour.py \
  /data/SPECULOOSPipeline/PipelineOutput/v3/Io/output/20260128/1234567890123456789/ \
  --filt I+z \
  --ap 6 \
  --targname TRAPPIST-1 \
  --targ_teff 2600 \
  --sp_id SP1234567 \
  --date 20260128
```

**Requirements:**
- Target Teff between 2000-3200K
- ESO PWV data available for observation date
- Pre-computed atmospheric model grid for filter

**Effect:** Typically reduces RMS by 5-20% for I+z observations under variable PWV conditions.

---

#### Monitoring Pipeline Execution

**Real-Time Monitoring:**

```bash
tail -f /data/SPECULOOSPipeline/PipelineOutput/v3/Io/logs/20260128_1_v3.log
```

**Check Plate Solve Success Rate:**

```bash
grep "SUCCESS\|FAILED" /data/SPECULOOSPipeline/PipelineOutput/v3/Io/logs/20260128_1_v3.log | \
  awk '{print $NF}' | sort | uniq -c
```

**Identify Problem Frames:**

```bash
grep "FAILED" /data/SPECULOOSPipeline/PipelineOutput/v3/Io/logs/20260128_1_v3.log
```

---

### Environment Variables

Key variables defined in `docker-compose.server.yml`:

| Variable | Purpose | Default |
|----------|---------|---------|
| `ORCHARD_PATH` | Pipeline source code location | `/opt/orchard/src` |
| `PYTHONPATH` | Python module search path | `/opt/orchard/src` |
| `N_CORES` | Default parallelisation level | `20` |
| `GAIADATABASEPATH` | Local Gaia catalogue database | `/gaia_database/gaia_tmass_16_jm_cut.db` |

Override defaults via `--cores` flag or by modifying docker-compose file.

---

## Pipeline Stages

The pipeline consists of 12 stages (T1-T12) that can be controlled individually via command-line flags.

### Stage T1: Image Classification

**Module:** `src/main/createlists.py`

Sorts raw nightly observations into science targets and calibration frames based on FITS header keywords (`IMAGETYP`, `EXPTIME`). Creates organised file lists that drive subsequent processing.

**Output:**
```
{run}_image_{target}.list    # Science frames per target
{run}_bias.list              # Bias frames
{run}_dark.list              # Dark frames
{run}_flat_{filter}.list     # Flat fields per filter
```

**Handles:** Multiple image extensions (.fits, .fts, .fz, .bz2), target name sanitisation, filter-specific organisation.

---

### Stage T2: Master Bias Creation

**Module:** `src/calibration/pipebias.py`

Creates master bias frames through sigma-clipped median combination of bias exposures. Computes readout noise from consecutive bias frames and measures overscan levels.

**Process:**
1. Load all bias frames
2. Apply overscan correction to each frame
3. Sigma-clipped median combination (removes cosmic rays)
4. Calculate readout noise: `RON = σ(B1-B2) × gain / √2`
5. Write master bias and metrics

**Output:**
- `{run}_MasterBias.fits`
- `readoutnoise.dat` (in electrons)
- `overscan.dat` (median level in ADU)

**Quality Metrics:** Typical readout noise: 5-15 e⁻ for ANDOR, 15-25 e⁻ for SPIRIT

---

### Stage T3: Master Dark Creation

**Module:** `src/calibration/pipedark.py`

Creates master dark frames to correct for thermal current. Supports two modes depending on instrument configuration:

**Combined Mode (ANDOR):**
- Combines all dark exposures into single master
- Scales to ADU/s for flexible exposure time matching
- Assumes linear dark current

**Exposure-Matched Mode (SPIRIT):**
- Creates separate master darks for each exposure time
- Better handles non-linear dark behaviour in NIR detectors
- Automatically selects appropriate dark during reduction

**Process:**
1. Load all dark frames
2. Apply overscan and bias correction
3. Normalise by exposure time (or group by exposure)
4. Sigma-clipped median combination
5. Calculate dark current: `DC = gain × median(dark)`

**Output:**
- `{run}_MasterDark.fits` (combined mode)
- `{run}_MasterDark_{exptime}s.fits` (matched mode)
- `darkcurrent.dat` (in e⁻/s)

**Quality Metrics:** Typical dark current: 0.1-0.5 e⁻/s for ANDOR, 5-50 e⁻/s for SPIRIT

---

### Stage T4: Master Flat Creation

**Module:** `src/calibration/pipeflat.py`

Creates normalised master flat fields from dome or twilight flats. Processes flats separately per filter, applying bias and dark corrections before combining.

**Process:**
1. Apply overscan, bias, and dark corrections to each flat
2. Normalise each flat by its median flux
3. Sigma-clipped median combination (rejects cosmic rays)
4. Normalise master to unity mean
5. Generate variance and standard deviation maps

**Output:**
- `{run}_MasterFlat_{filter}.fits` (one per filter)
- `variance.fts` (quality assessment)
- `std.fts` (pixel-to-pixel variation)

**Purpose:** Corrects for pixel sensitivity variations, vignetting, and dust shadows.

---

### Stage T5: Bad Pixel Map Creation

**Module:** `src/calibration/pipebadpixel.py`

Identifies and maps bad pixels for interpolation during reduction. Currently only used for SPIRIT NIR camera where hot pixels significantly impact source detection.

**Detection Methods:**

1. **Hot Pixels:** Dark > median + 5σ
2. **Cold Pixels:** Dark < median - 5σ  
3. **Flatbad Pixels:** Flat < 0.1 (low sensitivity)

**Process:**
1. Load shortest-exposure master dark
2. Detect hot/cold pixels via sigma thresholding
3. Detect flatbad pixels from normalised master flat
4. Combine all maps (logical OR)
5. Calculate cleaned dark current

**Output:**
- `{run}_BadPixelMap.fits` (boolean mask)
- `hot_pixel_map.fits`, `cold_pixel_map.fits`, `flatbad_pixel_map.fits`
- Updated `darkcurrent.dat` (with bad pixels masked)

**Statistics:** Typical bad pixel fractions: <0.1% for ANDOR, 1-3% for SPIRIT

---

### Stage T6: Science Image Reduction

**Module:** `src/calibration/pipered.py`

Applies all calibration corrections to science frames, producing reduced images ready for astrometry and photometry.

**Corrections Applied:**
1. Overscan subtraction (sigma-clipped median)
2. Bias subtraction
3. Dark subtraction (exposure-time scaled or matched)
4. Flat field division (filter-specific)
5. Bad pixel interpolation (SPIRIT only)
6. CRPIX adjustment for trimming

**Process:**
```python
reduced = (raw - overscan - bias - dark × exptime) / flat
reduced = interpolate_bad_pixels(reduced, bad_pixel_map)
```

**Output:**
- `proc{image}.fits` for each science frame
- `{filter}_processed.dat` lists of successfully reduced images
- FITS headers updated with calibration history

**Parallel Processing:** Supports multiprocessing for faster execution on multi-core systems.

---

### Stage T7: Astrometric Calibration

**Modules:** `src/astrom/astrometry.py`, `src/astrom/pointer_wcs.py`

Performs plate solving to establish World Coordinate System (WCS) solutions for each reduced science frame. Uses local Gaia DR3 database for fast, offline astrometric calibration.

**Why After Reduction for SPIRIT:**  
NIR images are noisy with many hot pixels. Calibration removes these false sources, dramatically improving plate solve success rates.

**Process:**

1. **Source Detection** (multiscale algorithm):
   - Cleans image via background subtraction and filtering
   - Detects sources at multiple smoothing scales (1, 2, 3 pixels)
   - Removes edge artefacts and duplicates
   - Selects brightest 16 sources for matching

2. **Catalogue Query:**
   - Extracts target RA/Dec from FITS header
   - Queries local Gaia database with 1.2× FOV margin
   - Applies proper motion corrections to observation date
   - Returns brightest ~32 Gaia sources in field

3. **Pattern Matching:**
   - Uses `twirl` library for fast triangle pattern matching
   - Matches detected sources to Gaia positions
   - Computes WCS transformation (6-parameter plate solution)
   - Validates with minimum 4 matched stars

4. **WCS Propagation:**
   - Writes WCS keywords to reduced image header
   - Copies WCS back to raw image header (with CRPIX offsets)
   - Enables both reduced and raw images for subsequent processing

**Output:**
- Updated FITS headers with WCS keywords (CRVAL, CRPIX, CD matrix)
- Plate solve success/failure logged to console
- Diagnostic plots (if verbose mode enabled)

**Parameters:**
- `--force-platesolve`: Re-solve images that already have WCS
- `--timeout`: Maximum time per image (default: 60s)
- `--cores`: Number of parallel processes

**Performance:** Typical plate solve time: 2-5s per image with local database, ~60s timeout handles difficult fields.

**Success Rates:** >95% for ANDOR (optical), >85% for SPIRIT (NIR after reduction).

---

### Stage T8: Stack Image & Catalogue Creation

**Module:** `src/photometry/ZLP_create_cat.py`

Creates deep stacked images from multiple science frames and generates source catalogues for photometric reference. Essential for establishing the field star catalogue used in differential photometry.

**Process:**

1. **Frame Selection:**
   - Selects subset of images from middle of night (lowest background)
   - Default: 50 frames (configurable via `--nfiles`)
   - Filters out frames with failed WCS solutions

2. **Image Stacking:**
   - Uses `casutools.imstack` for WCS-based stacking
   - Weighted median combination with confidence maps
   - Sigma-clipped outlier rejection (5σ threshold)
   - Exposure time normalisation

3. **Source Detection:**
   - Runs `casutools.imcore` on stacked image
   - Lower detection threshold than individual frames (2σ vs 20σ)
   - Detects ~100-1000 sources per field depending on crowding

4. **Gaia Cross-Matching:**
   - Matches detected sources to Gaia DR3 via CDS VizieR
   - Retrieves proper motions, parallaxes, G-band magnitudes
   - Stores cross-match in `Gaia_Crossmatch` FITS extension
   - Reports match percentage for quality assessment

5. **Backup System:**
   - Preserves best stack catalogue when subsequent runs fail
   - Automatically restores backup if current run produces poor results
   - Useful for nights with degraded data quality

**Output:**
- `{GAIA_ID}_outstack_{filter}.fits` - Deep stacked image (~50× deeper than single frame)
- `{GAIA_ID}_stack_catalogue_{filter}.fits` - Source catalogue with:
  - Image coordinates (X, Y pixels)
  - Equatorial coordinates (RA, Dec)
  - Aperture photometry (13 apertures)
  - Gaia DR3 cross-match (proper motions, G-mags, source IDs)
- `{filter}_stacked.dat` - List of images included in stack

**Quality Metrics:**
- Number of detected sources
- Gaia cross-match percentage (typically >80% for good fields)
- Stack depth (number of contributing frames)

**Reuse:** Stack catalogues are stored in `StackImages/` and reused across nights for the same target, avoiding repeated stacking.

---

### Stage T9: Aperture Photometry

**Modules:** `src/photometry/ZLP_app_photom.py`, `src/photometry/wcs_photom.py`, `src/astrom/proper_motion.py`

Performs forced aperture photometry at catalogue positions across all science frames. Applies proper motion corrections and measures image quality metrics.

**Process:**

1. **Proper Motion Correction:**
   - Loads stack catalogue with Gaia cross-match
   - Calculates new RA/Dec for each star from Gaia epoch (J2016.0) to observation date
   - Formula: `new_ra = ra + (pmra/cos(dec)) × Δt`, where Δt in years
   - Accounts for parallax effects on transverse velocities
   - Creates PM-corrected catalogue for photometry

   **Why Critical:** High proper motion stars (>50 mas/yr) shift by >0.5 pixels over 3 years. Without correction, apertures miss star centres causing flux loss and systematic errors.

2. **Forced Photometry:**
   - Uses `casutools.imcore_list` for aperture photometry at fixed positions
   - Places apertures at PM-corrected catalogue coordinates
   - Extracts flux in circular apertures (default: 4 pixel radius)
   - Measures local background and applies aperture corrections
   - Parallel processing supported (default: 20 cores)

3. **Image Quality Assessment:**
   
   **PSF Measurements:**
   - Divides image into 3×3 grid
   - Fits 2D Gaussian PSF to bright stars in each region
   - Extracts FWHM semi-major axis (PSF_a), semi-minor axis (PSF_b), rotation angle (PSF_t)
   - Detects focus gradients and field aberrations

   **Seeing Calculation:**
   - Median FWHM across all detected stars
   - Converted to arcseconds using plate scale (0.35 "/pixel for SPIRIT, 0.64 "/pixel for ANDOR)
   - Written to FITS header as `SEEING` keyword

   **Cloud Detection:**
   - Bulk image structure metric (S/N of large-scale features)
   - Flags high-variance intervals for quality filtering
   - Stored as `CLOUD_S` header keyword

4. **Frame Tracking Monitoring:**
   - Calculates RA/Dec shift from previous frame
   - Identifies telescope tracking errors or guiding issues
   - Headers: `RA_MOVE`, `DEC_MOVE` (arcseconds), `SKY_MOVE` (total shift)

**Output:**
- `{image}.phot` files for each science frame (FITS tables)
  - Extension 1: Photometry table with aperture fluxes (13 apertures per source)
  - Headers: SEEING, FWHM, CLOUD_S, PSF parameters, frame shifts
- `{filelist}_phot` - List of successfully processed frames
- `{GAIA_ID}_stack_catalogue_pm.fits` - PM-corrected reference catalogue

**Parameters:**
- `--apsize`: Aperture radius in pixels (default: 4, typically matches FWHM)
- `--nproc`: Parallel processes (default: 1 for debugging, 8-20 for production)
- `--norunwcs`: Skip WCS refinement (use existing solutions)

**Quality Validation:**
- Removes frames with failed WCS
- Logs frames with unusual PSF or tracking anomalies
- All metrics available for downstream filtering

---

### Stage T10: Differential Photometry

**Module:** `src/diff_photometry/SPlightcurve.py`

Generates final calibrated light curves through differential photometry. Selects optimal comparison stars, performs iterative weighting, and applies atmospheric corrections.

**Process:**

1. **Target Identification:**
   - Matches target to Gaia DR3 source ID
   - Retrieves effective temperature from Filippazzo catalogue or FITS header
   - Locates target in photometry catalogues

2. **Comparison Star Selection:**
   - Filters candidates by:
     - Brightness (within 2 mag of target)
     - Low intrinsic variability (RMS < threshold)
     - Angular separation (avoid close companions)
     - No saturation or flags
   - Typical selection: 10-50 comparison stars per field

3. **Iterative Differential Photometry:**
   - Creates weighted artificial comparison: `C(t) = Σ wᵢ × Fᵢ(t)`
   - Computes differential flux: `flux_diff = F_target / C(t)`
   - Iteratively adjusts weights to minimise RMS
   - Algorithm converges in 5-10 iterations
   - Poor comparison stars automatically down-weighted

4. **PWV Atmospheric Correction** (optional):
   - Only applied for cool stars (2000K < Teff < 3200K) at Paranal
   - Loads pre-computed atmospheric model grid
   - Interpolates correction using:
     - Real-time PWV from ESO monitor
     - Observation airmass
     - Target/comparison Teff differential
   - Applies correction: `flux_corr = flux_diff × correction_factor`
   - Most significant in I+z filter where H₂O bands overlap
   - Effects: 0.1-1% depending on PWV conditions

5. **Bad Weather Masking:**
   - Identifies high-variance intervals (clouds, seeing variations)
   - Automatically flags suspect data points
   - Quality flags propagated to output

6. **Time Binning:**
   - Bins light curve to specified interval (default: 5 minutes)
   - Uses median for robustness against outliers
   - Calculates binned uncertainties

**Output:**

**FITS Files:**
- `{GAIA_ID}_{filter}_{aperture}_diff.fits` containing:
  - Target differential light curve
  - PWV-corrected light curve (if applicable)
  - Artificial comparison light curve
  - Individual comparison star light curves
  - Stellar catalogue (Gaia parameters, Teff, colours)
  - Image metadata (airmass, FWHM, seeing, etc.)
  - Quality flags and weights

**ASCII Files:**
- `{GAIA_ID}_{filter}_{aperture}_MCMC` - Light curve for transit modelling:
  ```
  # JD          BJD         Flux    Error
  2459123.456  2459123.459  0.9995  0.0003
  2459123.460  2459123.463  1.0001  0.0003
  ...
  ```

**Plots:**
- Unbinned and binned light curves
- PWV-corrected light curves
- Multi-night stitched light curves
- Diagnostic plots: airmass vs flux, FWHM vs flux, background trends
- PWV comparison (original/processed/interpolated)
- RMS vs magnitude for comparison star quality assessment

**Parameters:**
- `--ap`: Aperture number (3-8, corresponding to different radii)
- `--binning`: Time binning in minutes (default: 5)
- `--teff`: Target effective temperature (auto-retrieved if not provided)

**Quality Metrics:**
- Point-to-point scatter (unbinned RMS)
- Binned RMS (typically 0.1-0.5% for bright targets)
- Comparison star weights
- PWV correction effectiveness (RMS before/after)

---

### Stage T11: PDF Report Generation

**Module:** `src/reporting/pdf_report_catriona.py`

Generates comprehensive PDF reports summarising nightly processing results, data quality, and light curve characteristics.

**Contents:**
- Pipeline execution summary
- Master calibration frame statistics
- Per-target photometry quality metrics
- Light curve thumbnails
- Data completeness assessment

**Output:**
- `{TELESCOPE}_{DATE}_pipeline_report.pdf`

---

### Stage T12: Version Migration

**Implementation:** Handled in `ZLP_pipeline.sh`

Manages migration of successful processing results from development version (v3) to production version (v2). Only migrates nights where differential photometry completed successfully.

**Process:**
1. Checks for presence of `*_diff.fits` files (indicates successful pipeline completion)
2. If found:
   - Moves v2 data to v2_old (backup)
   - Copies v3 results to v2 (production)
   - Moves calibration frames, reduced data, reports, logs
3. If not found:
   - Leaves v2 unchanged
   - v3 results remain in development tree

**Purpose:** Allows testing pipeline improvements in v3 while maintaining stable v2 for science operations.

---

## Key Concepts

### Differential Photometry

Dividing target flux by weighted combination of comparison stars cancels systematic effects:

```
flux_diff(t) = F_target(t) / Σ[wᵢ × Fᵢ(t)]
```

**Removes:**
- Atmospheric extinction variations
- Telescope throughput changes
- Detector sensitivity drifts
- Guiding errors (if field rotates)

**Requirements:**
- Stable comparison stars (low intrinsic variability)
- Similar brightness to target (±2 mag)
- Spatially distributed across field

**Typical Performance:**
- Bright targets (I < 12 mag): 0.1-0.3% precision (5-min bins)
- Faint targets (I > 14 mag): 0.5-1.0% precision

---

### Proper Motion Correction

Stars move across the sky due to transverse velocities. Without correction:

| PM (mas/yr) | Δt (years) | Pixel shift (0.35"/pix) | Impact |
|-------------|------------|------------------------|--------|
| 10 | 3 | 0.09 | Negligible |
| 50 | 3 | 0.43 | Centroiding error |
| 100 | 3 | 0.86 | Significant flux loss |
| 500 | 3 | 4.3 | Star leaves aperture |

The pipeline automatically corrects catalogue positions from Gaia epoch (J2016.0) to observation date using Gaia-measured proper motions.

---

### Local Gaia Database

The pipeline uses a local copy of Gaia DR3 cross-matched with 2MASS for fast astrometric calibration:

**Database Structure:**
- Sharded by declination (1° bins: -90_-89, -89_-88, ..., 89_90)
- Contains: RA, Dec, pmra, pmdec, parallax, G-mag, J-mag
- Size: ~50 GB
- Location: `/gaia_database/gaia_tmass_16_jm_cut.db`

**Query Performance:**
- Local: 0.1-0.5s per field
- VizieR: 5-30s per field (network dependent)

**Advantage:** Enables rapid plate solving of 1000+ images per night without network dependency.

---

### Bad Pixel Handling

**ANDOR (Optical):**  
Hot pixels are rare (<0.1%) and don't significantly impact source detection. No correction applied.

**SPIRIT (NIR):**  
Hot pixels are common (1-3%) and comparable in brightness to faint stars. Three-step correction:

1. **Detection:** Sigma-based thresholding in darks and flats
2. **Masking:** Set bad pixels to NaN
3. **Interpolation:** Replace with median of surrounding 8 pixels

**Impact on Plate Solving:**  
Without correction, hot pixels are detected as "stars" by source finder, causing false pattern matches and WCS failures. With correction, plate solve success rate improves from ~60% to >85% for SPIRIT.

---

### PWV Atmospheric Correction

Water vapour in Earth's atmosphere absorbs starlight, particularly in near-infrared bands. The absorption varies with:
- **PWV amount** (mm of precipitable water)
- **Airmass** (path length through atmosphere)
- **Stellar temperature** (spectral shape affects overlap with H₂O bands)

**Effect in I+z filter:**
- Cool stars (2600K): Strong absorption in H₂O bands
- Hot stars (5000K): Minimal absorption (different spectral shape)
- Differential effect: 0.1-1% depending on PWV

**Correction Method:**
1. Pre-compute atmospheric transmission models for grid of (PWV, airmass, Teff)
2. Query real-time PWV from ESO monitor at Paranal
3. Interpolate correction factor for target/comparison temperature difference
4. Apply to differential light curve

**Applicability:**
- Only Paranal observations (ESO PWV data available)
- Only cool targets (2000-3200K) where effect is significant
- Only optical filters with H₂O band overlap (I+z, z')

---

## Quality Monitoring

The pipeline generates comprehensive quality metrics tracked over time:

### Per-Night Metrics

**Calibration Quality:**
- Readout noise (e⁻) - should be stable over time
- Dark current (e⁻/s) - temperature dependent
- Overscan level (ADU) - bias offset stability

**Astrometric Quality:**
- Plate solve success rate (%) - should be >90%
- WCS RMS (arcsec) - typical: 0.1-0.3"
- Number of matched Gaia sources

**Photometric Quality:**
- Seeing (arcsec) - weather dependent
- Zero point (mag) - transparency indicator
- Gaia catalogue match percentage

### Per-Target Metrics

**Light Curve Quality:**
- Point-to-point RMS (%) - unbinned scatter
- Binned RMS (%) - N-minute bins
- Comparison star weights - stability indicator

**Data Completeness:**
- Number of epochs
- Time baseline
- Duty cycle

### Accessing QC Data

Quality control metrics are stored in `/data/SPECULOOSPipeline/PipelineOutput/{VERSION}/{TELESCOPE}/reports/QC/`.

View time-series plots:

```bash
ls /data/SPECULOOSPipeline/PipelineOutput/v3/Io/reports/QC/*.png
```

---

## Troubleshooting

### Common Issues

#### 1. Plate Solving Failures

**Symptoms:** Many images with "FAILED" in plate solve logs.

**Causes:**
- Poor image quality (clouds, out of focus)
- Incorrect telescope/filter in FITS headers
- Hot pixels (SPIRIT without bad pixel correction)
- Wrong plate scale in configuration

**Solutions:**
```bash
# Check plate solve logs
grep "FAILED" logs/20260128_1_v3.log

# Force re-solve with higher timeout
./main/ZLP_pipeline.sh --force-platesolve --only_T7 1 /data/SPECULOOSPipeline 20260128 2 20 Io

# Verify WCS keywords in failed images
fitsheader {failed_image}.fits | grep -E "RA|DEC|CRVAL|CRPIX"
```

---

#### 2. No Comparison Stars Found

**Symptoms:** Light curve generation fails with "insufficient comparison stars" error.

**Causes:**
- Sparse field (few stars)
- Target very bright or very faint
- Poor Gaia cross-match (<50%)

**Solutions:**
```bash
# Check stack catalogue quality
fitsheader {stack_catalogue}.fits | grep NAXIS2  # Number of sources

# Examine Gaia cross-match
fitsheader {stack_catalogue}.fits | grep -i gaia

# Adjust comparison star selection criteria in SPlightcurve.py
# - Relax magnitude range
# - Accept more distant comparison stars
```

---

#### 3. Memory Errors During Stacking

**Symptoms:** Pipeline crashes at T8 with memory allocation errors.

**Causes:**
- Too many input frames (>100)
- Large image dimensions
- Insufficient RAM

**Solutions:**
```bash
# Reduce number of stacked frames
./main/ZLP_pipeline.sh 1 /data/SPECULOOSPipeline 20260128 2 20 Io
# Edit ZLP_create_cat.py: nfiles=30 (default: 50)

# Process serially to reduce memory pressure
./main/ZLP_pipeline.sh --cores 1 --only_T8 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

---

#### 4. Catastrophic WCS Failures

**Symptoms:** All images in a night fail plate solving.

**Causes:**
- Telescope not pointing at header coordinates
- Wrong telescope in FITS header
- Corrupted Gaia database

**Solutions:**
```bash
# Verify FITS headers
fitsheader {image}.fits | grep -E "TELESCOP|RA|DEC"

# Test Gaia database access
python -c "import sqlite3; conn = sqlite3.connect('/gaia_database/gaia_tmass_16_jm_cut.db'); print('DB OK')"

# Run plate solving in verbose mode
# Edit astrometry.py: verbose=True
./main/ZLP_pipeline.sh --only_T7 --cores 1 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

---

#### 5. PWV Correction Failures

**Symptoms:** `water_vapour.py` crashes or produces corrupted light curves.

**Causes:**
- PWV data unavailable for observation date
- Target temperature outside grid range
- Corrupted atmospheric model files

**Solutions:**
```bash
# Check PWV data availability
ls /data/SPECULOOSPipeline/water_grid_I+z_*.npy

# Verify target temperature
fitsheader {diff_file}.fits | grep TEFF

# Skip PWV correction if not critical
# Comment out PWV correction call in SPlightcurve.py
```

---

### Log File Locations

| Log Type | Path | Contents |
|----------|------|----------|
| Pipeline execution | `PipelineOutput/{VERSION}/{TELESCOPE}/logs/{DATE}_*.log` | Complete stdout/stderr |
| Plate solve summary | Embedded in pipeline log | Success/failure per image |
| Photometry errors | Embedded in pipeline log | Failed frames, quality issues |
| Docker logs | `docker logs orchard-server` | Container startup/errors |

---

### Getting Help

**Check Logs First:**

```bash
# View most recent log
ls -lt /data/SPECULOOSPipeline/PipelineOutput/v3/*/logs/*.log | head -1 | xargs tail -100

# Search for errors
grep -i "error\|exception\|failed" /data/SPECULOOSPipeline/PipelineOutput/v3/Io/logs/20260128_*.log
```

**Validate Setup:**

```bash
# Check Docker container
docker ps | grep orchard

# Verify mounts
docker exec orchard-server ls /data/SPECULOOSPipeline
docker exec orchard-server ls /gaia_database

# Test Python environment
docker exec orchard-server python -c "import astropy, photutils, fitsio; print('Environment OK')"
```

---

## Development

### Code Structure

```
src/
├── main/                      # Pipeline orchestration
│   ├── ZLP_pipeline.sh        # Main pipeline script
│   └── createlists.py         # Image classification
│
├── calibration/               # Calibration stages (T2-T6)
│   ├── instrument_config.json # Telescope/instrument parameters
│   ├── pipebias.py           # Master bias creation
│   ├── pipedark.py           # Master dark creation
│   ├── pipeflat.py           # Master flat creation
│   ├── pipebadpixel.py       # Bad pixel detection
│   ├── pipered.py            # Science reduction
│   └── pipeutils.py          # Common calibration utilities
│
├── astrom/                    # Astrometry (T7)
│   ├── astrometry.py         # WCS solving orchestration
│   ├── pointer_wcs.py        # Plate solving implementation
│   ├── proper_motion.py      # PM corrections
│   └── twirl_speculoos.py    # Pattern matching wrapper
│
├── photometry/                # Photometry (T8-T9)
│   ├── ZLP_create_cat.py     # Stack & catalogue creation
│   ├── ZLP_app_photom.py     # Aperture photometry orchestration
│   ├── wcs_photom.py         # Forced photometry implementation
│   ├── casutools.py          # CASU pipeline wrappers
│   ├── gaia_dr2_test.py      # Gaia cross-matching
│   └── photom_quality_checks.py  # Image quality metrics
│
├── diff_photometry/           # Differential photometry (T10)
│   └── SPlightcurve.py       # Light curve generation
│
├── water_vapour/              # PWV correction
│   ├── water_vapour.py       # Atmospheric correction
│   └── pwvGrid.py            # Model grid generation
│
├── reporting/                 # Reports & QC (T11)
│   ├── pdf_report_catriona.py
│   ├── QC.py                 # Quality metrics tracking
│   └── email_eon.py          # Status notifications
│
├── condense/                  # Data aggregation
│   └── zlp_condense.py       # Photometry file merging
│
└── utils/                     # Helper utilities
    ├── correct_target_names.py
    ├── gaia_id_from_schedule.py
    └── already_run_targs.py
```

---

### Adding a New Telescope

1. **Add instrument configuration** to `src/calibration/instrument_config.json`:

```json
{
  "telescopes": {
    "NEW-TELESCOPE": {
      "camera_type": {
        "gain": 1.0,
        "saturation_threshold": 60000,
        "naxis2": 2048,
        "dark_matching_exptime": false,
        "overscan": {
          "top_row": null,
          "bottom_row": null,
          "left_col": null,
          "right_col": null
        },
        "trim": {
          "top_row": 1,
          "bottom_row": -1,
          "left_col": 1,
          "right_col": -1
        }
      }
    }
  }
}
```

2. **Create data directories:**

```bash
mkdir -p /data/SPECULOOSPipeline/Observations/NEW-TELESCOPE/images
mkdir -p /data/SPECULOOSPipeline/PipelineOutput/v3/NEW-TELESCOPE/output
mkdir -p /data/SPECULOOSPipeline/PipelineOutput/v3/NEW-TELESCOPE/reports
mkdir -p /data/SPECULOOSPipeline/PipelineOutput/v3/NEW-TELESCOPE/logs
```

3. **Test with sample data:**

```bash
./main/ZLP_pipeline.sh --cores 1 1 /data/SPECULOOSPipeline {DATE} 2 20 NEW-TELESCOPE
```

4. **Calibrate plate scale** using field with known astrometry.

---

### Modifying Pipeline Stages

**Example: Adjust Detection Threshold**

Edit `src/photometry/ZLP_create_cat.py`:

```python
# Line ~200: Change detection threshold for stack catalogue
casutools.imcore(outstack_name, outcatname,
                 threshold=2,  # Change from 2 to 1.5 for deeper detection
                 confidence_map=outstackconf_name,
                 verbose=verbose,
                 ellfile=True)
```

**Testing:**

```bash
# Test on single night with modified code
./main/ZLP_pipeline.sh --only_T8 1 /data/SPECULOOSPipeline 20260128 2 20 Io

# Compare source counts
fitsheader {old_catalogue}.fits | grep NAXIS2
fitsheader {new_catalogue}.fits | grep NAXIS2
```

---

### Contributing

The pipeline is currently maintained for SPECULOOS operations at Cambridge. For modifications:

1. **Create feature branch** from `main`
2. **Test thoroughly** on representative data
3. **Document changes** in code comments and commit messages
4. **Submit for review** before merging to main

**Development Workflow:**

```bash
# Create feature branch
git checkout -b feature/improved-wcs-fitting

# Make changes and test
./main/ZLP_pipeline.sh --cores 1 1 /data/SPECULOOSPipeline 20260128 2 20 Io

# Commit with descriptive message
git add src/astrom/pointer_wcs.py
git commit -m "Improve WCS fitting robustness for sparse fields

- Add fallback to linear interpolation when cubic fails
- Increase Gaia query radius for edge cases
- Add diagnostic output for failed matches"

# Push and create PR
git push origin feature/improved-wcs-fitting
```

---

## Performance

### Typical Processing Times

| Stage | Single Core | 20 Cores | Notes |
|-------|-------------|----------|-------|
| T1: Classification | 10s | 10s | I/O bound |
| T2: Master Bias | 30s | 30s | Memory bound |
| T3: Master Dark | 45s | 45s | Memory bound |
| T4: Master Flat | 60s | 60s | Memory bound |
| T5: Bad Pixel Map | 15s | 15s | Computation light |
| T6: Reduction | 20 min | 3 min | **Highly parallel** |
| T7: Plate Solving | 30 min | 5 min | **Highly parallel** |
| T8: Stack Creation | 5 min | 2 min | Mixed |
| T9: Photometry | 15 min | 3 min | **Highly parallel** |
| T10: Diff Photom | 2 min | 2 min | Serial per target |
| T11: PDF Report | 30s | 30s | Light computation |
| **Total** | **~75 min** | **~15 min** | For typical night |

**Assumptions:** 1000 science frames, 5 targets, Io telescope, typical observing conditions.

---

### Bottlenecks

**Calibration Stages (T2-T5):**  
Memory-intensive operations (median of large arrays). Limited parallelisation benefit. Rarely the bottleneck for single-night processing.

**Plate Solving (T7):**  
CPU-intensive pattern matching. Highly parallel. Bottleneck for large datasets. Consider running overnight for multi-night reprocessing.

**Stack Creation (T8):**  
I/O and memory intensive. WCS solving during stack creation can be slow. Consider using pre-solved images or reducing number of stacked frames.

---

### Optimisation Strategies

**For Regular Operations:**
```bash
# Use available cores
./main/ZLP_pipeline.sh --cores 20 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

**For Reprocessing:**
```bash
# Skip unnecessary stages
./main/ZLP_pipeline.sh --no_T2 --no_T3 --no_T4 --no_T5 --cores 20 \
  1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

**For Debugging:**
```bash
# Serial processing for clear error messages
./main/ZLP_pipeline.sh --cores 1 --only_T7 1 /data/SPECULOOSPipeline 20260128 2 20 Io
```

---

## Acknowledgements

Orchard builds upon work by numerous contributors to the SPECULOOS project and relies on several open-source astronomical software packages:

**Core Dependencies:**
- **CASU Tools** - Cambridge Astronomical Survey Unit pipeline software
- **Astropy** - Community Python library for astronomy
- **twirl** - Fast pattern matching for astrometric calibration
- **Gaia DR3** - ESA mission providing astrometric catalogue
- **2MASS** - Two Micron All Sky Survey photometric catalogue

**Original Pipeline Development:**  
The SPECULOOS photometry pipeline was originally developed by Catriona Murray and others as part of their PhD research. This modernised version (Orchard) maintains the core scientific algorithms while improving robustness, configurability, and deployment.

**SPECULOOS Team:**  
The SPECULOOS project is led by Michaël Gillon (University of Liège) with contributions from institutions worldwide.

---

## License

Copyright © 2025 University of Cambridge

Licensed under the MIT License. See LICENSE file for details.

---

## Version History

**v3.0** (Current - February 2025)
- Complete Python 2.7 → 3.9 migration
- Dockerised deployment
- JSON-based instrument configuration
- Local Gaia DR3 database integration
- NIR camera support (SPIRIT)
- Improved error handling and logging
- Modular architecture for easier maintenance

**v2.0** (Historical)
- Production pipeline for SPECULOOS
- Python 2.7 codebase
- Hardcoded instrument parameters
- Online catalogue queries

---

## Contact

For questions about the pipeline or SPECULOOS data:

**Pipeline Maintainer:**  
Matthew Hooton  
Cavendish Laboratory, University of Cambridge  
mh2143@cam.ac.uk

**SPECULOOS Project:**  
https://www.speculoos.uliege.be/

**Data Access:**  
Light curves and observing logs available through the SPECULOOS data portal.