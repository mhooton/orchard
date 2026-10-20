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

### Telescope Name Aliases

Settings are looked up by the frame's `TELESCOP` header. A frame that names its
telescope differently from the configuration is mapped on by `telescope_aliases`:

```json
"telescope_aliases": {
  "ACP->Artemis": "Artemis"
}
```

Artemis frames say `ACP->Artemis` from the start of the survey until 2022-06-15
(ACP's telescope hub was in the chain) and `Artemis` from 2022-06-16. It is the same
camera either side: the flat-field pixel response correlates at 0.92 across the change,
as it does between neighbouring nights. Without the alias, older frames fall back to
generic defaults (no trim or overscan, no plate scale) and T8 stops at the WCS check.

### Gain, Read Noise and Saturation per Night

For the Andors, `gain` and `saturation_threshold` above are only the last fallback. No raw
header records the readout mode, and the cameras have been read out in several (measured
October 2026):

| Mode | Gain (e⁻/ADU) | Bias noise (ADU) | Saturation used (raw ADU) | Nights |
|------|---------------|------------------|---------------------------|--------|
| `NORMAL` | 0.97–1.13, by telescope and date | 5.7–6.9 | 62,000 | most |
| `GAIN3.8` | 3.6–3.9 | 2.4–2.6 | Europa 34,000; Artemis 37,300 | Europa 2019-07-25 to 10-26; Artemis 2020-08-06/07 and from 2025-10-29 |
| `FAST` | 1.19 | 9.0 | 59,500 | Ganymede 2026-02-27 to 03-01 |
| `PED7500` | 1.32 | 17.6, bias level 7,546 | 39,500 | Ganymede 2026-03-02/03 |
| `PED2900` | ~6.5, estimated | 6.3, bias level ~2,900 | 22,800 | Artemis 2026-08-06/07 |

`src/calibration/readout.py` gives each calibrated frame its own values:

1. its mode, from the noise and level of its own overscan strip (the mode can change
   during a night);
2. that mode's read noise and gain, measured from the night's own bias pairs and
   twilight-flat pairs (frames borrowed from other nights are not used);
3. failing that, the dated table `REGIMES` in `readout.py`, and only then this configuration.

T2 writes the night's measurements to `reduction/readout.json`.
`python -m calibration.readout <reduction dir>` prints them for any night without writing
anything. Instruments without an overscan (SPIRIT, MOANA) keep the configured values.

To do: nothing yet collects these values over time or flags a change of mode
([#31](https://github.com/mhooton/orchard/issues/31)).

### Instrument-Specific Features

**ANDOR (Optical):**
- Combined dark frame (scaled by exposure time)
- Overscan correction with sigma clipping
- 2044×2048 pixel trimmed output
- No bad pixel correction
- Readout-smear removal for shutterless frames (SSO ANDORs, `"smear"` block: readout `axis`, trail `direction`, `readout_time` in s, `n_shifts`)

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
| `--force-platesolve` | T7: ignore the night's WCS cache, solve every frame and refresh its cache rows ([docs/wcs-cache.md](docs/wcs-cache.md)) | Disabled |
| `--save-flat-products` | T4: also write each calibrated flat (`reduction/flats/proc*`), `std.fts` and `variance.fts`. Nothing in the pipeline reads them | Disabled |
| `--frame-format F` | After T11, store the night's calibrated frames as `int16` (whole ADU, half the size), `rice` (int16 plus lossless Rice compression, `proc*.fits.fz`) or `float32` (as T6 wrote them) | `int16` |
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

#### Job Queue on appct

`src/jobqueue/` is a priority queue for pipeline runs on appct, with an ESO watcher that queues each night as
soon as its frames are at ESO and a daily 30-night look-back. It runs inside `orchard-server`:

```bash
docker exec orchard-server python -m jobqueue status
```

Design, shadow mode, the proposed crontab and the switch-over are in [docs/job-queue.md](docs/job-queue.md).

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
4. Measure the night's read noise and gain per readout mode from its own bias and flat pairs
   (`readout.py`, see [Gain, Read Noise and Saturation per Night](#gain-read-noise-and-saturation-per-night))
5. Write master bias and metrics

**Output:**
- `{run}_MasterBias.fits`
- `readout.json` (the night's read noise and gain per readout mode, and how they were found)
- `readoutnoise.dat` (in electrons, for the mode of the night's bias frames)
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
5. Calculate dark current: `DC = gain × median(dark)`, with the night's gain

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
2. De-smear flats taken without the shutter (see [Shutterless frames](#shutterless-frames-readout-smear)); the master records how many in `NSMEARED`
3. Normalise each flat by its median flux
4. Sigma-clipped median combination (rejects cosmic rays)
5. Normalise master to unity mean
6. With `--save-flat-products` only: write each calibrated flat and the per-pixel standard deviation and variance across the flats

**Output:**
- `{run}_MasterFlat_{filter}.fits` (one per filter)
- With `--save-flat-products`: `flats/proc*` (each calibrated flat), `std.fts` and `variance.fts`. They are off by default because nothing reads them; together they took 3.8 TB of PipelineOutput by October 2026

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
4. Readout-smear removal for frames taken without the shutter (SSO ANDORs only)
5. Flat field division (filter-specific)
6. Bad pixel interpolation (SPIRIT only)
7. CRPIX adjustment for trimming

**Process:**
```python
signal = raw - overscan - bias - dark × exptime
signal = desmear(signal) if shutterless else signal   # before the flat
reduced = signal / flat
reduced = interpolate_bad_pixels(reduced, bad_pixel_map)
```

**Output:**
- `proc{image}.fits` for each science frame
- `{filter}_processed.dat` lists of successfully reduced images
- FITS headers updated with calibration history, and with the frame's own readout values:
  `GAIN` (e⁻/ADU), `RDNOISE` and `RON` (e⁻), `SATURATE` (ADU in the calibrated frame, i.e. raw
  saturation minus the overscan and bias level), `READMODE`, `GAINSRC` (`night`, `table`,
  `estimate` or `config`) and `OSNOISE` (overscan noise of the raw frame, ADU).
  CASU reads `GAIN` for every flux error; condense copies the values into `output.fits`
  and T10 uses them for the noise model, the I+z magnitudes and the saturation flags.

**Parallel Processing:** Supports multiprocessing for faster execution on multi-core systems.

#### Shutterless frames (readout smear)

**Module:** `src/calibration/smear.py`

Without a mechanical shutter the CCD keeps integrating while it is read out, so every source leaves a one-sided trail along the readout axis: each pixel downstream of a source receives `t_read / (N × t_exp)` of that source's flux in the same row. On the SSO ANDORs the trail runs along +x (NAXIS1), `t_read` = 4.33 s and N = 2048, so the trail is 1.1×10⁻⁴ of the source per pixel at 20 s and the sky picks up a ramp reaching ~22% of the sky level at the far edge. `desmear()` inverts this exactly, row by row, after bias/dark subtraction and **before** the flat (the trail charge was collected by the upstream pixels, with their sensitivity).

Whether a frame was shutterless is decided per frame:
1. a row in `src/calibration/shutter_overrides.csv` for that `TELESCOP` whose UTC window holds `DATE-OBS`, else
2. the `CAMSHUTT` header (written by Astra on all ANDORs since 2026-07-17, but set by hand, so it can be wrong; a value that is not a recognisable true/false is reported and taken as shutter used), else
3. shutter used.

The readout time of a shutterless frame comes from, in order: the override window, if the row sets `readout_time`; Astra's `READOUTM` card, written on every Andor frame since 2026-10-18 (Astra 1.3), an index into 5 MHz / 3 MHz / 1 MHz / 0.05 MHz looked up in `"readout_rate_time"` under `"smear"` in `instrument_config.json` (1 MHz 4.33 s, 3 MHz 1.5 s); else the readout mode T6 infers from the overscan, `READMODE`, looked up in `"readout_time"` (`NORMAL` 4.33 s, `FAST` 1.5 s, the 3 MHz readout of Ganymede 2026-02-27 to 03-01). A frame whose `READOUTM` names 5 or 0.05 MHz (never measured; the pixel-clock arithmetic predicts ~0.94 s and ~86 s), or a shutterless frame in any other `READMODE` (`PED7500`, `GAIN3.8`, `UNKNOWN`, …), is **not** de-smeared, because a wrong `t_read` would add a ramp of its own: it gets `SMEARCOR = F` and the reason in `SMEARERR`, and T4/T6 print a warning. Add a measured readout time to the config (or to the override window) before reprocessing such a night. T4 does the same for flats, and a master flat counts the shutterless inputs it could not de-smear in `NSMRLEFT`.

Every T6 product of an instrument with smear parameters records the decision: `SMEARCOR` (T/F), and when corrected `SMEARTRD` (t_read used; the comment says whether from the override or `READMODE`), `SMEARTAU`, `SMEARSRC` (override/CAMSHUTT), `SMSATLEV` and `SMSATROW` (rows containing a clipped pixel). A pixel counts as clipped at the frame's own saturation level, `SATURATE` plus the pedestal T6 subtracted (overscan + master-bias median), which is far below 62,000 ADU in Ganymede's fast and 7,500-pedestal modes; frames without `SATURATE` use the config's `saturation_threshold`. Shuttered frames gain only `SMEARCOR = F`; their pixels are unchanged.

**Known limit, flagged rather than corrected:** a saturated pixel's lost charge is not recorded, so its trail is under-corrected in every pixel downstream of it in that row. T6 records those rows in the frame header rather than storing a mask: `SMTRAXIS`/`SMTRDIR` give the geometry and `SMTR0001`, `SMTR0002`, … hold runs of consecutive rows as `'L:P,P,…'` (first row of the run, then the x of the first clipped pixel in each row, 1-based). That is a few cards per saturated star. T9 (`condense_photometry.py`) turns these into `SMEAR_AP` in `output.fits`, one int16 per star and frame: the smallest of the 13 apertures that overlaps a trail pixel. It is 0 if none does and −1 if not assessed, so a star is affected at aperture *j* when `0 < SMEAR_AP <= j`. T10 never uses a star that is affected in any frame as a comparison star (`SMEARED` in the weights table). It marks the target's affected frames in `SMEAR_FLAG` (the `FLAGS` table and the last column of the `_MCMC` file) and keeps them in the light curve. When choosing the best aperture, it never picks one where the target is on a trail if another aperture is clear; `*_bestap.txt` gains a third field, 1 when the chosen aperture is on a trail. How big the residual is, measured on Ganymede 20250408: behind a G = 8.8 star it is ~300 ADU per pixel across its 12–14 clipped rows, 17× the sky noise. Nothing shows upstream, or in the same rows on the shuttered night before. That puts 2–10% into TOI-1227's apertures 7–8.

**Checking a night:** `src/calibration/smear_check.py` measures the trail behind bright isolated stars (or the ramp across a flat) and says whether frames were shutterless and what `t_read` they imply. Use it before adding an override window — a wrong window de-smears shuttered frames, and on a 7 s flat that puts a ~60% ramp into the master flat:

```bash
python calibration/smear_check.py @frames.txt --desmear
```

An asymmetric PSF fakes a small trail (+0.3 s on Ganymede at 20 s, −1.5 s on Io at 39 s with the default windows), which is negligible for telling shutterless from shuttered but matters when measuring `t_read` itself.

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

4. **WCS Output:**
   - Writes the WCS keywords to the calibrated frame's header (never to the raw frame)
   - Keeps the night's solutions in the WCS cache, `<BASEDIR>/PipelineOutput/wcs_cache/<TEL>/<night>.ecsv`

**WCS cache** ([docs/wcs-cache.md](docs/wcs-cache.md)): before solving, T7 applies the solution an earlier run kept
for a frame, header only, when the raw frame, the calibrated frame's shape, the solver (`SOLVER_EPOCH`), the Gaia
database and astropy/twirl are unchanged and a code release or the running commit wrote it. Only the other frames
are solved; their solutions are added to the table afterwards. A WCS already in a frame's header (Astra's,
PinPoint's, ESO's) never counts as solved. Failed frames are not kept, and any fault in the cache only means the
frame is solved.

**Output:**
- Updated FITS headers with WCS keywords (CRVAL, CRPIX, PC matrix), and HISTORY `T7: WCS solved` or
  `T7: WCS from cache`
- The night's WCS cache table
- Plate solve success/failure logged to console, with a `WCS cache:` line for frames taken from the cache
- Diagnostic plots (if verbose mode enabled)

**Parameters:**
- `--force-platesolve`: ignore the WCS cache, solve every frame and refresh its rows (SSO pipeline jobs pass it
  until the cache is validated: `force_platesolve` in `jobqueue/config.py`)
- `--timeout`: Maximum time per image (default: 60s)
- `--cores`: Number of parallel processes
- `WCS_CACHE_DIR` (environment): keep the cache tables under this directory instead

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
   - Matches detected sources to Gaia DR3 from the local database (`GAIADATABASEPATH`),
     propagating every catalogue position from J2016.0 to the night with its proper motion
   - Retrieves proper motions, parallaxes, G-band magnitudes
   - Stores cross-match in `Gaia_Crossmatch` FITS extension (both `GAIA_DR2_ID` and `GAIA_DR3_ID`)
   - Also returns hand-added objects from the database's `custom_sources` table (see below)
   - Reports match percentage for quality assessment

5. **Target Identification** (`utils/target_management.py`):
   - Reads the primary target's Gaia ID from the night's schedule plan file, then applies
     `calibration/target_id_aliases.csv` for plan files known to carry a wrong ID
   - Matches that ID against **either** `GAIA_DR2_ID` or `GAIA_DR3_ID` of the crossmatch
   - If nothing matches, looks the ID up directly in the local database and logs one of two
     distinct statuses: `PRIMARY_NOT_IN_DB` (the star is missing from our database, a build
     problem that will fail on every night) or `PRIMARY_NOT_IN_FOV` (the star is in the
     database but not in tonight's catalogue)
   - Coordinate fallback: propagates the star's position to the night with its proper motion
     (from the database, else the plan-file J2000 position without proper motion) and accepts
     the nearest detection within 2 arcsec, provided it carries no conflicting Gaia ID. The
     schedule ID and any database astrometry/photometry are then injected into that catalogue
     row so every later stage keys on it as usual
   - Falls back to the TOI table for `TOI-*` names, and always adds 40 pc target-list stars
     in the field as secondary targets
   - Records `PRIM_HOW`, `PRIM_ST`, `PRIM_ID` and `SCHED_ID` in the `Gaia_Crossmatch` header
     so the method used on each night is auditable

6. **Backup System:**
   - Preserves best stack catalogue when subsequent runs fail
   - Automatically restores backup if current run produces poor results
   - A backup is validated by content (the target must be in it, by DR2 or DR3 ID) and is
     rejected when more than 30 days old for a target without proper motion, because catalogue
     positions are propagated to the night and a row without proper motion stays put
   - Useful for nights with degraded data quality

**Targets missing from Gaia DR2, from the local database, or from Gaia altogether**

Wolf 359 (Sp1056+0700) is the cautionary example: it has no Gaia DR2 entry, so the scheduler
wrote the ID of a G=16 neighbour into its 2021 to 2022 plan files, and the database build
dropped it because it required a DR2 counterpart. Every product from four years of data was
either the neighbour's light curve or no photometry at all. Three cases are handled:

- *In DR3 but not DR2* (bright, very high proper-motion nearby stars): keep the database
  built with `gaia-tmass-sqlite/db_maker.py` at or after the `g_cut_dr3_only` change, or run
  `augment_dr3_only.py` on an existing database. Nothing else is needed.
- *Name resolves to no ID* (missing from schedule, TOI table and target list): the coordinate
  fallback above finds it. 59 rows of `ml_40pc.txt` carry a zero ID;
  `utils/resolve_master_list_gaia_ids.py` identifies them in DR3 and its report feeds the
  generated table (see [The 40 pc Target Table](#the-40-pc-target-table)). Its 12 *weak*
  matches — neither a proper motion nor a 2MASS counterpart, so a neighbour is as likely as
  the target — are recorded but deliberately left unidentified; `--allow-weak` overrides that
  once a human has reviewed them.
- *Not in Gaia at all* (faint T/Y dwarfs): add a row to
  `calibration/supplementary_sources.csv` with a synthetic numeric ID above 9e18 and load it
  with `python -m utils.supplementary_sources load`. Estimate G and colours from the spectral
  type; Teff feeds the PWV correction.
- *In Gaia DR2 and DR3 but left out by the database build* (no DR2 2MASS match, or a DR2/DR3
  strip boundary; T8 logs `PRIMARY_NOT_IN_DB`): add the row the build would have written to
  `calibration/gaia_db_manual_rows.csv` and run `python -m utils.gaia_db_manual_rows apply`
  (see [Local Gaia Database](#local-gaia-database)). The CSV's header defines a code for each
  non-standard choice (proper motion from DR2 for a 2-parameter DR3 solution, from 2MASS when
  Gaia has none, and so on), and every row lists its codes and sources.

For a coordinate fallback to work on a fast-moving star the proper motion must be known, so
put the star in the database or the supplementary table rather than relying on plan-file
coordinates alone.

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
- `{GAIA_ID}_{filter}_{date}_output.fits` (condensed by `condense_photometry.py`): per-star, per-frame fluxes in 13 apertures, positions, peaks, and `SMEAR_AP` (stars on a saturated smear trail; see [Shutterless frames](#shutterless-frames-readout-smear))

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
   - Matches the target's Gaia DR2 or DR3 source ID (from the output filename) to the catalogue
   - Retrieves effective temperature from Filippazzo catalogue or FITS header
   - Locates target in photometry catalogues

2. **Comparison Star Selection** (`stable_weights`, the default; `--alc legacy` restores the pre-2026-10 chain):
   - Candidates: every star that is not saturated, not low-flux, not on a saturated star's smear trail in any frame
     (shutterless frames; see [Shutterless frames](#shutterless-frames-readout-smear)) and no more than 3 mag fainter
     than the target
   - `preclean` drops stars that disagree with the rest of the field in more than 20% of frames
   - A smooth distance kernel, `1 / (1 + (sep / STABLE_KERNEL_ARCMIN)²)`, favours stars near the target
   - There is no magnitude-threshold search, no top-N cut and no reset, so a small change in the photometry gives a
     small change in the weights. The legacy chain could swap its whole ensemble on a 0.3% flux change.

3. **Iterative Differential Photometry:**
   - Creates weighted artificial comparison: `C(t) = Σ wᵢ × Fᵢ(t)`
   - Computes differential flux: `flux_diff = F_target / C(t)`
   - Weights are inverse variances, `wᵢ ∝ kernelᵢ / σᵢ²`. Each star's σᵢ is measured against the ALC of all the
     *other* stars (leave-one-out), so no star can lower its own σ by dominating the ALC.
   - The weights are iterated to convergence, usually in a few iterations
   - `WEIGHTS_<ap>`: W1 = 1/σ², W2 = distance kernel, W3 = eligibility (0/1), W_TOTAL = their normalised product

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
  - Quality flags and weights (`FLAGS`: `BW_FLAG` bad weather, read by the promotion check; `FLAG` the same bad-weather flag, read by the portal ingest; `SMEAR_FLAG` target aperture on a saturated smear trail. `WEIGHTS_{ap}`: `SATURATED`, `LOWFLUX`, `SMEARED`)

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

A night with no reduced targets - no science images at all, or every scheduled
target failing before T10 - still produces a report; the target overview page
carries a note in place of the per-target columns.

**Resource use:**
- `--nproc` sets the worker count for the image alignment behind the stack image.
  `ZLP_pipeline.sh` passes the pipeline's `CORES` budget; the default is 1.
- `--align-timeout` bounds one target's alignment (default 1800 s). On expiry the
  stack panel is dropped and the rest of the report is still produced.
- The module pins `OMP_NUM_THREADS` and friends to 1 unless the caller sets them,
  so the alignment workers do not each spawn a thread per host core.
- `ZLP_pipeline.sh` runs the whole stage under `timeout ${PDF_REPORT_TIMEOUT}`
  (default 3600 s) and treats a failure as a warning, so the report can never
  hold up T12.

---

### After T11: Frame Packing

**Implementation:** `pack_frames` in `ZLP_pipeline.sh`, calling `utils/framepack.py`

Calibrated frames (`<target>/1/proc*.fits`) are 85% of PipelineOutput. T6 writes
them as float32; nothing in a normal run reads them after T11, so before T12 moves
them to v2 they are rewritten according to `--frame-format`:

- `int16` (default): whole ADU, BITPIX 16, BSCALE 1 and a per-frame BZERO. No
  pixel moves by more than 0.5 ADU except saturated cores beyond the frame's
  65,534 ADU window, and non-finite pixels are stored as BLANK (read back as NaN).
  Readers that apply the scaling (astropy, fitsio, CASU) see float values as before.
- `rice`: int16, then lossless RICE_1 tile compression, renamed `proc*.fits.fz`.
  The pipeline itself never reads `.fz` frames: anything rerun on them must be
  decompressed first (`framepack.py decompress`).
- `float32`: leave the frames as T6 wrote them.

Every file is written to a hidden temporary file, read back and checked, then
renamed into place with the original's mode and modification time. A frame that
fails a check is left untouched, and a failure never stops the night. Each run
appends one row per frame to `logs/framepack_<DATE>_<RUNNAME>.csv`. The stage is
skipped when the run did not make or touch frames (`--only_T12`, for instance).

`framepack.py` also converts frames outside a pipeline run, for example:

```bash
python utils/framepack.py int16 <night>/<target>/1 --workers 8 --manifest framepack.csv
python utils/framepack.py decompress <night>/<target>/1
```

Never run it on a night the pipeline is working on.

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

### The 40 pc Target Table

The pipeline identifies a target by Gaia ID and by nothing else, so the target table is
what decides whether a star can be observed at all.

**Two files, one of them generated**

| File | Role |
|------|------|
| `ml_40pc.txt` | The sample as circulated. Source of truth for every astrophysical parameter. Never edited — nine analysis scripts outside this pipeline read it, including the pipeline and flare papers. |
| `ml_40pc_v2.csv` | What the pipeline reads. A generated artefact, rebuilt by `utils/build_target_table.py`; never hand-edited, so a new 40 pc release is a regeneration rather than a merge. |

`ml_40pc.txt` has a comma-separated header over whitespace-separated data rows. Readers use
`ascii.read(path, delimiter=' ', header_start=0, data_start=1)`, so the commas end up inside
the column names — the pipeline used to look up a column literally called `Gaia_ID,`. Only
`Program` escapes, because nothing follows it. Appending a column to that header silently
renames `Program` to `Program,`. The generated table simply uses one delimiter throughout.

**Identifiers**

The legacy `Gaia_ID` column holds a Gaia **DR2** source ID under an unqualified name. It
becomes two explicit columns, and either may be empty:

- `Gaia_DR2_ID` — empty for the four stars with no DR2 entry. Wolf 359 (Sp1056+0700) is the
  reason this matters: DR2's transit matching could not follow a star moving 4.7"/yr, so it
  has no DR2 row, and four years of SPECULOOS data were attributed to a G=16 background star.
- `Gaia_DR3_ID` — empty where DR3 has no counterpart, or where the DR2 source split into
  several DR3 sources. An identifier that cannot be pinned down is left out on purpose: the
  pipeline would otherwise adopt it at face value.

A DR2 source ID is **not** a DR3 source ID. Looking one up directly in `gaiadr3.gaia_source`
fails twice over — it returns nothing for a changed identifier, and where that number was
reused it returns a different star, silently. `gaiadr3.dr2_neighbourhood` is the table that
relates the two releases, and is what `utils/resolve_gaia_dr3_map.py` uses.

**Distances**

`parallax`, `DR3_dist_pc` and their errors come from Gaia DR3 via the resolved DR3
identifier. `Dis`/`e_Dis` are deliberately left alone: `M`, `R` and `T_eff` were derived from
them and have *not* been recomputed, so overwriting `Dis` would leave each row internally
inconsistent. Read `Dis` as "the distance the derived parameters assume" and `DR3_dist_pc` as
"the best distance we have".

**No row is ever dropped.** `low_conf_dist`, `dist_gt_40pc` and `flagged` (the criteria
Ben Rackham applied to the November 2024 revision, kept so the flags stay comparable) are
advisory columns, not a filter. The table's purpose is to find as many M and L dwarfs as
possible, so a star that has drifted beyond the 40 pc selection boundary stays in with the
flag set.

**Regenerating**

```bash
# offline and deterministic: a pure function of the three committed inputs
python3 -m utils.build_target_table /path/to/ml_40pc.txt -o ml_40pc_v2.csv

# only when the Gaia resolution itself needs rebuilding (this is the one
# part that touches the network)
python3 -m utils.resolve_gaia_dr3_map /path/to/ml_40pc.txt \
    --outdir calibration --report resolved.csv
```

The generator warns if the source row count changes, which is how a new sample release
announces itself.

**Which file the pipeline reads:** `TARGET_LIST` in the environment wins; otherwise
`ml_40pc_v2.csv` in the data root if present, falling back to `ml_40pc.txt`. Both formats are
read through `utils/target_list.py`, so either works.

---

### Local Gaia Database

The pipeline uses a local copy of Gaia DR3 cross-matched with 2MASS for fast astrometric calibration:

**Database Structure:**
- Sharded by declination (1° bins: -90_-89, -89_-88, ..., 89_90)
- Contains: RA, Dec (J2016), pmra, pmdec, parallax, G, G-RP, BP-RP, Teff, DR3 and DR2 IDs, J
- Size: ~41 GB
- Location: `$GAIADATABASEPATH`, in production `/gaia_database/gaia_dr3_unified_16jcut.db`
  (`orchard-server-spirit` still reads the older `gaia_tmass_16_jm_cut.db`)

**How it is built, and what it leaves out.** `gaia-tmass-sqlite/db_maker.py` keeps a DR3 source only
if its closest DR2 counterpart has 2MASS J < 16 through Gaia DR2's 2MASS cross-match, looked up in the
same 0.01° declination strip. That drops real stars DR2 never matched to 2MASS (mostly poor astrometry:
binaries) and some fast movers whose DR2 and DR3 positions fall in different strips: 470 targets of
`ml_40pc` on 2026-10-17. A star missing from the database cannot be identified as primary or secondary
on any night. Two hand-curated additions fill the gaps, each from a CSV in this repo:

- `custom_sources` table, from `calibration/supplementary_sources.csv`: objects in no Gaia release
  (synthetic IDs above 9e18); `python -m utils.supplementary_sources load`.
- rows in the band tables themselves, from `calibration/gaia_db_manual_rows.csv`: real Gaia stars
  the build left out, each with codes saying how every value was chosen;
  `python -m utils.gaia_db_manual_rows check` (read-only) and `apply`. `apply` records each row in the
  database's `manual_rows` table, so the database says which of its rows were made by hand and how
  (`python -m utils.gaia_db_manual_rows list`). It is idempotent: re-run it after any rebuild. It writes
  in place with a short lock, safe while jobs read the database; run it on appct inside orchard-server.

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
│   ├── wcs_cache.py          # T7's per-night solution cache
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
│   └── QC.py                 # Quality metrics tracking
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

Versions are git tags `vMAJOR.MINOR.PATCH` on `main`. Every FITS file the pipeline creates records the version
that made it (`ORCHVER`), and every later in-place edit records itself (`ORCHMOD` and a HISTORY card). See
[docs/versioning.md](docs/versioning.md) for when each part of the number goes up and what carries it.

**v4.0.0** (October 2026): the first tagged version
- Every FITS product, calibration frame and catalogue records the orchard version that made it, and later in-place
  edits record theirs
- Code releases: queue jobs run from frozen, versioned copies of the code
- From here on, "v2" and "v3" name only the output folders, which control what the portal ingests

**v3.0** (February 2025)
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