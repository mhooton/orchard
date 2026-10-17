# Target-nights that can never give a light curve

[`unrecoverable_nights.csv`](unrecoverable_nights.csv) lists target-nights whose raw data cannot give a light curve
of that target, whatever we change in the pipeline. Completeness counts (the ledger's gap, the "observed but not on
the portal" hours) should leave them out, and nobody should spend a rerun on them.

## What qualifies

Only a reason in the data itself:

- **TARGET_NOT_IN_FIELD**: the telescope pointed somewhere else, so the target is in none of the frames.
- **NO_ASTROMETRY**: no frame can be plate solved reliably (too few catalogue stars in the field), so no aperture can
  be placed on the target.

Add a new code when a new kind of reason turns up, and describe it here.

## What does not qualify

Anything a rerun, a download or a fix could recover: pipeline failures, a target missing from the local Gaia
database, a wrong ID in a target list or plan file, frames missing from disk but held at ESO, calibration frames in the
wrong orientation, a night that was never processed. Those belong to the completeness ledger's work packages.

## Adding an entry

One row per telescope, night and target, with the evidence that rules out every route to a light curve, and the date it
was checked. `frames` and `exposure_h` count the target's light frames on disk. `night` is the pipeline's night
directory (the evening date).
