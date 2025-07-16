#!/usr/bin/env python3
"""
Test script to verify proper motion correction behavior
"""

import fitsio
import numpy as np
import datetime
import astropy.time
from astropy.coordinates import SkyCoord
import astropy.units as u


def test_proper_motion_correction(original_cat, pm_corrected_cat, obs_date_string):
    """
    Test whether proper motion correction uses observation date or 2015 epoch

    Args:
        original_cat: Path to original stack catalogue
        pm_corrected_cat: Path to proper motion corrected catalogue
        obs_date_string: Observation date in YYYYMMDD format
    """

    print("=== PROPER MOTION CORRECTION TEST ===\n")

    # Load both catalogues
    try:
        with fitsio.FITS(original_cat) as f:
            orig_data = f[1].read()
            # Try different extension names for Gaia data
            gaia_data = None
            for ext_name in ['Gaia_Crossmatch', 'GAIA_CROSSMATCH', 'gaia']:
                try:
                    gaia_data = f[ext_name].read()
                    print(f"Found Gaia data in extension: {ext_name}")
                    break
                except:
                    continue

            if gaia_data is None:
                print("ERROR: No Gaia crossmatch data found in original catalogue")
                return

        with fitsio.FITS(pm_corrected_cat) as f:
            pm_data = f[1].read()

    except Exception as e:
        print(f"ERROR loading catalogues: {e}")
        return

    # Convert observation date
    obs_dt = datetime.datetime.strptime(obs_date_string, "%Y%m%d")
    obs_jd = astropy.time.Time(obs_dt).jd

    # Gaia DR2 reference epoch (2015.5)
    gaia_epoch_jd = astropy.time.Time("2015-07-02").jd

    print(f"Observation date: {obs_date_string} (JD: {obs_jd:.2f})")
    print(f"Gaia DR2 epoch: 2015-07-02 (JD: {gaia_epoch_jd:.2f})")
    print(f"Time difference: {(obs_jd - gaia_epoch_jd) / 365.25:.3f} years\n")

    # Find stars with significant proper motion (>50 mas/yr)
    # Try both uppercase and lowercase column names
    pmra_col = gaia_data['PMRA'] if 'PMRA' in gaia_data.dtype.names else gaia_data['pmra']
    pmdec_col = gaia_data['PMDEC'] if 'PMDEC' in gaia_data.dtype.names else gaia_data['pmdec']

    high_pm_mask = np.sqrt(pmra_col ** 2 + pmdec_col ** 2) > 50
    high_pm_indices = np.where(high_pm_mask & ~np.isnan(pmra_col))[0]

    if len(high_pm_indices) == 0:
        print("WARNING: No stars with high proper motion found. Test may be inconclusive.")
        # Use any star with PM data
        high_pm_indices = np.where(~np.isnan(pmra_col))[0][:5]

    print(f"Testing {len(high_pm_indices)} stars with proper motion data:\n")

    # Test each high PM star
    for i, idx in enumerate(high_pm_indices[:5]):  # Test first 5

        print(f"--- Star {i + 1} (Index {idx}) ---")

        # Original coordinates (try both RA/ra and DEC/dec)
        orig_ra = orig_data['RA'][idx] if 'RA' in orig_data.dtype.names else orig_data['ra'][idx]
        orig_dec = orig_data['DEC'][idx] if 'DEC' in orig_data.dtype.names else orig_data['dec'][idx]

        # PM-corrected coordinates from pipeline
        corrected_ra = pm_data['RA'][idx] if 'RA' in pm_data.dtype.names else pm_data['ra'][idx]
        corrected_dec = pm_data['DEC'][idx] if 'DEC' in pm_data.dtype.names else pm_data['dec'][idx]

        # Gaia proper motion data
        pmra = pmra_col[idx]  # mas/yr
        pmdec = pmdec_col[idx]  # mas/yr

        print(f"Original RA/Dec: {orig_ra:.6f}, {orig_dec:.6f}")
        print(f"PM RA/Dec: {pmra:.2f}, {pmdec:.2f} mas/yr")
        print(f"Pipeline corrected RA/Dec: {corrected_ra:.6f}, {corrected_dec:.6f}")

        # Calculate expected positions for both scenarios

        # Scenario 1: Using observation date (CORRECT)
        time_years_obs = (obs_jd - gaia_epoch_jd) / 365.25
        expected_ra_obs = orig_ra + (pmra / 3600000.0) * time_years_obs / np.cos(np.radians(orig_dec))
        expected_dec_obs = orig_dec + (pmdec / 3600000.0) * time_years_obs

        # Scenario 2: Using 2015 epoch (BUG - no change)
        expected_ra_bug = orig_ra  # No change
        expected_dec_bug = orig_dec  # No change

        print(f"Expected if obs date used: {expected_ra_obs:.6f}, {expected_dec_obs:.6f}")
        print(f"Expected if 2015 used: {expected_ra_bug:.6f}, {expected_dec_bug:.6f}")

        # Calculate differences
        diff_obs_ra = abs(corrected_ra - expected_ra_obs) * 3600 * 1000  # mas
        diff_obs_dec = abs(corrected_dec - expected_dec_obs) * 3600 * 1000  # mas

        diff_bug_ra = abs(corrected_ra - expected_ra_bug) * 3600 * 1000  # mas
        diff_bug_dec = abs(corrected_dec - expected_dec_bug) * 3600 * 1000  # mas

        print(f"Difference from obs-date prediction: {diff_obs_ra:.1f}, {diff_obs_dec:.1f} mas")
        print(f"Difference from 2015-epoch prediction: {diff_bug_ra:.1f}, {diff_bug_dec:.1f} mas")

        # Determine which scenario matches better
        total_diff_obs = np.sqrt(diff_obs_ra ** 2 + diff_obs_dec ** 2)
        total_diff_bug = np.sqrt(diff_bug_ra ** 2 + diff_bug_dec ** 2)

        if total_diff_obs < total_diff_bug:
            verdict = "✓ OBSERVATION DATE"
        else:
            verdict = "✗ 2015 EPOCH (BUG)"

        print(f"VERDICT: {verdict}")
        print()

    # Summary analysis
    print("=== SUMMARY ===")
    orig_ra_col = orig_data['RA'] if 'RA' in orig_data.dtype.names else orig_data['ra']
    orig_dec_col = orig_data['DEC'] if 'DEC' in orig_data.dtype.names else orig_data['dec']
    corr_ra_col = pm_data['RA'] if 'RA' in pm_data.dtype.names else pm_data['ra']
    corr_dec_col = pm_data['DEC'] if 'DEC' in pm_data.dtype.names else pm_data['dec']

    total_pm_change = np.sqrt(np.sum((corr_ra_col - orig_ra_col) ** 2 +
                                     (corr_dec_col - orig_dec_col) ** 2))

    if total_pm_change < 1e-8:  # Essentially no change
        print("❌ LIKELY BUG: No significant coordinate changes detected")
        print("   Pipeline appears to be using 2015 epoch instead of observation date")
    else:
        print("✅ LIKELY CORRECT: Significant coordinate changes detected")
        print("   Pipeline appears to be using observation date")

    print(f"\nTotal coordinate change magnitude: {total_pm_change:.2e} degrees")


def quick_debug_test():
    """
    Quick test by adding debug prints to the proper motion code
    """
    print("=== QUICK DEBUG TEST ===\n")
    print("Add this debug line to calculate_pm() function in proper_motion.py:")
    print()
    print("def calculate_pm(stars, cat_jd, new_jd):")
    print("    print(f'DEBUG: cat_jd={cat_jd}, new_jd={new_jd}, diff={(cat_jd-new_jd)/365:.3f} years')")
    print("    for s in range(len(stars['ra'])):")
    print("        ...")
    print()
    print("Then run your pipeline and check the output:")
    print("- If diff ≈ 0 years → using 2015 epoch (BUG)")
    print("- If diff = years since 2015 → using observation date (CORRECT)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test proper motion correction")
    parser.add_argument("--debug-only", action="store_true",
                        help="Just show debug instructions")
    parser.add_argument("original_cat", nargs="?", help="Path to original stack catalogue")
    parser.add_argument("pm_cat", nargs="?", help="Path to PM-corrected catalogue")
    parser.add_argument("obs_date", nargs="?", help="Observation date (YYYYMMDD)")

    args = parser.parse_args()

    if args.debug_only:
        quick_debug_test()
    elif args.original_cat and args.pm_cat and args.obs_date:
        test_proper_motion_correction(args.original_cat, args.pm_cat, args.obs_date)
    else:
        parser.error("Either use --debug-only or provide all three arguments: original_cat pm_cat obs_date")