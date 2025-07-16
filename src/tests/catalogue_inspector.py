#!/usr/bin/env python3
"""
Quick script to inspect catalogue structure and diagnose PM issues
"""

import fitsio
import sys


def inspect_catalogue(catfile):
    """Inspect the structure of a FITS catalogue"""

    print(f"=== INSPECTING: {catfile} ===\n")

    try:
        with fitsio.FITS(catfile) as f:
            print(f"Number of extensions: {len(f)}")

            # List all extensions
            for i in range(len(f)):
                try:
                    extname = f[i].get_extname() if hasattr(f[i], 'get_extname') else f"Extension_{i}"
                    exttype = type(f[i]).__name__
                    print(f"  [{i}] {extname} ({exttype})")
                except:
                    print(f"  [{i}] Unknown extension")
            print()

            # Check main catalogue (usually extension 1)
            if len(f) > 1:
                print("Main catalogue columns:")
                try:
                    cols = f[1].get_colnames()
                    for col in cols:
                        print(f"  {col}")
                    print()
                except Exception as e:
                    print(f"  Error reading columns: {e}")

                # Check for Gaia extension by trying different names
                gaia_ext = None
                for i in range(len(f)):
                    try:
                        extname = f[i].get_extname() if hasattr(f[i], 'get_extname') else ""
                        if 'gaia' in extname.lower():
                            gaia_ext = i
                            break
                    except:
                        pass

                if gaia_ext is not None:
                    print(f"\nGaia extension found at index {gaia_ext}!")
                    try:
                        gaia_cols = f[gaia_ext].get_colnames()
                        print("Gaia columns:")
                        for col in gaia_cols:
                            print(f"  {col}")
                    except Exception as e:
                        print(f"Error reading Gaia columns: {e}")
                else:
                    print("\n❌ NO Gaia extension found")

                    # Check if Gaia data is in main table
                    try:
                        cols = f[1].get_colnames()
                        gaia_like_cols = [col for col in cols if
                                          'gaia' in col.lower() or 'pm' in col.lower() or 'parallax' in col.lower()]
                        if gaia_like_cols:
                            print(f"Possible Gaia columns in main table: {gaia_like_cols}")
                        else:
                            print("No obvious Gaia/PM columns found in main table")
                    except:
                        print("Could not check main table columns")

            else:
                print("❌ No data extensions found")

    except Exception as e:
        print(f"❌ ERROR reading catalogue: {e}")


def check_pm_file_exists(pm_catfile):
    """Check if PM-corrected file exists and why it might not"""

    print(f"\n=== CHECKING PM FILE: {pm_catfile} ===")

    import os
    if os.path.exists(pm_catfile):
        print("✅ PM-corrected catalogue exists")
        inspect_catalogue(pm_catfile)
    else:
        print("❌ PM-corrected catalogue does NOT exist")

        # Check directory for similar files
        dirname = os.path.dirname(pm_catfile)
        basename = os.path.basename(pm_catfile)

        print(f"\nFiles in directory {dirname}:")
        try:
            files = os.listdir(dirname)
            stack_files = [f for f in files if 'stack' in f]
            for f in sorted(stack_files):
                print(f"  {f}")
        except:
            print("  Could not list directory")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python catalogue_inspector.py <catalogue_file> [pm_catalogue_file]")
        sys.exit(1)

    catfile = sys.argv[1]
    inspect_catalogue(catfile)

    if len(sys.argv) > 2:
        pm_catfile = sys.argv[2]
        check_pm_file_exists(pm_catfile)