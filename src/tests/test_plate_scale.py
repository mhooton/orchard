"""
Plate-scale resolution tests.

The header values here are the real ones measured on 2026-09-27 from frames
of Callisto 20251122 (SPIRIT) and Artemis 20260810 (ANDOR). The SPIRIT case
is the regression that mattered: FOCALLEN is 8.0 *metres* but its comment
says '[mm]', and reading the unit off the comment made the plate scale 1000x
too large.

Run with pytest, or directly:  python src/tests/test_plate_scale.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from astropy.io import fits

from calibration.pipeutils import (PLATE_SCALE_SANE_RANGE, plate_scale_from_config,
                                   plate_scale_from_header, resolve_plate_scale)

# arcsec/pixel, from instrument_config.json
SPIRIT_CONFIG = {'arcsec_per_pixel': 0.31}
ANDOR_CONFIG = {'arcsec_per_pixel': 0.34}


def spirit_header():
    """Callisto 20251122, procSPECU4.2025-11-22T23:53:00.934.fits (1024x1280)."""
    h = fits.Header()
    h['XPIXSZ'] = (12.0, 'Pixel Width in microns (after binning)')
    h['FOCALLEN'] = (8.0, '[mm] Focal length of telescope')
    return h


def andor_header():
    """Artemis 20260810, procSp1924+7533-S001-R001-C001-r.fits (2046x2044)."""
    h = fits.Header()
    h['XPIXSZ'] = (13.5, 'Pixel Width in microns (after binning)')
    h['FOCALLEN'] = (8000.0, 'Focal length of telescope in mm')
    return h


def check(label, condition, detail=''):
    print(('PASS  ' if condition else 'FAIL  ') + label + (f'   {detail}' if detail else ''))
    return condition


def main():
    ok = True

    # --- SPIRIT: metres, mis-commented as [mm] ------------------------------
    ps = plate_scale_from_header(spirit_header())
    ok &= check('SPIRIT header plate scale ~0.309 arcsec/px',
                abs(ps - 0.3094) < 1e-3, f'got {ps:.4f}')
    ok &= check('SPIRIT plate scale agrees with config 0.31 to 1%',
                abs(ps - 0.31) / 0.31 < 0.01, f'got {ps:.4f}')
    ps_r = resolve_plate_scale(spirit_header(), SPIRIT_CONFIG)
    ok &= check('SPIRIT resolved plate scale is sane',
                PLATE_SCALE_SANE_RANGE[0] <= ps_r <= PLATE_SCALE_SANE_RANGE[1],
                f'got {ps_r:.4f}')

    # The old comment-sniffing behaviour, for contrast: 309 arcsec/px, which
    # inflated the Gaia query to ~110 deg and timed every frame out.
    old = 12.0e-6 / (8.0 * 1e-3)
    import math
    old_arcsec = math.degrees(math.atan(old)) * 3600.0
    ok &= check('the old comment-based reading really was 1000x out',
                abs(old_arcsec / ps - 1000.0) < 5.0, f'{old_arcsec:.1f} vs {ps:.4f}')

    # --- ANDOR: genuinely millimetres --------------------------------------
    ps_a = plate_scale_from_header(andor_header())
    ok &= check('ANDOR header plate scale ~0.348 arcsec/px',
                abs(ps_a - 0.3481) < 1e-3, f'got {ps_a:.4f}')
    ps_ar = resolve_plate_scale(andor_header(), ANDOR_CONFIG)
    ok &= check('ANDOR resolved plate scale matches header',
                abs(ps_ar - ps_a) < 1e-9, f'got {ps_ar:.4f}')

    # --- config fallbacks ---------------------------------------------------
    ok &= check('placeholder arcsec_per_pixel ("FILL") is rejected',
                plate_scale_from_config({'arcsec_per_pixel': 'FILL'}) is None)
    ok &= check('missing params is rejected',
                plate_scale_from_config(None) is None)
    ok &= check('no header keywords falls back to config',
                resolve_plate_scale(fits.Header(), ANDOR_CONFIG) == 0.34)
    ok &= check('no header and no config yields None',
                resolve_plate_scale(fits.Header(), None) is None)

    # --- an absurd header is overridden by the config -----------------------
    bad = fits.Header()
    bad['XPIXSZ'] = 13.5
    bad['FOCALLEN'] = 0.0001          # 0.1 mm: implausible, out of sane range
    ok &= check('implausible header plate scale falls back to config',
                resolve_plate_scale(bad, ANDOR_CONFIG) == 0.34,
                f'got {resolve_plate_scale(bad, ANDOR_CONFIG)}')

    # --- a header that disagrees wildly with config loses ------------------
    disagree = fits.Header()
    disagree['XPIXSZ'] = 13.5
    disagree['FOCALLEN'] = 1.0        # 1 m -> 2.8 arcsec/px, 8x the config
    got = resolve_plate_scale(disagree, ANDOR_CONFIG)
    ok &= check('header disagreeing with config by >4x falls back to config',
                got == 0.34, f'got {got}')

    # --- binning must still be honoured by the header path -----------------
    binned = fits.Header()
    binned['XPIXSZ'] = (27.0, 'Pixel Width in microns (after binning)')
    binned['FOCALLEN'] = (8000.0, 'Focal length of telescope in mm')
    got = resolve_plate_scale(binned, ANDOR_CONFIG)
    ok &= check('2x2 binned header (27 um) is trusted over the config',
                abs(got - 0.6963) < 1e-3, f'got {got:.4f}')

    print()
    print('ALL PASS' if ok else 'FAILURES PRESENT')
    return 0 if ok else 1


# pytest entry point
def test_plate_scale():
    assert main() == 0


if __name__ == '__main__':
    sys.exit(main())
