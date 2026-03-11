import fitsio

__all__ = ['wcs_succeeded', 'set_wcs_status']

HEADER_KEY = 'wcscompl'


def wcs_succeeded(fname, hdu=0):
    '''
    Return true if the wcs has succeeded for a file
    '''
    from astropy.io import fits
    from calibration.pipeutils import get_instrument_parameters
    import math

    with fits.open(fname) as hdul:
        header = hdul[hdu].header

        # Check the explicit success flag first
        if HEADER_KEY in header and not header[HEADER_KEY]:
            return False

        # Check that CD matrix keywords are present and non-zero
        cd_vals = {}
        for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
            val = header.get(key, None)
            if val is None or val == 0.0:
                return False
            cd_vals[key] = val

        # Check plate scale against known physical value
        params = get_instrument_parameters(hdul, trimmed=True)
        expected_arcsec_per_pixel = params['arcsec_per_pixel']
        expected_deg_per_pixel = expected_arcsec_per_pixel / 3600.0

        determinant = cd_vals['CD1_1'] * cd_vals['CD2_2'] - cd_vals['CD1_2'] * cd_vals['CD2_1']
        measured_deg_per_pixel = math.sqrt(abs(determinant))

        tolerance = 0.1  # 10%
        if abs(measured_deg_per_pixel - expected_deg_per_pixel) / expected_deg_per_pixel > tolerance:
            return False

    return True


def set_wcs_status(fname, succeeded, hdu=0):
    '''
    Set the corresponding header key to `succeeded`
    '''
    with fitsio.FITS(fname, 'rw') as outfile:
        outfile[hdu].write_key(HEADER_KEY, succeeded, comment='WCS succeeded?')
