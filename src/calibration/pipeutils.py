from contextlib import contextmanager
import bz2
from astropy.io import fits

def detect_instrument(hdul):
    """
    Centralized instrument detection from FITS HDUList

    Args:
        hdul: Opened FITS HDUList

    Returns:
        str: Instrument name ('andor', 'spirit', 'moana')
    """
    try:
        hdr = hdul[0].header

        if 'SPECULOOS' in hdr['TELESCOP'] or hdr['TELESCOP'] == 'Artemis':
            if hdr['NAXIS2'] == 1280:
                return 'spirit'  # Fixed: spirit is 1280
            else:
                return 'andor'  # Fixed: andor is the larger format
        elif 'MOANA' in hdr['TELESCOP']:
            return 'moana'
        else:
            # Fallback - changed to andor as requested
            print(f"Warning: Unknown telescope, defaulting to 'andor'")
            return 'andor'

    except Exception as e:
        print(f"Warning: Could not detect instrument from header: {e}")
        return 'andor'  # Default fallback changed to andor


def extract_overscan(hdul):
    """
    Extract overscan region based on instrument type

    Args:
        hdul: FITS HDU list
    """
    instrument = detect_instrument(hdul)

    # Instrument-specific overscan extraction
    if instrument == 'spirit':
        return 0.
    elif instrument == 'andor':
        return hdul[0].data[-15:, 4:]
    elif instrument == 'moana':
        return 0.
    else:
        return 0.


def image_trim(hdul):
    """
    Trim image based on instrument type

    Args:
        hdul: FITS HDU list
    """
    instrument = detect_instrument(hdul)

    # Instrument-specific image trimming
    if instrument == 'spirit':
        return hdul[0].data
    elif instrument == 'andor':
        return hdul[0].data[22:2066, 2:2048]
    elif instrument == 'moana':
        return hdul[0].data[1:-1, 1:-1]
    else:
        return hdul[0].data


@contextmanager
def open_fits_file(filename):
    if filename.endswith('.bz2'):
        with bz2.BZ2File(filename) as uncompressed:
            with fits.open(uncompressed) as infile:
                yield infile
    else:
        with fits.open(filename) as infile:
            yield infile

# class NullPool(object):
#     """A dummy pool class that executes map() sequentially instead of in parallel."""
#
#     def __init__(self, *args, **kwargs):
#         pass
#
#     def map(self, fn, args):
#         return list(map(fn, args))  # Consistent parameter name + list conversion
#
#     def close(self):
#         pass
#
#     def join(self):
#         pass
