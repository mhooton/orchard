from contextlib import contextmanager
import bz2
from astropy.io import fits
import numpy as np
from astropy.convolution import interpolate_replace_nans
import json
import os


def load_instrument_config():
    """Load instrument configuration from JSON file"""
    config_path = os.path.join(os.path.dirname(__file__), 'instrument_config.json')
    with open(config_path, 'r') as f:
        return json.load(f)

def detect_instrument(hdul, trimmed=False):
    """
    Centralized instrument detection from FITS HDUList using configuration file

    Args:
        hdul: Opened FITS HDUList
        trimmed: Boolean indicating if the image has been trimmed (default: False)

    Returns:
        str: Instrument name ('andor', 'spirit', etc.)
    """
    try:
        hdr = hdul[0].header
        telescop = hdr.get('TELESCOP', '').strip()
        naxis2 = hdr.get('NAXIS2', 0)

        # Load configuration
        config = load_instrument_config()

        # First, try to match telescope name exactly
        if telescop in config['telescopes']:
            telescope_config = config['telescopes'][telescop]

            # Loop through instruments for this telescope and match by naxis2
            for instrument_name, instrument_params in telescope_config.items():
                config_naxis2 = instrument_params['naxis2']

                # If trimmed=True, calculate expected trimmed naxis2
                if trimmed and 'trim' in instrument_params:
                    trim_config = instrument_params['trim']
                    top_row = trim_config.get('top_row', 0) or 0
                    bottom_row = trim_config.get('bottom_row', config_naxis2)

                    # Handle negative indexing
                    if bottom_row and bottom_row < 0:
                        bottom_row = config_naxis2 + bottom_row + 1
                    elif bottom_row is None:
                        bottom_row = config_naxis2

                    expected_naxis2 = bottom_row - top_row
                else:
                    expected_naxis2 = config_naxis2

                if expected_naxis2 == naxis2:
                    return instrument_name

            # If no naxis2 match, return the first (and likely only) instrument
            instrument_names = list(telescope_config.keys())
            if instrument_names:
                print(f"Warning: No NAXIS2 match for {telescop}, using {instrument_names[0]}")
                return instrument_names[0]

        # Fallback: search all telescopes for naxis2 match
        for telescope_name, telescope_config in config['telescopes'].items():
            for instrument_name, instrument_params in telescope_config.items():
                config_naxis2 = instrument_params['naxis2']

                if trimmed and 'trim' in instrument_params:
                    trim_config = instrument_params['trim']
                    top_row = trim_config.get('top_row', 0) or 0
                    bottom_row = trim_config.get('bottom_row', config_naxis2)
                    if bottom_row and bottom_row < 0:
                        bottom_row = config_naxis2 + bottom_row + 1
                    elif bottom_row is None:
                        bottom_row = config_naxis2
                    expected_naxis2 = bottom_row - top_row
                else:
                    expected_naxis2 = config_naxis2

                if expected_naxis2 == naxis2:
                    print(f"Warning: Telescope name '{telescop}' not recognized, but matched instrument by NAXIS2")
                    return instrument_name

        # Final fallback
        print(f"Warning: Unknown telescope '{telescop}' with NAXIS2={naxis2}, defaulting to 'andor'")
        return 'andor'

    except Exception as e:
        print(f"Warning: Could not detect instrument from header: {e}")
        return 'andor'

def get_instrument_parameters(hdul, trimmed=False):
    """
    Get instrument parameters from config file given a FITS HDUList

    Args:
        hdul: Opened FITS HDUList
        trimmed: Boolean indicating if the image has been trimmed

    Returns:
        dict: Dictionary containing instrument parameters
    """
    try:
        hdr = hdul[0].header
        telescop = hdr.get('TELESCOP', '')

        # Load configuration and detect instrument
        config = load_instrument_config()
        instrument = detect_instrument(hdul, trimmed=trimmed)

        # Get telescope config
        if telescop in config['telescopes']:
            telescope_config = config['telescopes'][telescop]

            # Get instrument parameters
            if instrument in telescope_config:
                return telescope_config[instrument]
            else:
                raise ValueError(f"Instrument '{instrument}' not found for telescope '{telescop}'")
        else:
            raise ValueError(f"Telescope '{telescop}' not found in configuration")

    except Exception as e:
        print(f"Warning: Could not load instrument parameters: {e}")
        # Return default fallback parameters
        return {
            'gain': 1.0029,
            'saturation_threshold': 62000,
            'naxis2': 2088,
            'bad_pixel_correction': False
        }

def get_instrument_parameter(hdul, parameter_name):
    """
    Get a specific instrument parameter from config file

    Args:
        hdul: Opened FITS HDUList
        parameter_name: Name of the parameter to retrieve

    Returns:
        The parameter value, or None if not found
    """
    params = get_instrument_parameters(hdul)
    return params.get(parameter_name)


def get_acquisition_software(hdul):
    """Determine acquisition software from observer header"""
    observer = hdul[0].header.get('observer', '')

    if observer == 'Astra':
        return 'astra'
    elif observer == 'Gavin Boyle [NDA #1]':
        return 'moana'
    else:
        return 'acp'  # Default for anything else

def get_jd_header_name(hdul):
    """Get the correct JD header name based on acquisition software"""
    config = load_instrument_config()
    acquisition = get_acquisition_software(hdul)
    return config['acquisition'][acquisition]['jd_header_name']

def apply_image_slice(hdul, config_key, fallback_value):
    """
    Apply image slicing based on instrument configuration

    Args:
        hdul: FITS HDU list
        config_key: Configuration key to look for ('overscan', 'trim', etc.)
        fallback_value: Value to return if no configuration found

    Returns:
        Sliced image data or fallback value if no config
    """
    try:
        # Get instrument parameters
        params = get_instrument_parameters(hdul)

        # Check if slice configuration exists for this instrument
        if config_key not in params:
            return fallback_value

        slice_config = params[config_key]
        data = hdul[0].data

        # Extract slice region using configured values
        top_row = slice_config.get('top_row')
        bottom_row = slice_config.get('bottom_row')
        left_col = slice_config.get('left_col')
        right_col = slice_config.get('right_col')

        # Handle None values (equivalent to no slicing on that axis)
        row_slice = slice(top_row, bottom_row)
        col_slice = slice(left_col, right_col)

        return data[row_slice, col_slice]

    except Exception as e:
        print(f"Warning: Could not apply {config_key} slice: {e}")
        return fallback_value


def extract_overscan(hdul):
    """
    Extract overscan region based on instrument configuration

    Args:
        hdul: FITS HDU list

    Returns:
        Overscan data array or 0 if no overscan configured
    """
    return apply_image_slice(hdul, 'overscan', 0.)


def image_trim(hdul):
    """
    Trim image based on instrument configuration

    Args:
        hdul: FITS HDU list

    Returns:
        Trimmed image data or original data if no trim configured
    """
    return apply_image_slice(hdul, 'trim', hdul[0].data)


def interpolate_pixels(frame):
    """
    Interpolate bad pixels in the frame using a kernel

    Args:
        frame: 2D numpy array representing the image frame
    """
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])
    result_frame = interpolate_replace_nans(frame, kernel)

    return result_frame

def clean_bad_pixels(image, bad_pixel_map):
    """
    Clean bad pixels in image data

    Args:
        image: 2D image array to clean
        dark: Dark frame data
        flat: Flat frame data

    Returns:
        Cleaned image array
    """
    # try:
    image = image.copy()

    # bad_pixel_map = open_fits_file(bpm_name)
    image[bad_pixel_map] = np.nan
    image = interpolate_pixels(image)
    return image
    # except Exception as e:
    #     print(f"Warning: Could not clean bad pixels: {e}")
    #     return image

def create_dark_dict(darknames, bp_correct=False):
    dark_dict = {}
    bpm_name_dict = {}
    for d in darknames:
        try:
            # Parse exposure time from filename
            basename = os.path.basename(d)
            if '_MasterDark.fits' in basename and basename.endswith('_MasterDark.fits'):
                key = 'combined'
            elif '_MasterDark_' in basename and 's.fits' in basename:
                key = int(basename.split('_MasterDark_')[1].replace('s.fits', ''))
            else:
                print(f"Warning: Could not parse dark filename: {basename}")

            # Combined dark file
            with open_fits_file(d) as hdulist:
                dark_dict[key] = hdulist[0].data
            print(f"Loaded master dark: {basename}")

            if bp_correct:
                bpm_name = f"{d.split('MasterDark')[0]}BadPixelMap{d.split('MasterDark')[1]}"
                bpm_name_dict[key] = bpm_name
                print(f"Generated bad pixel map {bpm_name}")

        except Exception as e:
            print(f"Warning: Could not load dark file {d}: {e}")
            continue

    print(list(dark_dict.keys()))
    if bp_correct:
        return dark_dict, bpm_name_dict
    else:
        return dark_dict

def create_flat_dict(flatnames):
    flat_dict = {}
    for f in flatnames:
        filt = f.split('.fits')[0].split('_')[-1]
        with open_fits_file(f) as hdulist:
            flat_dict[filt] = hdulist[0].data

    print(list(flat_dict.keys()))
    return flat_dict


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
