import subprocess as sp
import os
from threading import Lock
from calibration.pipeutils import open_fits_file, detect_instrument


def get_confmap_path(instrument_name):
    """
    Get the path to the confmap for a given instrument
    """
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up to parent directory and into photometry/confmaps
    confmap_dir = os.path.join(os.path.dirname(script_dir), 'photometry', 'confmaps')
    confmap_file = f"confmap_{instrument_name}.fits"
    confmap_path = os.path.join(confmap_dir, confmap_file)

    if not os.path.exists(confmap_path):
        raise FileNotFoundError(f"Confmap not found: {confmap_path}")

    return confmap_path


# In casutools.py, modify auto_detect_confmap() to add debug info:
def auto_detect_confmap(filelist_path):
    print(f"DEBUG: Auto-detecting confmap from filelist: {filelist_path}")
    try:
        with open(filelist_path, 'r') as f:
            first_file = f.readline().strip()

        print(f"DEBUG: First file in list: {first_file}")

        # Open the first file to detect instrument
        with open_fits_file(first_file) as hdul:
            instrument = detect_instrument(hdul)

        confmap_path = get_confmap_path(instrument)

        print(f"Auto-detected instrument: {instrument}")
        print(f"Using confmap: {confmap_path}")
        print(f"DEBUG: Confmap file exists: {os.path.exists(confmap_path)}")

        return confmap_path

    except Exception as e:
        print(f"Error auto-detecting confmap: {e}")
        raise


def find_imstack():
    '''
    Function to find imstack, as it has been renamed on ngtshead
    '''
    names = ['casu_imstack', 'imstack']
    for name in names:
        try:
            sp.check_call(['which', name], stderr=sp.PIPE, stdout=sp.PIPE)
        except sp.CalledProcessError:
            pass
        else:
            return name


def construct_filelist_argument(filelist):
    '''
    Wrapper around constructing a filelist
    '''
    return '@{0}'.format(filelist)


lock = Lock()


def run_command(cmd, verbose=False):
    '''
    Wraps subprocess to run the command
    '''
    str_cmd = list(map(str, cmd))
    if verbose == 'True':
        with lock:
            print(' '.join(str_cmd))

    # run command
    try:
        sp.check_call(str_cmd)  # , shell=True)
    except sp.CalledProcessError as casuexc:
        print("error code", casuexc.returncode, casuexc.output)
        with lock:
            print(' '.join(str_cmd))


def imstack(filelist, confidence_map=None,
            outstack='outstack.fits',
            outconf='outconf.fits',
            verbose=False,
            catalogues=""):
    '''
    Runs the casu task `imstack`

    Args:
        filelist: Path to filelist or the filelist content
        confidence_map: Path to confidence map. If None, auto-detects from first file in filelist
        outstack: Output stack filename
        outconf: Output confidence filename
        verbose: Verbose output
        catalogues: Catalogue parameter for imstack
    '''
    verbose = True

    # Auto-detect confmap if not provided
    if confidence_map is None:
        print("No confmap specified, auto-detecting...")
        confidence_map = auto_detect_confmap(filelist)
    elif confidence_map == 'auto':
        # Explicit auto-detection request
        confidence_map = auto_detect_confmap(filelist)

    cmd = ['imstack', construct_filelist_argument(filelist),
           confidence_map, catalogues, outstack, outconf]

    run_command(cmd, verbose=verbose)


def imcore(input_file, output_table,
           ipix=6,
           threshold=2.0,
           confidence_map='noconf',
           rcore=4,
           filtfwhm=1,
           ellfile=False,
           casu_verbose=False,
           verbose=False):
    '''
    Runs the casu task `imcore`
    '''
    cmd = ['imcore', input_file, confidence_map, output_table, ipix, threshold,
           '--filtfwhm', filtfwhm, '--rcore', rcore]

    if casu_verbose:
        cmd.append('--verbose')

    if not ellfile:
        cmd.append('--noell')

    run_command(cmd, verbose=verbose)


def imcore_list(input_file, listfile, output_file,
                threshold=2.0,
                confidence_map='noconf',
                rcore=4,
                cattype=6,
                casu_verbose=False,
                noell=True,
                verbose=False):
    '''
    Runs the casu task `imcore_list`
    '''
    cmd = ['imcore_list', input_file, confidence_map, listfile, output_file,
           threshold, '--rcore', rcore, '--cattype', cattype]

    if noell:
        cmd.append('--noell')

    if casu_verbose:
        cmd.append('--verbose')

    run_command(cmd, verbose=verbose)


def wcsfit(infile, incat, catsrc, verbose=True, site='cds'):
    '''
    Runs the casu task `wcsfit`. Local fits reference is required
    '''
    cmd = ['wcsfit', infile, incat, '--site', site, '--catsrc', catsrc]
    run_command(cmd, verbose=verbose)
