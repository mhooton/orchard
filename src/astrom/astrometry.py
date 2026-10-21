# This module makes a proper WCS-header of a FITS file.

# One BLAS thread per process, set before numpy is first imported. T7's parallelism is its process pool, and BLAS
# thread pools started inside the forked workers deadlock: on appct (2026-10-20/21) workers whose OpenBLAS had started
# its 64 threads sat in futex_wait for hours, out of reach of the 60 s SIGALRM, until the queue's job timeout. numpy
# and scipy each bring their own OpenBLAS; both read these.
import os
for _var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_var] = '1'

from astropy.io import fits
from functools import partial
import argparse
import timeit
import astrom.find_pa
from calibration.pipeutils import detect_instrument, get_instrument_parameters, open_fits_file, proc_name
import fnmatch
import shutil
from astrom.twirl_speculoos import twirl_wcs
from astrom.pointer_wcs import pointer_wcs
from astrom import frame_pool, wcs_cache
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)

# The plate solver's generation, recorded in every WCS cache row (astrom/wcs_cache.py, docs/wcs-cache.md). Bump it
# in any change that alters the solutions T7 writes (pointer_wcs's detection, catalogue query, matching or
# acceptance, the twirl call), so cached solutions from the old solver are no longer used.
# 2: the second pass from the nearest solved frame (neighbour_jobs, known_scale.solve_from_neighbour).
SOLVER_EPOCH = 2

# The second pass. A frame that failed with at least this many sources is tried again from the solution of the frame
# nearest to it in time that has one, solved this run or taken from the WCS cache, if that frame is no further away
# than NEIGHBOUR_MAX_GAP. Trailed and timed-out frames are not retried.
NEIGHBOUR_MIN_SOURCES = 5
NEIGHBOUR_MAX_GAP = 7200.0      # seconds

# A pool worker that has had a frame this long past the per-frame timeout is stuck where SIGALRM cannot reach it (in C,
# like the workers spinning in OpenBLAS's sched_yield on 2026-10-21) and is killed (astrom/frame_pool.py). Generous:
# a long C call, a Gaia query or a slow disk, delays the alarm's exception until it returns.
HANG_GRACE = 240                # seconds


def add_astrometry_with_timeout(infile, timeout, ext, db_path, total_files, file_num=None, neighbour=None):
    """
    Wrapper that runs add_astrometry with a signal-based timeout.

    Parameters
    ----------
    infile : str
        Path to image file
    timeout : int
        Timeout in seconds
    ext : str
        File extension
    db_path : str
        Path to Gaia database
    total_files : int
        Total number of files being processed
    file_num : int, optional
        File number for progress display
    neighbour : dict, optional
        WCS cards of the nearest frame with a solution, for the second pass

    Returns
    -------
    dict
        Result dictionary from add_astrometry or timeout error
    """
    import signal

    class TimeoutException(Exception):
        pass

    def timeout_handler(signum, frame):
        raise TimeoutException(f"Plate solving timed out after {timeout}s")

    # `progress` is filled in place by pointer_wcs as it works, so when the
    # alarm fires we can say how far the solve got and which stage was
    # running. Reporting bare zeros instead used to read as "the Gaia query
    # returned nothing", which is a different fault entirely.
    progress = {}

    # Set up the timeout
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)

    try:
        result = add_astrometry(infile, ext, db_path,
                                file_num=file_num, total_files=total_files,
                                progress=progress, neighbour=neighbour)
        signal.alarm(0)  # Cancel the alarm
        return result
    except TimeoutException as e:
        signal.alarm(0)  # Cancel the alarm
        print(f"[TIMEOUT] {infile}: {str(e)}")
        return _partial_result(progress, str(e))
    except Exception as e:
        signal.alarm(0)  # Cancel the alarm
        print(f"[ERROR] {infile}: {str(e)}")
        return _partial_result(progress, str(e))

# The stages pointer_wcs moves through, in order. Used to work out which
# counters a killed solve had actually established.
_STAGE_ORDER = ['starting', 'opening_file', 'clean_image', 'plate_scale',
                'star_detection', 'gaia_query', 'star_matching', 'complete']

# A counter is only meaningful once the stage that fills it has finished, so
# it is trustworthy only if the solve died in a *later* stage than this.
_COUNTER_FILLED_IN = {
    'sources_detected': 'star_detection',
    'sources_used': 'star_detection',
    'gaia_queried': 'gaia_query',
    'gaia_used': 'gaia_query',
    'matches': 'star_matching',
}


def _partial_result(progress, error):
    """
    Build a failure result from however far pointer_wcs got.

    pointer_wcs seeds every counter at zero, so a counter's value alone
    cannot tell us whether its stage ran. We therefore keep only the counters
    whose stage completed before the one we died in, and report the rest as
    None so the log prints '?'. Reporting them all as 0 is what made a solve
    killed during the Gaia query look like a Gaia query that found nothing.
    """
    stage = progress.get('stage', 'unknown')
    try:
        reached = _STAGE_ORDER.index(stage)
    except ValueError:
        reached = -1

    result = {
        'success': False,
        'error': f"{error} during {stage}",
        'crpix': None,
        'crval': None,
        'stage': stage,
    }
    for counter, filled_in in _COUNTER_FILLED_IN.items():
        if reached > _STAGE_ORDER.index(filled_in):
            result[counter] = progress.get(counter)
        else:
            result[counter] = None
    return result


def _read_list(path):
    """The non-blank lines of a file list."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def raw_frames_by_proc(raw_list):
    """{calibrated frame name: raw frame path} from T1's image list, by T6's naming rule (pipeutils.proc_name).
    Read only: T7 never writes a raw frame."""
    if not raw_list or not os.path.exists(raw_list):
        return {}
    by_proc = {}
    for raw in _read_list(raw_list):
        by_proc.setdefault(proc_name(raw), raw)
    return by_proc


def apply_cached_solutions(infiles, cache_path, raw_by_proc, run, force):
    """Apply the night's cached solutions to the frames that have a usable one.

    Returns (frames to solve, number applied, {frame: (raw name, validity fields)}). The fields are taken before
    anything is solved and are what a new solve's cache row is keyed on. A failure of any kind costs only the
    cache: the frame is solved instead."""
    code = wcs_cache.code_identity()
    entries = {}
    if not force:
        try:
            entries = wcs_cache.load(cache_path)
        except Exception as e:
            print(f"WARNING: WCS cache not used: {e}", flush=True)
    print(f"WCS cache: {cache_path} ({len(entries)} rows" +
          (", ignored: --force-platesolve re-solves every frame and refreshes its row)" if force else ")"),
          flush=True)

    to_solve, frames, misses, applied = [], {}, {}, 0
    for f in infiles:
        try:
            raw = raw_by_proc.get(os.path.basename(f))
            if raw is None:
                why = 'no raw frame'
            else:
                raw_name = os.path.basename(raw)
                fields = dict(run, **wcs_cache.frame_fields(raw, fits.getheader(f)))
                frames[f] = (raw_name, fields)
                why = 'forced' if force else wcs_cache.check(entries.get(raw_name), fields, code)
                if why is None:
                    wcs_cache.apply(f, entries[raw_name][wcs_cache.CARDS], WCS_KEYWORDS)
                    applied += 1
                    continue
        except Exception as e:
            print(f"WARNING: WCS cache not used for {os.path.basename(f)}: {e}", flush=True)
            why = 'unusable'
        misses[why] = misses.get(why, 0) + 1
        to_solve.append(f)

    tally = ', '.join(f"{why} {n}" for why, n in sorted(misses.items(), key=lambda kv: -kv[1]))
    print(f"WCS cache: {applied} frames from the cache, {len(to_solve)} to solve" +
          (f" ({tally})" if tally else ""), flush=True)
    return to_solve, applied, frames


def store_new_solutions(cache_path, solved, frames, results):
    """Write a cache row for every frame this run solved, from its header as now on disk: one read-change-write
    of the night's table. Failed and timed-out frames get no row. Never raises."""
    new = {}
    for f, result in zip(solved, results):
        if not (result and result.get('success', False)) or f not in frames:
            continue
        raw_name, fields = frames[f]
        try:
            cards = wcs_cache.cards_from_header(fits.getheader(f), WCS_KEYWORDS)
            new[raw_name] = wcs_cache.make_entry(raw_name, fields, cards, matches=result.get('matches'),
                                                 proc=os.path.abspath(f))
        except Exception as e:
            print(f"WARNING: no WCS cache row for {os.path.basename(f)}: {e}", flush=True)
    if not new:
        return
    try:
        written = wcs_cache.store(cache_path, new)
        print(f"WCS cache: {written} rows written to {cache_path}", flush=True)
    except Exception as e:
        print(f"WARNING: WCS cache not updated ({cache_path}): {e}", flush=True)


def _obs_seconds(header):
    """A frame's DATE-OBS in seconds (any fixed origin), or None."""
    from datetime import datetime
    value = str(header.get('DATE-OBS', '')).strip()
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return (datetime.strptime(value, fmt) - datetime(2000, 1, 1)).total_seconds()
        except ValueError:
            pass
    return None


def worth_a_second_pass(result):
    """A failed frame with enough round sources: not trailed, not cut short by the timeout."""
    if not result or result.get('success', False):
        return False
    error = result.get('error') or ''
    if error.startswith('Sources trailed') or 'timed out' in error:
        return False
    return (result.get('sources_detected') or 0) >= NEIGHBOUR_MIN_SOURCES


def neighbour_jobs(to_solve, results, infiles, read_header=fits.getheader):
    """[(index into to_solve, WCS cards of its neighbour)] for the second pass: each frame worth one, with the frame
    nearest to it in time that has a solution, solved this run or taken from the WCS cache (every listed frame not
    in to_solve), and no more than NEIGHBOUR_MAX_GAP away."""
    retry = [i for i, r in enumerate(results) if worth_a_second_pass(r)]
    if not retry:
        return []
    solving = set(to_solve)
    solved = {f for f, r in zip(to_solve, results) if r and r.get('success', False)}
    headers, times = {}, {}
    for f in [f for f in infiles if f in solved or f not in solving] + [to_solve[i] for i in retry]:
        try:
            headers[f] = read_header(f)
            times[f] = _obs_seconds(headers[f])
        except Exception:
            times[f] = None
    have = [f for f in infiles if (f in solved or f not in solving) and times.get(f) is not None
            and headers[f].get('CTYPE1')]
    jobs = []
    for i in retry:
        t = times.get(to_solve[i])
        if t is None or not have:
            continue
        nearest = min(have, key=lambda f: abs(times[f] - t))
        if abs(times[nearest] - t) > NEIGHBOUR_MAX_GAP:
            continue
        jobs.append((i, {k: headers[nearest][k] for k in WCS_KEYWORDS if k in headers[nearest]}))
    return jobs


def _solve_job(job, timeout, ext, db_path, total_files):
    """One frame in a pool worker. job = (frame, None) in the first pass, (frame, the neighbour's cards) in the
    second."""
    infile, neighbour = job
    return add_astrometry_with_timeout(infile, timeout=timeout, ext=ext, db_path=db_path, total_files=total_files,
                                       neighbour=neighbour)


def _lost_frame(job, why):
    """The result of a frame the pool lost with its worker (or whose solve raised): failed, with the reason, logged
    as any failed frame is."""
    result = {'success': False, 'error': why, 'crpix': None, 'crval': None, 'stage': 'unknown'}
    result.update((counter, None) for counter in _COUNTER_FILLED_IN)
    print(_summary_line('', os.path.basename(job[0]), result), flush=True)
    return result


def _solve_in_pool(work, args, total_files):
    """_solve_job on every job of `work`, in args.nproc worker processes; the results in order (astrom/frame_pool.py:
    a dead or hung worker costs only its frame)."""
    worker = partial(_solve_job, timeout=args.timeout, ext=args.ext, db_path=args.db_path, total_files=total_files)
    return frame_pool.run(worker, work, int(args.nproc), hang_after=args.timeout + HANG_GRACE, lost=_lost_frame,
                          label=lambda job: os.path.basename(job[0]))


def main(args):
    import sys
    sys.stdout.flush()
    print("*************************", flush=True)
    print("Starting plate solving", flush=True)
    print("*************************", flush=True)
    start_time = timeit.default_timer()
    print(f"Input file: {args.filelist}", flush=True)
    print(f"Extension: {args.ext}", flush=True)
    print(f"Force plate solve: {args.force_platesolve}", flush=True)

    newly_solved = 0
    failed_solve = 0

    # A WCS already in a frame's header (Astra's, PinPoint's, ESO's, or orchard's copied into the raw frame by
    # older code) never counts as solved: every frame is solved, or takes its own solution from the WCS cache.
    try:
        infiles = _read_list(args.filelist)
    except Exception as e:
        print("Importing Astrometry failed: " + str(e), flush=True)
        infiles = []
    total_files = len(infiles)
    print(f"Total files: {total_files}", flush=True)
    if not infiles:
        print("No frames to plate solve", flush=True)
        return

    # Get parameters once from the first file in the list
    with open_fits_file(infiles[0]) as hdul:
        params = get_instrument_parameters(hdul, trimmed=True)  # <-- Add trimmed=True

    # Extract trim offsets once
    trim_offsets = (0, 0)  # default
    if 'trim' in params:
        trim_config = params['trim']
        trim_offsets = (trim_config.get('left_col', 0), trim_config.get('top_row', 0))

    # Solutions from earlier runs on this night (astrom/wcs_cache.py)
    to_solve, from_cache, frames = infiles, 0, {}
    use_cache = False
    if args.wcs_cache:
        raw_by_proc = raw_frames_by_proc(args.raw_list)
        try:
            run = wcs_cache.run_fields(args.db_path, trim_offsets, SOLVER_EPOCH)
        except Exception as e:
            print(f"WARNING: WCS cache not used: {e}", flush=True)
        else:
            if raw_by_proc:
                use_cache = True
                cache_start = timeit.default_timer()
                to_solve, from_cache, frames = apply_cached_solutions(infiles, args.wcs_cache, raw_by_proc, run,
                                                                      args.force_platesolve)
                print(f"WCS cache: applied in {timeit.default_timer() - cache_start:.1f} s", flush=True)
            else:
                print(f"WARNING: WCS cache not used: no raw frames listed in {args.raw_list}", flush=True)
    else:
        print("WCS cache: none given", flush=True)

    files_to_process = len(to_solve)
    print(f"Files to process: {files_to_process}", flush=True)

    if not to_solve:
        results = []
    elif int(args.nproc) > 1:
        # Parallel processing with timeout
        print(f"Running in parallel mode with {args.nproc} processes", flush=True)
        results = _solve_in_pool([(f, None) for f in to_solve], args, files_to_process)
    else:
        # Sequential processing with timeout
        print(f"Running in serial mode (nproc={args.nproc})", flush=True)
        results = []
        for i, infile in enumerate(to_solve, 1):
            results.append(add_astrometry_with_timeout(infile,
                                                       timeout=args.timeout,
                                                       ext=args.ext,
                                                       db_path=args.db_path,
                                                       total_files=files_to_process,
                                                       file_num=i))

    # The second pass: failed frames tried again from the nearest frame with a solution (neighbour_jobs)
    jobs = neighbour_jobs(to_solve, results, infiles)
    if jobs:
        print(f"Second pass: {len(jobs)} frames tried again from the nearest frame with a solution", flush=True)
        work = [(to_solve[i], cards) for i, cards in jobs]
        if int(args.nproc) > 1:
            again = _solve_in_pool(work, args, files_to_process)
        else:
            again = [_solve_job(w, timeout=args.timeout, ext=args.ext, db_path=args.db_path,
                                total_files=files_to_process) for w in work]
        for (i, _), result in zip(jobs, again):
            results[i] = result
        print(f"Second pass: {sum(1 for r in again if r and r.get('success', False))} of {len(jobs)} solved",
              flush=True)

    # Count successes and failures
    for result in results:
        if result and result.get('success', False):
            newly_solved += 1
        else:
            failed_solve += 1

    if use_cache:
        store_new_solutions(args.wcs_cache, to_solve, frames, results)

    try:
        for f in to_solve:
            new_f = f.replace(args.ext, 'new')
            f1 = f.replace(args.ext, 'axy')
            f2 = f.replace(args.ext, 'wcs')
            f3 = f.replace(args.ext, 'solved')
            f4 = f.replace(args.ext, 'corr')
            f5 = f.replace("." + args.ext, '-indx.xyls')
            f6 = f.replace(args.ext, 'match')
            f7 = f.replace(args.ext, 'rdls')
            if os.path.isfile(new_f):

                hdulist_old = fits.open(f)
                hdr_old = hdulist_old[0].header
                data = hdulist_old[0].data
                hdulist_new = fits.open(new_f)
                hdr_new = hdulist_new[0].header
                pa, jd, solve = find_pa.pa(new_f)

                hdr_old.set('PA', pa)
                # hdr_old.set('CTYPE1',hdr_new['CTYPE1']) #'RA---TAN')
                # hdr_old.set('CTYPE2', hdr_new['CTYPE2']) #'DEC--TAN')
                # hdr_old.set('CRVAL1',hdr_new['CRVAL1'])
                # hdr_old.set('CRVAL2',hdr_new['CRVAL2'])
                # hdr_old.set('CRPIX1',hdr_new['CRPIX1'])
                # hdr_old.set('CRPIX2',hdr_new['CRPIX2'])
                # hdr_old.set('CD1_1',hdr_new['CD1_1'])
                # hdr_old.set('CD2_1',hdr_new['CD2_1'])
                # hdr_old.set('CD1_2',hdr_new['CD1_2'])
                # hdr_old.set('CD2_2',hdr_new['CD2_2'])
                for k in WCS_KEYWORDS:
                    if hdr_old.get(k) == None:
                        hdr_old.set(k, hdr_new[k])
                # if hdr_old.get('RA') == None:
                #     hdr_old.set('RA', '21 31 48.113')
                #     hdr_old.set('DEC', '48 26 29.330')
                hdr_old.set('ASTSOLVE', 'T')
                fits.writeto(f, data, hdr_old, overwrite=True)

                for fname in [new_f, f1, f2, f3, f4, f5, f6, f7]:
                    os.remove(fname)

        elapsed = timeit.default_timer() - start_time
        print('Total time taken for astrom: ' + str(elapsed / 60.) + ' minutes')

    except Exception as e:
        print("Replacing Astrometry failed: " + str(e))
        print("on file: " + f)

    # Print summary
    print("\n" + "=" * 50)
    print("\n" + "=" * 50)
    print("PLATE SOLVING SUMMARY")
    print("=" * 50)
    print(f"Total images in list: {total_files}")
    print(f"Images with a WCS from the cache: {from_cache}")
    print(f"Images processed: {files_to_process}")
    print(f"Successfully plate solved: {newly_solved}")
    print(f"Failed plate solving: {failed_solve}")
    if files_to_process > 0:
        print(f"Success rate: {newly_solved / files_to_process * 100:.1f}%")
    print("=" * 50)


def get_trim_offsets(instrument):
    """
    Get pixel offsets for trimming based on instrument type.

    Parameters
    ----------
    instrument : str
        Instrument name ('spirit', 'andor', 'moana')

    Returns
    -------
    tuple
        (x_offset, y_offset) - pixels removed from left and top edges
    """
    if instrument == 'spirit':
        return (0, 0)  # No trimming
    elif instrument == 'andor':
        return (2, 22)  # Left=2, Top=22 based on [22:2066, 2:2048]
    elif instrument == 'moana':
        return (1, 1)  # Based on [1:-1, 1:-1]
    else:
        return (0, 0)  # Unknown instrument, no adjustment


# Define WCS keywords once at module level
WCS_KEYWORDS = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CUNIT1', 'CUNIT2',
                'WCSAXES', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1',
                'CDELT2', 'LONPOLE', 'LATPOLE', 'RADESYS']


def _summary_line(prefix, filename, result):
    """A frame's line in the T7 log (parsed by utils/promote_policy.py's SUMMARY)."""
    # A counter of None means that stage never ran; print '?' so it is not
    # mistaken for a stage that ran and came back empty.
    def _n(key):
        value = result.get(key)
        return '?' if value is None else value

    if result['success']:
        crpix1, crpix2 = result['crpix']
        crval1, crval2 = result['crval']
        return (f"{prefix} {filename} SUCCESS - "
                f"Sources: {_n('sources_detected')}→{_n('sources_used')} "
                f"Gaia: {_n('gaia_queried')}→{_n('gaia_used')} "
                f"Matches: {_n('matches')} "
                f"Method: {result.get('method') or '?'} "
                f"CRPIX=({crpix1:.2f},{crpix2:.2f}) "
                f"CRVAL=({crval1:.6f},{crval2:.6f})")
    return (f"{prefix} {filename} FAILED - "
            f"Sources: {_n('sources_detected')}→{_n('sources_used')} "
            f"Gaia: {_n('gaia_queried')}→{_n('gaia_used')} "
            f"Matches: {_n('matches')} "
            f"Error: {result['error']}")


def add_astrometry(f, ext, db_path, file_num=None, total_files=None, progress=None, neighbour=None):
    if fnmatch.fnmatch(f, '*.' + ext):
        file = f
        # print("Add astrometry.net data to header of file: " + file)

        # Attempt WCS solving
        # result = twirl_wcs(str(file), verbose=False)
        result = pointer_wcs(str(file), db_path, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True,
                             verbose=False, progress=progress, neighbour=neighbour)

        # Create summary line
        filename = os.path.basename(file)
        if file_num is not None and total_files is not None:
            prefix = f"[{file_num}/{total_files}]"
        else:
            prefix = ""
        print(_summary_line(prefix, filename, result))
        return result

    return None  # Add this for non-matching files

        #
        # hdulist = fits.open(file)
        # prihdr = hdulist[0].header
        #
        # # Extract RA and DEC
        # ra_val = prihdr.get('RA')
        # dec_val = prihdr.get('DEC')
        #
        # if ra_val is not None and dec_val is not None:
        #     # Process RA
        #     if " " in str(ra_val):
        #         strra = str(ra_val.replace(" ", ":"))
        #     else:
        #         strra = str(ra_val)
        #
        #     # Process DEC
        #     if " " in str(dec_val):
        #         strdec = str(dec_val.replace(" ", ":"))
        #     else:
        #         strdec = str(dec_val)
        #
        #     if detect_instrument(hdulist) == 'spirit':
        #         print("Running solve-field with spirit astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/spirit_astrometry.cfg " \
        #                      f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                      f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                      f"--scale-low 0.32 --scale-high 0.38 --scale-units arcsecperpix {file}"
        #     elif detect_instrument(hdulist) == 'andor':
        #         print("Running solve-field with andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"
        #     else:
        #         print("Instrument not spirit or andor, defaulting to andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"

        # print(f"Running: {cmd}")
        # os.system(cmd)
        # else:
        #     print(f"Warning: Missing RA and/or DEC in {file}, skipping astrom")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filelist')
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help='Number of parallel processes (default: 1)')
    parser.add_argument('-e', '--ext')
    parser.add_argument('-d', '--db_path')
    parser.add_argument('-r', '--raw-list',
                        help="T1's raw image list: matches each frame to its raw frame for the WCS cache (read only)")
    parser.add_argument('-w', '--wcs-cache',
                        help="this night's WCS cache table (ZLP: <BASEDIR>/PipelineOutput/wcs_cache/<TEL>/<night>.ecsv)")
    parser.add_argument('-f', '--force-platesolve', action='store_true',
                        help='Ignore the WCS cache: solve every frame and refresh its cache row')
    parser.add_argument('-t', '--timeout', type=int, default=60,
                        help='Timeout per image in seconds (default: 60)')
    main(parser.parse_args())