#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPECULOOS Observation Time Calculator (Python 3.9+)

Calculates total observation time for each target across the telescope network.
Handles timestamp-based duration calculation with fallback to exposure time method.
"""

import os
import re
import csv
import argparse
import sys
import numpy as np
from datetime import datetime
from collections import defaultdict
from pathlib import Path

try:
    from astropy.io import fits

    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False
    print("Warning: astropy not available. Fallback method will be disabled.")

try:
    import matplotlib

    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update({'font.size': 16})
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available. Histogram plots will be disabled.")

try:
    from astropy.table import Table
    from astropy.io import ascii

    ASTROPY_TABLE_AVAILABLE = True
except ImportError:
    ASTROPY_TABLE_AVAILABLE = False
    print("Warning: astropy.table not available. ByTelescope files will use simple format.")


class ObservationCalculator:
    def __init__(self):
        self.telescopes = ['Artemis', 'Callisto', 'Europa', 'Ganymede', 'Io']
        self.base_path = '/data/SPECULOOSPipeline'
        self.observations_path = f'{self.base_path}/Observations'
        # self.pipeline_path = f'{self.base_path}/PipelineOutput/v2'

        # Logging
        self.log_file = None

        # Data storage
        self.target_data = defaultdict(lambda: {
            'total_hours': 0.0,
            'sessions': [],
            'total_images': 0,
            'telescopes': set(),
            'dates': []
        })

        # Statistics
        self.stats = {
            'total_targets': 0,
            'total_sessions': 0,
            'timestamp_successes': 0,
            'exposure_fallbacks': 0,
            'failures': 0,
            'errors': [],
            'failed_calculations': []
        }

    def setup_logging(self, csv_filename):
        """Setup logging to file with same base name as CSV output."""
        log_filename = csv_filename.replace('.csv', '.log')
        self.log_file = open(log_filename, 'w')
        self.log_info(f"Log file created: {log_filename}")

    def close_logging(self):
        """Close the log file."""
        if self.log_file:
            self.log_file.close()
            self.log_file = None

    def _write_log(self, level, message):
        """Write message to both console and log file."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"[{timestamp}] {level} - {message}"

        # Write to console
        print(f"[{level}] {message}")

        # Write to log file if available
        if self.log_file:
            self.log_file.write(log_message + '\n')
            self.log_file.flush()

    def log_info(self, message):
        """Simple logging replacement for Python 3"""
        self._write_log("INFO", message)

    def log_warning(self, message):
        """Simple logging replacement for Python 3"""
        self._write_log("WARNING", message)

    def log_error(self, message):
        """Simple logging replacement for Python 3"""
        self._write_log("ERROR", message)

    def log_debug(self, message):
        """Simple logging replacement for Python 3"""
        # Uncomment next line for debug output
        # self._write_log("DEBUG", message)
        pass

    def normalize_filepath(self, original_path, telescope, date):
        """
        Normalize file paths to current directory structure.
        Extract filename and reconstruct path using current structure.
        """
        # Extract just the filename from the original path
        filename = os.path.basename(original_path)

        # Reconstruct the correct current path
        normalized_path = f"{self.observations_path}/{telescope}/images/{date}/{filename}"

        return normalized_path

    def validate_target_name(self, target_name, first_filepath, telescope, date):
        """Validate target name against coordinates and log discrepancies."""
        if not ASTROPY_AVAILABLE:
            return

        # Get coordinates from first FITS file
        normalized_path = self.normalize_filepath(first_filepath, telescope, date)

        # Check for compressed version first, then uncompressed
        if os.path.exists(normalized_path + '.fz'):
            file_to_open = normalized_path + '.fz'
        elif os.path.exists(normalized_path):
            file_to_open = normalized_path
        else:
            return  # Can't validate without file

        try:
            with fits.open(file_to_open) as hdul:
                # For .fz files, data is usually in extension 1, for .fits in extension 0
                header_ext = 1 if file_to_open.endswith('.fz') else 0

                # Check if the required extension exists
                if len(hdul) <= header_ext:
                    self.log_warning(f"FITS file {file_to_open} missing expected extension {header_ext}")
                    return

                header = hdul[header_ext].header

                # Get coordinates with error checking
                if 'RA' in header and 'DEC' in header:
                    try:
                        ra = float(header['RA'])
                        dec = float(header['DEC'])
                    except (ValueError, TypeError) as e:
                        self.log_debug(f"Invalid coordinate values in {file_to_open}: {e}")
                        return

                    # Generate expected coordinate-based name
                    expected_name = self.generate_coordinate_name(ra, dec)

                    # Handle special targets
                    if "TRAPPIST" in target_name:
                        # Don't validate TRAPPIST targets against coordinates
                        return

                    # Only validate SPECULOOS coordinate-based targets
                    if (len(target_name) >= 6 and
                            target_name.startswith("SP") and
                            target_name[2:6].isdigit()):

                        # Check if coordinate-based name matches (case-insensitive)
                        if target_name.upper() != expected_name.upper():
                            # Add to wrong names tracking
                            if not hasattr(self, 'wrong_names'):
                                self.wrong_names = []
                                self.corrected_names = []

                            if target_name.ljust(30) not in self.wrong_names:
                                self.wrong_names.append(target_name.ljust(30))
                                self.corrected_names.append(expected_name)
                                self.log_debug(f"Name mismatch: {target_name} -> {expected_name}")

                    # Don't validate non-SPECULOOS targets (WASP, TESS, etc.)

                else:
                    self.log_debug(f"Missing RA/DEC headers in {file_to_open}")

        except Exception as e:
            self.log_debug(f"Could not validate target name for {target_name} (file: {file_to_open}): {e}")

    def generate_coordinate_name(self, ra, dec):
        """Generate SPECULOOS target name from RA/Dec coordinates."""
        try:
            from astropy.coordinates import SkyCoord
            from astropy import units as u

            if " " not in str(ra):
                c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
                ra_hms = c.ra.hms
                ra_str = str(int(ra_hms[0])).zfill(2) + str(int(ra_hms[1])).zfill(2)

                dec_dms = c.dec.dms
                if int(dec_dms[0]) > 0:
                    sign = '+'
                    deg = str(int(dec_dms[0])).zfill(2)
                    arcmin = str(int(dec_dms[1])).zfill(2)
                else:
                    sign = "-"
                    deg = str(int(dec_dms[0]))[1:].zfill(2)
                    arcmin = str(int(dec_dms[1]))[1:].zfill(2)
                dec_str = sign + deg + arcmin

                return "SP" + ra_str + dec_str
            else:
                ra_str = str(ra).replace(" ", "")
                dec_str = str(dec).replace(" ", "")
                return "SP" + ra_str[:4] + dec_str[:5]
        except Exception:
            return "UNKNOWN"

    def parse_timestamp_from_filename(self, filepath):
        """Extract timestamp from FITS filename."""
        filename = os.path.basename(filepath)

        # Pattern for SPECU1.2022-04-06T01:43:15.640.fits or similar
        # Also handle .fits.fz files and underscore format (T00_48_18 instead of T00:48:18)
        patterns = [
            # Standard colon format
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3})\.fits(?:\.fz)?$',  # .sss format
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{2})\.fits(?:\.fz)?$',  # .ss format
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})\.fits(?:\.fz)?$',  # no milliseconds
            # Underscore format (T00_48_18 instead of T00:48:18)
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}\.\d{3})\.fits(?:\.fz)?$',  # .sss format with underscores
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}\.\d{2})\.fits(?:\.fz)?$',  # .ss format with underscores
            r'\.(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2})\.fits(?:\.fz)?$'  # no milliseconds with underscores
        ]

        for pattern in patterns:
            match = re.search(pattern, filename)
            if match:
                timestamp_str = match.group(1)
                try:
                    # Convert underscores back to colons for datetime parsing
                    if '_' in timestamp_str:
                        timestamp_str = timestamp_str.replace('_', ':')

                    # Handle different timestamp formats
                    if '.' in timestamp_str:
                        # With milliseconds - Python 3.9+ handles this directly
                        # But need to pad/truncate to 6 digits for microseconds
                        base_time, fractional = timestamp_str.split('.')
                        # Pad or truncate to 6 digits (microseconds)
                        if len(fractional) == 3:  # milliseconds
                            fractional = fractional + '000'
                        elif len(fractional) == 2:  # centiseconds
                            fractional = fractional + '0000'
                        elif len(fractional) < 6:
                            fractional = fractional.ljust(6, '0')
                        else:
                            fractional = fractional[:6]

                        timestamp_str = f"{base_time}.{fractional}"
                        return datetime.strptime(timestamp_str, '%Y-%m-%dT%H:%M:%S.%f')
                    else:
                        # Without milliseconds
                        return datetime.strptime(timestamp_str, '%Y-%m-%dT%H:%M:%S')
                except ValueError:
                    continue

        return None

    def get_exposure_time(self, filepath, telescope, date):
        """Get exposure time from FITS header using normalized path."""
        if not ASTROPY_AVAILABLE:
            return None

        # Normalize the filepath to current directory structure
        normalized_path = self.normalize_filepath(filepath, telescope, date)

        # Check for compressed version first, then uncompressed
        if os.path.exists(normalized_path + '.fz'):
            file_to_open = normalized_path + '.fz'
            self.log_debug(f"Using compressed file: {file_to_open}")
        elif os.path.exists(normalized_path):
            file_to_open = normalized_path
        else:
            self.log_warning(f"Neither compressed nor uncompressed file exists: {normalized_path}")
            return None

        try:
            with fits.open(file_to_open) as hdul:
                # For .fz files, data is usually in extension 1, for .fits in extension 0
                header_ext = 1 if file_to_open.endswith('.fz') else 0

                # Check if the required extension exists
                if len(hdul) <= header_ext:
                    self.log_warning(f"FITS file {file_to_open} missing expected extension {header_ext}")
                    return None

                header = hdul[header_ext].header

                # Try different header keywords
                for keyword in ['EXPTIME', 'EXPOSURE']:
                    if keyword in header:
                        try:
                            return float(header[keyword])
                        except (ValueError, TypeError) as e:
                            self.log_warning(f"Invalid {keyword} value in {file_to_open}: {e}")
                            continue
        except Exception as e:
            self.log_warning(f"Could not read exposure time from {file_to_open}: {e}")

        return None

    def calculate_duration_from_timestamps(self, filepaths):
        """Calculate observation duration from timestamps with gap detection."""
        timestamps = []
        failed_parses = 0

        for filepath in filepaths:
            timestamp = self.parse_timestamp_from_filename(filepath)
            if timestamp:
                timestamps.append(timestamp)
            else:
                failed_parses += 1

        # Log if no timestamps could be parsed
        if len(timestamps) == 0:
            self.log_warning(f"No filenames matched timestamp format for {len(filepaths)} files")
            return None
        elif failed_parses > 0:
            self.log_warning(f"Could not parse timestamps from {failed_parses} of {len(filepaths)} files")

        if len(timestamps) < 2:
            return None

        timestamps.sort()

        # Calculate time differences between consecutive images
        time_diffs = []
        for i in range(1, len(timestamps)):
            diff_seconds = (timestamps[i] - timestamps[i - 1]).total_seconds()
            time_diffs.append(diff_seconds)

        if len(time_diffs) == 0:
            return None

        # Calculate median cadence (normal time between images)
        time_diffs_sorted = sorted(time_diffs)
        median_cadence = time_diffs_sorted[len(time_diffs_sorted) // 2]

        # Identify gaps (differences > 3x median cadence)
        gap_threshold = median_cadence * 3.0
        total_gap_time = 0.0
        gap_count = 0

        for diff in time_diffs:
            if diff > gap_threshold:
                gap_time = diff - median_cadence  # Subtract normal cadence, count excess as gap
                total_gap_time += gap_time
                gap_count += 1

        # Total duration minus detected gaps
        total_duration = (timestamps[-1] - timestamps[0]).total_seconds()
        observation_duration = total_duration - total_gap_time

        # Log gap detection results
        if gap_count > 0:
            self.log_info(
                f"Detected {gap_count} gaps totalling {total_gap_time / 60.0:.1f}min (median cadence: {median_cadence:.1f}s)")

        return observation_duration / 3600.0  # Convert to hours

    def calculate_duration_from_exposure(self, filepaths, telescope, date):
        """Calculate duration using cadence from first 10 files, scaled to total count."""
        if not filepaths:
            return None

        if not ASTROPY_AVAILABLE:
            self.log_warning("Astropy not available for exposure time method")
            return None

        # Sample first 10 files (or all if fewer than 10)
        sample_size = min(10, len(filepaths))
        sample_files = filepaths[:sample_size]

        jd_times = []
        exposure_time = None
        failed_reads = 0

        for filepath in sample_files:
            normalized_path = self.normalize_filepath(filepath, telescope, date)

            # Check for compressed version first, then uncompressed
            if os.path.exists(normalized_path + '.fz'):
                file_to_open = normalized_path + '.fz'
            elif os.path.exists(normalized_path):
                file_to_open = normalized_path
            else:
                failed_reads += 1
                self.log_debug(f"File not found: {normalized_path}")
                continue

            try:
                with fits.open(file_to_open) as hdul:
                    # For .fz files, data is usually in extension 1, for .fits in extension 0
                    header_ext = 1 if file_to_open.endswith('.fz') else 0

                    # Check if the required extension exists
                    if len(hdul) <= header_ext:
                        self.log_warning(f"FITS file {file_to_open} missing expected extension {header_ext}")
                        failed_reads += 1
                        continue

                    header = hdul[header_ext].header

                    # Get JD time (try multiple header keywords)
                    jd_time = None
                    for jd_keyword in ['JD', 'JD-OBS', 'MJD-OBS']:
                        if jd_keyword in header:
                            try:
                                jd_val = float(header[jd_keyword])
                                # Convert MJD to JD if needed
                                if jd_keyword.startswith('MJD'):
                                    jd_val += 2400000.5
                                jd_time = jd_val
                                break
                            except (ValueError, TypeError) as e:
                                self.log_debug(f"Invalid {jd_keyword} value in {file_to_open}: {e}")
                                continue

                    if jd_time:
                        jd_times.append(jd_time)
                    else:
                        failed_reads += 1
                        self.log_debug(f"No valid JD time found in {file_to_open}")

                    # Get exposure time from first file
                    if exposure_time is None:
                        for exp_keyword in ['EXPTIME', 'EXPOSURE']:
                            if exp_keyword in header:
                                try:
                                    exposure_time = float(header[exp_keyword])
                                    break
                                except (ValueError, TypeError) as e:
                                    self.log_debug(f"Invalid {exp_keyword} value in {file_to_open}: {e}")
                                    continue

            except Exception as e:
                failed_reads += 1
                self.log_debug(f"Could not read FITS file {file_to_open}: {e}")
                continue

        # Log if no headers could be read
        if len(jd_times) == 0:
            self.log_warning(f"No JD headers found in {sample_size} sampled FITS files")
            if exposure_time is not None:
                # Fallback to simple exposure time calculation with estimated readout
                estimated_readout = 3.0  # seconds
                cycle_time = exposure_time + estimated_readout
                total_seconds = len(filepaths) * cycle_time
                self.log_info(
                    f"Using fallback: exposure + estimated readout ({exposure_time:.1f}s + {estimated_readout:.1f}s) x {len(filepaths)} images")
                return total_seconds / 3600.0
            return None
        elif failed_reads > 0:
            self.log_warning(f"Could not read JD headers from {failed_reads} of {sample_size} sampled files")

        if len(jd_times) < 2:
            return None

        # Calculate median cadence from sample
        jd_times.sort()
        time_diffs = []
        for i in range(1, len(jd_times)):
            diff_seconds = (jd_times[i] - jd_times[i - 1]) * 86400.0  # JD difference to seconds
            time_diffs.append(diff_seconds)

        time_diffs_sorted = sorted(time_diffs)
        median_cadence = time_diffs_sorted[len(time_diffs_sorted) // 2]

        # Scale to total number of images
        total_duration = median_cadence * (len(filepaths) - 1)  # N-1 intervals for N images

        self.log_info(
            f"Sampled {len(jd_times)} files: median cadence {median_cadence:.1f}s, scaled to {len(filepaths)} images")

        return total_duration / 3600.0  # Convert to hours

    def process_target_list(self, list_filepath, telescope, date):
        """Process a single target list file."""
        if not os.path.exists(list_filepath):
            return None

        try:
            with open(list_filepath, 'r') as f:
                filepaths = [line.strip() for line in f if line.strip()]

            if not filepaths:
                return None

            # Extract target name from list filename
            list_filename = os.path.basename(list_filepath)
            target_match = re.match(r'1_image_(.+)\.list$', list_filename)
            if not target_match:
                self.log_warning(f"Could not extract target name from {list_filename}")
                return None

            target_name = target_match.group(1).upper()

            # Validate target name against coordinates (for wrong_target_names.txt)
            self.validate_target_name(target_name, filepaths[0], telescope, date)

            # Try timestamp-based calculation first
            duration = self.calculate_duration_from_timestamps(filepaths)
            method = 'timestamp'

            if duration is None:
                # Fallback to exposure time method
                duration = self.calculate_duration_from_exposure(filepaths, telescope, date)
                method = 'exposure'

                if duration is None:
                    failure_info = {
                        'target': target_name,
                        'telescope': telescope,
                        'date': date,
                        'num_files': len(filepaths)
                    }
                    self.stats['failed_calculations'].append(failure_info)
                    self.log_warning(f"Could not calculate duration for {target_name} on {telescope} {date}")
                    self.stats['failures'] += 1
                    return None

            # Record the observation
            session_data = {
                'telescope': telescope,
                'date': date,
                'duration_hours': duration,
                'num_images': len(filepaths),
                'method': method
            }

            self.target_data[target_name]['total_hours'] += duration
            self.target_data[target_name]['sessions'].append(session_data)
            self.target_data[target_name]['total_images'] += len(filepaths)
            self.target_data[target_name]['telescopes'].add(telescope)
            self.target_data[target_name]['dates'].append(date)

            # Update statistics
            if method == 'timestamp':
                self.stats['timestamp_successes'] += 1
            else:
                self.stats['exposure_fallbacks'] += 1

            self.stats['total_sessions'] += 1

            self.log_info(f"Processed {target_name}: {duration:.3f}h using {method} method ({len(filepaths)} images)")
            return session_data

        except Exception as e:
            error_msg = f"Error processing {list_filepath}: {e}"
            self.log_error(error_msg)
            self.stats['errors'].append(error_msg)
            return None

    def process_telescope_date(self, telescope, date):
        """Process all targets for a specific telescope and date."""
        # Try v2 first, then fallback to v1
        pipeline_dir_v2 = f"{self.base_path}/PipelineOutput/v2/{telescope}/output/{date}/reduction"
        pipeline_dir_v1 = f"{self.base_path}/PipelineOutput/v1/{telescope}/output/{date}/reduction"

        pipeline_dir = None
        version = None

        # Check v2 first
        if os.path.exists(pipeline_dir_v2):
            v2_lists = [f for f in os.listdir(pipeline_dir_v2) if f.startswith('1_image_') and f.endswith('.list')]
            if v2_lists:
                pipeline_dir = pipeline_dir_v2
                version = "v2"

        # Fallback to v1 if v2 not available or empty
        if pipeline_dir is None and os.path.exists(pipeline_dir_v1):
            v1_lists = [f for f in os.listdir(pipeline_dir_v1) if f.startswith('1_image_') and f.endswith('.list')]
            if v1_lists:
                pipeline_dir = pipeline_dir_v1
                version = "v1"

        # If neither version has data
        if pipeline_dir is None:
            self.log_debug(f"No pipeline data found for {telescope}/{date}")
            return 0

        self.log_debug(f"Using {version} for {telescope}/{date}")
        processed_count = 0

        # Find all list files
        for filename in os.listdir(pipeline_dir):
            if filename.startswith('1_image_') and filename.endswith('.list'):
                list_filepath = os.path.join(pipeline_dir, filename)
                if self.process_target_list(list_filepath, telescope, date):
                    processed_count += 1

        return processed_count

    def collect_telescope_date_data(self, telescope, date, all_observations):
        """Collect observation data for a specific telescope and date without aggregating."""
        # Try v2 first, then fallback to v1
        pipeline_dir_v2 = f"{self.base_path}/PipelineOutput/v2/{telescope}/output/{date}/reduction"
        pipeline_dir_v1 = f"{self.base_path}/PipelineOutput/v1/{telescope}/output/{date}/reduction"

        pipeline_dir = None
        version = None

        # Check v2 first
        if os.path.exists(pipeline_dir_v2):
            v2_lists = [f for f in os.listdir(pipeline_dir_v2) if f.startswith('1_image_') and f.endswith('.list')]
            if v2_lists:
                pipeline_dir = pipeline_dir_v2
                version = "v2"

        # Fallback to v1 if v2 not available or empty
        if pipeline_dir is None and os.path.exists(pipeline_dir_v1):
            v1_lists = [f for f in os.listdir(pipeline_dir_v1) if f.startswith('1_image_') and f.endswith('.list')]
            if v1_lists:
                pipeline_dir = pipeline_dir_v1
                version = "v1"

        # If neither version has data
        if pipeline_dir is None:
            self.log_debug(f"No pipeline data found for {telescope}/{date}")
            return 0

        self.log_debug(f"Using {version} for {telescope}/{date}")
        processed_count = 0

        # Process all list files and store data without aggregating
        for filename in os.listdir(pipeline_dir):
            if filename.startswith('1_image_') and filename.endswith('.list'):
                list_filepath = os.path.join(pipeline_dir, filename)
                session_data = self.process_target_list_no_aggregate(list_filepath, telescope, date)
                if session_data:
                    target_name = session_data['target_name']
                    all_observations[date][target_name][telescope] = session_data
                    processed_count += 1

        return processed_count

    def process_target_list_no_aggregate(self, list_filepath, telescope, date):
        """Process a target list file but return data instead of aggregating."""
        if not os.path.exists(list_filepath):
            return None

        try:
            with open(list_filepath, 'r') as f:
                filepaths = [line.strip() for line in f if line.strip()]

            if not filepaths:
                return None

            # Extract target name from list filename
            list_filename = os.path.basename(list_filepath)
            target_match = re.match(r'1_image_(.+)\.list$', list_filename)
            if not target_match:
                self.log_warning(f"Could not extract target name from {list_filename}")
                return None

            target_name = target_match.group(1).upper()

            # Validate target name against coordinates (for wrong_target_names.txt)
            self.validate_target_name(target_name, filepaths[0], telescope, date)

            # Try timestamp-based calculation first
            duration = self.calculate_duration_from_timestamps(filepaths)
            method = 'timestamp'

            if duration is None:
                # Fallback to exposure time method
                duration = self.calculate_duration_from_exposure(filepaths, telescope, date)
                method = 'exposure'

                if duration is None:
                    failure_info = {
                        'target': target_name,
                        'telescope': telescope,
                        'date': date,
                        'num_files': len(filepaths)
                    }
                    self.stats['failed_calculations'].append(failure_info)
                    self.log_warning(f"Could not calculate duration for {target_name} on {telescope} {date}")
                    self.stats['failures'] += 1
                    return None

            # Return session data instead of aggregating
            session_data = {
                'target_name': target_name,
                'telescope': telescope,
                'date': date,
                'duration_hours': duration,
                'num_images': len(filepaths),
                'method': method
            }

            # Update statistics
            if method == 'timestamp':
                self.stats['timestamp_successes'] += 1
            else:
                self.stats['exposure_fallbacks'] += 1

            self.stats['total_sessions'] += 1

            self.log_info(f"Collected {target_name}: {duration:.3f}h using {method} method ({len(filepaths)} images)")
            return session_data

        except Exception as e:
            error_msg = f"Error processing {list_filepath}: {e}"
            self.log_error(error_msg)
            self.stats['errors'].append(error_msg)
            return None

    def process_date_observations(self, date, date_observations):
        """Process all observations for a specific date, handling multi-telescope duplicates."""
        for target_name, telescope_data in date_observations.items():
            if len(telescope_data) > 1:
                # Multiple telescopes observed this target on this date
                telescope_names = list(telescope_data.keys())
                durations = [data['duration_hours'] for data in telescope_data.values()]
                max_duration = max(durations)
                max_telescope = telescope_names[durations.index(max_duration)]

                self.log_info(
                    f"Multi-telescope observation {date}/{target_name}: {len(telescope_data)} telescopes, using max {max_duration:.3f}h from {max_telescope}")

                # Use the session data from the telescope with maximum duration
                best_session = telescope_data[max_telescope]
            else:
                # Only one telescope observed this target
                best_session = list(telescope_data.values())[0]

            # Add to target_data (this is where aggregation finally happens)
            target_name = best_session['target_name']
            duration = best_session['duration_hours']

            if target_name not in self.target_data:
                self.target_data[target_name] = {
                    'total_hours': 0.0,
                    'sessions': [],
                    'total_images': 0,
                    'telescopes': set(),
                    'dates': []
                }

            self.target_data[target_name]['total_hours'] += duration
            self.target_data[target_name]['sessions'].append(best_session)
            self.target_data[target_name]['total_images'] += best_session['num_images']
            self.target_data[target_name]['telescopes'].add(best_session['telescope'])
            self.target_data[target_name]['dates'].append(date)

    def apply_minimum_threshold_filter(self):
        """Remove targets with less than 2 hours of observation time."""
        min_threshold = 2.0  # hours
        targets_to_remove = []

        # Identify targets below threshold
        for target_name, data in self.target_data.items():
            if data['total_hours'] < min_threshold:
                targets_to_remove.append(target_name)

        # Remove targets below threshold
        removed_count = 0
        for target_name in targets_to_remove:
            total_hours = self.target_data[target_name]['total_hours']
            self.target_data.pop(target_name, None)
            removed_count += 1
            self.log_debug(
                f"Removed target {target_name} with {total_hours:.3f} hours (below {min_threshold:.1f}h threshold)")

        if removed_count > 0:
            self.log_info(f"Removed {removed_count} targets with less than {min_threshold:.1f} hours observation time")
        else:
            self.log_info(f"No targets removed by {min_threshold:.1f}h minimum threshold")

    def get_date_range(self, telescope, start_date=None, end_date=None):
        """Get list of dates to process for a telescope."""
        obs_dir = f"{self.observations_path}/{telescope}/images"

        if not os.path.exists(obs_dir):
            self.log_warning(f"Observations directory not found: {obs_dir}")
            return []

        # Get all date directories
        all_dates = []
        for item in os.listdir(obs_dir):
            if re.match(r'^\d{8}$', item):  # YYYYMMDD format
                all_dates.append(item)

        all_dates.sort()

        # Filter by date range if specified
        if start_date or end_date:
            filtered_dates = []
            for date in all_dates:
                if start_date and date < start_date:
                    continue
                if end_date and date > end_date:
                    continue
                filtered_dates.append(date)
            return filtered_dates

        return all_dates

    def run_analysis(self, telescope_filter=None, start_date=None, end_date=None):
        """Run the main analysis with multi-telescope duplicate handling."""
        telescopes_to_process = [telescope_filter] if telescope_filter else self.telescopes

        self.log_info(f"Starting analysis for telescopes: {telescopes_to_process}")
        if start_date or end_date:
            self.log_info(f"Date range: {start_date or 'beginning'} to {end_date or 'end'}")

        # Phase 1: Collect all observations
        self.log_info("Phase 1: Collecting all observations...")
        all_observations = defaultdict(lambda: defaultdict(dict))  # {date: {target: {telescope: session_data}}}
        all_dates = set()
        total_processed = 0

        for telescope in telescopes_to_process:
            self.log_info(f"Collecting data from telescope: {telescope}")

            dates = self.get_date_range(telescope, start_date, end_date)
            self.log_info(f"Found {len(dates)} dates for {telescope}")
            all_dates.update(dates)

            for date in dates:
                processed = self.collect_telescope_date_data(telescope, date, all_observations)
                total_processed += processed

                if processed > 0:
                    self.log_info(f"  {date}: {processed} targets collected")

        # Phase 2: Process chronologically to handle multi-telescope observations
        self.log_info(f"Phase 2: Processing {len(all_dates)} dates chronologically...")
        for date in sorted(all_dates):
            self.process_date_observations(date, all_observations[date])

        # Apply minimum threshold filter (same as nights_per_field.py)
        self.apply_minimum_threshold_filter()

        self.stats['total_targets'] = len(self.target_data)
        self.log_info(f"Analysis complete. Processed {total_processed} target observations.")

    def generate_csv_output(self, output_file='observation_summary.csv'):
        """Generate simple table output with target names and hours."""
        if not self.target_data:
            self.log_warning("No data to write to output")
            return

        # Create reports directory if it doesn't exist
        reports_dir = '/data/SPECULOOSPipeline/reports'
        Path(reports_dir).mkdir(parents=True, exist_ok=True)

        output_path = os.path.join(reports_dir, 'SurveyTotal')

        with open(output_path, 'w') as f:
            # Write table header with exact spacing
            f.write("|                    Target |  Hours |\n")
            f.write("| ------------------------- | ------ |\n")

            # Sort targets by name and write data with exact spacing
            for target_name in sorted(self.target_data.keys()):
                data = self.target_data[target_name]
                # Convert target name to uppercase
                target_upper = target_name.upper()
                # Format with exact spacing: 25 chars right-aligned for target, 6 chars for hours
                # Calculate digits before decimal point
                digits_before_decimal = len(str(int(data['total_hours'])))
                decimal_places = 4 - digits_before_decimal

                # Ensure we don't go negative on decimal places
                decimal_places = max(0, decimal_places)

                # Format with calculated decimal places
                hours_str = f"{data['total_hours']:.{decimal_places}f}"

                f.write(f"| {target_upper:>25} | {hours_str:>6} |\n")

        self.log_info(f"Survey total output written to {output_path}")

        # Generate ObservedVsProcessed output
        self.generate_observed_vs_processed_output()

        # Generate ByTelescope files
        self.generate_by_telescope_files()

        # Generate histogram plots
        self.generate_histogram_plots()

        # Write wrong target names file if any found
        if hasattr(self, 'wrong_names') and self.wrong_names:
            wrong_names_file = '/data/SPECULOOSPipeline/reports/wrong_target_names.txt'
            with open(wrong_names_file, 'w') as f:
                for i in range(len(self.wrong_names)):
                    f.write(f"{self.wrong_names[i]}\t{self.corrected_names[i]}\n")
            self.log_info(f"Wrong target names written to {wrong_names_file}")

    def generate_by_telescope_files(self):
        """Generate ByTelescope files in the same format as nights_per_field.py."""
        if not self.target_data:
            self.log_warning("No data to write to ByTelescope files")
            return

        reports_dir = '/data/SPECULOOSPipeline/reports'

        # Prepare data for ByTelescope format
        targets = []
        telescopes = []
        hours = []

        # Extract all sessions data
        for target_name, data in sorted(self.target_data.items()):
            for session in data['sessions']:
                targets.append(target_name)
                telescopes.append(session['telescope'])
                hours.append(session['duration_hours'])

        if ASTROPY_TABLE_AVAILABLE:
            # Use astropy.Table format (same as nights_per_field.py)
            self._generate_astropy_by_telescope(reports_dir, targets, telescopes, hours)
        else:
            # Fallback to simple text format
            self._generate_simple_by_telescope(reports_dir, targets, telescopes, hours)

    def _generate_astropy_by_telescope(self, reports_dir, targets, telescopes, hours):
        """Generate ByTelescope files using astropy.Table format."""
        # Create test row to set column widths
        test_row = ['-------------------------', '----------', '------']
        hdr = ['Target', 'Telescope', 'Hours']
        dtypes = ['U25', 'U10', 'U6']  # Changed from 'S' (bytes) to 'U' (Unicode) for Python 3

        # Create the table
        report = Table(names=hdr, dtype=dtypes)
        report.add_row(test_row)

        # Add data rows
        for i in range(len(targets)):
            row = [targets[i], telescopes[i], f"{hours[i]:.4f}"]
            report.add_row(row)

        # Write ByTelescope file
        output_path = os.path.join(reports_dir, 'ObservationsByTelescope')
        ascii.write(report, output_path, format='fixed_width', overwrite=True)
        self.log_info(f"ByTelescope file written to {output_path}")

    def generate_observed_vs_processed_output(self):
        """Generate ObservedVsProcessed output comparing observed vs processed hours."""
        if not self.target_data:
            self.log_warning("No data to write to ObservedVsProcessed output")
            return

        # Create reports directory if it doesn't exist
        reports_dir = '/data/SPECULOOSPipeline/reports'
        Path(reports_dir).mkdir(parents=True, exist_ok=True)

        output_path = os.path.join(reports_dir, 'ObservedVsProcessed')

        # Debug counters
        total_sessions_checked = 0
        dirs_found = 0
        dirs_with_diff_files = 0

        with open(output_path, 'w') as f:
            # Write table header
            f.write("|                    Target | Observed Hours | Processed Hours | Difference |\n")
            f.write("| ------------------------- | -------------- | --------------- | ---------- |\n")

            # Process each target
            for target_name in sorted(self.target_data.keys()):
                data = self.target_data[target_name]
                target_upper = target_name.upper()

                observed_hours = data['total_hours']
                processed_hours = 0.0

                # Check each session to see if it has been processed
                for session in data['sessions']:
                    telescope = session['telescope']
                    date = session['date']
                    total_sessions_checked += 1

                    # Find the actual directory name (case-insensitive search)
                    base_dir = f"/data/SPECULOOSPipeline/PipelineOutput/v2/{telescope}/output/{date}"
                    processed_dir = None

                    # Debug: Log first few attempts
                    if total_sessions_checked <= 5:
                        self.log_info(f"DEBUG: Checking base directory: {base_dir}")

                    if os.path.exists(base_dir):
                        # Get all directories in the base directory
                        try:
                            all_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

                            # Find directory that matches target name (case-insensitive)
                            for dir_name in all_dirs:
                                if dir_name.upper() == target_name.upper():
                                    processed_dir = os.path.join(base_dir, dir_name)
                                    break

                            if total_sessions_checked <= 5:
                                self.log_info(
                                    f"DEBUG: Found {len(all_dirs)} target directories, looking for match with {target_name}")
                                if processed_dir:
                                    self.log_info(f"DEBUG: Matched directory: {processed_dir}")
                                else:
                                    self.log_info(f"DEBUG: No matching directory found. Available: {all_dirs[:5]}")

                        except OSError:
                            if total_sessions_checked <= 5:
                                self.log_info("DEBUG: Could not list directory contents")

                    if processed_dir and os.path.exists(processed_dir):
                        dirs_found += 1
                        if total_sessions_checked <= 5:
                            self.log_info("DEBUG: Directory exists, listing contents...")

                        # Look for files ending in "5_diff.fits"
                        diff_files = [filename for filename in os.listdir(processed_dir) if
                                      filename.endswith('5_diff.fits')]

                        if total_sessions_checked <= 5:
                            all_files = os.listdir(processed_dir)
                            self.log_info(
                                f"DEBUG: Found {len(all_files)} total files, {len(diff_files)} ending in 5_diff.fits")
                            if len(all_files) > 0:
                                self.log_info(f"DEBUG: Sample files: {all_files[:3]}")

                        if diff_files:
                            dirs_with_diff_files += 1
                            # Add this session's hours to processed total
                            processed_hours += session['duration_hours']
                            if total_sessions_checked <= 5:
                                self.log_info(f"DEBUG: Found diff files, adding {session['duration_hours']:.3f}h")
                    else:
                        if total_sessions_checked <= 5:
                            self.log_info("DEBUG: No matching directory found")

                # Calculate difference
                difference = observed_hours - processed_hours

                # Format hours with appropriate decimal places (same logic as SurveyTotal)
                def format_hours(hours):
                    digits_before_decimal = len(str(int(hours)))
                    decimal_places = max(0, 4 - digits_before_decimal)
                    return f"{hours:.{decimal_places}f}"

                observed_str = format_hours(observed_hours)
                processed_str = format_hours(processed_hours)
                difference_str = format_hours(difference)

                # Write the row with proper spacing
                f.write(f"| {target_upper:>25} | {observed_str:>14} | {processed_str:>15} | {difference_str:>10} |\n")

        # Print debug summary
        self.log_info("DEBUG SUMMARY:")
        self.log_info(f"Total sessions checked: {total_sessions_checked}")
        self.log_info(f"Directories found: {dirs_found}")
        self.log_info(f"Directories with 5_diff.fits files: {dirs_with_diff_files}")
        self.log_info(f"ObservedVsProcessed output written to {output_path}")

    def _generate_simple_by_telescope(self, reports_dir, targets, telescopes, hours):
        """Generate ByTelescope files using simple text format."""
        output_path = os.path.join(reports_dir, 'ObservationsByTelescope')

        with open(output_path, 'w') as f:
            f.write(f"| {'Target':>25} | {'Telescope':>10} | {'Hours':>6} |\n")
            f.write(f"| {'-' * 25:>25} | {'-' * 10:>10} | {'-' * 6:>6} |\n")

            for i in range(len(targets)):
                f.write(f"| {targets[i]:>25} | {telescopes[i]:>10} | {hours[i]:>6.4f} |\n")

        self.log_info(f"ByTelescope file (simple format) written to {output_path}")

    def generate_histogram_plots(self):
        """Generate cumulative histogram plots in the same format as nights_per_field.py."""
        if not MATPLOTLIB_AVAILABLE:
            self.log_warning("Matplotlib not available. Skipping histogram generation.")
            return

        if not self.target_data:
            self.log_warning("No data for histogram generation")
            return

        # Collect total hours for all targets
        all_hours = [data['total_hours'] for data in self.target_data.values()]

        # Filter for targets with >= 20 hours (same as nights_per_field.py)
        hours_hist = [h for h in all_hours if h >= 20.0]

        if not hours_hist:
            self.log_info("No targets with >= 20 hours for histogram")
            return

        # Create bins (same as nights_per_field.py)
        bins = np.arange(20, np.ceil(np.nanmax(hours_hist)) + 10, 10)

        # Create the plot
        plt.figure(figsize=(16, 16))
        plt.hist(hours_hist, bins=bins, cumulative=-1, facecolor='green', edgecolor='black')
        plt.grid(visible=True, axis='y')  # Changed 'b' to 'visible' for matplotlib 3.5+
        plt.title(f"Cumulative Histogram of Number of Observed Hours for {len(hours_hist)} Targets")
        plt.xlabel("Number of Hours Observed")
        plt.ylabel("Number of Targets (with >20 hours Observation)")

        # Save the plot
        output_path = '/data/SPECULOOSPipeline/reports/Observations_hist.png'
        plt.savefig(output_path)
        plt.close()

        self.log_info(f"Histogram plot saved to {output_path}")
        self.log_info(f"Histogram includes {len(hours_hist)} targets with >= 20 hours observation")

    def generate_detailed_csv(self, output_file='observation_detailed.csv'):
        """Generate detailed CSV with individual session data."""
        if not self.target_data:
            self.log_warning("No data to write to detailed CSV")
            return

        fieldnames = [
            'Target_Name',
            'Telescope',
            'Date',
            'Duration_Hours',
            'Number_of_Images',
            'Calculation_Method'
        ]

        with open(output_file, 'w', newline='') as csvfile:  # Python 3 uses newline='' instead of 'wb'
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for target_name in sorted(self.target_data.keys()):
                data = self.target_data[target_name]
                for session in data['sessions']:
                    row = {
                        'Target_Name': target_name,
                        'Telescope': session['telescope'],
                        'Date': session['date'],
                        'Duration_Hours': f"{session['duration_hours']:.3f}",
                        'Number_of_Images': session['num_images'],
                        'Calculation_Method': session['method']
                    }
                    writer.writerow(row)

        self.log_info(f"Detailed CSV output written to {output_file}")

    def print_failed_calculations(self):
        """Print details of failed calculations."""
        if not self.stats['failed_calculations']:
            return

        print("\nFAILED CALCULATIONS DETAILS:")
        print("-" * 60)
        for failure in self.stats['failed_calculations']:
            print(f"  {failure['telescope']}/{failure['date']}: {failure['target']} ({failure['num_files']} files)")

    def print_summary_statistics(self):
        """Print summary statistics."""
        print("\n" + "=" * 60)
        print("OBSERVATION ANALYSIS SUMMARY")
        print("=" * 60)
        print(f"Total targets found: {self.stats['total_targets']}")
        print(f"Total observation sessions: {self.stats['total_sessions']}")
        print(f"Successful timestamp calculations: {self.stats['timestamp_successes']}")
        print(f"Exposure time fallback calculations: {self.stats['exposure_fallbacks']}")
        print(f"Failed calculations: {self.stats['failures']}")

        if self.target_data:
            total_hours = sum(data['total_hours'] for data in self.target_data.values())
            total_images = sum(data['total_images'] for data in self.target_data.values())
            print(f"Total observation time: {total_hours:.1f} hours")
            print(f"Total images processed: {total_images}")

            # Top targets by observation time
            target_items = list(self.target_data.items())
            top_targets = sorted(target_items, key=lambda x: x[1]['total_hours'], reverse=True)[:5]
            print("\nTop 5 targets by observation time:")
            for target, data in top_targets:
                print(f"  {target}: {data['total_hours']:.1f}h")

        if self.stats['errors']:
            print(f"\nErrors encountered: {len(self.stats['errors'])}")
            for error in self.stats['errors'][:5]:  # Show first 5 errors
                print(f"  {error}")
            if len(self.stats['errors']) > 5:
                print(f"  ... and {len(self.stats['errors']) - 5} more")

        self.print_failed_calculations()
        print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description='Calculate observation times for SPECULOS targets')
    parser.add_argument('--telescope', choices=['Artemis', 'Callisto', 'Europa', 'Ganymede', 'Io'],
                        help='Process only this telescope')
    parser.add_argument('--start-date', help='Start date (YYYYMMDD format)')
    parser.add_argument('--end-date', help='End date (YYYYMMDD format)')
    parser.add_argument('--output', default='observation_summary.csv',
                        help='Output CSV filename (default: observation_summary.csv)')
    parser.add_argument('--detailed', action='store_true',
                        help='Also output detailed session data')

    args = parser.parse_args()

    # Validate date formats
    for date_arg, date_name in [(args.start_date, 'start-date'), (args.end_date, 'end-date')]:
        if date_arg and not re.match(r'^\d{8}$', date_arg):
            print(f"Error: {date_name} must be in YYYYMMDD format")
            return 1

    # Create calculator and run analysis
    calc = ObservationCalculator()

    # Setup logging based on output filename
    calc.setup_logging(args.output)

    try:
        calc.run_analysis(args.telescope, args.start_date, args.end_date)

        # Generate outputs
        calc.generate_csv_output()
        calc.print_summary_statistics()

        # Generate detailed session data if requested
        if args.detailed:
            detailed_file = args.output.replace('.csv', '_detailed.csv')
            calc.generate_detailed_csv(detailed_file)

    finally:
        # Always close the log file
        calc.close_logging()

    return 0


if __name__ == '__main__':
    sys.exit(main())