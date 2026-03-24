"""
ESO Archive downloader for SPECULOOS data
Python 3 version using modern ESO APIs
"""

import sys
print(sys.path)
print(sys.version)
import smtplib
import time
import datetime as dt
import os
import glob
import timeit
import numpy as np
from download.request_eso import ESODownloader
from download.astra_transform import astra_transform, transformation_check
from dotenv import load_dotenv
# Load .env file from src/reporting directory
env_path = os.path.join(os.path.dirname(__file__), '..', 'reporting', '.env')
load_dotenv(dotenv_path=env_path)

# Load credentials from environment variables
ESO_username = os.getenv('ESO_USERNAME')
ESO_password = os.getenv('ESO_PASSWORD')

def get_email_config():
    """Load email configuration from environment variables (same as email_eon.py)."""
    config = {
        'from_addr': os.getenv('EMAIL_ADDRESS'),
        'password': os.getenv('EMAIL_PASSWORD'),
        'smtp_server': os.getenv('SMTP_SERVER', 'smtp.gmail.com:587'),
        'to_addr_list': os.getenv('TO_EMAIL_LIST', '').split(',') if os.getenv('TO_EMAIL_LIST') else [],
        'cc_addr_list': os.getenv('CC_EMAIL_LIST', '').split(',') if os.getenv('CC_EMAIL_LIST') else []
    }

    # Clean up email lists (remove empty strings and whitespace)
    config['to_addr_list'] = [email.strip() for email in config['to_addr_list'] if email.strip()]
    config['cc_addr_list'] = [email.strip() for email in config['cc_addr_list'] if email.strip()]

    return config

# Debug: Check if variables are loaded
if ESO_username and ESO_password:
    print(f"✓ ESO credentials loaded for user: {ESO_username}")
else:
    print("✗ ESO credentials not found in environment variables")
    sys.exit(1)

# Test email config
try:
    config = get_email_config()
    print(f"✓ Email config loaded - sending from: {config['from_addr']}")
except:
    print("✗ Email configuration incomplete")


# AFTER:
def parse_timestamp_from_dp_id(dp_id):
    """
    Extract ISO timestamp from ESO dp_id
    Example: SPECU4.2025-12-16T04:49:08.869 -> 2025-12-16T04:49:08.869
    """
    try:
        # dp_id format: INSTRUMENT.YYYY-MM-DDTHH:MM:SS.mmm (no .fits extension)
        # Split on first dot to separate instrument from timestamp
        if '.' not in dp_id:
            return None

        # Everything after the first dot is the timestamp
        timestamp = dp_id.split('.', 1)[1]

        # Validate it looks like a timestamp (basic check)
        if 'T' in timestamp and len(timestamp) > 10:
            return timestamp

        return None
    except Exception as e:
        print(f"Warning: Could not parse timestamp from {dp_id}: {e}")
        return None

def subtract_time_buffer(timestamp_str, buffer_minutes=2):
    """
    Subtract buffer minutes from ISO timestamp string
    Returns new timestamp string in format expected by ESO TAP: 'YYYY-MM-DD HH:MM:SS'
    """
    try:
        from datetime import datetime, timedelta
        # Parse timestamp (format: 2025-12-16T04:49:08.869)
        dt_obj = datetime.fromisoformat(timestamp_str)
        # Subtract buffer
        dt_buffered = dt_obj - timedelta(minutes=buffer_minutes)
        # Return in 'YYYY-MM-DD HH:MM:SS' format (no milliseconds, space separator)
        return dt_buffered.isoformat()
    except Exception as e:
        print(f"Warning: Could not subtract buffer from {timestamp_str}: {e}")
        return timestamp_str

def get_images(tel, telname, stime, etime, date, imgdir, down_im, test_mode=False, verbose=False, query_only=False, skip_transform=False):
    """
    Download images from ESO archive using modern API

    Args:
        tel: Program ID (e.g., '60.A-9009(C)')
        telname: Telescope name (e.g., 'Ganymede')
        stime: Start time in YYYY-MM-DD format
        etime: End time in YYYY-MM-DD format
        date: Date string for logging
        imgdir: Destination directory
        down_im: List of already downloaded images
        test_mode: If True, only download the first image found

    Returns:
        tuple: (num_downloaded, total_eso_files)
    """
    print("Getting images from ESO archive")
    if test_mode:
        print("*** TEST MODE: Will only download the first image found ***")

    try:
        # Initialize modern ESO downloader
        downloader = ESODownloader(ESO_username, ESO_password)

        total_eso_files = 0
        downloaded_files = []

        for data_type in ['SCIENCE', 'CALIB']:
            t0 = timeit.default_timer()

            # Pagination loop for handling >10k results
            all_files = []
            seen_dp_ids = set()  # Track unique dp_ids to avoid duplicates from pagination overlap
            query_start_time = stime
            pagination_attempt = 0
            max_pagination_attempts = 50  # Safety limit

            while pagination_attempt < max_pagination_attempts:
                pagination_attempt += 1

                # Query for files using modern TAP service
                if pagination_attempt == 1:
                    print(f"Querying for {data_type} images...")
                else:
                    print(f"Pagination query {pagination_attempt} starting from {query_start_time}...")

                # Use MJD if available (for pagination), otherwise use date strings
                if pagination_attempt == 1:
                    files = downloader.query_files(tel, query_start_time, etime, data_type)
                    # Store end MJD for pagination
                    end_mjd = downloader._date_to_mjd_night_end(etime)
                else:
                    # Pagination query with MJD
                    files = downloader.query_files(tel, stime, etime, data_type,
                                                   start_mjd=query_start_mjd,
                                                   end_mjd=end_mjd)

                if not files:
                    if pagination_attempt == 1:
                        print(f"WARNING: no {data_type} results found for {telname} on {date}")
                    break

                # Filter out duplicates from pagination overlap
                unique_files_in_batch = []
                duplicates_in_batch = 0
                for file_row in files:
                    dp_id = file_row[0]  # dp_id is first column
                    if dp_id not in seen_dp_ids:
                        seen_dp_ids.add(dp_id)
                        unique_files_in_batch.append(file_row)
                    else:
                        duplicates_in_batch += 1

                num_files = len(files)
                num_unique = len(unique_files_in_batch)
                print(
                    f"Retrieved {num_files} files in this batch ({num_unique} unique, {duplicates_in_batch} duplicates)")
                all_files.extend(unique_files_in_batch)

                # Check if we hit the 10k limit
                if num_files >= 10000:
                    print(f"Hit 10,000 result limit - paginating to get remaining files...")

                    # Extract timestamp from last file
                    last_dp_id = files[-1][0]  # dp_id is first column
                    last_timestamp = parse_timestamp_from_dp_id(last_dp_id)

                    if last_timestamp:
                        # Subtract buffer and convert to MJD for next query
                        buffered_timestamp = subtract_time_buffer(last_timestamp, buffer_minutes=2)
                        query_start_mjd = downloader._timestamp_to_mjd(buffered_timestamp)

                        if query_start_mjd:
                            print(f"Last file timestamp: {last_timestamp}")
                            print(f"Next query will start from: {buffered_timestamp} (MJD: {query_start_mjd:.6f})")
                            # Store MJD for next iteration
                            query_start_time = None  # Signal to use MJD instead
                        else:
                            print(f"WARNING: Could not convert timestamp to MJD - stopping pagination")
                            break
                    else:
                        print(f"WARNING: Could not parse timestamp from {last_dp_id} - stopping pagination")
                        break
                else:
                    # Got fewer than 10k results - we're done
                    break

            if pagination_attempt >= max_pagination_attempts:
                print(f"WARNING: Hit maximum pagination attempts ({max_pagination_attempts}) - may be missing files")

            files = all_files  # Use the combined results

            t1 = timeit.default_timer()
            elapsed = t1 - t0
            print(f"Time taken to query archive for {data_type} images: {elapsed / 60:.3f} minutes")

            if not files:
                continue

            num_files = len(files)
            print(f"Total {num_files} {data_type} images for {telname} on {date} (after pagination)")

            if num_files >= 10000 and pagination_attempt == 1:
                print(f"WARNING: Exactly 10000 files but pagination not triggered - check logic")

            total_eso_files += num_files

            if num_files == 0:
                continue

            # Filter out already downloaded files
            files_to_download = []
            already_downloaded_count = 0

            for file_row in files:
                dp_id = file_row[0]  # dp_id is first column
                if dp_id not in down_im:
                    files_to_download.append(file_row)
                    # In test mode, only download the first image
                    if test_mode:
                        print(f"*** TEST MODE: Found first image {dp_id}, stopping search ***")
                        break
                else:
                    already_downloaded_count += 1

            # Print summary statistics
            total_found = len(files)
            need_downloading = len(files_to_download)
            print(f"Archive summary for {data_type}:")
            print(f"  Total images found on archive: {total_found}")
            print(f"  Already downloaded: {already_downloaded_count}")
            print(f"  Need downloading: {need_downloading}")

            # In test mode, if we found a file to download, skip other data types
            if test_mode and files_to_download:
                print("*** TEST MODE: Skipping remaining data types ***")
                total_eso_files += num_files  # Still count all available files
                # Download the single test file
                if files_to_download:
                    print(f"Downloading 1 test file from {data_type}...")
                    print(dt.datetime.now())

                    try:
                        batch_downloaded = downloader.download_files(
                            files_to_download,
                            imgdir,
                            verbose=verbose,
                            query_only=query_only
                        )

                        # Skip astra_transform in query_only or skip_transform mode
                        if not query_only and not skip_transform:
                            # Run astra_transform on downloaded files
                            dp_ids = [file_row[0] for file_row in files_to_download]
                            astra_transform(dp_ids, imgdir, 1)

                        downloaded_files.extend(batch_downloaded)
                        print(dt.datetime.now())
                        print(f"*** TEST MODE: Successfully downloaded test file ***")

                    except Exception as e:
                        print(f"Test download failed: {e}")
                        print(dt.datetime.now())

                break  # Exit the data_type loop

            # Download files if any need downloading
            if files_to_download:
                print(f"Downloading {len(files_to_download)} new {data_type} files...")
                print(dt.datetime.now())

                try:
                    # Download files using modern API
                    batch_downloaded = downloader.download_files(
                        files_to_download,
                        imgdir,
                        verbose=verbose,
                        query_only=query_only,
                        parallel=True  # Use parallel downloads for efficiency
                    )

                    # Skip astra_transform in query_only or skip_transform mode
                    if not query_only and not skip_transform:
                        # Run astra_transform on downloaded files
                        dp_ids = [file_row[0] for file_row in files_to_download]
                        astra_transform(dp_ids, imgdir, 1)

                    downloaded_files.extend(batch_downloaded)
                    print(dt.datetime.now())

                except Exception as e:
                    print(f"Download failed for {data_type} batch: {e}")
                    print(dt.datetime.now())
                    continue

        num_downloaded = len(downloaded_files)
        print(f"Successfully downloaded {num_downloaded} files out of {total_eso_files} available")

        return num_downloaded, total_eso_files

    except Exception as e:
        print(f"Error in get_images: {e}")
        return 0, 0


def check_num_images(imgdir, num_im_eso):
    """Check if the correct number of images were downloaded"""
    ok = True
    num_im_dir = 0

    try:
        xlist = os.listdir(imgdir)
        for x in xlist:
            if x[-5:] == ".fits" or x[-4:] == ".fts":
                num_im_dir = num_im_dir + 1

        print(num_im_dir, num_im_eso)
    except Exception as e:
        print(e)
        num_im_dir = 0

    if num_im_eso != num_im_dir:
        print(f"WARNING: Not all images downloaded in {imgdir}")
        ok = False

    return ok, num_im_dir


def num_im_transfer(tel, date, dir):
    """Get number of images transferred from Chile"""
    found = False

    transferlog = f"{dir}/{tel}/transfer_log.txt"
    try:
        with open(transferlog, 'r') as reader:
            for line in reader:
                line = line.strip()
                d = line.split(" ")[0]
                num_transfer = int(line.split(" ")[1])

                if d == date:
                    found = True
                    break
    except FileNotFoundError:
        print(f"Transfer log not found: {transferlog}")
        return 0, False

    if found == False:
        print("Transfer Log not updated yet!")
        num_transfer = 0

    return num_transfer, found


def check_download_images(num_i_transfer, num_i):
    """Check if downloaded images match transferred images"""
    ok = True
    if num_i_transfer != num_i:
        ok = False
        print("Number of images transferred to ESO Archive is different to number of images on archive")
    return ok


def write_to_log(im, i_transfer, i_eso, dir, tel, telroute, date):
    """Write download results to log file using simple CSV format"""
    logfile = f"{dir}/{telroute}/download_log.csv"

    # Create header if file doesn't exist
    if not os.path.exists(logfile):
        with open(logfile, 'w') as f:
            f.write("Night,Telescope,ESO_Archive,Transferred,Downloaded\n")

    # Read existing data
    existing_data = {}
    if os.path.exists(logfile):
        with open(logfile, 'r') as f:
            lines = f.readlines()[1:]  # Skip header
            for line in lines:
                if line.strip():
                    parts = line.strip().split(',')
                    if len(parts) >= 5:
                        key = (parts[0], parts[1])  # (date, telescope)
                        existing_data[key] = parts

    # Update or add this entry
    key = (date, tel)
    existing_data[key] = [date, tel, str(i_eso), str(i_transfer), str(im)]

    # Write all data back (sorted by date)
    with open(logfile, 'w') as f:
        f.write("Night,Telescope,ESO_Archive,Transferred,Downloaded\n")
        for entry in sorted(existing_data.values()):
            f.write(','.join(entry) + '\n')

    print(f"Written to log file: {logfile}")


def clean_cache():
    """Clean astropy cache (legacy function, may not be needed)"""
    cache = '/home/cam217/.astropy/cache/astroquery/Eso/'
    if os.path.exists(cache):
        for f in os.listdir(cache):
            file_path = os.path.join(cache, f)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)


def send_final_email_summary(telescope_states, date, telnames, test_mode, is_recent,
                            grace_days, total_attempts):
    """Send comprehensive email summary after all retry attempts complete"""

    # Build summary for each telescope
    telescope_summaries = []
    any_failures = False
    any_success = False

    for idx, state in telescope_states.items():
        telname = state['telname']
        success = state['success']
        attempts = state['attempts']

        if not attempts:
            continue

        last_attempt = attempts[-1]

        # Build attempt history
        history = []
        for att in attempts:
            history.append(
                f"  Attempt {att['attempt_num']}: "
                f"Archive={att['archive_count']}, "
                f"Downloaded={att['downloaded_count']}, "
                f"Expected={att['expected_count'] if att['expected_count'] is not None else 'N/A'}"
            )
        history_str = "\n".join(history)

        # Determine status message
        if success:
            any_success = True
            if last_attempt['expected_count'] is not None:
                status = f"✓ SUCCESS: {last_attempt['downloaded_count']}/{last_attempt['expected_count']} images"
            else:
                status = f"✓ SUCCESS: {last_attempt['downloaded_count']} images downloaded"
        else:
            any_failures = True
            if last_attempt['expected_count'] is not None and last_attempt['expected_count'] > 0:
                status = f"✗ INCOMPLETE: {last_attempt['downloaded_count']}/{last_attempt['expected_count']} images"
            else:
                status = f"✗ INCOMPLETE: Transfer log not updated (0 images expected)"

        telescope_summaries.append(f"""
{telname}:
{status}
Attempts: {len(attempts)}/{total_attempts}
{history_str}
        """.strip())

    # Construct overall subject and message
    if test_mode:
        subject = "ESO Test Download Complete"
        status_line = "Test mode - downloaded first image only"
    elif any_failures:
        if is_recent:
            subject = f"ESO Download INCOMPLETE after {total_attempts} attempts"
            status_line = f"Some telescopes incomplete after {total_attempts} attempts over ~{sum([15, 30] + [60] * (total_attempts - 3))} minutes"
        else:
            subject = "ESO Download Complete (old data)"
            status_line = "Downloaded all available archive data (no retry for old observations)"
    else:
        subject = "ESO Download SUCCESS"
        status_line = f"All telescopes complete (took {len(telescope_states[0]['attempts'])} attempt(s))" if telescope_states else "Complete"

    message = f"""
ESO Archive Download Summary
Date: {date}
Data age: {"Recent (within grace period)" if is_recent else f"Old (>{grace_days} days)"}
Status: {status_line}

{'='*60}
TELESCOPE RESULTS:
{'='*60}

{chr(10).join(telescope_summaries)}

{'='*60}
    """.strip()

    # Add recommendations for failures
    if any_failures and is_recent:
        message += """

RECOMMENDED ACTIONS:
- Images may still be transferring from Chile to ESO archive
- Transfer log may not be updated yet
- Run the same command again to retry
- Check transfer_log.txt for expected image counts
        """

    # Send email
    try:
        config = get_email_config()
        if config['to_addr_list']:
            sendemail(
                from_addr=config['from_addr'],
                to_addr_list=config['to_addr_list'],
                cc_addr_list=config['cc_addr_list'],
                subject=subject,
                message=message,
                login=config['from_addr'],
                password=config['password'],
                smtpserver=config['smtp_server']
            )
            print("Email notification sent")
        else:
            print("No email recipients configured, skipping email notification")
    except Exception as e:
        print(f"Email sending failed: {e}")


def eso_download(dir, telname, sdate, edate, wait, test_mode=False, query_only=False, verbose=False,
                 transfer_log_grace_days=5, max_retries=13, skip_transform=False):
    """
    Main download function with retry logic

    Args:
        transfer_log_grace_days: Number of days after observation to enforce transfer log checking (default 5)
        max_retries: Maximum number of download attempts (default 13, covering ~12 hours)
    """

    # Generate retry intervals: [15, 30, 60, 60, 60, ...] up to max_retries-1 intervals
    retry_intervals = [15, 30] + [60] * (max_retries - 3) if max_retries > 2 else [15, 30][:max_retries-1]

    if telname == None:
        tel = ['60.A-9009(A)', '60.A-9009(B)', '60.A-9009(C)', '60.A-9009(D)']
        telname = ['Io', 'Europa', 'Ganymede', 'Callisto']
        telroute = ['Observations/' + s for s in telname]
    else:
        telname_str = str(telname)
        if telname_str == 'Io':
            tel = ['60.A-9009(A)']
        elif telname_str == 'Europa':
            tel = ['60.A-9009(B)']
        elif telname_str == 'Ganymede':
            tel = ['60.A-9009(C)']
        elif telname_str == 'Callisto':
            tel = ['60.A-9009(D)']
        telname = [telname_str]
        telroute = ['Observations/' + telname[0]]


    s = dt.datetime.strptime(sdate, "%Y%m%d")
    e = dt.datetime.strptime(edate, "%Y%m%d")
    days = (e - s).days

    print(f"Downloading {days} days of data from ESO archive")
    if test_mode:
        print("*** TEST MODE: Will only download the first image found ***")

    # Process each date
    for d in (s + dt.timedelta(n) for n in range(days)):
        os.system(f"rm -rf {dir}/.astropy/cache/astroquery/Eso/*")

        date = d.strftime("%Y%m%d")
        stime = d.strftime("%Y-%m-%d")
        etime = (d + dt.timedelta(days=1)).strftime("%Y-%m-%d")

        # Check if this is recent data (within grace period)
        days_since_observation = (dt.datetime.now() - d).days
        is_recent = days_since_observation <= transfer_log_grace_days

        if is_recent:
            print(f"Data is {days_since_observation} days old (within {transfer_log_grace_days} day grace period)")
            print("Will check against transfer log and retry if needed")
        else:
            print(f"Data is {days_since_observation} days old (beyond {transfer_log_grace_days} day grace period)")
            print("Will download whatever is on archive without retry logic")

        # State tracking for retry logic (per telescope)
        telescope_states = {}
        for idx in range(len(tel)):
            telescope_states[idx] = {
                'success': False,
                'attempts': [],  # List of (attempt_num, archive_count, downloaded_count, expected_count)
                'telname': None,
                'final_message': None
            }

        # Retry loop - skip entirely for old data
        max_attempts = max_retries if is_recent else 1

        for attempt in range(max_attempts):
            attempt_start_time = dt.datetime.now()

            if attempt > 0:
                print(f"\n{'='*80}")
                print(f"RETRY ATTEMPT {attempt + 1}/{max_attempts}")
                print(f"Started at {attempt_start_time.strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"{'='*80}\n")
            else:
                print(f"\n{'='*80}")
                print(f"INITIAL ATTEMPT")
                print(f"Started at {attempt_start_time.strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"{'='*80}\n")

            all_telescopes_successful = True

            # Process each telescope
            for t in range(len(tel)):
                # Skip if this telescope already succeeded
                if telescope_states[t]['success']:
                    print(f"Skipping {telname[t]} - already successfully downloaded")
                    continue

                num_i_eso, num_i_transfer, num_i, count = 0, 0, 0, 0
                ok, found = False, False

                imgdir = f"{dir}/{telroute[t]}/images/{date}"
                telescope_states[t]['telname'] = telname[t]

                print(f"Processing {telname[t]} on {date}")

                # Get expected count from transfer log
                num_i_transfer, found = num_im_transfer(telroute[t], date, dir)
                print(f"Transfer log: {num_i_transfer} images (found={found})")

                # Create directory if needed
                if not os.path.exists(imgdir):
                    os.makedirs(imgdir)
                else:
                    print(f"{imgdir} already exists")

                # Get already downloaded files
                downloaded_files = []
                if os.path.exists(imgdir) and os.listdir(imgdir):
                    print(f"Images already exist for {telname[t]} on {date}")
                    downloaded_files = [a.replace(".fits", "") for a in os.listdir(imgdir)
                                        if a.endswith(('.fits', '.fts'))]

                    # Apply telescope-specific filtering
                    if telname[t] == 'Callisto':
                        if (int(date) > 20220509 and int(date) < 20230317) or int(date) > 20250225:
                            downloaded_files = transformation_check(imgdir, downloaded_files, "SPIRIT")
                        else:
                            downloaded_files = transformation_check(imgdir, downloaded_files, "ANDOR")
                    elif telname[t] == 'Ganymede' and int(date) > 20250303:
                        downloaded_files = transformation_check(imgdir, downloaded_files, "ANDOR")

                print(f"Already downloaded: {len(downloaded_files)} files")

                # Download images
                num_i, num_i_eso = get_images(tel[t], telname[t], stime, etime, date, imgdir,
                                              downloaded_files, test_mode, verbose, query_only, skip_transform)

                # Calculate totals
                total_downloaded = len(downloaded_files) + num_i

                # Validation
                ok_transfer = check_download_images(num_i_eso, num_i)
                ok_num, num_i_dir = check_num_images(imgdir, num_i_eso)
                ok = ok_transfer and ok_num

                if not ok:
                    num_i = num_i_dir
                    total_downloaded = num_i_dir
                    if not ok_transfer:
                        print("WARNING: Issue with number of images on ESO Archive")
                    if not ok_num:
                        print("WARNING: Issue with number of downloaded images on server")

                # Record this attempt
                telescope_states[t]['attempts'].append({
                    'attempt_num': attempt + 1,
                    'archive_count': num_i_eso,
                    'downloaded_count': total_downloaded,
                    'expected_count': num_i_transfer if found else None,
                    'newly_downloaded': num_i
                })

                # Determine success for this telescope
                if is_recent:
                    # For recent data, check against transfer log
                    if not found or num_i_transfer == 0:
                        # Transfer log not updated or shows 0 images - need to retry
                        print(f"Transfer log not ready for {telname[t]} - will retry")
                        telescope_states[t]['success'] = False
                        all_telescopes_successful = False
                    elif total_downloaded == num_i_transfer:
                        # Perfect match!
                        print(f"SUCCESS: {telname[t]} has all {num_i_transfer} expected images")
                        telescope_states[t]['success'] = True
                    elif num_i_eso < num_i_transfer:
                        # Archive doesn't have everything yet
                        print(f"Archive incomplete: {num_i_eso}/{num_i_transfer} for {telname[t]} - will retry")
                        telescope_states[t]['success'] = False
                        all_telescopes_successful = False
                    elif total_downloaded < num_i_transfer:
                        # We're missing some downloads
                        print(f"Download incomplete: {total_downloaded}/{num_i_transfer} for {telname[t]} - will retry")
                        telescope_states[t]['success'] = False
                        all_telescopes_successful = False
                    else:
                        # Got everything we expected
                        print(f"SUCCESS: {telname[t]} complete")
                        telescope_states[t]['success'] = True
                else:
                    # For old data, just accept what we got
                    if total_downloaded == num_i_eso:
                        print(f"Downloaded all {num_i_eso} images available on archive for {telname[t]}")
                        telescope_states[t]['success'] = True
                    else:
                        print(f"Partial download: {total_downloaded}/{num_i_eso} for {telname[t]}")
                        # Still mark as success for old data (no retries)
                        telescope_states[t]['success'] = True

                # Write to log
                write_to_log(total_downloaded, num_i_transfer, num_i_eso, dir, telname[t], telroute[t], str(date))

                # Cleanup zipped files
                print("Deleting remaining zipped files")
                zip_files = glob.glob(f"{imgdir}/*ts.Z")
                for zip_file in zip_files:
                    try:
                        os.remove(zip_file)
                    except:
                        pass

            # Check if we should retry
            if all_telescopes_successful or not is_recent:
                print(f"\n{'='*80}")
                print("All telescopes complete - exiting retry loop")
                print(f"{'='*80}\n")
                break

            # If not the last attempt and not all successful, wait before retry
            if attempt < max_attempts - 1 and not all_telescopes_successful:
                wait_minutes = retry_intervals[attempt] if attempt < len(retry_intervals) else 60
                print(f"\n{'='*80}")
                print(f"Not all telescopes complete. Waiting {wait_minutes} minutes before retry...")
                print(f"Next attempt will be {attempt + 2}/{max_attempts}")
                print(f"{'='*80}\n")
                time.sleep(wait_minutes * 60)

        # After all attempts, send email summary (skip in query-only mode)
        if not query_only:
            send_final_email_summary(telescope_states, date, telname, test_mode, is_recent,
                                    transfer_log_grace_days, max_attempts)


def sendemail(from_addr, to_addr_list, cc_addr_list, subject, message, login, password, smtpserver):
    """Send email using SMTP with proper error handling (same as email_eon.py)."""
    import logging
    logger = logging.getLogger(__name__)

    try:
        # Build email header
        header = f'From: {from_addr}\n'
        header += f'To: {",".join(to_addr_list)}\n'
        if cc_addr_list:
            header += f'Cc: {",".join(cc_addr_list)}\n'
        header += f'Subject: {subject}\n'

        full_message = header + "\n" + str(message)

        # Combine to and cc lists for actual sending
        all_recipients = to_addr_list + cc_addr_list

        # Connect and send
        server = smtplib.SMTP(smtpserver)
        server.starttls()
        server.login(login, password)

        problems = server.sendmail(from_addr, all_recipients, full_message)
        server.quit()

        if problems:
            print(f"Some email recipients failed: {problems}")
        else:
            print("Email sent successfully")

    except smtplib.SMTPAuthenticationError:
        print("SMTP Authentication failed. Check your email address and password.")
        print("For Gmail, you may need to use an App Password instead of your regular password.")
        raise
    except smtplib.SMTPRecipientsRefused as e:
        print(f"All recipients were refused: {e}")
        raise
    except smtplib.SMTPException as e:
        print(f"SMTP error occurred: {e}")
        raise
    except Exception as e:
        print(f"Unexpected error sending email: {e}")
        raise


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Download SPECULOOS data from ESO archive',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download one night for Io
  python SSO_download.py --dir /data/SPECULOOSPipeline --telescope Io --sdate 20250728 --edate 20250729

  # Test mode (download only first image)
  python SSO_download.py --dir /data/SPECULOOSPipeline --telescope Io --sdate 20250728 --edate 20250729 --test

  # Query only (list files without downloading)
  python SSO_download.py --dir /data/SPECULOOSPipeline --telescope Io --sdate 20250728 --edate 20250729 --query-only

  # Verbose mode
  python SSO_download.py --dir /data/SPECULOOSPipeline --telescope Io --sdate 20250728 --edate 20250729 --verbose
  
  # Custom retry settings
  python SSO_download.py --dir /data/SPECULOOSPipeline --telescope Io --sdate 20250728 --edate 20250729 --transfer-log-grace-days 7 --max-retries 10
        """
    )

    # Required arguments
    parser.add_argument('--dir', required=True,
                        help='Base directory for SPECULOOS data (e.g., /data/SPECULOOSPipeline)')
    parser.add_argument('--telescope', required=True,
                        choices=['Io', 'Europa', 'Ganymede', 'Callisto', 'all'],
                        help='Telescope name or "all" for all telescopes')
    parser.add_argument('--sdate', required=True,
                        help='Start date in YYYYMMDD format (e.g., 20250728)')
    parser.add_argument('--edate', required=True,
                        help='End date in YYYYMMDD format (e.g., 20250729)')

    # Optional arguments
    parser.add_argument('--wait', action='store_true',
                        help='Wait for images to be transferred from Chile (for recent data)')
    parser.add_argument('--test', action='store_true',
                        help='Test mode: download only the first image found')
    parser.add_argument('--query-only', action='store_true',
                        help='Query only: list matching files without downloading')
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose mode: show detailed download progress')
    parser.add_argument('--skip-transform', action='store_true',
                        help='Skip astra_transform after downloading (download raw files only)')
    parser.add_argument('--transfer-log-grace-days', type=int, default=5,
                        help='Number of days after observation to check transfer log and retry (default: 5)')
    parser.add_argument('--max-retries', type=int, default=13,
                        help='Maximum number of download attempts for recent data (default: 13, ~12 hours)')

    args = parser.parse_args()

    # Print all arguments for debugging
    print("Arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")
    print()

    # Convert telescope argument
    telescope = None if args.telescope == 'all' else args.telescope

    # Mode indicators
    if args.test:
        print("*** TEST MODE ENABLED - Will only download first image ***")
    if args.query_only:
        print("*** QUERY ONLY MODE - Will list files but not download ***")
    if args.verbose:
        print("*** VERBOSE MODE ENABLED ***")

    if args.skip_transform:
        print("*** SKIP TRANSFORM MODE - astra_transform will not be run ***")

    print(dt.datetime.now())
    eso_download(
        dir=args.dir,
        telname=telescope,
        sdate=args.sdate,
        edate=args.edate,
        wait=str(args.wait).lower(),
        test_mode=args.test,
        query_only=args.query_only,
        verbose=args.verbose,
        transfer_log_grace_days=args.transfer_log_grace_days,
        max_retries=args.max_retries,
        skip_transform=args.skip_transform
    )