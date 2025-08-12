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

def get_images(tel, telname, stime, etime, date, imgdir, down_im, test_mode=False, verbose=False, query_only=False):
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

            # Query for files using modern TAP service
            print(f"Querying for {data_type} images...")
            files = downloader.query_files(tel, stime, etime, data_type)

            t1 = timeit.default_timer()
            elapsed = t1 - t0
            print(f"Time taken to query archive for {data_type} images: {elapsed / 60:.3f} minutes")

            if not files:
                print(f"WARNING: no {data_type} results found for {telname} on {date}")
                continue
            elif len(files) >= 10000:
                print(f"WARNING: returned maximum number (10000) of {data_type} images for {telname} on {date}")

            num_files = len(files)
            print(f"{num_files} {data_type} images for {telname} on {date}")
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

                        # Skip astra_transform in query_only mode
                        if not query_only:
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

                    # Skip astra_transform in query_only mode
                    if not query_only:
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


def eso_download(dir, telname, sdate, edate, wait, test_mode=False, query_only=False, verbose=False):
    """Main download function"""

    if telname == None:
        tel = ['60.A-9009(A)', '60.A-9009(B)', '60.A-9009(C)', '60.A-9009(D)']
        telname = ['Io', 'Europa', 'Ganymede', 'Callisto']
        telroute = ['Observations/' + s for s in telname]
    else:
        telname = str(telname)
        if telname == 'Io':
            tel = ['60.A-9009(A)']
        elif telname == 'Europa':
            tel = ['60.A-9009(B)']
        elif telname == 'Ganymede':
            tel = ['60.A-9009(C)']
        elif telname == 'Callisto':
            tel = ['60.A-9009(D)']
        telname = [telname]
        telroute = ['Observations/' + telname[0]]

    s = dt.datetime.strptime(sdate, "%Y%m%d")
    e = dt.datetime.strptime(edate, "%Y%m%d")
    days = (e - s).days

    print(f"Downloading {days} days of data from ESO archive")
    if test_mode:
        print("*** TEST MODE: Will only download the first image found ***")

    # for date in range(int(sdate), int(edate) + 1):
    for d in (s + dt.timedelta(n) for n in range(days)):
        os.system(f"rm -rf {dir}/.astropy/cache/astroquery/Eso/*")

        date = d.strftime("%Y%m%d")
        stime = d.strftime("%Y-%m-%d")
        etime = (d + dt.timedelta(days=1)).strftime("%Y-%m-%d")

        # Only check against transfer log for data from within the past 5 days
        new_data = s > (dt.datetime.now() - dt.timedelta(days=5))
        if new_data:
            if wait.lower() == "true":
                print("Set to wait, skipped if no images.")
            else:
                print("Images less than 5 days old. Compare against transfer log.")

        for t in range(len(tel)):
            num_i_eso, num_i_transfer, num_i, count = 0, 0, 0, 0
            ok, found = False, False

            imgdir = f"{dir}/{telroute[t]}/images/{date}"

            print(f"Downloading images for {telname[t]} on {date}")
            num_i_transfer, found = num_im_transfer(telroute[t], date, dir)
            print(f"{num_i_transfer} images transferred from Chile.")

            while (found == False and count < 20 and new_data) and wait.lower() == "true":
                now = dt.datetime.now()
                current_time = now.strftime("%H:%M:%S")
                print(f"Current Time = {current_time}")
                print("Images have not yet been transferred. Waiting 30 minutes...")
                # wait for 30 minutes and then check if the images have been transferred
                time.sleep(1800)
                count = count + 1
                print("Checking number of images transferred to ESO archive...")
                num_i_transfer, found = num_im_transfer(telroute[t], date, dir)

            if found == False and wait.lower() == "true":
                print("ESO download timed out.")
            else:
                count = 0
                print(tel, stime, etime)
                # Note: num_i_eso will be determined by get_images function
                print("Querying ESO archive for available images...")

                while count < 20 and new_data and wait.lower() == "true":
                    # For now, we'll get the count from get_images
                    break

                if count < 20:
                    # Create directory if needed
                    if not os.path.exists(imgdir):
                        os.makedirs(imgdir)
                    else:
                        print(f"{imgdir} already exists!")

                    # Get already downloaded files (works for empty directories too)
                    downloaded_files = []
                    if os.path.exists(imgdir) and os.listdir(imgdir):
                        print(f"WARNING: Images already exist for {telname[t]} on {date}.")
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
                    else:
                        print("Folder exists, but there are no images!" if os.path.exists(
                            imgdir) else "Created new directory")

                    print("Already downloaded files:", len(downloaded_files))

                    # Single call to get_images
                    num_i, num_i_eso = get_images(tel[t], telname[t], stime, etime, date, imgdir,
                                                  downloaded_files, test_mode, verbose, query_only)

                    # Single validation block
                    ok_transfer = check_download_images(num_i_eso, num_i)
                    ok_num, num_i_dir = check_num_images(imgdir, num_i_eso)
                    ok = ok_transfer and ok_num

                    if not ok:
                        num_i = num_i_dir
                        if not ok_transfer:
                            print("WARNING: Issue with number of images on ESO Archive.")
                        if not ok_num:
                            print("WARNING: Issue with number of downloaded images on server.")
                else:
                    print("WARNING: Different number of images transferred as on the ESO archive.")

            print(telname)
            write_to_log(num_i, num_i_transfer, num_i_eso, dir, telname[t], telroute[t], str(date))

            # DELETE ALL ZIPPED FILES (legacy cleanup)
            print("Delete all remaining zipped files.")
            zip_files = glob.glob(f"{imgdir}/*ts.Z")
            for zip_file in zip_files:
                try:
                    os.remove(zip_file)
                except:
                    pass

            # Don't send emails in query-only mode
            if query_only:
                print("Query-only mode: skipping email notification")
            else:
                # Determine email content based on actual status
                total_on_archive = num_i_eso
                newly_downloaded = num_i  # Use the return value from get_images
                already_had = len(downloaded_files) if downloaded_files else 0
                total_local = newly_downloaded + already_had

                if test_mode:
                    # Test mode specific messaging
                    remaining = total_on_archive - already_had - newly_downloaded
                    if newly_downloaded > 0:
                        subject = f"ESO Test Download SUCCESS: {newly_downloaded} test image downloaded"
                        message = f"""
                        Test download completed for {', '.join(telname)} on {date}:
                        - Total images available on ESO archive: {total_on_archive}
                        - Already downloaded: {already_had}
                        - Test download: {newly_downloaded} image
                        - Remaining to download: {remaining}

                        Status: TEST SUCCESSFUL - Ready for full download

                        To download all remaining images, run without --test flag:
                        python download/SSO_download.py --dir [DIR] --telescope {telname[0] if telname else 'TARGET'} --sdate {date[:4]}{date[4:6]}{date[6:8]} --edate {date[:4]}{date[4:6]}{date[6:8]}
                        """
                    else:
                        subject = f"ESO Test Download: No new images needed"
                        message = f"""
                        Test download check for {', '.join(telname)} on {date}:
                        - Total images available on ESO archive: {total_on_archive}
                        - Already downloaded: {already_had}
                        - Test result: All images already present locally

                        Status: TEST COMPLETE - No download needed
                        """

                elif newly_downloaded == 0 and total_local == total_on_archive:
                    # All files were already downloaded
                    subject = f"ESO Archive Status: All {total_on_archive} images already present"
                    message = f"""
                    Archive check for {', '.join(telname)} on {date}:
                    - Total images on ESO archive: {total_on_archive}
                    - Already downloaded: {already_had}
                    - Newly downloaded: 0
                    Status: Complete (no action needed)
                    """

                elif newly_downloaded > 0 and total_local == total_on_archive:
                    # Successfully downloaded missing files
                    subject = f"ESO Download Complete: {newly_downloaded} new images"
                    message = f"""
                    Download completed for {', '.join(telname)} on {date}:
                    - Total images on ESO archive: {total_on_archive}
                    - Already had: {already_had}
                    - Newly downloaded: {newly_downloaded}
                    Status: Complete
                    """

                elif total_local < total_on_archive:
                    # Partial download - something went wrong
                    missing = total_on_archive - total_local
                    subject = f"ESO Download INCOMPLETE: {missing} images missing"
                    message = f"""
                    Download INCOMPLETE for {', '.join(telname)} on {date}:
                    - Total images on ESO archive: {total_on_archive}
                    - Already had: {already_had}
                    - Newly downloaded: {newly_downloaded}
                    - Still missing: {missing}
                    Status: REQUIRES ATTENTION

                    Possible causes:
                    - Network timeouts during download
                    - Authentication issues  
                    - Disk space problems
                    - File corruption during transfer

                    To retry missing files, run the same command again.
                    """

                else:
                    # Edge case - more local files than on archive?
                    subject = f"ESO Archive Anomaly: Local file count mismatch"
                    message = f"""
                    Unusual situation for {', '.join(telname)} on {date}:
                    - Total images on ESO archive: {total_on_archive}
                    - Total local files: {total_local}
                    Status: REVIEW NEEDED
                    """

                # Send email for meaningful events
                should_send_email = (
                        test_mode or  # Always send test results
                        newly_downloaded > 0 or  # New downloads happened
                        total_local != total_on_archive  # Mismatch that needs attention
                )

                if should_send_email:
                    try:
                        config = get_email_config()
                        if config['to_addr_list']:
                            sendemail(
                                from_addr=config['from_addr'],
                                to_addr_list=config['to_addr_list'],
                                cc_addr_list=config['cc_addr_list'],
                                subject=subject,
                                message=message.strip(),
                                login=config['from_addr'],
                                password=config['password'],
                                smtpserver=config['smtp_server']
                            )
                            print("Email notification sent")
                        else:
                            print("No email recipients configured, skipping email notification")
                    except Exception as e:
                        print(f"Email sending failed: {e}")
                else:
                    print("No significant changes to report, skipping email notification")


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

    print(dt.datetime.now())
    eso_download(
        dir=args.dir,
        telname=telescope,
        sdate=args.sdate,
        edate=args.edate,
        wait=str(args.wait).lower(),
        test_mode=args.test,
        query_only=args.query_only,
        verbose=args.verbose
    )