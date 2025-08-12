#!/usr/bin/env python3
"""
Modern ESO Archive downloader using Python 3, OAuth2.0 authentication and DataLink protocol.
Replaces the old Python 2 downloadportalapi/store method with proper VO standards.
"""

import requests
import json
import time
import os
from datetime import datetime, timedelta


class ESODownloader:
    """
    Modern ESO Archive downloader using OAuth2.0 authentication and DataLink protocol.
    Replaces the old downloadportalapi/store method with proper VO standards.
    """

    def __init__(self, username, password):
        self.username = username
        self.password = password
        self.token = None
        self.token_expires = 0
        self.session = requests.Session()

        # Base URLs for ESO services
        self.oauth_url = "https://www.eso.org/sso/oidc/token"
        self.tap_url = "https://archive.eso.org/tap_obs"
        self.datalink_base = "https://archive.eso.org/datalink/links"

    def ensure_valid_token(self):
        """Get a fresh token if needed (valid for 8 hours)"""
        current_time = time.time()

        # Get new token if we don't have one or it expires in < 30 minutes
        if not self.token or (current_time + 1800) > self.token_expires:
            print("Getting fresh OAuth2.0 token...")
            self.get_new_token()
        else:
            remaining_hours = (self.token_expires - current_time) / 3600
            print("Using existing token (expires in {:.1f} hours)".format(remaining_hours))

    def get_new_token(self):
        """Get fresh OAuth2.0 token from ESO using their official method"""
        try:
            # Based on ESO's official eso_programmatic.py getToken() function
            response = requests.get(
                self.oauth_url,
                params={
                    "response_type": "id_token token",
                    "grant_type": "password",
                    "client_id": "clientid",
                    "username": self.username,
                    "password": self.password
                },
                verify=False,  # For compatibility with your existing SSL setup
                timeout=30
            )

            if response.status_code == 200:
                token_data = response.json()
                # ESO uses 'id_token' not 'access_token'
                self.token = token_data['id_token']
                # ESO tokens are valid for 8 hours
                self.token_expires = time.time() + 28800  # 8 hours
                print("Successfully obtained OAuth2.0 token")
                return True
            else:
                print(f"Token request failed with status {response.status_code}: {response.text}")
                return False

        except Exception as e:
            print(f"Error getting OAuth2.0 token: {str(e)}")
            return False

    def query_files(self, prog_id, start_date, end_date, dp_cat='SCIENCE'):
        """
        Query ESO TAP service for files matching criteria

        Args:
            prog_id: Program ID (e.g., '60.A-9009(C)')
            start_date: Start date in YYYY-MM-DD format (astronomical night begins 15:00 UTC this day)
            end_date: End date in YYYY-MM-DD format (astronomical night ends 15:00 UTC this day)
            dp_cat: Data category ('SCIENCE' or 'CALIB')

        Returns:
            List of file metadata dictionaries
        """
        self.ensure_valid_token()

        # Convert to astronomical night boundaries (15:00 UTC to 15:00 UTC next day)
        start_mjd = self._date_to_mjd_night_start(start_date)
        end_mjd = self._date_to_mjd_night_end(end_date)

        # ADQL query for raw data - using only basic columns that definitely exist
        query = f"""
        SELECT TOP 10000
            dp_id, 
            object, 
            ra, 
            dec, 
            mjd_obs,
            instrument,
            dp_cat,
            datalink_url
        FROM dbo.raw 
        WHERE prog_id = '{prog_id}' 
            AND dp_cat = '{dp_cat}' 
            AND dp_id LIKE 'SPECU%'
            AND mjd_obs BETWEEN {start_mjd} AND {end_mjd}
        ORDER BY mjd_obs
        """

        try:
            # Submit synchronous TAP query
            tap_params = {
                'REQUEST': 'doQuery',
                'LANG': 'ADQL',
                'FORMAT': 'json',
                'QUERY': query
            }

            headers = {
                'Authorization': f'Bearer {self.token}'
            }

            response = requests.post(
                self.tap_url + '/sync',
                data=tap_params,
                headers=headers,
                verify=False,
                timeout=300
            )

            if response.status_code == 200:
                data = response.json()
                files = data.get('data', [])
                print(
                    f"Found {len(files)} {dp_cat} files for program {prog_id} between {start_date} 15:00 UTC and {end_date} 15:00 UTC")
                return files
            else:
                print(f"TAP query failed with status {response.status_code}: {response.text}")
                return []

        except Exception as e:
            print(f"Error querying TAP service: {str(e)}")
            return []

    def get_download_urls(self, dp_ids):
        """
        Get download URLs for files using direct dataportal URLs

        Args:
            dp_ids: List of ESO data product IDs

        Returns:
            Dictionary mapping dp_id to download URL
        """
        self.ensure_valid_token()

        download_urls = {}

        for dp_id in dp_ids:
            try:
                # Use direct dataportal URL (simpler and more reliable than DataLink)
                download_url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"
                download_urls[dp_id] = download_url
                print(f"Using direct download URL for {dp_id}")

            except Exception as e:
                print(f"Error getting download URL for {dp_id}: {str(e)}")

        return download_urls

    def download_files(self, file_metadata, destination, parallel=False, verbose=False, query_only=False):
        """
        Download files to destination directory

        Args:
            file_metadata: List of file metadata from TAP query
            destination: Directory to save files
            parallel: Whether to download files in parallel (default: False)
            verbose: Whether to show detailed output (default: False)
            query_only: If True, only print file info without downloading (default: False)

        Returns:
            List of successfully downloaded file paths (empty if query_only=True)
        """
        self.ensure_valid_token()

        # Extract dp_ids from metadata
        dp_ids = [row[0] for row in file_metadata]  # dp_id is first column

        if query_only:
            print(f"\n=== QUERY ONLY MODE - Found {len(dp_ids)} files ===")
            for i, row in enumerate(file_metadata, 1):
                dp_id = row[0]
                # Show additional metadata if available
                if len(row) > 1:
                    obj_name = row[1] if len(row) > 1 else "N/A"
                    ra = row[2] if len(row) > 2 else "N/A"
                    dec = row[3] if len(row) > 3 else "N/A"
                    mjd = row[4] if len(row) > 4 else "N/A"
                    instrument = row[5] if len(row) > 5 else "N/A"
                    print(
                        f"{i:3d}. {dp_id} | Object: {obj_name} | RA: {ra} | Dec: {dec} | MJD: {mjd} | Instrument: {instrument}")
                else:
                    print(f"{i:3d}. {dp_id}")
            print(f"=== End of {len(dp_ids)} files ===\n")
            return []  # Return empty list since no files were downloaded

        # Rest of the existing download logic...
        if not os.path.exists(destination):
            os.makedirs(destination)

        # Extract dp_ids from metadata
        dp_ids = [row[0] for row in file_metadata]  # dp_id is first column

        # Get download URLs
        download_urls = self.get_download_urls(dp_ids)

        if not download_urls:
            print("No download URLs obtained")
            return []

        successful_downloads = []

        headers = {
            'Authorization': 'Bearer {}'.format(self.token)
        }

        print("Starting download of {} files...".format(len(download_urls)))

        for i, (dp_id, download_url) in enumerate(download_urls.items()):
            try:
                if verbose:
                    print(f"Downloading {i + 1}/{len(download_urls)}: {dp_id}")

                response = requests.get(
                    download_url,
                    headers=headers,
                    verify=True,  # Enable proper SSL verification
                    timeout=300,
                    stream=True
                )

                if verbose:
                    print(f"Response status: {response.status_code}")
                    print(f"Content-Type: {response.headers.get('Content-Type', 'Unknown')}")
                    print(f"Content-Length: {response.headers.get('Content-Length', 'Unknown')}")

                if response.status_code == 200:
                    content_type = response.headers.get('Content-Type', '')

                    # Check if we got HTML instead of a file (authentication failure)
                    if 'text/html' in content_type:
                        print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Got HTML response instead of file")
                        if verbose:
                            print("First 200 characters of response:")
                            print(response.text[:200])
                        continue

                    # Check content disposition for filename hints
                    content_disposition = response.headers.get('Content-Disposition', '')
                    if verbose:
                        print(f"Content-Disposition: {content_disposition}")

                    # Determine appropriate file extension
                    if 'application/gzip' in content_type or '.gz' in content_disposition:
                        filename = f"{dp_id}.fits.gz"
                    elif 'application/x-compress' in content_type or '.Z' in content_disposition:
                        filename = f"{dp_id}.fits.Z"
                    else:
                        filename = f"{dp_id}.fits"

                    filepath = os.path.join(destination, filename)

                    with open(filepath, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)

                    file_size = os.path.getsize(filepath)
                    if verbose:
                        print(f"Downloaded {filename} ({file_size / 1024 / 1024:.1f} MB)")

                    # Check if file starts with FITS header or compression magic
                    with open(filepath, 'rb') as f:
                        first_bytes = f.read(10)
                        if verbose:
                            print(f"First 10 bytes: {first_bytes}")

                        # Check for compression magic numbers
                        if first_bytes.startswith(b'\x1f\x8b'):  # gzip magic
                            if verbose:
                                print("✓ Gzip compressed file detected")
                            # Rename file to .gz if not already
                            if not filename.endswith('.gz'):
                                new_filepath = filepath + '.gz'
                                os.rename(filepath, new_filepath)
                                filepath = new_filepath
                                filename = filename + '.gz'
                        elif first_bytes.startswith(b'\x1f\x9d'):  # compress (.Z) magic
                            if verbose:
                                print("✓ Compress (.Z) file detected")
                            # Rename file to .Z if not already
                            if not filename.endswith('.Z'):
                                new_filepath = filepath + '.Z'
                                os.rename(filepath, new_filepath)
                                filepath = new_filepath
                                filename = filename + '.Z'
                        elif first_bytes.startswith(b'SIMPLE'):
                            if verbose:
                                print("✓ Valid FITS file detected")
                        else:
                            if verbose:
                                print("⚠ File doesn't appear to be FITS or compressed")
                                print(f"Content preview: {first_bytes}")

                    # Decompress if needed
                    decompression_status = "downloaded"
                    if filename.endswith('.gz'):
                        if verbose:
                            print(f"Decompressing {filename}...")
                        import gzip
                        try:
                            with gzip.open(filepath, 'rb') as f_in:
                                with open(filepath[:-3], 'wb') as f_out:  # Remove .gz extension
                                    f_out.write(f_in.read())
                            os.remove(filepath)  # Remove compressed file
                            final_filepath = filepath[:-3]
                            decompression_status = "decompressed"
                            if verbose:
                                print(f"✓ Decompressed to {os.path.basename(final_filepath)}")
                        except Exception as e:
                            print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Decompression failed - {e}")
                            continue
                    elif filename.endswith('.Z'):
                        if verbose:
                            print(f"Decompressing {filename}...")
                        try:
                            import subprocess
                            # Use system uncompress command for .Z files
                            result = subprocess.run(['uncompress', filepath], capture_output=True, text=True)
                            if result.returncode == 0:
                                final_filepath = filepath[:-2]  # Remove .Z extension
                                decompression_status = "decompressed"
                                if verbose:
                                    print(f"✓ Decompressed to {os.path.basename(final_filepath)}")
                            else:
                                print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Uncompress failed - {result.stderr}")
                                continue
                        except Exception as e:
                            print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Decompression failed - {e}")
                            continue
                    else:
                        final_filepath = filepath

                    # Success message (non-verbose)
                    if not verbose:
                        print(f"✓ {i + 1}/{len(download_urls)} {dp_id}: {decompression_status}")

                    successful_downloads.append(final_filepath)

                else:
                    print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Download failed (status {response.status_code})")
                    if verbose:
                        if response.status_code == 401:
                            print("Authentication may have failed")
                        elif response.status_code == 403:
                            print("Access forbidden - check permissions")

            except Exception as e:
                print(f"✗ {i + 1}/{len(download_urls)} {dp_id}: Error - {str(e)}")

        return successful_downloads

    def _date_to_mjd_night_start(self, date_str):
        """Convert YYYY-MM-DD date to MJD for astronomical night start (15:00 UTC)"""
        dt = datetime.strptime(date_str + ' 15:00:00', '%Y-%m-%d %H:%M:%S')
        # Calculate MJD: days since 1858-11-17 00:00 UTC
        mjd = (dt - datetime(1858, 11, 17)).total_seconds() / 86400.0
        return mjd

    def _date_to_mjd_night_end(self, date_str):
        """Convert YYYY-MM-DD date to MJD for astronomical night end (15:00 UTC same day)"""
        dt = datetime.strptime(date_str + ' 15:00:00', '%Y-%m-%d %H:%M:%S')
        # Don't add a day - the date_str is already the end date
        mjd = (dt - datetime(1858, 11, 17)).total_seconds() / 86400.0
        return mjd


def store(dp_ids, dirname, username, password):
    """
    Compatibility function to replace the old store() function.
    Now uses modern ESO APIs instead of downloadportalapi/store.

    Args:
        dp_ids: List of ESO data product IDs
        dirname: Destination directory
        username: ESO username
        password: ESO password

    Returns:
        List of downloaded file paths
    """
    print("Using modern ESO API for download...")
    print(f"Number of files to download: {len(dp_ids)}")

    try:
        downloader = ESODownloader(username, password)

        # Create fake metadata structure for compatibility
        file_metadata = [[dp_id] for dp_id in dp_ids]

        # Download files
        downloaded_files = downloader.download_files(file_metadata, dirname)

        print(f"Modern ESO download completed: {len(downloaded_files)} files")
        return downloaded_files

    except Exception as e:
        print(f"Modern ESO download failed: {str(e)}")
        print("This may require updating your ESO_archive.py to use the new ESODownloader class directly")
        return []


# Example usage for testing
if __name__ == "__main__":
    import sys

    # Test the new downloader
    if len(sys.argv) > 2:
        username = sys.argv[1]
        password = sys.argv[2]
    else:
        username = input("ESO Username: ")
        password = input("ESO Password: ")

    downloader = ESODownloader(username, password)

    # Query for files
    files = downloader.query_files(
        prog_id="60.A-9009(C)",
        start_date="2025-05-02",
        end_date="2025-05-03",
        dp_cat="SCIENCE"
    )

    # Download files
    if files:
        downloaded = downloader.download_files(files, "/tmp/eso_test")
        print(f"Downloaded {len(downloaded)} files")
    else:
        print("No files found to download")