"""The first step is to create an SMTP object, each object is used for connection
with one server."""
import smtplib
import argparse
import os
import sys
import logging
from astropy.io import ascii
from dotenv import load_dotenv
load_dotenv()  # This loads .env automatically

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_email_config():
    """Load email configuration from environment variables."""
    config = {
        'from_addr': os.getenv('EMAIL_ADDRESS'),
        'password': os.getenv('EMAIL_PASSWORD'),
        'smtp_server': os.getenv('SMTP_SERVER', 'smtp.gmail.com:587'),
        'to_addr_list': os.getenv('TO_EMAIL_LIST', '').split(',') if os.getenv('TO_EMAIL_LIST') else [],
        'cc_addr_list': os.getenv('CC_EMAIL_LIST', '').split(',') if os.getenv('CC_EMAIL_LIST') else []
    }

    # Validate required configuration
    missing_vars = []
    if not config['from_addr']:
        missing_vars.append('EMAIL_ADDRESS')
    if not config['password']:
        missing_vars.append('EMAIL_PASSWORD')
    if not config['to_addr_list'] or config['to_addr_list'] == ['']:
        missing_vars.append('TO_EMAIL_LIST')

    if missing_vars:
        logger.error(f"Missing required environment variables: {', '.join(missing_vars)}")
        logger.error("Please set the following environment variables:")
        logger.error("  EMAIL_ADDRESS - Your email address")
        logger.error("  EMAIL_PASSWORD - Your email password/app password")
        logger.error("  TO_EMAIL_LIST - Comma-separated list of recipient emails")
        logger.error("Optional:")
        logger.error("  SMTP_SERVER - SMTP server (default: smtp.gmail.com:587)")
        logger.error("  CC_EMAIL_LIST - Comma-separated list of CC emails")
        sys.exit(1)

    # Clean up email lists (remove empty strings and whitespace)
    config['to_addr_list'] = [email.strip() for email in config['to_addr_list'] if email.strip()]
    config['cc_addr_list'] = [email.strip() for email in config['cc_addr_list'] if email.strip()]

    return config

def sendemail(from_addr, to_addr_list, cc_addr_list, subject, message, login, password, smtpserver):
    """Send email using SMTP with proper error handling."""
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
        logger.info(f"Connecting to SMTP server: {smtpserver}")
        server = smtplib.SMTP(smtpserver)
        server.starttls()
        server.login(login, password)

        logger.info(f"Sending email to: {', '.join(all_recipients)}")
        problems = server.sendmail(from_addr, all_recipients, full_message)
        server.quit()

        if problems:
            logger.warning(f"Some recipients failed: {problems}")
        else:
            logger.info("Email sent successfully")

    except smtplib.SMTPAuthenticationError:
        logger.error("SMTP Authentication failed. Check your email address and password.")
        logger.error("For Gmail, you may need to use an App Password instead of your regular password.")
        raise
    except smtplib.SMTPRecipientsRefused as e:
        logger.error(f"All recipients were refused: {e}")
        raise
    except smtplib.SMTPException as e:
        logger.error(f"SMTP error occurred: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error sending email: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Send pipeline report emails')
    parser.add_argument('date', help='Date of the observation')
    parser.add_argument('-t', '--target', required=True, help='Target name')
    parser.add_argument('-r', '--runname', required=True, help='Run name')
    parser.add_argument('-tel', '--telescope', required=True, help='Telescope name')
    parser.add_argument('-rep', '--report', required=True, help='Path to report file')
    parser.add_argument('-e', '--end', required=True, choices=['0', '1'],
                       help='End flag: 1 for pipeline report, 0 for new stack notification')
    args = parser.parse_args()

    # Load email configuration
    config = get_email_config()

    # Extract arguments
    date = args.date
    targ = args.target
    report = args.report
    run = args.runname
    tel = args.telescope
    end = args.end

    try:
        report_table = ascii.read(report, format='fixed_width')
        logger.info(f"Successfully loaded report file: {report}")
    except Exception as e:
        logger.error(f"Failed to read report file {report}: {e}")
        sys.exit(1)

    message = ''

    if end == "1":
        # Pipeline completion report
        subject = f'Pipeline Report: {date} {targ}'

        # Search for matching entry in report table
        for i in range(len(report_table)):
            if (str(report_table[i]['Night']) == date and
                str(report_table[i]['Target']) == targ and
                str(report_table[i]['Telescope']) == tel and
                str(report_table[i]['Runname']) == run):
                message = str(report_table[i])
                break

        if message == '':
            message = f"WARNING: no result found in {tel} report for run {run} for {targ} on {date}"
            logger.warning(message)
    else:
        # New stack notification
        message = f"New Stack for {targ} on {date}"
        subject = 'New Stack'

    try:
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
    except Exception as e:
        logger.error(f"Failed to send email: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()