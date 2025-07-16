import sys
import os
from astropy.io import ascii
from astropy.table import Table,Column
from tabulate import tabulate

def update_file(fname,date,target,tel,run,t):
    update = False
    # read in the existing report file
    report = ascii.read(fname,format='fixed_width')

    # for each line in the report file:
    for i in range(len(report)):
        # if this date/target/run has already been run then update this line to the report file
        if str(report[i]['Night'])==date and str(report[i]['Target'])==target and \
            str(report[i]['Telescope'])==tel and str(report[i]['Runname'])==run:
            # print "Update Report"
            update = True

            # reset the line to T*=[o]
            for j in range(1,11):
                report[i]['T'+str(j)] = '[o]'
            # for all tasks up to and including this one replace [o] with [x]
            for j in range(1,int(t)+1):
                report[i]['T'+str(j)] = '[x]'

            ascii.write(report, fname, format='fixed_width',overwrite=True)

    # if this is a new date/target/run then write a new line to the report file
    if update==False:
        # print "Write new line to report"
        row = [date, target, tel, run, '[o]', '[o]', '[o]', '[o]', '[o]', '[o]', '[o]', '[o]', '[o]', '[o]']
        # for all tasks up to and including this one replace [o] with [x]
        for k in range(int(t)):
            row[k + 4] = '[x]'
        report.add_row(row)
        ascii.write(report, fname, format='fixed_width', overwrite=True)

def create_file(fname,date,target,tel,run,t):
    # if the report file does not exist then initialise it
    row = [date,target,tel,run,'[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]']
    # create a test row to set the width of the columns when reading in the table
    test_row = ['--------','--------------------','---------','--------------------','---','---','---','---','---','---','---','---','---','---']
    # for all tasks up to and including this one replace [o] with [x]
    for i in range(int(t)):
        row[i+4]='[x]'
    # column headers:
    hdr = ['Night','Target','Telescope','Runname','T1','T2','T3','T4','T5','T6','T7','T8','T9','T10']
    # column data types:
    dtypes = ['S8','S20','S9','S20','S3','S3','S3','S3','S3','S3','S3','S3','S3','S3']

    # create the table, add the test row and first row to the table and write to report file
    report = Table(names=hdr,dtype=dtypes)
    report.add_row(test_row)
    report.add_row(row)
    ascii.write(report,fname,format='fixed_width')


def main(args):
    fname = args[1]
    date = args[2]
    target = args[3]
    tel = args[4]
    run = args[5]
    t = args[6]

    # if the report file does not exist then create it, otherwise update the file
    if os.path.exists(fname):
        update_file(fname,date,target,tel,run,t)
    else:
        create_file(fname,date,target,tel,run,t)

if __name__ == '__main__':
    main(sys.argv)