import sys
import os
from astropy.io import ascii
from astropy.table import Table,Column
# from termcolor import colored
import numpy as np
from tabulate import tabulate

def update_file(fname,date,target,tel,run,col,val):
    update = False
    # read in the existing report file
    report = ascii.read(fname,format='fixed_width')
    # report.show_in_browser(jsviewer=True)
    hdr = report.colnames
    test_row = ['--------', '--------------------', '---------', '--------------------', '---', '---', '---', '---',
                '---', '---', '---', '---', '---', '----', '----', '----', '----', '----']

    # test_row = ['--------', '--------------------', '---------', '--------------------', '----', '----', '----', '----',
    #             '--------', '--------', '------', '------', '------', '------', '------', '------', '------', '------',
    #             '------']
    qual = thresholds(hdr.index(col)-4, val)
    # colnum = hdr.index(col)
    if qual == 'g':
        q = "[x]"
    elif qual == 'y':
        q = "[-]"
    elif qual == 'r':
        q = "[o]"

    # for each line in the report file:
    for i in range(len(report)):
        # if this date/target/run has already been run then update this line to the report file
        if str(report[i]['Night'])==date and str(report[i]['Target'])==target and \
            str(report[i]['Telescope'])==tel and str(report[i]['Runname'])==run:
            update = True
            report[i][col] = q
            row = report[i]
            # for j in range(colnum+1,len(hdr)):
            #     report[i][hdr[j]] = test_row[j]

            ascii.write(report, fname, format='fixed_width',overwrite=True)

    # if this is a new date/target/run then write a new line to the report file
    if update==False:
        row = [date, target, tel, run, '---','---','---','---','---','---','---','---','---','----',
           '----','----','----','----']
        colnum = hdr.index(col)
        row[colnum] = q
        report.add_row(row)
        ascii.write(report, fname, format='fixed_width', overwrite=True)


def thresholds(col_i,val):
    hdr = ['Night', 'Target', 'Telescope', 'Runname', "QC1", 'QC2', "QC3", "QC4", "QC5", "QC6", "QC7", "QC8", "QC9",
           "QC10", "QC11", "QC12", "QC13", "QC14"]

    # column headers:
    # QC1 = '#science'
    # QC2='#bias'
    # QC3= '#dark'
    # QC4='#flat'
    # QC5='flat flux'
    # QC6='RON'
    # QC7='Dark Current'
    # QC8 ='Overscan'
    # QC9='astrom/wcs?'
    # QC10='#sources Cat'
    # QC11='% Gaia Xmatch'
    # QC12='PSF'
    # QC13='FWHM'
    # QC14='Bkg'

    thresh_op = ['LT','LT','LT','LT','GT','GT','GT','GT','LT','LT','LT','LT','GT','GT']
    thresh_yellow = [10,3,3,3,1000000,8,1,305,10,50,30,0,8,500]
    thresh_red = [1,1,1,1,1000000,100,100,1000,1,5,1,0,100,2000]

    q = 'g'

    if thresh_op[col_i] == 'LT':
        if int(val) < thresh_yellow[col_i]:
            q='y'
        if int(val) < thresh_red[col_i]:
            q='r'
    if thresh_op[col_i] == 'GT':
        if int(val) > thresh_yellow[col_i]:
            q='y'
        if int(val) > thresh_red[col_i]:
            q='r'

    return q

# def print_row(hdr,testr,r):
#
#     for j in range(len(hdr)):
#         hdr[j] = hdr[j].ljust(len(testr[j]))
#
#     print "\t".join(hdr)
#     # print "\t".join(testr)
#     # sys.stdout.write("\t".join(r[0:3]) + "\t")
#     for i in range(len(hdr)):
#         if r[i] == '[x]':
#             sys.stdout.write(colored(r[i],'green') + "\t")
#         elif r[i] == '[-]':
#             sys.stdout.write(colored(r[i],'yellow') + "\t")
#         elif r[i] == '[o]':
#             sys.stdout.write(colored(r[i],'red') + "\t")
#         else:
#             sys.stdout.write(r[i].ljust(len(testr[i])) + "\t")
#
#     sys.stdout.flush()

def create_file(fname,date,target,tel,run,col,val):
    # if the report file does not exist then initialise it
    row = [date,target,tel,run, '---','---','---','---','---','---','---','---','---','----',
           '----','----','----','----']
    # create a test row to set the width of the columns when reading in the table
    test_row = ['--------','--------------------','---------','--------------------','---','---','---','---',
                '---','---','---','---','---','----','----','----','----','----']
    # for all tasks up to and including this one replace [o] with [x]

    hdr = ['Night','Target','Telescope','Runname',"QC1",'QC2',"QC3","QC4","QC5","QC6","QC7","QC8","QC9","QC10","QC11",
           "QC12","QC13","QC14"]

    colnum = hdr.index(col)

    qual = thresholds(colnum-4, val)
    if qual == 'g':
        q = "[x]"
    elif qual == 'y':
        q = "[-]"
    elif qual == 'r':
        q = "[o]"

    row[colnum]=q
    # column data types:
    dtypes = ['S8','S20','S9','S20','S3','S3','S3','S3','S3','S3','S3','S3','S3','S4','S4','S4','S4','S4']

    # create the table, add the test row and first row to the table and write to report file
    report = Table(names=hdr,dtype=dtypes)
    report.add_row(test_row)
    report.add_row(row)
    ascii.write(report,fname,format='fixed_width')

    # print_row(hdr,test_row,row)


def main(args):
    fname = args[1]
    date = args[2]
    target = args[3]
    tel = args[4]
    run = args[5]
    col = args[6]
    val = args[7]

    # if the report file does not exist then create it, otherwise update the file
    if os.path.exists(fname):
        update_file(fname,date,target,tel,run,col,val)
    else:
        create_file(fname,date,target,tel,run,col,val)

if __name__ == '__main__':
    main(sys.argv)