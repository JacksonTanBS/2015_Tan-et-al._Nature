#!/usr/bin/python3

# This script automates the FTP download of ISCCP data.

import os
import sys
from ftplib import FTP

outpath = '/media/Sentinel/data/ISCCP/D1/'
year = int(sys.argv[1])

# establish a connection each year because server seems prone to timeouts
if not os.path.exists('%s%4d' % (outpath, year)):
    os.system('mkdir -p %s%4d' % (outpath, year))

if year == 1983: months = range(7, 13)
else: months = range(1, 13)

for month in months:

    while True:   # this will repeat the month if TimeoutError or OSError (both from connection issues)

        try:

            ftp = FTP('eclipse.ncdc.noaa.gov')
            ftp.login()
            ftp.cwd('/pub/isccp/d1/%4d/' % year)
            files = ftp.nlst('ISCCP.D1.*.GLOBAL.%4d.%02d.*.GPC' % (year, month))

            # SFTP in the file
            for ff in files:
                with open('%s%4d/%s' % (outpath, year, ff), 'wb') as f:
                    ftp.retrbinary('RETR %s' % ff, f.write)

            ftp.quit()

        except (TimeoutError, OSError):

            print('Timeout error for %4d %02d. Repeating...' % (year, month))
            ftp.quit()
            continue

        break
