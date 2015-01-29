#!/usr/bin/python3

# This script automates the FTP download of TRMM 3A25 data.

import os
import gzip
from ftplib import FTP

outpath = '/media/Sentinel/data/TRMM/3A25/'
pps_login = ''    # put your login here
pps_psswd = ''    # put your password here

if not os.path.exists(outpath): os.system('mkdir -p %s' % outpath)

ftp = FTP('arthurhou.pps.eosdis.nasa.gov')
ftp.login(pps_login, pps_psswd)

for year in range(1998, 2010):
    for month in range(1, 13):

        ftp.cwd('/trmmdata/ByDate/V07/%4d/%02d/01/' % (year, month))
        filename = '3A25.%4d%02d01.7' % (year, month)

        # SFTP in the file
        with open('%s%s.HDF.gz' % (outpath, filename), 'wb') as f:
            ftp.retrbinary('RETR %s.HDF.gz' % filename, f.write)

        # decompress the file
        with open('%s%s.HDF' % (outpath, filename), 'wb') as f:
            with gzip.open('%s%s.HDF.gz' % (outpath, filename), 'r') as g:
                f.write(g.read())
        os.remove('%s%s.HDF.gz' % (outpath, filename))

ftp.quit()
