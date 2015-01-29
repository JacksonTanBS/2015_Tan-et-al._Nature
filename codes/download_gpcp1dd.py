#!/usr/bin/python3

# This script downloads and/or converts the GPCP 1DD V1.2 data from monthly Fortran binary files to yearly netCDF format. You may choose whether to download or to convert the files by changing the flags in the settings below.

# DISCLAIMER: Use this script and its outputs at your own risk. I cannot guarantee that the conversion is flawless.

import os
import sys

# settings
download = True
convert = True
path = '/media/Sentinel/data/GPCP/1DD/'

path_raw = '%sraw/' % path
path_nc = '%snetcdf/' % path
nlat, nlon = 180, 360


if download:

    from urllib.request import urlretrieve

    gpcp_url = 'ftp://rsd.gsfc.nasa.gov/pub/1dd-v1.2/'

    # make directory to downloaded raw files if nonexistent
    if not os.path.exists(path_raw): os.makedirs(path_raw)

    for year in range(1997, 2013):
        for month in range(1, 13):

            file_raw = 'gpcp_1dd_v1.2_p1d.%4d%02d.gz' % (year, month)
            urlretrieve('%s%s' % (gpcp_url, file_raw), '%s%s' % (path_raw, file_raw))

if convert:

    import gzip
    import numpy as np
    from struct import unpack
    from calendar import isleap
    from datetime import datetime, timedelta
    from netCDF4 import Dataset, date2num

    for year in range(1997, 2013):

        nday = 365 + isleap(year)

        P = np.ma.masked_all([nday, nlat, nlon])
        day = 0

        for month in range(1, 13):

            file_raw = 'gpcp_1dd_v1.2_p1d.%4d%02d.gz' % (year, month)
            with gzip.open('%s%s' % (path_raw, file_raw), 'r') as g: raw = g.read()

            # separate the header (and record for troubleshooting purposes)
            head = raw[:1440]
            data = raw[1440:]

            if len(data) % 4: sys.exit('Error: binary data length does not match!')

            converted = np.array(unpack('>' + 'f' * (len(data) // 4), data))
            mday = len(converted) / (nlat * nlon)

            # extract binary file as 32-bit float (REAL*8) in big-endian byte order
            P[day : day + mday] = converted.reshape(-1, nlat, nlon)

            day += mday

        if np.ma.count_masked(P): sys.exit('Error: unknown missing data in %d' % year)

        #--- Commencing the annoyingly verbose setting up of the netCDF file ---#

        if not os.path.exists(path_nc): os.makedirs(path_nc)
        file_nc = 'gpcp_1dd_v1.2_p1d.%4d.nc' % year
        nc = Dataset('%s%s' % (path_nc, file_nc), 'w', format = 'NETCDF3_CLASSIC')

        # create dimensions
        lat = nc.createDimension('latitude', nlat)
        lon = nc.createDimension('longitude', nlon)
        time = nc.createDimension('time', nday)

        # create variables
        lats = nc.createVariable('latitude', 'f4', ('latitude',))
        lons = nc.createVariable('longitude', 'f4', ('longitude',))
        times = nc.createVariable('time', 'f4', ('time',))
        precip = nc.createVariable('precip', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%d' % -99999.)

        # create global attributes
        nc.title = 'GPCP 1DD Precipitation V1.2'
        nc.history = 'Converted to netCDF format by Jackson Tan (jackson.tan@monash.edu) on ' \
            + datetime.utcnow().strftime('%Y-%m-%d %H%M') + ' UTC. See reference for original data.'
        nc.references = 'http://precip.gsfc.nasa.gov/'

        # create variable attributes
        lats.standard_name = 'latitude'
        lats.long_name = 'latitude of grid box mid-point'
        lats.units = 'degree_north'
        lons.standard_name = 'longitude'
        lons.long_name = 'longitude of grid box mid-point'
        lons.units = 'degree_east'
        times.units = 'days since 1996-10-01'
        times.calendar = 'gregorian'
        precip.long_name = 'daily precipitation'
        precip.units = 'mm / day'

        # write time-independent data
        lats[:] = np.linspace(89.5, -89.5, nlat)
        lons[:] = np.linspace(0.5, 359.5, nlon)

        # write time info
        if year == 1996:
            dates = [datetime(year, 10, 1) + n * timedelta(days = 1) for n in range(nday)]
            times[:] = date2num(dates, units = times.units, calendar = times.calendar)
        else:
            dates = [datetime(year, 1, 1) + n * timedelta(days = 1) for n in range(nday)]
            times[:] = date2num(dates, units = times.units, calendar = times.calendar)

        precip[:] = P

        nc.close()
