#!/usr/bin/python3

# This script downloads and concatenates TRMM 3B42 precipitation dataset, producing a monthly file.
# DISCLAIMER: Use this script and its outputs at your own risk. I cannot guarantee that the conversion is flawless.
# Written by Jackson Tan (jackson.tan@monash.edu), 06/06/14

import os
import sys
import gzip
from calendar import monthrange
from datetime import datetime, timedelta
from urllib.request import urlretrieve
from netCDF4 import Dataset, MFDataset, date2num
from glob import glob

#--- Select variables to copy over ---#

# Note: a month of each variable has an uncompressed size of about 600 MB
# and a compressed size of about 50 MB.

err_cat        = False  # TMPA random error estimate
hqp_cat        = False  # Pre-gauge-adjusted microwave precipitation estimate
irp_cat        = False  # Pre-gauge-adjusted infrared precipitation estimate
pcp_cat        = True   # TMPA precipitation estimate
satobstime_cat = False  # Satellite observation time minus time at timestep
source_cat     = False  # Flag to show source of data

#--- Set year, month and directories ---#

year, month = int(sys.argv[1]), int(sys.argv[2])
path = '/media/Sentinel/data/TRMM/3B42/'
path_nc  = '%snetcdf/' % path
path_tmp = '%stmp/' % path
if not os.path.exists(path_nc ): os.makedirs(path_nc )
if not os.path.exists(path_tmp): os.makedirs(path_tmp)

#--- Download the netCDF files ---#

for day in range(1, monthrange(year, month)[1] + 1):
    for hour in range(0, 24, 3):

        # defining the variations in the URL
        d = datetime(year, month, day, hour) - timedelta(hours = 3)
        yd = d.timetuple().tm_yday
        yr = d.timetuple().tm_year
        if year < 2000: fmt = '7'
        else: fmt = '7A'
        date = '{0:4d}{1:02d}{2:02d}.{3:02d}'.format(year, month, day, hour)

        # guess the URL (see end of script)
        url1 = 'http://disc2.nascom.nasa.gov/daac-bin/OTF/HTTP_services.cgi?'
        url2 = 'FILENAME=%2Fftp%2Fdata%2Fs4pa%2FTRMM_L3%2FTRMM_3B42%2F{0}%2F{1:03d}%'.format(yr, yd)
        url3 = '2F3B42.{0}.{1}.HDF.Z&FORMAT=L2d6aXA&LABEL=3B42.{0}.{1}.nc.gz'.format(date, fmt)
        url4 = '&SHORTNAME=TRMM_3B42&SERVICE=HDF_TO_NetCDF&VERSION=1.02&DATASET_VERSION=007'
        url = url1 + url2 + url3 + url4

        # download all the files into a temporary folder first
        filename = '3B42.{0}.nc'.format(date)
        urlretrieve(url, '%s%s.gz' % (path_tmp, filename))

        # use file size check to ensure that downloaded zip file is not empty
        if os.path.getsize('%s%s.gz' % (path_tmp, filename)) < 1024:
            sys.exit('Error: downloaded file is less than 1 kB. Defined URL is likely wrong.')

        # decompress the files
        with open('%s%s' % (path_tmp, filename), 'wb') as f:
            with gzip.open('%s%s.gz' % (path_tmp, filename), 'r') as g:
                f.write(g.read())
        os.remove('%s%s.gz' % (path_tmp, filename))

#--- Concatenate the netCDF files ---#

# set up the output netCDF file
date = '{0:4d}{1:02d}'.format(year, month)
filename = '3B42.{0}'.format(date)
nc = Dataset('%s%s.nc' % (path_nc, filename), 'w', format = 'NETCDF3_CLASSIC')

# multi-read the netCDF files and extract basic data
f = MFDataset('%s%s*.nc' % (path_tmp, filename))
nlat, nlon = len(f.variables['latitude']), len(f.variables['longitude'])
nday = monthrange(year, month)[1] * 8

# create global attributes
nc.title = 'TRMM 3B42 Precipitation Product'
nc.Conventions = 'COARDS'
nc.comments = 'file created by grads using lats4d available from http://dao.gsfc.nasa.gov/software/grads/lats4d/'
nc.model = 'geos/das'
nc.center = 'gsfc'
nc.history = 'Concatenated to monthly files by Jackson Tan (jackson.tan@monash.edu) on ' \
    + datetime.utcnow().strftime('%Y-%m-%d %H%M') + ' UTC.'

# create dimensions
latitude = nc.createDimension('latitude', nlat)
longitude = nc.createDimension('longitude', nlon)
time = nc.createDimension('time', nday)

# write the compulsory variables
latitude = nc.createVariable('latitude', 'f8', ('latitude',))
latitude.units = f.variables['latitude'].units
latitude.long_name = f.variables['latitude'].long_name
latitude[:] = f.variables['latitude'][:]
longitude = nc.createVariable('longitude', 'f8', ('longitude',))
longitude.units = f.variables['longitude'].units
longitude.long_name = f.variables['longitude'].long_name
longitude[:] = f.variables['longitude'][:]
time = nc.createVariable('time', 'f8', ('time',))
time.units = 'hours since 1998-10-01 00'
time.calendar = 'standard'
dates = [datetime(year, 1, 1, 0) + nn * timedelta(hours = 3) for nn in range(nday)]
time[:] = date2num(dates, units = time.units, calendar = time.calendar)

# write the optional variables
if err_cat:
    err = nc.createVariable('err', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    err.units = 'mm / hr'
    err.long_name = 'relative error'
    err.grid_name = f.variables['err'].grid_name
    err.grid_type = f.variables['err'].grid_type
    err.level_description = f.variables['err'].level_description
    err.time_statistic = f.variables['err'].time_statistic
    err[:] = f.variables['err'][:]
if hqp_cat:
    hqp = nc.createVariable('hqp', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    hqp.units = 'mm / hr'
    hqp.long_name = 'high quality precipitation'
    hqp.grid_name = f.variables['hqp'].grid_name
    hqp.grid_type = f.variables['hqp'].grid_type
    hqp.level_description = f.variables['hqp'].level_description
    hqp.time_statistic = f.variables['hqp'].time_statistic
    hqp[:] = f.variables['hqp'][:]
if irp_cat:
    irp = nc.createVariable('irp', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    irp.units = 'mm / hr'
    irp.long_name = 'IR precipitation'
    irp.grid_name = f.variables['irp'].grid_name
    irp.grid_type = f.variables['irp'].grid_type
    irp.level_description = f.variables['irp'].level_description
    irp.time_statistic = f.variables['irp'].time_statistic
    irp[:] = f.variables['irp'][:]
if pcp_cat:
    pcp = nc.createVariable('pcp', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    pcp.units = 'mm / hr'
    pcp.long_name = 'precipitation'
    pcp.grid_name = f.variables['pcp'].grid_name
    pcp.grid_type = f.variables['pcp'].grid_type
    pcp.level_description = f.variables['pcp'].level_description
    pcp.time_statistic = f.variables['pcp'].time_statistic
    pcp[:] = f.variables['pcp'][:]
if satobstime_cat:
    satobstime = nc.createVariable('satobstime', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    satobstime.units = 'minutes'
    satobstime.long_name = 'observation time'
    satobstime.grid_name = f.variables['satobstime'].grid_name
    satobstime.grid_type = f.variables['satobstime'].grid_type
    satobstime.level_description = f.variables['satobstime'].level_description
    satobstime.time_statistic = f.variables['satobstime'].time_statistic
    satobstime.description = 'satellite observation time - time at timestep'
    satobstime[:] = f.variables['satobstime'][:]
if source_cat:
    source = nc.createVariable('source', 'f4', ('time', 'latitude', 'longitude'), fill_value = '%f' % -9999.9)
    source.units = ''
    source.long_name = 'source'
    source.grid_name = f.variables['source'].grid_name
    source.grid_type = f.variables['source'].grid_type
    source.level_description = f.variables['source'].level_description
    source.time_statistic = f.variables['source'].time_statistic
    source[:] = f.variables['source'][:]

nc.close()
f.close()

#--- Remove the temporary files ---#

[os.remove(ff) for ff in glob('%s3B42.%s*.nc' % (path_tmp, date))]

#--- Guessing the URL ---#

# To determine the URL of the file, I requested for 3B42 through Mirador,
# proceeding until I reach the point where I am given a list of URLs and
# instructed to use wget to retrieve the files. The URL used here is inferred
# from this list of URLs, and should apply for all the files. If the URL is
# incorrect, the urlretrieve function will download a gzipped shell of very
# small file size, which will then be identified by the subsequent line. In
# such a situation, you will need to check the URL to ensure that it is
# correct.
#
# Mirador: http://mirador.gsfc.nasa.gov/
