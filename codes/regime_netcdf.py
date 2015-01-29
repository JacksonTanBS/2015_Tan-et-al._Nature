#!/usr/local/bin/python

# This Python script reads the raw ASCII data of the ISCCP regime information and converts it to netCDF format. This script does not interpolate the regime field.

import numpy as np
import calendar as cal
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num

# preliminaries
res = 'daily'    # 'daily' (i.e. average over daytime) or '3hr'
year_start, year_end = 1983, 2009
inpath = '/home/btan1/Monash/Data/trend/raw_%s/' % res    # path where the uninterpolated regime files reside
outpath = '/home/btan1/Monash/Data/trend/netcdf_%s/' % res    # path to output the interpolated regime files
nbox = 3784    # number of entries for each day

lat_new, lon_new = np.mgrid[-33.75 : 34 : 2.5, 1.25 : 360 : 2.5]    # coordinates of new grid
format = ('%4d', '%02d', '%02d', '%6.2f', '%6.2f', '%d', '%d', '%d')    # output format

for year in range(year_start, year_end + 1):

    # define the months for the year
    if year == 1983: months = list(range(7, 13))
    else: months = list(range(1, 13))

    print('Processing for year %d...' % year)

    #--- SET UP THE NETCDF FILE HERE #---

    f = Dataset('%sregimes_%4d.nc' % (outpath, year), 'w', format = 'NETCDF3_CLASSIC')

    # create dimensions
    box = f.createDimension('box', nbox)
    time = f.createDimension('time', None)

    # create variables
    times = f.createVariable('time', 'f4', ('time',))
    box = f.createVariable('box', 'i2', ('box',))
    box_lat = f.createVariable('box_lat', 'f4', ('box',))
    box_lon = f.createVariable('box_lon', 'f4', ('box',))
    R9 = f.createVariable('Regime9', 'i2', ('time', 'box'), fill_value = '%d' % 0)
    R8 = f.createVariable('Regime8', 'i2', ('time', 'box'), fill_value = '%d' % 0)
    surf = f.createVariable('Surface', 'i2', ('time', 'box'), fill_value = '%d' % 9)
    satid = f.createVariable('SatelliteID', 'i2', ('time', 'box'), fill_value = '%d' % 0)
    angle = f.createVariable('SatelliteAngle', 'f4', ('time', 'box'), fill_value = '%f' % 0)

    # create global attributes
    f.description = 'ISCCP Tropical Cloud Regimes'
    f.author = 'Jackson Tan, jackson.tan@monash.edu'
    f.history = 'Created ' + datetime.utcnow().strftime('%Y-%m-%d %H%M') + ' UTC'

    # create variable attributes
    times.units = 'hours since 1983-01-01 00:00:00'
    times.calendar = 'gregorian'
    box.name = 'ISCCP box (35 N/S)'
    box_lat.name = 'Latitude of mid-point of box'
    box_lon.name = 'Longitude of mid-point of box'
    R9.name = 'ISCCP tropical cloud regimes (2 CD sub-regimes + 7)'
    R9.key = '1: CD1, 2: CD2, 3: CC, 4: IM, 5: IC, 6: ST, 7: SS1, 8: SS2, 9: SS3, 0: missing'
    R8.name = 'ISCCP tropical cloud regimes (8 original regimes)'
    R8.key = '1: CD, 2: CC, 3: IM, 4: IC, 5: ST, 6: SS1, 7: SS2, 8: SS3, 0: missing'
    surf.name = 'Surface type designated by ISCCP'
    surf.key = '0: water; 1: land; 2: coast; 9: missing'
    satid.name = 'Satellite ID code'
    satid.key = 'See http://isccp.giss.nasa.gov/docs/sathist.html.'
    angle.name = 'Cosine of daily-mean satellite zenith angle'

    # write time info
    if res == 'daily':
        nday = (np.sum(cal.mdays[months[0] : months[-1] + 1]) + cal.isleap(year))
        dates = [datetime(year, months[0], 1) + n * timedelta(days = 1) for n in range(nday)]
    elif res == '3hr':
        nday = (np.sum(cal.mdays[months[0] : months[-1] + 1]) + cal.isleap(year)) * 8
        dates = [datetime(year, months[0], 1) + n * timedelta(hours = 3) for n in range(nday)]
    times[:] = date2num(dates, units = times.units, calendar = times.calendar)

    #--- FINISH SETTING UP NETCDF FILE ---#

    #--- DO THE INTERPOLATION AND OUTPUT DATA TO NETCDF FILE ---#

    dd = 0    # counter for day of year
    for month in months:

        # read the raw data
        filename = 'regimes_' + str(year) + str(month).zfill(2) + '.asc'
        g = open(inpath + filename, 'r')
        data = np.array([[np.float(jj) for jj in line.split()] for line in g.readlines()])
        g.close()

        # extract required information
        lats, lons = data[:nbox, 3 + (res == '3hr')], data[:nbox, 4 + (res == '3hr')]
        S = data[:, 5 + (res == '3hr')].reshape(-1, nbox)
        r9 = data[:, 6 + (res == '3hr')].reshape(-1, nbox)
        r8 = data[:, 7 + (res == '3hr')].reshape(-1, nbox)
        sat = data[:, 8 + (res == '3hr')].reshape(-1, nbox)
        ang = data[:, 9 + (res == '3hr')].reshape(-1, nbox)

        box_lat[:] = lats
        box_lon[:] = lons

        # write time-dependent data
        R9[dd : dd + len(S)] = r9
        R8[dd : dd + len(S)] = r8
        surf[dd : dd + len(S)] = S
        satid[dd : dd + len(S)] = sat
        angle[dd : dd + len(S)] = ang

        dd += len(S)

    f.close()
