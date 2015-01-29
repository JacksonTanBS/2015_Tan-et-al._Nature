#!/usr/bin/python3

# A script to interpolate ERA-Interim omega at 500 hPa to match ISCCP regular area boxes.

import numpy as np
import calendar as cal
from scipy.interpolate import griddata
from netCDF4 import Dataset, date2num
from trend_funcs import read_netcdf

datapath = '/home/btan1/Monash/Data/'

# open datasets downloaded from ERA-Interim
f = Dataset('%serai/w500_1990-1999.nc' % datapath, 'r')
g = Dataset('%serai/w500_2000-2009.nc' % datapath, 'r')

time = f.variables['time']
lats = f.variables['latitude']
lons = f.variables['longitude']
lat0 = np.where(lats[:] ==  35.)[0][0]
lat1 = np.where(lats[:] == -35.)[0][0] + 1
w1 = f.variables['w']
w2 = g.variables['w']

# set up omega grid to our domain of interest
t1_start = np.where(time[:] == date2num(cal.datetime.datetime(1998, 1, 1), time.units))[0][0]
w1_daily = np.array([np.ma.mean(w1[tt : tt + 4, lat0 : lat1], 0) for tt in range(t1_start, len(w1), 4)])
w1_daily_ext = np.ma.concatenate([w1_daily, w1_daily[:, :, 0 : 1]], 2)
w2_daily = np.array([np.ma.mean(w2[tt : tt + 4, lat0 : lat1], 0) for tt in range(0, len(w2), 4)])
w2_daily_ext = np.ma.concatenate([w2_daily, w2_daily[:, :, 0 : 1]], 2)
w_daily_ext = np.ma.concatenate([w1_daily_ext, w2_daily_ext], 0)
lons_ext = np.ma.concatenate([lons, [360,]])
pts = np.array([[lat, lon] for lat in lats[lat0 : lat1] for lon in lons_ext.flatten()])

# define ISCCP grid
lat30 = np.where(np.abs(read_netcdf(datapath, 2007, variable = 'box_lat')) < 30)[0]
box_lat = read_netcdf(datapath, 2007, 1, variable = 'box_lat')[lat30]
box_lon = read_netcdf(datapath, 2007, 1, variable = 'box_lon')[lat30]

t0 = 0
for year in range(1998, 2010):

    w500 = np.zeros([365 + cal.isleap(year), len(lat30)], dtype = np.float32)
    w500[:] = np.nan
    t1 = t0 + 365 + cal.isleap(year)

    for dd, day in enumerate(range(t0, t1)):
        w500[dd] = griddata(pts, w_daily_ext[day].flatten(), (box_lat, box_lon), method = 'cubic')

    np.save('%strend/precip_data/w500_%d.npy' % (datapath, year), w500)

    t0 = t1
