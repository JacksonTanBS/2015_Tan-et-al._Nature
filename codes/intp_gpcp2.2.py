#!/usr/bin/python3

# A script to interpolate the GPCP 2.2 (monthly) data to match ISCCP regular area boxes.

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from trend_funcs import read_netcdf
from scipy.interpolate import griddata
from netCDF4 import Dataset, date2num

datapath = '/home/btan1/Monash/Data/'
gpcppath = '/media/Sentinel/data/GPCP/2.2/'

f = Dataset('%sprecip.mon.mean.nc' % (gpcppath,))

lats = f.variables['lat']
lons = f.variables['lon']
gpcp = f.variables['precip']
time = f.variables['time']

# define the cutoffs
t0 = np.where(time[:] == date2num(datetime(1983, 7, 1), time.units))[0][0]
t1 = np.where(time[:] == date2num(datetime(2009, 12, 1), time.units))[0][0] + 1
lat0 = np.where(lats[:] ==  33.75)[0][0]
lat1 = np.where(lats[:] == -33.75)[0][0] + 1
nlat = lat1 - lat0
nlon = len(lons)

box_lat = read_netcdf(datapath, 2007, 1, variable = 'box_lat')
box_lon = read_netcdf(datapath, 2007, 1, variable = 'box_lon')

P = np.ma.zeros([t1 - t0, len(box_lat)])
tt = 0

for tt in range(t1 - t0):

    gpcp_ext = np.ma.concatenate([gpcp[t0 + tt, lat0 : lat1], gpcp[t0 + tt, lat0 : lat1],
        gpcp[t0 + tt, lat0 : lat1]], 1).flatten()
    lons_ext = np.hstack([lons[:] - 360, lons[:], lons[:] + 360]).flatten()
    pts = np.array([[lat, lon] for lat in lats[lat0 : lat1] for lon in lons_ext])
    P[tt] = griddata(pts, gpcp_ext, (box_lat, box_lon), method = 'linear')

with open('%strend/precip_data/gpcp2.2.npy' % datapath, 'wb') as f:
    np.ma.dump(P, f)
