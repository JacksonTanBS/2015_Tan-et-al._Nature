#!/usr/bin/python3

# A script to interpolate the GPCP 1DD data to match ISCCP regular area boxes.

import numpy as np
import matplotlib.pyplot as plt
import sys
from trend_funcs import read_netcdf
from scipy.interpolate import griddata
from netCDF4 import Dataset

year = int(sys.argv[1])
datapath = '/home/btan1/Monash/Data/'
gpcppath = '/media/Sentinel/data/GPCP/1DD/netcdf/'

# ninja-load the lat/lon info
with Dataset('%sgpcp_1dd_v1.2_p1d.2007.nc' % (gpcppath,)) as f:
    lats = f.variables['latitude'][:].copy()
    lons = f.variables['longitude'][:].copy()

lat0 = np.where(lats ==  34.5)[0][0]
lat1 = np.where(lats == -34.5)[0][0] + 1
nlat = lat1 - lat0
nlon = len(lons)

box_lat = read_netcdf(datapath, 2007, 1, variable = 'box_lat')
box_lon = read_netcdf(datapath, 2007, 1, variable = 'box_lon')

# load data
with Dataset('%sgpcp_1dd_v1.2_p1d.%d.nc' % (gpcppath, year), 'r') as f:
    gpcp = f.variables['precip'][:, lat0 : lat1].copy()

P = np.zeros([len(gpcp), len(box_lat)])
P[:] = np.nan

for day in range(len(gpcp)):

    gpcp_ext = np.concatenate([gpcp[day], gpcp[day], gpcp[day]], 1).flatten()
    lons_ext = np.hstack([lons - 360, lons, lons + 360]).flatten()
    pts = np.array([[lat, lon] for lat in lats[lat0 : lat1] for lon in lons_ext])
    P[day] = griddata(pts, gpcp_ext, (box_lat, box_lon), method = 'linear')

np.save('%strend/precip_data/gpcp1dd_%d.npy' % (datapath, year), P)
