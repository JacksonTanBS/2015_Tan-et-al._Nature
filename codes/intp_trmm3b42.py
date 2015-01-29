#!/usr/bin/python3

# A script to interpolate the TRMM data to match ISCCP regular area boxes. The data is coarsened to 2.5 deg resolution first.

import numpy as np
import sys
from calendar import isleap
from trend_funcs import read_netcdf
from scipy.interpolate import griddata
from netCDF4 import Dataset

year = int(sys.argv[1])
datapath = '/home/btan1/Monash/Data/'
trmmpath = '/media/Sentinel/data/TRMM/3B42/netcdf/'

# ninja-load the lat/lon info
with Dataset('%s3B42.199801.nc' % (trmmpath,)) as f:
    lats = f.variables['latitude'][:].copy()
    lons = f.variables['longitude'][:].copy()

lat0 = np.where(lats == -34.875)[0][0]
lat1 = np.where(lats ==  34.875)[0][0] + 1
nlat = lat1 - lat0
nlon = len(lons)

# extract lat/lon coordinates of the TRMM grid
lats_avg = np.array([np.mean(lats[lat0 + ii : lat0 + ii + 10]) for ii in range(0, nlat, 10)])
lons_avg = np.array([np.mean(lons[jj : jj + 10]) for jj in range(0, nlon, 10)])
lons_ext = np.hstack([lons_avg, lons_avg + 360]).flatten()
pts = np.array([[lat, lon] for lat in lats_avg for lon in lons_ext])

# extract lat/lon coordinates of ISCCP grid
box_lat = read_netcdf(datapath, 2007, 1, variable = 'box_lat')
box_lon = read_netcdf(datapath, 2007, 1, variable = 'box_lon')

# load data
P = np.ma.masked_all([(365 + isleap(year)) * 8, len(box_lat)], dtype = np.float32)
mday = 0

for month in range(1, 13):

    f = Dataset('%s3B42.%d%02d.nc' % (trmmpath, year, month), 'r')
    nday = len(f.variables['pcp'])

    for day in range(nday):

        # average the data to 1 deg boxes so values are comparable to GPCP prior to interpolation
        trmm = f.variables['pcp'][day, lat0 : lat1] * 24
        trmm_avg = np.ma.masked_invalid([[np.ma.mean(trmm[ii : ii + 10, jj : jj + 10])
            for jj in range(0, nlon, 10)] for ii in range(0, nlat, 10)])

        # perform the interpolation
        trmm_ext = np.ma.concatenate([trmm_avg, trmm_avg], 1).flatten()
        mask = ~np.ma.getmaskarray(trmm_ext)
        mask_intp = griddata(pts, mask, (box_lat, box_lon), method = 'linear')
        P[mday + day] = np.ma.masked_where(mask_intp != 1,
            griddata(pts[mask], trmm_ext[mask], (box_lat, box_lon), method = 'linear'))

    mday += nday
    f.close()

with open('%strend/precip_data/trmm_%d.npy' % (datapath, year), 'wb') as f:
    np.ma.dump(P, f)
