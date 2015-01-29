#!/usr/bin/python3

# A script to interpolate the TRMM 3A25 data to match ISCCP regular area boxes. The data is coarsened to 1 deg resolution prior to interpolation.

import numpy as np
import sys
import gzip
import h5py
from trend_funcs import read_netcdf
from scipy.interpolate import griddata

datapath = '/home/btan1/Monash/Data/'
trmmpath = '/media/Sentinel/data/TRMM/3A25/'
nmon = 144

# define the lat/lon info (as it went missing from HDF4 to HDF5 conversion)
lats = np.arange(-36.75, 36.76, 0.5)
lons = np.arange(-179.75, 179.76, 0.5)

# choosing lat cutoff such that it is divisible by 5 (for coarsening purposes)
lat0 = np.where(lats == -34.75)[0][0]
lat1 = np.where(lats ==  34.75)[0][0] + 1
nlat = lat1 - lat0
nlon = len(lons)

# read ISCCP lat/lon
# NOTE: if your HDF5 installation is version 1.8.5, use this:
box_lat = np.load(datapath + 'trend/hdf_issue_bypass/box_lat.npy')
box_lon = np.load(datapath + 'trend/hdf_issue_bypass/box_lon.npy')
# otherwise, use this:
#box_lat = read_netcdf(datapath, 2007, 1, variable = 'box_lat')
#box_lon = read_netcdf(datapath, 2007, 1, variable = 'box_lon')
# see Note at bottom for explanation.

# pre-defining arrays for convective/stratiform area/rain fractions
x = np.ma.masked_all([nmon, len(box_lat)], dtype = np.float32)
cap = x.copy()      # convective pixels
sap = x.copy()      # stratiform pixels
tap = x.copy()      # total pixels
caf = x.copy()      # convective area fraction
saf = x.copy()      # stratiform area fraction
srr = x.copy()      # stratiform rain rate
trr = x.copy()      # convective rain rate
crf = x.copy()      # convective rain fraction
srf = x.copy()      # stratiform rain fraction
msr = x.copy()      # mean stratiform rain rate
mcr = x.copy()      # mean convective rain rate
tt = 0

for year in range(1998, 2010):
    for month in range(1, 13):

        filename = '3A25.%4d%02d01.7.HDF5' % (year, month)

        f = h5py.File('%s%s' % (trmmpath, filename), 'r')

        # read convective/stratiform/all pixels and unconditional rain rates
        pixA = f['Grids/Grid2/e_surfRainPix2'][:, lat0 : lat1].T
        pixC = f['Grids/Grid2/e_surfRainConvPix2'][:, lat0 : lat1].T
        pixS = f['Grids/Grid2/e_surfRainStratPix2'][:, lat0 : lat1].T
        rainA = f['Grids/Grid2/e_surfRainMean2'][:, lat0 : lat1].T * pixA
        rainC = f['Grids/Grid2/e_surfRainConvMean2'][:, lat0 : lat1].T * pixC
        rainS = f['Grids/Grid2/e_surfRainStratMean2'][:, lat0 : lat1].T * pixS

        # average to 2.5 deg boxes (sum pixels, average rain rate)
        pixA_avg = np.array([[np.ma.sum(pixA[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        pixC_avg = np.array([[np.ma.sum(pixC[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        pixS_avg = np.array([[np.ma.sum(pixS[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        rainA_avg = np.array([[np.ma.mean(rainA[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        rainC_avg = np.array([[np.ma.mean(rainC[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        rainS_avg = np.array([[np.ma.mean(rainS[ii : ii + 5, jj : jj + 5])
            for jj in range(0, nlon, 5)] for ii in range(0, nlat, 5)])
        lats_avg = np.array([np.mean(lats[lat0 + ii : lat0 + ii + 5]) for ii in range(0, nlat, 5)])
        lons_avg = np.array([np.mean(lons[jj : jj + 5]) for jj in range(0, nlon, 5)])

        # perform the interpolation
        pixA_ext = np.ma.concatenate([pixA_avg, pixA_avg], 1).flatten()
        pixC_ext = np.ma.concatenate([pixC_avg, pixC_avg], 1).flatten()
        pixS_ext = np.ma.concatenate([pixS_avg, pixS_avg], 1).flatten()
        rainA_ext = np.ma.concatenate([rainA_avg, rainA_avg], 1).flatten()
        rainC_ext = np.ma.concatenate([rainC_avg, rainC_avg], 1).flatten()
        rainS_ext = np.ma.concatenate([rainS_avg, rainS_avg], 1).flatten()
        lons_ext = np.hstack([lons_avg, lons_avg + 360]).flatten()
        pts = np.array([[lat, lon] for lat in lats_avg for lon in lons_ext])
        pixA_intp = griddata(pts, pixA_ext, (box_lat, box_lon), method = 'linear')
        pixC_intp = griddata(pts, pixC_ext, (box_lat, box_lon), method = 'linear')
        pixS_intp = griddata(pts, pixS_ext, (box_lat, box_lon), method = 'linear')
        rainA_intp = griddata(pts, rainA_ext, (box_lat, box_lon), method = 'linear')
        rainC_intp = griddata(pts, rainC_ext, (box_lat, box_lon), method = 'linear')
        rainS_intp = griddata(pts, rainS_ext, (box_lat, box_lon), method = 'linear')

        # calculate convective/stratiform area/rain fraction, masking grids with only one pixel
        cap[tt] = pixC_intp
        sap[tt] = pixS_intp
        tap[tt] = pixA_intp
        caf[tt] = np.ma.masked_where(pixA_intp < 1, pixC_intp / (pixA_intp + 1e-9))
        saf[tt] = np.ma.masked_where(pixA_intp < 1, pixS_intp / (pixA_intp + 1e-9))
        srr[tt] = rainS_intp
        trr[tt] = rainA_intp
        crf[tt] = np.ma.masked_where(rainA_intp < 1e-9, rainC_intp / (rainA_intp + 1e-9))
        srf[tt] = np.ma.masked_where(rainA_intp < 1e-9, rainS_intp / (rainA_intp + 1e-9))
        tt += 1

        f.close()

with open('%strend/precip_data/trmm3a25cap.npy' % datapath, 'wb') as f:
    np.ma.dump(cap, f)
with open('%strend/precip_data/trmm3a25sap.npy' % datapath, 'wb') as f:
    np.ma.dump(sap, f)
with open('%strend/precip_data/trmm3a25tap.npy' % datapath, 'wb') as f:
    np.ma.dump(tap, f)
with open('%strend/precip_data/trmm3a25caf.npy' % datapath, 'wb') as f:
    np.ma.dump(caf, f)
with open('%strend/precip_data/trmm3a25saf.npy' % datapath, 'wb') as f:
    np.ma.dump(saf, f)
with open('%strend/precip_data/trmm3a25srr.npy' % datapath, 'wb') as f:
    np.ma.dump(srr, f)
with open('%strend/precip_data/trmm3a25trr.npy' % datapath, 'wb') as f:
    np.ma.dump(trr, f)
with open('%strend/precip_data/trmm3a25crf.npy' % datapath, 'wb') as f:
    np.ma.dump(crf, f)
with open('%strend/precip_data/trmm3a25srf.npy' % datapath, 'wb') as f:
    np.ma.dump(srf, f)


# Note regarding box_lat and box_lon:
# h5py will throw a fit if you are using HDF5 1.8.5 and then use netCDF4 functions. See https://groups.google.com/forum/#!topic/h5py/AZQ30pSy-RI for explanation. To work around this issue, load box_lat/box_lon using netCDF4, then save the file using numpy:
# np.save(datapath + 'trend/hdf_issue_bypass/box_lat.npy', box_lat)
# np.save(datapath + 'trend/hdf_issue_bypass/box_lon.npy', box_lon)
