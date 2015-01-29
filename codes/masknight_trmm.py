#!/usr/bin/python3

# This code is an efficient (i.e. inconsequentially lazy) way to mask nighttime values of TRMM 3B42 so as to take the daytime-averages (as opposed to daily-averages). This is achieved by extracting the missing values of the 3-hourly regimes and assuming that missing values are nighttime values. This is true in most cases, since the exception of missing values due to satellite issues is relatively infrequent.

import numpy as np
from trend_funcs import read_netcdf3hr

datapath = '/home/btan1/Monash/Data/'
trmmpath = '%strend/precip_data/' % datapath

years = list(range(1998, 2010))

for year in years:

    # load the desired variables
    Praw = np.load('%strmm_%d.npy' % (trmmpath, year))
    R = read_netcdf3hr(datapath, year)
    P = np.ma.masked_where(np.ma.getmask(R), Praw)

    # do a consistency check
    if np.ma.any(Praw != P):
        print('Warning: masked array does not match parent for %d' % year)
        
    Pday = np.ma.masked_array([np.ma.mean(P[ii : ii + 8], 0) for ii in range(0, len(P), 8)])

    with open('%strmmdaytime_%d.npy' % (trmmpath, year), 'wb') as f:
        np.ma.dump(Pday, f)
