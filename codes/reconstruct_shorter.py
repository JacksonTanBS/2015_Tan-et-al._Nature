#!/usr/bin/python3

# This script compares the reconstructions of rainfall trends from the regime trends to the observed changes. Choice of either TRMM (trmm or trmmdaytime) or GPCP 1DD. 'trmm' is TRMM 3B42 at every 3 hr (not used in paper); 'trmmdaytime' is daytime-averaged TRMM 3B42, i.e. averaged to daily with nighttime values masked (see masknight_trmm.py).

import numpy as np
import calendar as cal
import os
from trend_funcs import read_netcdf, read_netcdf3hr

pdata = 'trmmdaytime'    # 'gpcp1dd', 'trmm', or 'trmmdaytime'
combine = True    # if True, treats CR1-3 as one; else, separately
datapath = '/home/btan1/Monash/Data/'
if   pdata == 'gpcp1dd': years = list(range(1997, 2010))
elif pdata == 'trmm': years = list(range(1998, 2010))
elif pdata == 'trmmdaytime': years = list(range(1998, 2010))
nbox = 3784
if combine: nreg = 1
else: nreg = 8
nmon = len(years) * 12

if combine:
    file_F_r = '%strend/P_reconstruct/%s_data_F(1-3)_r.npy' % (datapath, pdata)
    file_P_r = '%strend/P_reconstruct/%s_data_P(1-3)_r.npy' % (datapath, pdata)
else:
    file_F_r = '%strend/P_reconstruct/%s_data_F_r.npy' % (datapath, pdata)
    file_P_r = '%strend/P_reconstruct/%s_data_P_r.npy' % (datapath, pdata)
file_P_sum ='%strend/P_reconstruct/%s_data_P_sum.npy' % (datapath, pdata)

# calculate the trends of regime FOC and precipitation
F_r = np.ma.masked_all([nreg, nbox, nmon])
P_r = np.ma.masked_all([nreg, nbox, nmon])
P_sum = np.ma.masked_all([nbox, nmon])

tt = 0
for year in years:

    P = np.load('%strend/precip_data/%s_%d.npy' % (datapath, pdata, year))

    t0 = 0
    for month in range(1, 13):

        if pdata == 'trmm': t1 = t0 + cal.monthrange(year, month)[1] * 8
        else: t1 = t0 + cal.monthrange(year, month)[1]

        P_sum[:, tt] = np.ma.masked_invalid([np.ma.mean(P[t0 : t1, box]) for box in range(nbox)])

        if pdata == 'trmm': R = read_netcdf3hr(datapath, year, month, variable = 'Regime8')
        else: R = read_netcdf(datapath, year, month, variable = 'Regime8')
        if len(R) != (t1 - t0): print('Warning: month range does not match regime length.')

        if combine:
            R[(R == 2) + (R == 3)] = 1

        for rr in range(nreg):
            window = R == rr + 1
            count = np.ma.sum(R == rr + 1, 0)
            total = np.float32(np.ma.count(R, 0))
            F_r[rr, :, tt] = np.ma.masked_where(total == 0, count / total)

            for box in range(nbox):
                if np.ma.sum(window[:, box]) != 0:
                    P_r[rr, box, tt] = np.ma.mean(P[t0 : t1, box][window[:, box]])

        t0 = t1
        tt += 1

with open(file_F_r, 'wb') as f: np.ma.dump(F_r, f)
with open(file_P_r, 'wb') as f: np.ma.dump(P_r, f)
if not combine:
    with open(file_P_sum, 'wb') as f: np.ma.dump(P_sum, f)
