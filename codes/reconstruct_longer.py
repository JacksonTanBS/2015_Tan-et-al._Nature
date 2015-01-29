#!/usr/bin/python3

# This script compares the reconstructions of rainfall trends from the regime trends to the observed changes in GPCP 2.2.

import numpy as np
import calendar as cal
from trend_funcs import read_netcdf

pdata = 'gpcp2.2'
datapath = '/home/btan1/Monash/Data/'
nbox = 3784
nreg = 8
nmonF = 318

# calculate the trends of regime FOC
file_F_r = '%strend/P_reconstruct/%s_data_F_r.npy' % (datapath, pdata)
F_r = np.ma.masked_all([nreg, nbox, nmonF])
tt = 0

for year in range(1983, 2010):

    # define the months for the year
    if year == 1983: months = list(range(7, 13))
    else: months = list(range(1, 13))

    for month in months:

        R = read_netcdf(datapath, year, month, variable = 'Regime8')
        for rr in range(nreg):
            count = np.sum(R == rr + 1, 0)
            total = np.float32(np.ma.count(R, 0))
            F_r[rr, :, tt] = np.ma.masked_where(total == 0, count / total)
        tt += 1

with open(file_F_r, 'wb') as f:
    np.ma.dump(F_r, f)
