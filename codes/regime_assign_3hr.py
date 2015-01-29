#!/usr/local/bin/python

# An alternative script to assigning the regime, but to all eight regimes, the CD sub-regimes + 7 regimes. The satellite ID code is included in the output to allow analysis within each domain and track satellite changes. This script produces 3-hourly regimes by using the native 3-hour ISCCP resolution.

import numpy as np
import sys
import os
from trend_funcs import equal_area
from calendar import monthrange

year_start, year_end = 1983, 2009
inpath = '/media/Sentinel/data/ISCCP/D1/'      # path where raw ISCCP binary files are kept
outpath = '/home/btan1/Monash/Data/trend/raw_3hr/'  # path where regime files are output to
cenpath = '/home/btan1/Monash/Data/trend/'            # path to the centroids file
fv = 255.    # fill value
box_start, box_end = 1406, 5190     # starting and ending boxes for 35S to 35N

# read the centroids file
C8 = np.loadtxt('%scentroids_8.dat' % cenpath)
C9 = np.loadtxt('%scentroids_subCD_3hr.txt' % cenpath)

# retrieve the lat/lon coordinates from the equal_area function
lat, lon = equal_area()
lat, lon = lat[box_start : box_end], lon[box_start : box_end]

def read_d1(year, month, day, hour, inpath = '/media/Sentinel/data/ISCCP/D1/'):
    '''Reads and returns the entire ISCCP D1 record for the specified 3-hour period.

    Arguments: year, month, day, hour, [inpath]
    Returns: array of D1 values in dims (box, variable)'''

    import os
    import sys

    cell = np.zeros([6596, 202])

    # load ISCCP binary file (accommodating the variation in versions in filename)
    if 'data' in locals(): del data
    for ver in (0, 1, 2, 3):
        filename = 'ISCCP.D1.%d.GLOBAL.%04d.%02d.%02d.%04d.GPC' % (ver, year, month, day + 1, hour * 100)
        if os.path.exists('%s%04d/%s' % (inpath, year, filename)):
            data = np.fromfile('%s%04d/%s' % (inpath, year, filename), dtype = 'B')
            break
    if 'data' not in locals():
        sys.exit('Input file for %04d %02d %02d %04d not found!' % (year, month, day + 1, hour * 100))

    # verify that year/month/day/hour matches with filename
    if data[2] != (year % 100) or data[3] != month or data[4] != (day + 1) or data[5] != hour:
        sys.exit('Meta data inconsistent with expectation for file %s.' % filename)

    # breaks the string of data into its respective cells (see Fig. 2.2 of documentation)
    for box in range(6596):
        index = 202 * (box + 1) + (box // 99 * 202)    # starting index of the cell (ignoring record prefix)
        cell[box] = data[index : index + 202]

    return cell

for year in range(year_start, year_end + 1):

    # define the months for the year
    if year == 1983: months = list(range(7, 13))
    else: months = list(range(1, 13))

    for month in months:

        output = []    # list to collect data

        for day in range(monthrange(year, month)[1]):
            for hour in range(0, 24, 3):

                cell = read_d1(year, month, day, hour)[box_start : box_end]

                for box in range(box_end - box_start):

                    if cell[box, 5] > 100:    # if night or masked
                        S = 9
                        R8, R9 = 0, 0
                        satID, satangle = 0, 0
                    else:
                        S = cell[box, 5]
                        hist = cell[box, 29 : 71] / cell[box, 10] * 100
                        satID = cell[box, 4]
                        satangle = cell[box, 7] / 100.

                        # assign the regimes (and sub-regimes)
                        R8 = np.sqrt(np.sum((hist - C8) ** 2, -1)).argmin() + 1
                        if R8 == 1:
                            R9 = np.sqrt(np.sum((hist - C9[:2]) ** 2, -1)).argmin() + 1
                        else:
                            R9 = R8 + 1

                    # record the desired data
                    output.append(( year, month, day + 1, hour * 100, lat[box], lon[box], S, R9, R8, satID, satangle))

        # save output to file
        format = ('%4d', '%02d', '%02d', '%04d', '%6.2f', '%6.2f', '%d', '%d', '%d', '%2d', '%5.3f')
        filename = 'regimes_%4d%02d.asc' % (year, month)
        np.savetxt(outpath + filename, output, fmt = format, delimiter = '  ')

'''
Notes on cell index (see Documentation).

joint-histogram = (cell[:, 29 : 71].T / cell[:, 10]).T
IR-histogram = (cell[:, 22 : 29].T / cell[:, 10]).T

cell[:, 22 : 29]    # CTP info only
cell[:, 29 : 71]    # CTP-tau info
cell[:, 10]    # total pixel count
cell[:, 11]    # cloudy pixel count
cell[:, 12]    # IR-cloudy pixel count
'''
