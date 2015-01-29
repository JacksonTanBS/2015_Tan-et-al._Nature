#!/usr/local/bin/python3

# This Python script provides supporting functions for the main plotting script.

import numpy as np
import calendar as cal
from netCDF4 import Dataset

def read_netcdf(datapath, year, month = 'all', variable = 'Regime9'):
    '''Extracts the variables from the modified regime netCDF files.

    Arguments: read_regime(datapath, year, [month, variable])
    Output: 2-dim array (dims: day, box) of type int16
    Requires: netCDF files of the modified regimes
    Remarks: 'variable' can be 'time', 'box', 'box_lat', 'box_lon', 'Regime9'
             'Regime8', 'Surface', 'SatelliteID', or 'SatelliteAngle'.
    '''

    f = Dataset('%strend/netcdf_daily/regimes_%d.nc' % (datapath, year), 'r')
    if variable in ('box', 'box_lat', 'box_lon'):
        data = f.variables[variable][:].copy()
    elif month == 'all':
        data = f.variables[variable][:].copy()
    elif year == 1983:
        day0 = (cal.datetime.datetime(year, month, 1).timetuple().tm_yday
            - cal.datetime.datetime(1983, 7, 1).timetuple().tm_yday)
        day1 = day0 + cal.monthrange(year, month)[1]
        data = f.variables[variable][day0 : day1].copy()
    else:
        day0 = cal.datetime.datetime(year, month, 1).timetuple().tm_yday - 1
        day1 = day0 + cal.monthrange(year, month)[1]
        data = f.variables[variable][day0 : day1].copy()
    f.close()

    return data


def isccp_latlon(datapath = '/home/btan1/Monash/Data/'):
    '''Returns the lat/lon boundaries of each ISCCP equal-area grid box.

    Arguments: isccp_latlon([datapath]):
    Output: a latitude boundary array, a longitude boundary list, and the box
            where a new latitude begins
    Requires: netCDF files of the modified regimes
    Remarks: Each 'row' of the array/list correspond to each latitude band.
    '''

    lat30 = np.where(np.abs(read_netcdf(datapath, 2007, variable = 'box_lat')) < 30)[0]
    lon = read_netcdf(datapath, 2007, variable = 'box_lon')[lat30]

    lats = np.vstack([np.arange(-30, 30, 2.5), np.arange(-27.5, 30.1, 2.5)]).T
    jumps = np.hstack([[0], np.where(np.diff(lon) < 0)[0] + 1, len(lon)])
    lons = []
    for jj in range(len(jumps) - 1):
        # take the slice of longitude for that latitude
        lon_midpt = lon[jumps[jj] : jumps[jj + 1]]
        lon_edges = lon_midpt[1:] - 0.5 * np.diff(lon_midpt)
        lons.append(np.hstack([0., lon_edges, 360.]))

    return lats, lons, jumps


def fitLine(x, y, alpha=0.05, printResidual=False):
    '''Fit a line with confidence intervals from bootstrapping method.
    Arguments: fitLine(x, y, [alpha = 0.05, residual=False])
    Output: gradient, intercept, (gradient C.I, intercept C.I,), (C.I. widths), 
            [residual])

    Modified from:
    Source: http://scipy-central.org/item/50/1/line-fit-with-confidence-intervals
    License: Creative Commons Zero (almost public domain) http://scpyce.org/cc0'''

    import scipy.stats as stats

    # Summary data
    n = len(x)			   # number of samples     
    
    Sxx = np.sum(x**2) - np.sum(x)**2/n
    Sxy = np.sum(x*y) - np.sum(x)*np.sum(y)/n    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    # Linefit
    b = Sxy/Sxx
    a = mean_y - b*mean_x
    
    # Residuals
    fit = lambda xx: a + b*xx    
    residuals = y - fit(x)
    
    var_res = np.sum(residuals**2)/(n-2)
    sd_res = np.sqrt(var_res)
    
    # Confidence intervals
    se_b = sd_res/np.sqrt(Sxx)
    se_a = sd_res*np.sqrt(np.sum(x**2)/(n*Sxx))
    
    df = n-2                            # degrees of freedom
    tval = stats.t.isf(alpha/2., df) 	# appropriate t value

    width_a=tval*se_a
    width_b=tval*se_b
    ci_a = a + width_a*np.array([-1,1])
    ci_b = b + width_b*np.array([-1,1])

    if printResidual:
        ri = {'residuals': residuals, 
            'var_res': var_res,
            'sd_res': sd_res,
            'alpha': alpha,
            'tval': tval,
            'df': df}
        return (b,a,(ci_b, ci_a),(width_b, width_a),ri)
    else:
        return (b,a,(ci_b, ci_a),(width_b, width_a))


def plot_isccp_map(Z, cmap = None, plot_map = True,
    fix_map_aspect = True, vmin = None, vmax = None, meridian = True):
    '''Perform a pcolor plot with the ISCCP irregular equal area grids.

    Arguments: plot_isccp_map(Z, [cmap, plot_map, fix_map_aspect, vmax])
    Output: pseudocolour plot
    Remarks: Z.shape must be (1, 3784)
    '''

    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    lats, lons, jumps = isccp_latlon()

    if vmax == None: vmax = max(np.max(Z), -np.min(Z))
    if cmap == None: cmap = plt.cm.RdBu
    if vmin == None: vmin = -vmax

    for jj in range(len(jumps) - 1):
        mp = plt.pcolor(lons[jj], lats[jj], Z[:, jumps[jj] : jumps[jj + 1]],
            vmin = vmin, vmax = vmax, cmap = cmap)
    if plot_map:
        m = Basemap(projection='cyl', llcrnrlat = -30, urcrnrlat = 30,
            llcrnrlon = 0, urcrnrlon = 360.1, fix_aspect = fix_map_aspect,
            resolution = 'c')
        m.drawcoastlines(color = 'k', linewidth = 0.5)
        m.drawparallels(np.arange(-30, 31, 15), color = 'k',
            linewidth = 0, labels = [1, 0, 0, 0])
        if meridian:
            m.drawmeridians(np.arange(0.,361.,60.), color = 'k',
                linewidth = 0, labels = [0, 0, 0, 1])

    return mp


def plot_isccp_sig(Z, nmc = 10000, a = 0.1, size = 2,
    datapath = '/home/btan1/Monash/Data/'):
    '''Plots stippling on top of the ISCCP maps.

    Arguments: Z (array to test for significance), [a, size, datapath]
    Output: None
    '''

    import matplotlib.pyplot as plt

    mask = np.ma.getmaskarray(delta(Z))
    sig = np.ma.masked_where(mask, permtest(Z, nmc) < a)
    lat30 = np.where(np.abs(read_netcdf(datapath, 2007, variable = 'box_lat')) < 30)[0]
    lat = read_netcdf(datapath, 2007, variable = 'box_lat')[lat30][sig]
    lon = read_netcdf(datapath, 2007, variable = 'box_lon')[lat30][sig]
    plt.scatter(lon, lat, s = size, c = '#ffff00', alpha = 0.9, edgecolor = 'none')

    return None


def delta(X):
    '''Calculates the difference between the second half from the first half
       for each grid box.

    Arguments: X (shape = (nbox, nmon))
    Output: delta(X) (shape = nbox)
    '''

    half = X.shape[1] / 2

    return np.ma.mean(X[:, half:], 1) - np.ma.mean(X[:, :half], 1)

def maxcorr(P, x, ntry = 10000):
    '''Calculates the maximum theoretically possible correlation between x and P,
       possible given the noise (standard deviation) in P.

    Arguments: P, x, [ntry]
    Output: float
    Remarks: P and x are not commutable
    '''

    mu_P = np.ma.mean(np.ma.abs(P))
    mu_x = np.ma.mean(np.ma.abs(x))
    return np.ma.mean([np.ma.corrcoef(x, x.copy() + np.random.normal(0, np.ma.std(P) * mu_x / mu_P, len(P)))[0, 1]
        for tt in range(ntry)])

def permtest(X0, nmc = 100):
    '''Calculates the p-value using a Monte Carlo permutation test. This test
       is non-parametric, i.e. it does not assume an underlying distribution.

    Arguments: X (shape = (nbox, nmon)), [nmc]
    Output: an array of length nbox
    Remarks: this function operates on the second dimension (i.e. nmon)

    Credit: Adapted from http://stackoverflow.com/a/24801874.'''

    X = X0.copy()    # making a copy in case shuffle() messes the original
    diff = np.ma.abs(delta(X))

    k = np.ma.zeros(len(X))
    for j in range(nmc):
        np.random.shuffle(X.T)
        k += diff < np.ma.abs(delta(X))
    return k / nmc



"""
def plot_isccp_sig(Z, nmc = 1000, a = 0.1, size = 2,
    datapath = '/home/btan1/Monash/Data/'):
    '''Plots stippling on top of the ISCCP maps.

    Arguments: Z (array to test for significance), [a, size, datapath]
    Output: None
    '''

    import matplotlib.pyplot as plt

    mask = np.ma.getmaskarray(delta(Z))
    sig = np.ma.masked_where(mask, permtest(Z, nmc) < a)
    lat30 = np.where(np.abs(read_netcdf(datapath, 2007, variable = 'box_lat')) < 30)[0]
    lat = read_netcdf(datapath, 2007, variable = 'box_lat')[lat30][sig]
    lon = read_netcdf(datapath, 2007, variable = 'box_lon')[lat30][sig]
    plt.scatter(lon, lat, s = size, c = '#ffff00', alpha = 0.9, edgecolor = 'none')
    #plt.scatter(lon, lat, s = size, c = 'k', alpha = 0.6, edgecolor = 'none')
    #plt.scatter(lon, lat, s = 5, c = 'k', marker = 'x', alpha = 0.5, linewidths = 0.5)
    #plt.scatter(lon, lat, s = 5, c = '#ffff00', marker = 'x', alpha = 0.8, linewidths = 0.5)

    return None
"""
