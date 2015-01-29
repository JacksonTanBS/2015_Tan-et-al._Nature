#!/usr/bin/python3

# This Python script contains some regularly-used functions for the trends project.

import numpy as np
import calendar as cal
from netCDF4 import Dataset

def read_netcdf(datapath, year, month = 'all', variable = 'Regime9'):
    '''Extracts the variables from the modified regime netCDF files (daily).

    Arguments: read_netcdf(datapath, year, [month, variable])
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


def read_netcdf3hr(datapath, year, month = 'all', variable = 'Regime9'):
    '''Extracts the variables from the modified regime netCDF files (3hr).

    Arguments: read_netcdf3hr(datapath, year, [month, variable])
    Output: 2-dim array (dims: day, box) of type int16
    Requires: netCDF files of the modified regimes
    Remarks: 'variable' can be 'time', 'box', 'box_lat', 'box_lon', 'Regime9'
             'Regime8', 'Surface', 'SatelliteID', or 'SatelliteAngle'.
    '''

    f = Dataset('%strend/netcdf_3hr/regimes_%d.nc' % (datapath, year), 'r')
    if month == 'all':
        data = f.variables[variable][:].copy()
    elif year == 1983:
        day0 = (cal.datetime.datetime(year, month, 1).timetuple().tm_yday
            - cal.datetime.datetime(1983, 7, 1).timetuple().tm_yday) * 8
        day1 = day0 + cal.monthrange(year, month)[1] * 8
        data = f.variables[variable][day0 : day1].copy()
    else:
        day0 = (cal.datetime.datetime(year, month, 1).timetuple().tm_yday - 1) * 8
        day1 = day0 + cal.monthrange(year, month)[1] * 8
        data = f.variables[variable][day0 : day1].copy()
    f.close()

    return data


def isccp_latlon(datapath = '/Volumes/Shared/Data/'):
    '''Returns the lat/lon boundaries of each ISCCP equal-area grid box.

    Arguments: isccp_latlon([datapath]):
    Output: a latitude boundary array, a longitude boundary list, and the box
            where a new latitude begins
    Requires: netCDF files of the modified regimes
    Remarks: Each 'row' of the array/list correspond to each latitude band.
    '''

    lon = read_netcdf(datapath, 2007, variable = 'box_lon')

    lats = np.vstack([np.arange(-35, 35, 2.5), np.arange(-32.5, 35.1, 2.5)]).T
    jumps = np.hstack([[0], np.where(np.diff(lon) < 0)[0] + 1, len(lon)])
    lons = []
    for jj in range(len(jumps) - 1):
        # take the slice of longitude for that latitude
        lon_midpt = lon[jumps[jj] : jumps[jj + 1]]
        lon_edges = lon_midpt[1:] - 0.5 * np.diff(lon_midpt)
        lons.append(np.hstack([0., lon_edges, 360.]))

    return lats, lons, jumps


def regular_angle(datapath = '/Volumes/Shared/Data/', year = 2006, month = 1,
        day = 0):
    '''Returns a map of the cosine of the zenith angle.

    Arguments: regular_angle([datapath, year, month, day]):
    Output: a latitude boundary array, a longitude boundary list, and the box
            where a new latitude begins
    Requires: netCDF files of the modified regimes
    Remarks: Each 'row' of the array/list correspond to each latitude band.
    '''

    from scipy.interpolate import griddata

    ang = read_netcdf(datapath, 2006, 1, variable = 'SatelliteAngle')[0]
    lat = read_netcdf(datapath, 2006, 1, variable = 'box_lat')
    lon = read_netcdf(datapath, 2006, 1, variable = 'box_lon')
    x, y = np.mgrid[1.25 : 360 : 2.5, -33.75 : 33.8 : 2.5]
    pts3 = np.vstack([np.hstack([lon - 360, lon, lon + 360]),
                np.hstack([lat, lat, lat])]).T
    ang3 = np.hstack([ang, ang, ang])

    return griddata(pts3, ang3, (x, y), method = 'linear')


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

def equal_area():
    '''Produces elat and elon, vectors of lat/lon corresponding to ISCCP
    280 km equal-area grids between 35°N/S. By Jackson Tan, based on
    equal_area.pro (IDL) by John Haynes, which is in turn based on
    read_i2_FD_prf.f by Y. C. Zhang.

    Arguments: none:
    Output: an array for latitude and an array for longitude
    '''

    maxbox = 6596
    dlon = 2.5
    dlat = 2.5

    dlontb = np.zeros(72)
    ncells = np.zeros(72, dtype = int)
    xlatb = np.zeros(72)
    xlate = np.zeros(72)
    loneqr = np.zeros(maxbox, dtype = int)
    lateqr = np.zeros(maxbox, dtype = int)
    elon = np.zeros(maxbox)
    elat = np.zeros(maxbox)  
      
    lonmax = int(360.0 / dlon)
    latmax = int(180.0 / dlat)
    lonmin = 1
    latmin = 1
      
    twopi = 2. * np.pi
    twopir = twopi / 360.
    re = 6371.2
    rdlat = dlat * twopir
    hezon = re * np.sin(rdlat)
    aezon = twopi * re * hezon
    aecell = (aezon * dlon) / 360.

    for lat in range(latmin, latmax + 1):
        xlatb[lat - 1] = dlat * (lat - 1) - 90.0
        xlate[lat - 1] = xlatb[lat - 1] + dlat
        if xlate[lat - 1] > 90.: xlate[lat - 1] = 90.
        if xlatb[lat - 1] > 90.: xlatb[lat - 1] = 90.
        rlatb = twopir * xlatb[lat - 1]
        rlate = twopir * xlate[lat - 1]
        htb = re * np.sin(rlatb)
        hte = re * np.sin(rlate)
        htzone = np.abs(hte - htb)
        azone = twopi * re * htzone
        cells = azone / aecell
        ncells[lat - 1] = int(cells + 0.5)

    for lat in range(latmin, latmax + 1):
        if ncells[lat - 1] > 0:
            dlontb[lat - 1] = 360. / ncells[lat - 1]
        else: dlontb[lat - 1] = 360.
      
    kount = -1
    for j in range(0, 72):
        for i in range(0, ncells[j]):
            kount = kount + 1
            loneqr[kount] = i + 1
            lateqr[kount] = j + 1

    for ibx in range(0, maxbox):
        lat = lateqr[ibx]
        elat[ibx] = (float(lat) - 0.5) * 2.5 - 90.
        elon[ibx] = (float(loneqr[ibx]) - 0.5) * dlontb[lat - 1]

    return elat, elon
