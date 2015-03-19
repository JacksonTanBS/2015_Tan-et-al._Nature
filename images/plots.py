#!/usr/local/bin/python3

''' This script plots the diagrams in this paper. To produce a single plot, 
run "python plots.py [option]" or "python3 plots.py [option]" (depending on
your Python installation), with [option] being Fig1, Fig2, ..., EDFig1, ...
etc. You can produce multiple plots by, e.g. "python plots.py Fig1 Fig2", 
but be warned that they are run simutaneously and may consume more memory. 
Note that EDFig5 is produced separately using output from this script.

Except for EDFig4 which needs some calculations, all plots may be produced 
very quickly. The map diagrams by default take longer because of the 
calculation of the stippling through Monte Carlo permutation test. You can 
speed it up by reducing the iterations for this process (default = 10,000); 
to do so, specify the variable "nmc" in the function "plot_isccp_sig", e.g.
"plot_isccp_sig(X, a = alpha, nmc = 1000)".

Jackson Tan (jackson.tan@monash.edu)
29 January 2015

'''
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import csv
from funcs import *

# preliminaries
option = sys.argv[1:]
datapath = '/home/btan1/Monash/Data/'
regimes = ('CR1', 'CR2', 'CR3', 'CR4', 'CR5', 'CR6', 'CR7', 'CR8')
lat30 = np.where(np.abs(read_netcdf(datapath, 2007, variable = 'box_lat')) < 30)[0]
nbox = len(lat30)
subpltlab = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k')

# quick function to create LaTeX labels for titles
lab1 = lambda x: '$\Delta f_{%s}$ $\overline{p_{%s}}$' % (x, x)
lab2 = lambda x: '$\overline{f_{%s}}$ $\Delta p_{%s}$' % (x, x)
lab3 = lambda x: '%s $+$ %s' % (lab1(x), lab2(x))

#----- SETTINGS -----#

# customise some plot settings
plt.rcParams['font.size'] = 8
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['axes.titlesize'] = 'medium'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.sans-serif'] = ['TeX Gyre Heros', 'Helvetica',
                                   'Bitstream Vera Sans']
plt.rcParams['pdf.fonttype'] = 42

# plot settings according to Nature guidelines
scol = 3.503    # single column (89 mm)
dcol = 7.204    # double column (183 mm)
flpg = 9.724    # full page length (247 mm)

#----- PLOTS FOR MAIN TEXT -----#

if 'Fig1' in option:

    from matplotlib.colors import BoundaryNorm

    C = np.loadtxt('%strend/centroids_8.txt' % datapath)

    bound = (0, 0.2, 1, 2, 3, 4, 6, 8, 10, 15, 20)
    x, y = np.meshgrid(np.arange(0.5, 7.5), np.arange(7.5, 0, -1))
    xlabels = ('0.02', '1.3', '3.6', '9.4', '23', '60', '1000')
    ylabels = ('1000', '800', '680', '560', '440', '310', '180', '10')

    plt.figure(figsize = (scol, 0.75 * scol))
    plt.subplots_adjust(bottom = 0.25, wspace = 0.15, hspace = 0.2)

    for rr in range(3):
        plt.subplot(221 + rr + (rr > 0))
        plt.pcolor(x, y, C[rr].reshape(7, 6), vmin = 0, vmax = 20, 
                   norm = BoundaryNorm(bound, 222), cmap = plt.cm.hot_r)
        if rr == 2: plt.text(-0.08, 1.05, subpltlab[rr], fontweight = 'bold', 
                             ha = 'center', va = 'center', transform = plt.gca().transAxes)
        else: plt.text(-0.4, 1.05, subpltlab[rr], fontweight = 'bold', 
                       ha = 'center', va = 'center', transform = plt.gca().transAxes)
        plt.setp(plt.gca(), 'yticks', np.arange(0.5, 8))
        plt.setp(plt.gca(), 'yticklabels', '')
        plt.axis('tight')
        if rr > 0:
            plt.setp(plt.gca(), 'xticks', np.arange(0.5, 7))
            plt.setp(plt.gca(), 'xticklabels', xlabels)
            for tick in plt.gca().xaxis.get_major_ticks():
                tick.label.set_rotation(45)
            plt.gca().xaxis.set_label_coords(0.5, -0.35)
            plt.xlabel('τ')
        else:
            plt.setp(plt.gca(), 'xticks', np.arange(0.5, 7))
            plt.setp(plt.gca(), 'xticklabels', '')
        if rr < 2:
            plt.setp(plt.gca(), 'yticks', np.arange(0.5, 8))
            plt.setp(plt.gca(), 'yticklabels', ylabels)
            plt.ylabel('CTP (hPa)')
        else:
            plt.setp(plt.gca(), 'yticks', np.arange(0.5, 8))
            plt.setp(plt.gca(), 'yticklabels', '')
        plt.axis('tight')
        plt.text(5.5, 1.5, regimes[rr], ha = 'center', va = 'center')

    cb = plt.colorbar(cax = plt.axes([0.125, 0.03, 0.775, 0.02]), 
                      orientation = 'horizontal', ticks = bound)
    cb.set_label('(%)')

    plt.savefig('Fig1.pdf', bbox_inches = 'tight')
    plt.close()

    np.savetxt('Fig1.txt', C[:3], delimiter = ',')

    with open('Fig1.txt', 'r') as f:
        with open('Fig1.csv', 'w') as g:
            g.write('"Source Data for Fig. 1"\n')
            g.write('"Centroids for CR1, CR2, CR3; values run from low τ to high τ, then low CTP to high CTP."\n')
            g.write(f.read())
    os.remove('Fig1.txt')

if 'Fig2' in option:

    import calendar as cal
    from statsmodels.graphics.boxplots import violinplot

    precippath = '%strend/precip_data/' % datapath
    years = list(range(1998, 2010))
    regimes = ('CR1', 'CR2', 'CR3')
    nreg = len(regimes)

    T = np.arange(years[0], years[-1] + 1, 1. / 12)
    Pall = [[] for rr in range(nreg + 1)]   # all precipitation values in regime
    P = np.zeros([nreg + 1, len(T)])        # monthly-mean precipitation
    F = np.zeros([nreg + 1, len(T)])        # monthly-mean FOC
    Ptotal = 0                              # total amount of precipitation

    tt = 0
    for yy, year in enumerate(years):

        precip = np.load('%strmmdaytime_%d.npy' % (precippath, year))[:, lat30]
        R = read_netcdf(datapath, year, variable = 'Regime8')[:, lat30]
        Ptotal += np.sum(precip)

        # calculate the monthly mean precipitation and regime FOC
        t0 = 0
        for month in range(1, 13):
            t1 = t0 + cal.monthrange(year, month)[1]
            for rr in range(nreg):
                P[rr, tt] = np.ma.mean([ii for ii in precip[t0 : t1][R[t0 : t1] == rr + 1]])
                F[rr, tt] = np.ma.sum(R[t0 : t1] == rr + 1) / np.float(np.ma.count(R[t0 : t1]))
                Pall[rr] += [ii for ii in precip[t0 : t1][R[t0 : t1] == rr + 1] if not np.ma.is_masked(ii)]
            P[nreg, tt] = np.ma.mean([ii for ii in precip[t0 : t1][R[t0 : t1] < 4]])
            F[nreg, tt] = np.ma.sum(R[t0 : t1] < 4) / np.float(np.ma.count(R[t0 : t1]))
            Pall[nreg] += [ii for ii in precip[t0 : t1][R[t0 : t1] < 4] if not np.ma.is_masked(ii)]

            t0 = t1
            tt += 1

    Pcon = [np.sum(Pall[rr]) / Ptotal for rr in range(nreg + 1)]    # fractional contribution

    x = np.arange(nreg + 1)

    plt.figure(figsize = (scol, 0.5 * scol))
    violinplot(Pall, ax = plt.gca(), positions = x, show_boxplot = False,
        plot_opts = {'violin_fc': '0.2', 'violin_ec': 'k', 'cutoff': True})
    bp = plt.boxplot(Pall, sym = '', widths = 0.4, positions = x, whis = [5, 95])
    plt.setp(bp['boxes'], color = 'black')
    plt.setp(bp['whiskers'], color = 'black', linestyle = '-')
    plt.setp(bp['medians'], color = 'black')
    plt.scatter(x, [np.mean(ii) for ii in Pall], s = 10, c = 'k')
    for rr in range(nreg + 1): plt.text(x[rr], 105, '{:2.0f}%'.format(Pcon[rr] * 100), ha = 'center')
    plt.ylim([0, 100])
    plt.setp(plt.gca(), 'xticks', x)
    plt.setp(plt.gca(), 'xticklabels', regimes + ('CR1-3',))
    plt.ylabel('$P$ (mm / day)')
    plt.grid(axis = 'y')
    plt.axvline(2.5, c = 'k', ls = ':', lw = 0.5)
    plt.savefig('Fig2.pdf')
    plt.close()

    # print to CSV file at %6.2f precision
    with open('Fig2.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['Source Data for Fig. 2'])
        w.writerow(['Precipitation values for CR1, CR2, CR3 and CR1-3'])
        w.writerow(['Note: the raw data here is sub-sampled by 50% to fit ' \
                    'into the 30 MB limit. Use the code to generate the full ' \
                    'set of raw data.'])
        w.writerows([['%6.2f' % p for p in P[::2]] for P in Pall])

if 'Fig3' in option:

    inpath1 = '%strend/P_reconstruct/' % datapath
    inpath2 = '%strend/precip_data/' % datapath
    nreg = len(regimes)
    alpha = 0.1

    Fr = np.load('%strmmdaytime_data_F_r.npy' % inpath1)[:, lat30]       # FOC over TRMM period
    Pr = np.load('%strmmdaytime_data_P_r.npy' % inpath1)[:, lat30]       # within-regime precip from TRMM
    Psum = np.load('%strmmdaytime_data_P_sum.npy' % inpath1)[lat30]      # monthly precip from TRMM
    Fsa = np.load('%strmm3a25saf.npy' % (inpath2,))[:, lat30].T          # stratiform area fraction
    Fsr = np.load('%strmm3a25srf.npy' % (inpath2,))[:, lat30].T          # stratiform rain fraction
    Fc = np.load('%strmm3a25caf.npy' % (inpath2,))[:, lat30].T           # convective area fraction
    F13 = np.load('%strmmdaytime_data_F(1-3)_r.npy' % inpath1)[:, lat30] # FOC of CR1-3
    P13 = np.load('%strmmdaytime_data_P(1-3)_r.npy' % inpath1)[:, lat30] # within-regime precip of CR1-3

    # calculate the time-means
    meanPr = np.ma.mean(Pr, 2)
    meanFr = np.ma.mean(Fr, 2)
    meanF13 = np.ma.mean(F13, 2)
    meanP13 = np.ma.mean(P13, 2)

    # calculate the terms
    term1 = Fr * meanPr[:, :, None]
    term2 = meanFr[:, :, None] * Pr
    term13 = np.ma.sum([F13 * meanP13[:, :, None], meanF13[:, :, None] * P13], 0)

    #vmax = max(np.percentile(DeltaPsum, 99), -np.percentile(DeltaPsum, 1))
    vmax = 1.6

    variables = (Psum, term1[0], term2[0], np.ma.sum(np.ma.concatenate([term1[1:3], term2[1:3]], 0), 0),
                 np.ma.sum(np.ma.concatenate([term1[:3], term2[:3]], 0), 0), Fsa, Fsr)
    titles = ('$\Delta P$', lab1('1'), lab2('1'), '%s $+$ %s' % (lab3('2'), lab3('3')),
              '$\Sigma_{i = 1}^3$ %s' % lab3('i'), '$\Delta f_{sa}$', '$\Delta f_{sr}$')
    nvar = len(variables)

    plt.figure(figsize = (dcol, nvar * 0.2 * dcol))
    plt.subplots_adjust(right = 0.85)
    for rr in range(nvar - 2):
        plt.subplot(nvar, 1, rr + 1)
        mp = plot_isccp_map(delta(variables[rr])[None], vmax = vmax, meridian = False)
        plot_isccp_sig(variables[rr], a = alpha)
        plt.title(titles[rr], loc = 'left')
        plt.text(-0.1, 1.1, subpltlab[rr], fontweight = 'bold', transform = plt.gca().transAxes)

    cb = plt.colorbar(mappable =  mp, cax = plt.axes([0.875, 0.34, 0.01, 0.555]),
        orientation = 'vertical', ticks = np.linspace(-vmax, vmax, 9), extend = 'both', extendfrac = 0.1 / 5)
    cb.set_label('mm / day')

    for rr in range(5, 7):
        plt.subplot(nvar, 1, rr + 1)
        mp = plot_isccp_map(delta(variables[rr])[None], vmax = 0.1, meridian = (rr == 6))
        plot_isccp_sig(variables[rr], a = alpha)
        plt.title(titles[rr], loc = 'left')
        plt.text(-0.1, 1.1, subpltlab[rr], fontweight = 'bold', transform = plt.gca().transAxes)
    plt.colorbar(mappable = mp, cax = plt.axes([0.875, 0.105, 0.01, 0.2]),
        orientation = 'vertical', ticks = np.linspace(-0.1, 0.1, 5), extend = 'both', extendfrac = 0.1 / 2)

    plt.savefig('Fig3.pdf')
    plt.close()

    # function to calculate RMSE
    rmse = lambda x, y : np.ma.sqrt(np.ma.mean((x - y) ** 2))

    # write to text file the correlations and RMSEs
    with open('tab_Fig3.txt', 'w') as f:
        f.write('Corr. (upper-right) and RMSE (lower-left) of plots in Fig. 3\n\n')
        for v1, var1 in enumerate(variables):
            if v1 == 0:
                f.write(''.join(['%7s' % ''] + ['%7.3s' % spl for spl in subpltlab[:len(variables)]]
                    + ['\n']))
            for v2, var2 in enumerate(variables):
                if v2 == 0:
                    f.write('%7.3s' % subpltlab[v1])
                if v1 > v2:
                    f.write('%7.3f' % rmse(delta(var1), delta(var2)))
                elif v1 < v2:
                    f.write('%7.3f' % np.ma.corrcoef(delta(var1), delta(var2))[0, 1])
                else:
                    f.write('%7s' % '-----')
            f.write('\n')

        f.write('\nMaximum theoretically possible correlations:\n\n')
        for vv, var in enumerate(variables):
            if vv == 0: f.write('%7.3s' % subpltlab[0])
            f.write('%7.3f' % maxcorr(variables[0], delta(var)))

    np.savetxt('Fig3.txt', [delta(var) for var in variables], delimiter = ',', fmt = '%6.3f')
    with open('Fig3.txt', 'r') as f:
        with open('Fig3.csv', 'w') as g:
            g.write('"Source Data for Fig. 3"\n')
            g.write('"Each line is one panel; each value is the change in each grid box ' \
                    'of the ISCCP equal-area grid between 30°N/S."\n')
            g.write(f.read())
    os.remove('Fig3.txt')
    

if 'Fig4' in option:

    inpath1 = '%strend/P_reconstruct/' % datapath
    inpath2 = '%strend/precip_data/' % datapath
    nreg = len(regimes)
    alpha = 0.1

    Fr = np.load('%sgpcp2.2_data_F_r.npy' % inpath1)[:, lat30]      # FOC over GPCP period
    Pr = np.load('%strmmdaytime_data_P_r.npy' % inpath1)[:, lat30]  # within-regime precip from TRMM
    Psum = np.load('%sgpcp2.2.npy' % inpath2).T[lat30]              # monthly precip from TRMM

    # calculate the terms
    meanPr = np.ma.mean(Pr, 2)
    term = Fr * meanPr[:, :, None]

    #vmax = max(np.percentile(DeltaPsum, 99), -np.percentile(DeltaPsum, 1))
    vmax = 1.2

    variables = (Psum, term[0], np.ma.sum(term[:3], 0))
    titles = ('$\Delta P$', lab1('1'), '$\Sigma_{i = 1}^3$ %s' % lab1('i'))
    nvar = len(variables)

    plt.figure(figsize = (dcol, 0.2 * nvar * dcol))
    plt.subplots_adjust(right = 0.85)
    for rr in range(nvar):
        plt.subplot(nvar, 1, rr + 1)
        mp = plot_isccp_map(delta(variables[rr])[None], vmax = vmax, meridian = (rr == 2))
        plot_isccp_sig(variables[rr], a = alpha)
        plt.title(titles[rr], loc = 'left')
        plt.text(-0.1, 1.1, subpltlab[rr], fontweight = 'bold', transform = plt.gca().transAxes)

    cb = plt.colorbar(mappable = mp, cax = plt.axes([0.875, 0.115, 0.01, 0.76]),
        orientation = 'vertical', ticks = np.linspace(-vmax, vmax, 5), extend = 'both', extendfrac = 0.1 / nvar)
    cb.set_label('mm / day')

    plt.savefig('Fig4.pdf')
    plt.close()

    # function to calculate RMSE
    rmse = lambda x, y : np.ma.sqrt(np.ma.mean((x - y) ** 2))

    # write to text file the correlations and RMSEs
    with open('tab_Fig4.txt', 'w') as f:
        f.write('Corr. (upper-right) and RMSE (lower-left) of plots in Fig. 4\n\n')
        for v1, var1 in enumerate(variables):
            if v1 == 0:
                f.write(''.join(['%7s' % ''] + ['%7.3s' % spl for spl in subpltlab[:len(variables)]]
                    + ['\n']))
            for v2, var2 in enumerate(variables):
                if v2 == 0:
                    f.write('%7.3s' % subpltlab[v1])
                if v1 > v2:
                    f.write('%7.3f' % rmse(delta(var1), delta(var2)))
                elif v1 < v2:
                    f.write('%7.3f' % np.ma.corrcoef(delta(var1), delta(var2))[0, 1])
                else:
                    f.write('%7s' % '-----')
            f.write('\n')

        f.write('\nMaximum theoretically possible correlations:\n\n')
        for vv, var in enumerate(variables):
            if vv == 0: f.write('%7.3s' % subpltlab[0])
            f.write('%7.3f' % maxcorr(variables[0], delta(var)))

    np.savetxt('Fig4.txt', [delta(var) for var in variables], delimiter = ',', fmt = '%6.3f')
    with open('Fig4.txt', 'r') as f:
        with open('Fig4.csv', 'w') as g:
            g.write('"Source Data for Fig. 4"\n')
            g.write('"Each line is one panel; each value is the change in each grid box ' \
                    'of the ISCCP equal-area grid between 30°N/S."\n')
            g.write(f.read())
    os.remove('Fig4.txt')


#----- SUPPLEMENTARY PLOTS -----#

if 'EDFig1' in option:

    from mpl_toolkits.basemap import Basemap
    from matplotlib.colors import BoundaryNorm, ListedColormap

    inpath = '%strend/P_reconstruct/' % datapath
    ntmon = 318

    F = np.load('%sgpcp2.2_data_F_r.npy' % inpath, 'r')[:, lat30]
    geog = np.ma.mean(F, 2)

    bounds = np.array([0] + list(np.logspace(-2, 0, 64)))
    lats, lons, jumps = isccp_latlon()

    # re-package the spectral colormap with a white at the start
    cmap = ListedColormap([(1.0, 1.0, 1.0, 1.0)] + [plt.cm.GnBu((ii + 1) /
        np.float(len(bounds))) for ii in range(len(bounds) - 1)])

    plt.figure(figsize = (scol, 0.75 * scol))
    plt.subplots_adjust(hspace = 0.3)
    for rr in range(3):
        plt.subplot(311 + rr)
        for jj in range(len(jumps) - 1):
            mp = plt.pcolor(lons[jj], lats[jj] + (-1, 1), geog[rr : rr + 1, jumps[jj] : jumps[jj + 1]],
                norm = BoundaryNorm(bounds, len(bounds)), cmap = cmap, rasterized = True)
        m = Basemap(projection='cyl', llcrnrlat = -30, urcrnrlat = 30,
            llcrnrlon = 0, urcrnrlon = 360, fix_aspect = True,
            resolution = 'c')
        m.drawcoastlines(color = 'k', linewidth = 0.5)
        m.drawparallels(np.arange(-30, 31, 15), color = 'k',
            linewidth = 0, labels = [1, 0, 0, 0])
        m.drawmeridians(np.arange(0.,361.,60.), color = 'k',
            linewidth = 0, labels = [0, 0, 0, 1])
        plt.text(-0.15, 1.2, subpltlab[rr], fontweight = 'bold', transform = plt.gca().transAxes)
        plt.title(regimes[rr], loc = 'left')
    plt.colorbar(cax = plt.axes([0.125, 0.02, 0.775, 0.015]), orientation = 'horizontal',
        spacing = 'proportional', ticks = np.linspace(0, 1, 11))
    plt.savefig('EDFig1.jpg', dpi = 300)
    plt.close()

    np.savetxt('EDFig1.txt', geog[:3], delimiter = ',', fmt = '%6.3f')
    with open('EDFig1.txt', 'r') as f:
        with open('EDFig1.csv', 'w') as g:
            g.write('"Source Data for Extended Data Fig. 1"\n')
            g.write('"Each line is one panel; each value is the average frequency ' \
                    'in each grid box of the ISCCP equal-area grid between 30°N/S."\n')
            g.write(f.read())
    os.remove('EDFig1.txt')

if 'EDFig2' in option:

    nreg = len(regimes)
    nmon = 318

    focs = np.ma.masked_all([nreg, nmon])
    foc_conv = np.ma.masked_all(nmon)
    tt = 0

    for year in range(1983, 2010):

        # define the months for the year
        if year == 1983: months = list(range(7, 13))
        else: months = list(range(1, 13))

        for month in months:

            R = read_netcdf(datapath, year, month, variable = 'Regime8')[:, lat30]

            if np.ma.count(R) > 5:
                count = [np.sum(R == rr + 1) for rr in range(nreg)]
                total = np.float(np.ma.count(R))
                focs[:, tt] = np.array(count) / total
                foc_conv[tt] = np.sum(count[:3]) / total

            tt += 1

    def deseason(foc):
        seasonmean = [np.ma.mean(foc[6 + mm : : 12]) for mm in range(12)]
        anom = np.array(seasonmean[6:] + seasonmean * 26) - np.mean(foc[6:])
        return foc - anom

    T = np.arange(1983 + 6 / 12., 2010., 1 / 12.)
    yticks = (np.linspace(0.03, 0.07, 3), np.linspace(0.06, 0.12, 3),
              np.linspace(0.08, 0.18, 3), np.linspace(0.22, 0.34, 3))

    plt.figure(figsize = (scol, 1.2 * scol))
    plt.subplots_adjust(hspace = 0.4)
    for rr, foc_raw in enumerate(focs[:3]):
        plt.subplot(411 + rr)
        foc = deseason(foc_raw)
        a, b, c, w = fitLine(T[~np.isnan(foc)], np.array(foc)[~np.isnan(foc)])
        print(a * len(T) / 12, w[0] * len(T) / 12)
        plt.plot(T, foc, c = 'k', lw = 0.5)
        plt.plot(T, a * T + b, 'k--')
        plt.xlim([1983, 2010])
        plt.setp(plt.gca(), 'xticklabels', '')
        plt.yticks(yticks[rr])
        plt.ylabel('$f_{%s}$' % regimes[rr][2:])
        plt.gca().yaxis.set_label_coords(-0.13, 0.5)
        plt.text(-0.22, 1.06, subpltlab[rr], fontweight = 'bold', 
            ha = 'center', va = 'center', transform = plt.gca().transAxes)
    plt.subplot(414)
    foc = deseason(foc_conv)
    a, b, c, w = fitLine(T[~np.isnan(foc)], np.array(foc)[~np.isnan(foc)])
    print(a * len(T) / 12, w[0] * len(T) / 12)
    plt.plot(T, foc, c = 'k', lw = 0.5)
    plt.plot(T, a * T + b, 'k--')
    plt.xlim([1983, 2010])
    plt.yticks(yticks[3])
    plt.ylabel('$f_{1-3}$')
    plt.gca().yaxis.set_label_coords(-0.13, 0.5)
    plt.text(-0.22, 1.06, 'd', fontweight = 'bold', 
        ha = 'center', va = 'center', transform = plt.gca().transAxes)

    plt.savefig('EDFig2.jpg', dpi = 300)
    plt.close()

    np.savetxt('EDFig2.txt', [deseason(foc) for foc in np.vstack([focs[:3], foc_conv])],
               delimiter = ',', fmt = '%6.3f')
    with open('EDFig2.txt', 'r') as f:
        with open('EDFig2.csv', 'w') as g:
            g.write('"Source Data for Extended Data Fig. 2"\n')
            g.write('"Each line is the FOC for each panel at each month, from Jul 1983 ' \
                    'to Dec 2009"\n')
            g.write(f.read())
    os.remove('EDFig2.txt')

    # calculate difference of means and its standard deviation (at ~95% C.I.)
    for foc_raw in focs[:3]:
        foc = deseason(foc_raw)
        diffmean = np.ma.mean(foc[(nmon / 2):]) - np.ma.mean(foc[:(nmon / 2)])
        diffstd = np.ma.sqrt((np.std(foc[(nmon / 2):]) ** 2 + np.std(foc[:(nmon / 2)]) ** 2) / (nmon / 2))
        print(diffmean, 2 * diffstd)
    foc = deseason(foc_conv)
    diffmean = np.ma.mean(foc[(nmon / 2):]) - np.ma.mean(foc[:(nmon / 2)])
    diffstd = np.ma.sqrt((np.std(foc[(nmon / 2):]) ** 2 + np.std(foc[:(nmon / 2)]) ** 2) / (nmon / 2))
    print(diffmean, 2 * diffstd)

if 'EDFig4' in option:

    inpath1 = '%strend/P_reconstruct/' % datapath
    inpath2 = '%strend/precip_data/' % datapath
    nreg = len(regimes)
    alpha = 0.1

    Fr = np.load('%sgpcp1dd_data_F_r.npy' % inpath1)[:, lat30]   # FOC over GPCP 1DD period
    Pr = np.load('%sgpcp1dd_data_P_r.npy' % inpath1)[:, lat30]   # within-regime precip from GPCP 1DD
    Psum = np.load('%sgpcp1dd_data_P_sum.npy' % inpath1)[lat30]  # monthly precip from GPCP 1DD

    # calculate the time-means
    meanPr = np.ma.mean(Pr, 2)
    meanFr = np.ma.mean(Fr, 2)

    # calculate the terms
    term1 = Fr * meanPr[:, :, None]
    term2 = meanFr[:, :, None] * Pr

    #vmax = max(np.percentile(DeltaPsum, 99), -np.percentile(DeltaPsum, 1))
    vmax = 1.6

    variables = (Psum, term1[0], term2[0], np.ma.sum(np.ma.concatenate([term1[1:3], term2[1:3]], 0), 0),
                 np.ma.sum(np.ma.concatenate([term1[:3], term2[:3]], 0), 0))
    titles = ('$\Delta P$', lab1('1'), lab2('1'), '%s $+$ %s' % (lab3('2'), lab3('3')),
              '$\Sigma_{i = 1}^3$ %s' % lab3('i'))
    nvar = len(variables)

    plt.figure(figsize = (dcol, 0.2 * nvar * dcol))
    plt.subplots_adjust(right = 0.85)
    for rr in range(nvar):
        plt.subplot(nvar, 1, rr + 1)
        mp = plot_isccp_map(delta(variables[rr])[None], vmax = vmax, meridian = (rr == 4))
        plot_isccp_sig(variables[rr], a = alpha)
        plt.title(titles[rr], loc = 'left')
        plt.text(-0.1, 1.1, subpltlab[rr], fontweight = 'bold', transform = plt.gca().transAxes)
    cb = plt.colorbar(mappable = mp, cax = plt.axes([0.875, 0.11, 0.01, 0.775]),
        orientation = 'vertical', ticks = np.linspace(-vmax, vmax, 5), extend = 'both', extendfrac = 0.1 / nvar)
    cb.set_label('mm / day')
    plt.savefig('EDFig4.jpg', dpi = 300)
    plt.close()

    # function to calculate RMSE
    rmse = lambda x, y : np.ma.sqrt(np.ma.mean((x - y) ** 2))

    # write to text file the correlations and RMSEs
    with open('tab_EDFig4.txt', 'w') as f:
        f.write('Corr. (upper-right) and RMSE (lower-left) of plots in ED Fig. 4\n\n')
        for v1, var1 in enumerate(variables):
            if v1 == 0:
                f.write(''.join(['%7s' % ''] + ['%7.3s' % spl for spl in subpltlab[:len(variables)]]
                    + ['\n']))
            for v2, var2 in enumerate(variables):
                if v2 == 0:
                    f.write('%7.3s' % subpltlab[v1])
                if v1 > v2:
                    f.write('%7.3f' % rmse(delta(var1), delta(var2)))
                elif v1 < v2:
                    f.write('%7.3f' % np.ma.corrcoef(delta(var1), delta(var2))[0, 1])
                else:
                    f.write('%7s' % '-----')
            f.write('\n')

    np.savetxt('EDFig4.txt', [delta(var) for var in variables], delimiter = ',', fmt = '%6.3f')
    with open('EDFig4.txt', 'r') as f:
        with open('EDFig4.csv', 'w') as g:
            g.write('"Source Data for Extended Data Fig. 4"\n')
            g.write('"Each line is one panel; each value is the change in each grid box ' \
                    'of the ISCCP equal-area grid between 30°N/S."\n')
            g.write(f.read())
    os.remove('EDFig4.txt')

if 'EDFig5' in option:

    import calendar as cal

    w_r = np.zeros([8, nbox, 144])    # within-regime omega
    f_r = np.zeros([8, nbox, 144])    # regime frequencies
    w_sum = np.zeros([nbox, 144])     # monthly omega

    tt = 0
    for year in range(1998, 2010):

        w500 = np.load('%strend/precip_data/w500_%d.npy' % (datapath, year))

        t0 = 0
        for month in range(1, 13):

            R = read_netcdf(datapath, year, month, variable = 'Regime8')[:, lat30]
            t1 = t0 + cal.monthrange(year, month)[1]
            w_sum[:, tt] = np.ma.masked_invalid([np.ma.mean(w500[t0 : t1, box]) for box in range(nbox)])

            for rr in range(8):

                window = R == rr + 1
                count = np.ma.sum(R == rr + 1, 0)
                total = np.float32(np.ma.count(R, 0))
                f_r[rr, :, tt] = np.ma.masked_where(total == 0, count / total)

                for box in range(nbox):
                    if np.ma.sum(window[:, box]) != 0:
                        w_r[rr, box, tt] = np.ma.mean(w500[t0 : t1, box][window[:, box]])

            t0 = t1
            tt += 1

    # calculate the terms
    meanWr = np.ma.mean(w_r, 2)
    meanFr = np.ma.mean(f_r, 2)
    term1 = f_r * meanWr[:, :, None]
    term2 = meanFr[:, :, None] * w_r

    alpha = 0.1

    labw1 = lambda x: '$\Delta f_{%s}$ $\overline{w_{%s}}$' % (x, x)
    labw2 = lambda x: '$\overline{f_{%s}}$ $\Delta w_{%s}$' % (x, x)
    labw3 = lambda x: '%s $+$ %s' % (labw1(x), labw2(x))

    variables = (w_sum, f_r[0])
    titles = ('$\Delta \omega$', '$\Delta f_1$')
    nvar = len(variables)

    plt.figure(figsize = (dcol, 0.2 * nvar * dcol))
    plt.subplots_adjust(right = 0.85)
    plt.subplot(nvar, 1, 1)
    mp = plot_isccp_map(delta(variables[0])[None], vmax = 0.02, meridian = False)
    plot_isccp_sig(variables[0], a = alpha)
    plt.title(titles[0], loc = 'left')
    plt.text(-0.1, 1.1, subpltlab[0], fontweight = 'bold', transform = plt.gca().transAxes)
    cb = plt.colorbar(mappable = mp, cax = plt.axes([0.875, 0.57, 0.01, 0.3]),
        orientation = 'vertical', ticks = np.linspace(-0.02, 0.02, 5), extend = 'both', extendfrac = 0.2 / nvar)
    cb.set_label('Pa / s')
    plt.subplot(nvar, 1, 2)
    mp = plot_isccp_map(delta(variables[1])[None], vmax = 0.05, meridian = False)
    plot_isccp_sig(variables[1], a = alpha)
    plt.title(titles[1], loc = 'left')
    plt.text(-0.1, 1.1, subpltlab[1], fontweight = 'bold', transform = plt.gca().transAxes)
    cb = plt.colorbar(mappable = mp, cax = plt.axes([0.875, 0.125, 0.01, 0.3]),
        orientation = 'vertical', ticks = np.linspace(-0.05, 0.05, 5), extend = 'both', extendfrac = 0.2 / nvar)

    plt.savefig('EDFig5.jpg', dpi = 300)
    plt.close()

    np.savetxt('EDFig5.txt', [delta(var) for var in variables], delimiter = ',', fmt = '%6.3f')
    with open('EDFig5.txt', 'r') as f:
        with open('EDFig5.csv', 'w') as g:
            g.write('"Source Data for Extended Data Fig. 5"\n')
            g.write('"Each line is one panel; each value is the change in each grid box ' \
                    'of the ISCCP equal-area grid between 30°N/S."\n')
            g.write(f.read())
    os.remove('EDFig5.txt')

if 'EDFig3' in option:

    '''Note: this option only produces the CSV files. To produce the table 
       image, use the spreadsheet to produce a PDF, then use the ImageMagick 
       command "convert -density 300 -trim EDFig3.pdf EDFig3.jpg" to convert
       it to jpeg.'''

    data = []

    with open('tab_Fig3.txt', 'r') as f:
        Fig3 = [ii.split() for ii in f.readlines()]

    with open('tab_Fig4.txt', 'r') as f:
        Fig4 = [ii.split() for ii in f.readlines()]

    with open('tab_EDFig4.txt', 'r') as f:
        EDFig4 = [ii.split() for ii in f.readlines()]

    data.append([float(ii) for ii in Fig3[3][2:]])
    data.append([float(ii[1]) for ii in Fig3[4 : 8]])
    data.append([float(ii) for ii in Fig4[3][2:]])
    data.append([float(ii[1]) for ii in Fig4[4 : 6]])
    data.append([float(ii) for ii in EDFig4[3][2:]])
    data.append([float(ii[1]) for ii in EDFig4[4 : 8]])

    # print to CSV file at %6.2f precision
    with open('EDFig3.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['"Source Data for Extended Data Fig. 3"'])
        w.writerow(['"Correlations and RMSEs for Figs. 3, 4 and Extended Fig. 4"'])
        w.writerows(data)
