#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2022-04-20

Copyright (c) 2022 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""

import numpy as np
import time
import os
import sys
import pdb
import glob
import datetime
import h5py  # @UnresolvedImport
from matplotlib import pyplot as plt
import cartopy.crs as ccrs

# from analyseTraczilla import getConvfiles, I_HIT, I_OLD, readCatalogFile, getCatalogFile

sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readidx107

dateWithNans = {'20200111': '20200113_003758', \
                '20220102': '20220101_080814', \
                '20220103': '20220101_080814', \
                '20220104': '20220107_092916', \
                '20220116': '20220115_123549', \
                '20220117': '20220115_123549'}


if __name__ == '__main__':
    mainDIr = '/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout'
    # mainDIr = '/home/ejohansson/Projects/LMD-Traczilla'
    mainDIr = mainDIr + '/Coord'
    
    outputHour = 1
    if outputHour != 3:
        mainDIr = '%s/%d' %(mainDIr, outputHour * 3600)
    
    h5Dir = os.path.join(mainDIr, 'H5')
    plotDir = os.path.join(mainDIr, 'Plots')
    
    if not os.path.isdir(h5Dir):
        os.makedirs(h5Dir)
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)
    
    dirsInOut = os.listdir(mainDIr)
    days = []
    for dio in dirsInOut:
        if not 'Coord-' in dio:
            continue
        days.append(int(dio.split('Coord-')[-1]))
    days.sort()
    # days = [20200113]#, 20200113]
    # days = [20191216]
    ticT = time.time()
    d = 0
    for day in days:
        
        if int(str(day)[0:4]) in [2020]:#, 2021]:
            continue
        d = d + 1
        print(day)
        print('%d - Day %d of %d' %(day, d, len(days)))
        outname = 'Coord-%s' %(day)
        partDir = os.path.join(mainDIr, outname)
        #: Read part_000 file
        f0 = os.path.join(partDir, 'part_000')
        part0 = readidx107(f0, quiet=True)
        inDays = [*set(part0['ir_start'])]
        nrInDays = len(inDays)
        lenInDays = int(part0['numpart'] / nrInDays)
        pfiles = glob.glob(partDir + '/part_*')
        pfi = []
        for pf in pfiles:
            pfi.append(int(pf.split('part_')[-1]))
        pind = np.argsort(pfi)
        pfiles = np.asarray(pfiles)[pind][1:]
        
        lons = np.zeros([part0['numpart'], len(pfiles)]) + np.nan
        lats = np.zeros([part0['numpart'], len(pfiles)]) + np.nan
        ps = np.zeros([part0['numpart'], len(pfiles)]) + np.nan
        ts = np.zeros([part0['numpart'], len(pfiles)]) + np.nan
        # tids = (np.zeros([part0['numpart'], len(pfiles)]) + np.nan).astype(str)
        outptime = []
        itime = []
        i = -1
        for pf in pfiles:
            partf = readidx107(pf, quiet=True)
            if len(partf['x']) == 0:
                continue
            i = i + 1
            itime.append(partf['itime'])
            outptime.append((datetime.datetime.strptime(str(part0['stamp_date']), "%Y%m%d%H%M%S") - datetime.timedelta(seconds=abs(partf['itime']))).strftime("%Y%m%d_%H%M%S"))
            # tidstr = [int((datums2D[i] - originDate).total_seconds()) for d in ]
            inds = partf['idx_back']-partf['idx_orgn']
            lons[inds, i] = partf['x']
            lats[inds, i] = partf['y']
            ps[inds, i] = partf['p']
            ts[inds, i] = partf['t']
        lons = lons[:, 0:i+1]
        lats = lats[:, 0:i+1]
        ps = ps[:, 0:i+1]
        ts = ts[:, 0:i+1]
        # outptime = np.asarray(outptime)
        # itime = np.asarray(itime)
        h5name = os.path.join(h5Dir, 'coord-%i' %day)
        print('Create H5')
        print(h5name)
        h5file = h5py.File(h5name + '.h5', 'w')
        grp = h5file.create_group('Initial_Values')
        grp.create_dataset('stamp_date', data = (part0['stamp_date'], ))
        grp.create_dataset('Internal-Step', data = (part0['step'], ))
        grp.create_dataset('Output-Step', data = ((datetime.datetime.strptime(outptime[0], "%Y%m%d_%H%M%S") - datetime.datetime.strptime(outptime[1], "%Y%m%d_%H%M%S")).seconds, ))
        if str(day) in dateWithNans.keys():
            dotmgr = dateWithNans[str(day)]
        else:
            dotmgr = str(day)
        grp.create_dataset('Day_of_temperature_mean_GPS-RO', data = (dotmgr, ))
        # part0['lhead']
        # part0['outnfmt']
        # part0['mode']
        # part0['itime']
        # part0['numpart']
        # part0['nact']
        # part0['idx_orgn']
        # part0['nact_lastO']
        # part0['nact_lastNM']
        # part0[['nact_lastNH']
        # part0['flag'] 
        # part0['idx_back']
        lons_st = np.where(part0['x']>180, part0['x']-360, part0['x'])
        for iday in inDays:
            idi = part0['ir_start'] == iday
            ids = (datetime.datetime.strptime(str(part0['stamp_date']), "%Y%m%d%H%M%S") - datetime.timedelta(seconds=int(iday*-1))).strftime("%Y%m%d_%H%M%S")
            grp = h5file.create_group(ids)
            grp_stv = grp.create_group('Start_Values')
            grp_stv.create_dataset('ir_start', data = (iday, ))
            
            grp_stv.create_dataset('lons (0 - 360)', data = part0['x'][idi])
            grp_stv.create_dataset('lons (-180 - 180)', data = lons_st[idi])
            grp_stv.create_dataset('lats', data = part0['y'][idi])
            grp_stv.create_dataset('pressure', data = part0['p'][idi])
            grp_stv.create_dataset('temperature', data = part0['t'][idi])

            grp_fpv = grp.create_group('Flexp_Values')
            grp_fpv.create_dataset('output_time', data = outptime)
            grp_fpv.create_dataset('itime', data = itime)
            grp_fpv.create_dataset('lons', data = lons[idi, :])
            grp_fpv.create_dataset('lats', data = lats[idi, :])
            grp_fpv.create_dataset('pressure', data = ps[idi, :])
            grp_fpv.create_dataset('temperature', data = ts[idi, :])
        
        h5file.close()
        
        
        
        lons_plot = np.where(lons>180, lons-360, lons)
        lons_plot = np.concatenate((np.zeros([lons_plot.shape[0], 1])+np.nan, lons_plot), axis=1)
        lats_plot = np.concatenate((np.zeros([lats.shape[0], 1])+np.nan, lats), axis=1)
        ms = []
        for n in range(lons_plot.shape[0]):
            for m in range(lons_plot.shape[1]-1):
            # ni  = np.where(np.isnan(lons_plot[n, :]))[0][-1]
                if np.isnan(lons_plot[n, m]) and not np.isnan(lons_plot[n, m + 1]):
                    lons_plot[n, m] = lons_st[n]
                    lats_plot[n, m] = part0['y'][n]
                    ms.append(m)
                    break
        
        

        print('Ploting')
        tic = time.time()
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
        if int(str(day)[0:4]) in [2020]: 
            ax.set_extent([-135, -45, -45, 45], crs=ccrs.PlateCarree())
        elif int(str(day)[0:4]) in [2021]: 
            ax.set_extent([-180, 180, -45, 45], crs=ccrs.PlateCarree())
        elif int(str(day)[0:4]) in [2022]: 
            ax.set_extent([-180, 180, -45, 45], crs=ccrs.PlateCarree())
        elif day in [20191216]:
            ax.set_extent([90, 180, -45, 45], crs=ccrs.PlateCarree())
        # elif day in [20200113]:
        # elif day in [20200111]:
        #     ax.set_extent([-135, -45, -45, 45], crs=ccrs.PlateCarree())
        # elif day in [20200119]:
        #     ax.set_extent([-135, -45, -45, 45], crs=ccrs.PlateCarree())
        else:
            ax.set_extent([-180, 180, -45, 45], crs=ccrs.PlateCarree())
        ax.coastlines()
        # [0::500]
        # step = int(lons[hits].shape[0]/250)
        # plt.plot([lons[hits][0::step], lons_0[hits][0::step]], [lats[hits][0::step], lats_0[hits][0::step]],
        #      linewidth=1, transform=ccrs.Geodetic())
        colours = ['r', 'b', 'g', 'c', 'm', 'y']
        tr =-1
        for iday in inDays:
            tr = tr + 1
            if tr > (len(colours) - 1):
                break
            idi = part0['ir_start'] == iday
            plt.plot([lons_plot[idi, :][tr, 0:-1], lons_plot[idi, :][tr, 1:]], [lats_plot[idi, :][tr, 0:-1], lats_plot[idi, :][tr, 1:]],
                 linewidth=1, transform=ccrs.Geodetic(), color=colours[tr])
            naninds = ~np.isnan(lons_plot[idi, :][tr, :])
            
            plt.plot(lons_plot[idi, :][tr, :][naninds][0], lats_plot[idi, :][tr, :][naninds][0], 'x', transform=ccrs.Geodetic(), color=colours[tr])
            plt.plot(lons_plot[idi, :][tr, :][naninds][-1], lats_plot[idi, :][tr, :][naninds][-1], 's', transform=ccrs.Geodetic(), color=colours[tr])
            # pdb.set_trace()
            
        # for tr in range(3):#lons_plot.shape[0]):
        #     plt.plot([lons_plot[tr, 0:-1], lons_plot[tr, 1:]], [lats[tr, 0:-1], lats[tr, 1:]],
        #          linewidth=1, transform=ccrs.Geodetic(), color=colours[tr])
        plt.tight_layout()
        figname = '%s/traj_%i' %(plotDir, day)
        fig.savefig(figname + '.png')
        toc = time.time()
        print(toc-tic)
        print(figname)

    pdb.set_trace()
        # h5file = h5py.File(os.path.join(h5Dir, 'coord-%i' %day), 'r')
        # pdb.set_trace()
        # h5file.close()
        # step = 3
        # hmax = 1800#3240 #1800 #: 75 days (75 * 24 = 1800)
        # dstep = datetime.timedelta(hours=step)
        # for phour in range(step, hmax + 1, step):