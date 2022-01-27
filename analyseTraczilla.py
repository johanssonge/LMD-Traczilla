#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2021-10-04

Copyright (c) 2021 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""

import numpy as np
import h5py  # @UnresolvedImport
import pdb
# import matplotlib
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt  # @UnresolvedImport
import cartopy.crs as ccrs  # @UnresolvedImport
import time
import os
import sys
import datetime
import pickle
import gzip

sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readpart107, readidx107

I_DEAD = 0x200000 #: Dead
I_HIT = 0x400000 #: Hit a cloud
I_CROSSED = 0x2000000 #: outside domain
I_DBORNE =  0x1000000 #: Lost after first step
I_OLD = 0x800000 #: Reached end of time without encounter a cloud
I_STOP = I_HIT + I_DEAD

magic = 0.84


#:----------------------------------------------------

def readParamFile(fn):
    objects = []
    # with (gzip.open(fn, "rb")) as openfile:
    with (open(fn, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    # if p0 is not None:
    #     retv = {'CM': np.zeros(p0['numpart']), 'TpH': np.zeros(p0['numpart']), \
    #             'SZA': np.zeros(p0['numpart']), 'SC': np.zeros(p0['numpart'])}
    #     temp = objects[0]
    #     sp = 0
    #     for da in temp.keys():
    #         for orbit in temp[da].keys():
    #             ep = sp + temp[da][orbit]['npart']
    #             retv['CM'][sp:ep] = temp[da][orbit]['Cloud_Mask']
    #             retv['TpH'][sp:ep] = temp[da][orbit]['Tropopause_Height']
    #             retv['SZA'][sp:ep] = temp[da][orbit]['SZA']
    #             retv['SC'][sp:ep] = temp[da][orbit]['Simplified_Categorization']
    #             #: Test to make sure the param korrespond to the part0
    #             if not temp[da][orbit]['longitudes'][0] == p0['x'][sp]:
    #                 print('Lon 0 test')
    #                 pdb.set_trace()
    #             #: Not sure why -1
    #             if not temp[da][orbit]['longitudes'][1] == p0['x'][ep-1]:
    #                 print('Lon 1 test')
    #                 pdb.set_trace()
    #             # pdb.set_trace()
    #             sp = ep
    # else:
    retv = objects[0]
    return retv


#:----------------------------------------------------

def readCatalogFile(fn, p0=None):
    objects = []
    with (gzip.open(fn, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    if p0 is not None:
        retv = {'CM': np.zeros(p0['numpart']), 'TpH': np.zeros(p0['numpart']), \
                'SZA': np.zeros(p0['numpart']), 'SC': np.zeros(p0['numpart']), \
                'lensel': np.zeros(p0['numpart']), 'npart': np.zeros(p0['numpart'])}
        temp = objects[0]
        sp = 0
        for da in temp.keys():
            for orbit in temp[da].keys():
                ep = sp + temp[da][orbit]['npart']
                retv['CM'][sp:ep] = temp[da][orbit]['Cloud_Mask']
                retv['TpH'][sp:ep] = temp[da][orbit]['Tropopause_Height']
                retv['SZA'][sp:ep] = temp[da][orbit]['SZA']
                retv['SC'][sp:ep] = temp[da][orbit]['Simplified_Categorization']
                retv['lensel'][sp:ep] = np.repeat(temp[da][orbit]['lensel'], (ep-sp))
                retv['npart'][sp:ep] = np.repeat(temp[da][orbit]['npart'], (ep-sp))
                #: Test to make sure the param korrespond to the part0
                if not temp[da][orbit]['longitudes'][0] == p0['x'][sp]:
                    print('Lon 0 test')
                    pdb.set_trace()
                #: Not sure why -1
                if not temp[da][orbit]['longitudes'][1] == p0['x'][ep-1]:
                    print('Lon 1 test')
                    pdb.set_trace()
                # pdb.set_trace()
                sp = ep
    else:
        retv = objects[0]
    return retv


def readConvFile(fn, usefk=False):
    """ 
    Read the resultfile from convsrcErikFullGridSatTropo
    Deafault is to use h5py but flamkuchen is an option, set usefk = True
    """
    if usefk:
        import flammkuchen as fk  # @UnresolvedImport
        h5f = fk.load(fn)
    else:
        h5f = h5py.File(fn, 'r')
    
    #: Some h5f starts with data. Do not know why
    if 'data' in h5f.keys():
        upperLev = 'data/'
    else:
        upperLev = ''
    
    rvs = h5f['%srvs' %upperLev][:]
    fls = h5f['%sflag_source' %upperLev][:]
    age = h5f['%ssrc' %upperLev]['age'][:]
    lons = h5f['%ssrc' %upperLev]['x'][:]
    lats = h5f['%ssrc' %upperLev]['y'][:]
    temp = h5f['%ssrc' %upperLev]['t'][:]
    pres = h5f['%ssrc' %upperLev]['p'][:]
    
    if not usefk:
        h5f.close()
    return rvs, fls, age, lons, lats, temp, pres


def getInitfiles(mD, years, months, dn, uD):
    if (not isinstance(years, list)) or (not isinstance(months, list)):
        print('Wrong type')
        pdb.set_trace()
    outnames = []
    i = -1
    for y in years:
        for m in months:
            i = i + 1
            #: Directories of the backward trajectories and name of the output file
            outname = 'CALIOP-EAD-%d%02d-%s' %(y, m, dn)
            if uD:
                outname = outname + '-DD'
            outnames.append(outname)
            trajDir = os.path.join(mD,outname) 
            initDir = os.path.join(trajDir, 'Initfiles')
            
            paramname = 'selDardar_Params-%s.pkl' %'-'.join(outname.split('-')[2:4])
            catalogname = paramname.replace('_Params-', '_Calalog-')
            paramFile = os.path.join(initDir, paramname)
            catalFile = os.path.join(initDir, catalogname)
            pf = readidx107(os.path.join(trajDir,'part_000'),quiet=False)
            cf = readCatalogFile(catalFile, pf)
            mf = readParamFile(paramFile)
            pf['x'] = np.where(pf['x'] > 180, pf['x'] - 360, pf['x'])
            if i == 0:
                retvP = pf.copy()
                retvC = cf.copy()
                retvM = mf.copy()
            else:
                for arname in retvP.keys():
                    if isinstance(retvP[arname], int):
                        retvP[arname] = [retvP[arname]]
                    if isinstance(pf[arname], int):
                        pf[arname] = [pf[arname]]
                    retvP[arname] = np.concatenate((np.asarray(retvP[arname]), np.asarray(pf[arname])))
                
                for arname in retvC.keys():
                    if isinstance(retvC[arname], int):
                        retvC[arname] = [retvC[arname]]
                    if isinstance(cf[arname], int):
                        cf[arname] = [cf[arname]]
                    retvC[arname] = np.concatenate((np.asarray(retvC[arname]), np.asarray(cf[arname])))
                
                for arname in retvM.keys():
                    if isinstance(retvM[arname], (int, np.int64, str, datetime.datetime)):
                        retvM[arname] = [retvM[arname]]
                    if isinstance(mf[arname], (int, np.int64, str, datetime.datetime)):
                        mf[arname] = [mf[arname]]
                    retvM[arname] = np.concatenate((np.asarray(retvM[arname]), np.asarray(mf[arname])))
    
    return retvP, retvC, retvM, outnames


def getConvfiles(oD, onames):
    i = -1
    for outname in onames:
        i = i + 1
        fname = os.path.join(oD, outname + '.h5')
        rvs, fls, age, lons, lats, temp, pres = readConvFile(fname)
        if i == 0:
            retvRvs = rvs.copy()
            retvFls = fls.copy()
            retvAge = age.copy()
            retvLons = lons.copy()
            retvLats = lats.copy()
            retvTemp = temp.copy()
            retvPres = pres.copy()
        else:
            retvRvs = np.concatenate((retvRvs, rvs))
            retvFls = np.concatenate((retvFls, fls))
            retvAge = np.concatenate((retvAge, age))
            retvLons = np.concatenate((retvLons, lons))
            retvLats = np.concatenate((retvLats, lats))
            retvTemp = np.concatenate((retvTemp, temp))
            retvPres = np.concatenate((retvPres, pres))
            
    return retvRvs, retvFls, retvAge, retvLons, retvLats, retvTemp, retvPres


def checkLons(lon, hit = None):
    """
    Make sure longitudes are in range
    
    Some dead pixel might have wondering outside the box 
    and have therefore goten some extra longitude values.
    
    This should not happend with pixels that hit an cloud
    """
    if hit is not None:
        if lon[hit].max() > 180:
            print('\n')
            print('Check lons')
            print('lons max = %f' %lon[hit].max())
            sys.exit()
        if lon[hit].min() < -180:
            print('\n')
            print('Check lons')
            print('lons min = %f' %lon[hit].max())
            sys.exit()
    lon = np.where(lon > 180, lon - 360, lon)
    lon = np.where(lon < -180, lon + 360, lon)
    return lon


def getDates(parf, inds=None):
    """
    Convert idx_back to datetime objects
    """
    stamp_year = int(str(parf['stamp_date'])[0:4])
    stamp_month = int(str(parf['stamp_date'])[4:6])
    stamp_day = int(str(parf['stamp_date'])[6:8])
    stamp_hour = int(str(parf['stamp_date'])[8:10])
    stamp_minute = int(str(parf['stamp_date'])[10:12])
    stamp_second = int(str(parf['stamp_date'])[12:14])
    
    originDate = datetime.datetime(stamp_year, stamp_month, stamp_day, stamp_hour, stamp_minute, stamp_second)
    if inds is None:
        inds = np.ones(parf['idx_back'].shape).astype(bool)
    tic = time.time()
    idxDates = [originDate + datetime.timedelta(seconds=float(t)) for t in parf['idx_back'][inds]]
    print(time.time() - tic)
    return originDate, idxDates


def conv2DHistToImage(h, x, y):
    xstep = x[3] - x[2]  # @UnusedVariable
    ystep = y[3] - y[2]
    yim = np.arange(-60, 60, ystep)
    yargmin = np.argmin(np.abs(yim - y[0]))
    yargmax = np.argmin(np.abs(yim - y[-1]))
    him = np.zeros([h.shape[0], yim.shape[0]])
    him[:, yargmin: yargmax] = h
    return him


if __name__ == '__main__':
    
    # if 'ciclad' in socket.gethostname():
        #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
    datPath = os.environ['HOME'].replace('/home/', '/data/')
    scrPath = os.environ['HOME'].replace('/home/', '/scratchu/')
    ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson'
    mainDir = '%s/flexout/STC/Calipso' %ekjDir
    outDir = '%s/flexout/STC/Calipso-OUT' %datPath
    plotDir = '%s/LMD-Traczilla/Plots' %scrPath
    
    compare = False
    year = [2008]
    month = [1]#[*range(1, 13)]#[1]
    dn = 'n'
    useDardar = True
    
    # outDir = '/data/ejohansson/flexout/STC/Calipso-OUT/Coldpoint'
    #: filename of outfiles
    # outnames = 'CALIOP-'+advect+'-'+date_end.strftime('%b%Y')+diffus+'-%s' %args.night 
    # if args.use_dardar:
    #     outnames = outnames + '-DD'
    print('Read Init-files')
    tic = time.time()
    part0, catalog, params, outnames = getInitfiles(mainDir, year, month, dn, useDardar)
    print('It took %d sec to read Init-files' %(time.time() - tic))
    #: Read the index file that contains the initial positions
    # print('numpart',part0['numpart'])
    #:
    
    lons_0 = part0['x']#np.where(part0['x'] > 180, part0['x'] - 360, part0['x'])
    lats_0 = part0['y']
    
    print('Read Conv-files')
    tic = time.time()
    rvs, fls, age, lons, lats, temp, pres = getConvfiles(outDir, outnames)
    print('It took %d sec to read Conv-files' %(time.time() - tic))
    
    # print(time.time() - tic)
    olds = (fls & I_OLD) == I_OLD
    hits = (fls & I_HIT) == I_HIT
    
    lons = checkLons(lons, hits)
    # originDate, idxDates = getDates(part0, hits)
    cloudy = (catalog['CM'] > 0)
    hits_cld = hits & cloudy
    hits_clr = hits & ~cloudy
    
    olds_cld = olds & cloudy
    olds_clr = olds & ~cloudy
    
    # step = 10
    # latarr = np.asarray(range(-90, 90, step))
    # lonarr = np.asarray(range(-180, 180, step))
    #
    # hitarr = np.zeros([latarr.shape[0], lonarr.shape[0]])
    # oldarr = np.zeros([latarr.shape[0], lonarr.shape[0]])
    # tic = time.time()
    # for i in range(latarr.shape[0]):
    #     latind = ((lats >= latarr[i]) & (lats < (latarr[i] + step)))
    #     for j in range(lonarr.shape[0]):
    #         lonind = ((lons >= lonarr[j]) & (lons < (lonarr[j] + step)))
    #         hitarr[i, j] = (latind & lonind & hits).sum()
    #         oldarr[i, j] = (latind & lonind & olds).sum()
    # toc = time.time()
    # print(toc-tic)
    if compare:
        outname = outnames[0]
        outDirCP = '%s/flexout/STC/Calipso-OUT/Coldpoint' %datPath
        outDirOrg = '%s/flexout/STC/Calipso-OUT/Org' %datPath
        coutname = outname.replace('-DD', '')
        fname2 = os.path.join(outDir, coutname + '.h5')
        fname3 = os.path.join(outDirCP, coutname + '.h5')
        fname4 = os.path.join(outDirOrg, coutname + '.h5')
        rvs2, fls2, age2, lons2, lats2, temp2, pres2 = readConvFile(fname2)
        rvs3, fls3, age3, lons3, lats3, temp3, pres3 = readConvFile(fname3)
        rvs4, fls4, age4, lons4, lats4, temp4, pres4 = readConvFile(fname4)
        hits2 = (fls2 & I_HIT) == I_HIT
        hits3 = (fls3 & I_HIT) == I_HIT
        hits4 = (fls4 & I_HIT) == I_HIT
        lons2 = checkLons(lons2, hits2)
        lons3 = checkLons(lons3, hits3)
        lons4 = checkLons(lons4, hits4)
    
    #: https://stackoverflow.com/questions/50611018/cartopy-heatmap-over-openstreetmap-background/50638691

    print('Plot')
    #: --- Plot ---
    fig = plt.figure()
    fig.suptitle('Age histogram')
    ax = fig.add_subplot(1,1,1)
    
    h = ax.hist(age[hits] / 86400, bins=400)#, density=True)
    # h = ax.hist(age[hits] / 86400, bins=np.logspace(np.log10(0.001),np.log10(42.0), 400))#bins=400)#, density=True)
    # ax.set_xscale('log')
    # fig.gca().set_xscale("log")
    ax.set_xlabel('Age [days]')
    ax.set_title('Hits pixels')
    ax.text(0.7, 0.9,'total = %d' %hits.sum(),
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
    plt.tight_layout()
    figname = '%s/hist_age_all' %plotDir
    fig.savefig(figname + '.png')
    # pdb.set_trace()
    
    #: --- Plot ---
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(age[hits_cld] / 86400, bins=400)#, density=True)
    ax.set_xlabel('Age [days]')
    ax.set_title('Hits and Cloudy pixels')
    ax.text(0.7, 0.9,'total = %d' %hits_cld.sum(),
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
    plt.tight_layout()
    figname = '%s/hist_age_cld' %plotDir
    fig.savefig(figname + '.png')
    
    #: --- Plot ---
    fig = plt.figure()
    # fig.suptitle('Age histogram')
    ax = fig.add_subplot(1,1,1)
    ax.hist(age[hits_clr] / 86400, bins=400)#, density=True)
    ax.set_xlabel('Age [days]')
    ax.set_title('Hits and Clear pixels')
    ax.text(0.7, 0.9,'total = %d' %hits_clr.sum(),
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
    plt.tight_layout()
    figname = '%s/hist_age_clr' %plotDir
    fig.savefig(figname + '.png')

    #: --- Plot ---
    if False:
        tic = time.time()
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 190, -90, 90], crs=ccrs.PlateCarree())
        ax.coastlines()
        # [0::500]
        step = int(lons[hits].shape[0]/250)
        plt.plot([lons[hits][0::step], lons_0[hits][0::step]], [lats[hits][0::step], lats_0[hits][0::step]],
             linewidth=1, transform=ccrs.Geodetic())
        plt.tight_layout()
        fig.savefig('%s/traj.png' %plotDir)
        toc = time.time()
        print(toc-tic)
    
    
    #: --- Plot ---
    # from matplotlib import ticker
    
    
    hh1, xedges1, yedges1 = np.histogram2d(lons[hits], lats[hits], bins=[180, 90])
    hh2, xedges2, yedges2 = np.histogram2d(lons_0[hits], lats_0[hits], bins=[180, 90])#, bins=400)#, density=True)
    hhim1 = conv2DHistToImage(hh1, xedges1, yedges1)
    hhim2 = conv2DHistToImage(hh2, xedges2, yedges2)

    fig = plt.figure()
    fig.suptitle('2D Hist - Hits pixels')
    ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    im1 = ax.imshow(hhim1.T, origin ='lower', extent = [-180, 180, -60, 60])#, vmin=vmin, vmax=vmax)
    ax.set_title('Trazilla')
    fig.subplots_adjust(right=0.89)
    cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
    cbar = fig.colorbar(im1, cbar_ax)
    ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    im2 = ax.imshow(hhim2.T, origin ='lower', extent = [-180, 180, -60, 60])#, vmin=vmin, vmax=vmax)
    ax.set_title('Dardar')
    fig.subplots_adjust(right=0.89)
    cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    cbar = fig.colorbar(im2, cax=cbar_ax)
    # plt.tight_layout()
    figname = '%s/2dhist_hits_all' %plotDir
    fig.savefig(figname + '.png')
    
#     fig = plt.figure()
#     fig.suptitle('2D Hist - Hits pixels')
#     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     # ax.gridlines(draw_labels=True)
#     # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, linewidth=2.001, color='k',alpha=0)
#     # gl.right_labels = gl.top_labels = False
#     # gl.ylocator = ticker.FixedLocator([60, -30, 0, 30, 60])
#     # gl.xlocator = ticker.FixedLocator([-180, -90, 0, 90, 180])
#     # do coordinate conversion of (x,y)
#     # xynps = ax.projection.transform_points(ccrs.Geodetic(), x, y)
#     # im1 = ax.imshow(hh1, vmin=vmin, vmax=vmax)
#     hh = ax.hist2d(lons[hits], lats[hits], bins=[180, 90])#, bins=400)#, density=True)
#     fig.subplots_adjust(right=0.9)
#     # cbar_ax = fig.add_axes([0.91, 0.54, 0.008, 0.24])
#     # fig.colorbar(hh[3], cax=cbar_ax)#, orientation='horizontal')
#     # import matplotlib as mpl
#     # cb1 = mpl.colorbar.ColorbarBase(ax, orientation='horizontal')
#     ax.set_xlabel('Longitude')
#     ax.set_ylabel('Latitude')
#     ax.set_title('Trazilla')
#     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     # im2 = ax.imshow(hh2, vmin=vmin, vmax=vmax, extent = [-180, 180, -60, 60])
#     hh = ax.hist2d(lons_0[hits], lats_0[hits], bins=[180, 90])#, bins=400)#, density=True)
#     # cbar_ax = fig.add_axes([pos_x, axpos.y0, 0.01, axpos.height])
#     # cbar = fig.colorbar(im, orientation='vertical', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
# #   cbar=fig.colorbar(im)#, cax=pos_cax)
#     # cbar.ax.tick_params(labelsize=fs)
#     # divider = make_axes_locatable(ax)
#     # cax = divider.append_axes("right", size="5%", pad=0.05)
#     fig.subplots_adjust(right=0.9)
#     # cbar_ax = fig.add_axes([0.91, 0.04, 0.008, 0.24])
#     # cbar = fig.colorbar(hh[3], cax=cbar_ax)#, orientation='horizontal')
#     ax.set_xlabel('Longitude')
#     ax.set_ylabel('Latitude')
#     ax.set_title('Dardar')
#
#     # plt.tight_layout()
#     fig.savefig('%s/2dhist_hits_all.png' %plotDir)
#

    
    
    fig = plt.figure()
    fig.suptitle('2D Hist - Hits and Cloudy pixels')
    ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.hist2d(lons[hits_cld], lats[hits_cld], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Trazilla')
    ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.hist2d(lons_0[hits_cld], lats_0[hits_cld], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Dardar')
    plt.tight_layout()
    fig.savefig('%s/2dhist_hits_cld.png' %plotDir)
    
    
    fig = plt.figure()
    fig.suptitle('2D Hist - Hits and Clear pixels')
    ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.hist2d(lons[hits_clr], lats[hits_clr], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Trazilla')
    ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.hist2d(lons_0[hits_clr], lats_0[hits_clr], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Dardar')
    plt.tight_layout()
    fig.savefig('%s/2dhist_hits_clr.png' %plotDir)
    
    pdb.set_trace()
    #: --- Plot ---
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist2d(lons[olds], lats[olds], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Old pixels')
    plt.tight_layout()
    fig.savefig('%s/2dhist_olds_all.png' %plotDir)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist2d(lons[olds_cld], lats[olds_cld], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Old and Cloudy pixels')
    plt.tight_layout()
    fig.savefig('%s/2dhist_olds_cld.png' %plotDir)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist2d(lons[olds_clr], lats[olds_clr], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Old and Clear pixels')
    plt.tight_layout()
    fig.savefig('%s/2dhist_olds_clr.png' %plotDir)
    
    #: --- Plot ---
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist2d(lons_0[hits], lats_0[hits], bins=[180, 90])#, bins=400)#, density=True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Hits pixels')
    plt.tight_layout()
    fig.savefig('%s/temp.png' %plotDir)
    
    
    #: rsv
    rvsH = rvs.copy()
    rvsH = np.where(hits, 1.6 * magic, np.nan)
    
    
    
    
    pdb.set_trace()
    
    
    
    
    
    
    
    
    
    
    #: --- Plot ---
    
    if compare:
        fig = plt.figure()
        fig.suptitle('2D histograms')
        ax = fig.add_subplot(2,2,1)
        ax.hist2d(age[hits] / 86400, pres[hits], bins=400)#, density=True)
        ax.invert_yaxis()
        ax.set_xlabel('Age [days]')
        ax.set_ylabel('Pressure')
        ax.set_title('WMO - Dardar - vshift 10')
        
        ax = fig.add_subplot(2,2,2)
        ax.hist2d(age2[hits2] / 86400, pres2[hits2], bins=400)#, density=True)
        ax.invert_yaxis()
        ax.set_xlabel('Age [days]')
        ax.set_ylabel('Pressure')
        ax.set_title('WMO - Calipso - vshift 10')
    
        ax = fig.add_subplot(2,2,3)
        ax.hist2d(age3[hits3] / 86400, pres3[hits3], bins=400)#, density=True)
        ax.invert_yaxis()
        ax.set_xlabel('Age [days]')
        ax.set_ylabel('Pressure')
        ax.set_title('Coldpoint - Calipso')
        
        ax = fig.add_subplot(2,2,4)
        ax.hist2d(age4[hits4] / 86400, pres4[hits4], bins=400)#, density=True)
        ax.invert_yaxis()
        ax.set_xlabel('Age [days]')
        ax.set_ylabel('Pressure')
        ax.set_title('WMO - Calipso')
        
    
        fig = plt.figure()
        fig.suptitle('Comparing age histograms')
        ax = fig.add_subplot(2,2,1)
        ax.hist(age[hits] / 86400, bins=200, density=True)
        # ax.set_ylim(0, 190000)
        ax.set_xlabel('Age [days]')
        # ax.set_ylabel('Number')
        ax.set_title('WMO - Dardar - vshift 10')
        ax.text(0.7, 0.9,'total = %d' %hits.sum(),
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        
        ax = fig.add_subplot(2,2,2)
        ax.hist(age2[hits2] / 86400, bins=200, density=True)
        # ax.set_ylim(0, 190000)
        ax.set_xlabel('Age [days]')
        # ax.set_yticks([])
        # ax.set_ylabel('Number')
        ax.set_title('WMO - Calipso - vshift 10')
        ax.text(0.7, 0.9,'total = %d' %hits2.sum(),
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        # ax.plot(age[hits], pres[hits])
    
        ax = fig.add_subplot(2,2,3)
        ax.hist(age3[hits3] / 86400, bins=200, density=True)
        # ax.set_ylim(0, 300000)
        ax.set_xlabel('Age [days]')
        # ax.set_ylabel('Number')
        ax.set_title('Coldpoint - Calipso')
        ax.text(0.7, 0.9,'total = %d' %hits3.sum(),
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        
        ax = fig.add_subplot(2,2,4)
        # ax.hist(dage[dhits] // 86400, bins=100)
        ax.hist(age4[hits4] / 86400, bins=200, density=True)
        # ax.set_ylim(0, 600000)
        ax.set_xlabel('Age [days]')
        # ax.set_ylabel('Number')
        ax.set_title('WMO - Calipso')
        ax.text(0.7, 0.9,'total = %d' %hits4.sum(),
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        
        plt.tight_layout()
        figname = '%s/comp_hist_age' %plotDir
        fig.savefig(figname + '.png')
    
    
    
    
    
    
    
    
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    # rotated_pole = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
    # ax.plot(0, 0, 'o', transform=rotated_pole)
    # ax.coastlines()
    # # ax.plot(age[argsort], rvs[argsort])
    #
    # fig.savefig('test.png')
    pdb.set_trace()
    