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
                'SZA': np.zeros(p0['numpart']), 'SC': np.zeros(p0['numpart'])}
        temp = objects[0]
        sp = 0
        for da in temp.keys():
            for orbit in temp[da].keys():
                ep = sp + temp[da][orbit]['npart']
                retv['CM'][sp:ep] = temp[da][orbit]['Cloud_Mask']
                retv['TpH'][sp:ep] = temp[da][orbit]['Tropopause_Height']
                retv['SZA'][sp:ep] = temp[da][orbit]['SZA']
                retv['SC'][sp:ep] = temp[da][orbit]['Simplified_Categorization']
                #: Test to make sure the param korrespond to the part0
                if not temp[da][orbit]['longitudes'][0] == p0['x'][sp]:
                    print('Lon 0 test')
                    pdb.set_trace()
                #: Not sure why -1
                if not temp[da][orbit]['longitudes'][1] == p0['x'][ep-1]:
                    print('Lon 1 test')
                    pdb.set_trace()
                
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
    retv = np.where(lon > 180, lon-360, lon)
    return retv


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



if __name__ == '__main__':
    
    # if 'ciclad' in socket.gethostname():
        #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
    datPath = os.environ['HOME'].replace('/home/', '/data/')
    scrPath = os.environ['HOME'].replace('/home/', '/scratchu/')
    mainDir = '%s/flexout/STC/Calipso' %datPath
    outDir = '%s/flexout/STC/Calipso-OUT' %datPath
    plotDir = '%s/LMD-Traczilla/Plots' %scrPath
    
    compare = False
    year = 2007
    month = 6
    dn = 'n'
    useDardar = True
    
    # outDir = '/data/ejohansson/flexout/STC/Calipso-OUT/Coldpoint'
    #: filename of outfiles
    # outnames = 'CALIOP-'+advect+'-'+date_end.strftime('%b%Y')+diffus+'-%s' %args.night 
    # if args.use_dardar:
    #     outnames = outnames + '-DD'
    
    #: Directories of the backward trajectories and name of the output file
    outname = 'CALIOP-EAD-%d%02d-%s' %(year, month, dn)
    if useDardar:
        outname = outname + '-DD'
    trajDir = os.path.join(mainDir,outname) 
    initDir = os.path.join(trajDir, 'Initfiles')
    
    paramname = 'selDardar_Params-%s.pkl' %'-'.join(outname.split('-')[2:4])
    catalogname = paramname.replace('_Params-', '_Calalog-')
    paramFile = os.path.join(initDir, paramname)
    catalFile = os.path.join(initDir, catalogname)
    

    #: Read the index file that contains the initial positions
    part0 = readidx107(os.path.join(trajDir,'part_000'),quiet=False)
    # print('numpart',part0['numpart'])
    #:
    
    
    fname = os.path.join(outDir, outname + '.h5')
    rvs, fls, age, lons, lats, temp, pres = readConvFile(fname)
    
    catalog = readCatalogFile(catalFile, part0)
    # print(time.time() - tic)
    old = (fls & I_OLD) == I_OLD
    hits = (fls & I_HIT) == I_HIT

    lons = checkLons(lons, hits)
    # originDate, idxDates = getDates(part0, hits)
    cloudy = (catalog['CM'] > 0)
    hits_cld = hits & cloudy
    hits_clr = hits & ~cloudy
    if compare:
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
    
    #: --- Plot ---
    fig = plt.figure()
    fig.suptitle('Age histogram')
    ax = fig.add_subplot(1,1,1)
    ax.hist(age[hits] / 86400, bins=400)#, density=True)
    ax.set_xlabel('Age [days]')
    ax.set_title('Hits pixels')
    ax.text(0.7, 0.9,'total = %d' %hits.sum(),
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
    plt.tight_layout()
    figname = '%s/hist_age_all' %plotDir
    fig.savefig(figname + '.png')
    
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
    