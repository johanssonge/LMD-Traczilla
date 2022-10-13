#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2022-06-28

Copyright (c) 2022 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""
import numpy as np
import socket
import os
import pdb
from matplotlib import pyplot as plt  # @UnresolvedImport
from matplotlib import colors as mcolors  # @UnresolvedImport
import cartopy.crs as ccrs  # @UnresolvedImport
from analyseTraczilla import getConvfiles, missing_months, readCatalogFile,\
    getCatalogFile, areas
from convsrcErikFullGridSatTropo import I_HIT, I_OLD, I_DEAD, I_CROSSED, I_DBORNE
import sys



def plotFlexpVar(outDir, plotDir, mainDir, dc=False):
    tot_cld_DEAD = 0
    tot_cld_HIT = 0
    tot_cld_CROSSED = 0
    tot_cld_DBORNE =  0
    tot_cld_OLD = 0
    tot_cld_TOT = 0
    
    tot_DEAD = 0
    tot_HIT = 0
    tot_CROSSED = 0
    tot_DBORNE =  0
    tot_OLD = 0
    tot_TOT = 0
    for year in range(2007, 2020):
        for mon in range(1, 13):
            #: Check for missing months
            if year in missing_months.keys():
                if mon in missing_months[year]:
                    continue
            print('Year = %d, Month = %02d' %(year, mon))
            
            outname = 'CALIOP-EAD-%d%02d-n-DD' %(year, mon)
            rvs_mon, fls_mon, age_mon, lons_monF, lats_monF, temp_monF, pres_monF = getConvfiles(outDir, ['', outname])  # @UnusedVariable
            
            tot_TOT = tot_TOT + fls_mon.shape[0]
            tot_HIT = tot_HIT + ((fls_mon & I_HIT) == I_HIT).sum()
            tot_OLD = tot_OLD + ((fls_mon & I_OLD) == I_OLD).sum()
            tot_CROSSED = tot_CROSSED + ((fls_mon & I_CROSSED) == I_CROSSED).sum()
            tot_DBORNE =  tot_DBORNE + ((fls_mon & I_DBORNE) == I_DBORNE).sum()
            tot_DEAD = tot_DEAD + ((fls_mon & I_DEAD) == I_DEAD).sum()
            
            if dc:
                partDir = os.path.join(mainDir, outname)
                initDir = os.path.join(partDir, 'Initfiles')
                catalogFile = os.path.join(initDir, 'selDardar_Calalog-%d%02d-n.pkl' %(year, mon))
                cf = getCatalogFile(catalogFile, checkForNewFile=False)
                cld = cf['iwc'] > 0     
                tot_cld_TOT = tot_cld_TOT + cld.sum()
                tot_cld_HIT = tot_cld_HIT + (((fls_mon & I_HIT) == I_HIT) & cld).sum()
                tot_cld_OLD = tot_cld_OLD + (((fls_mon & I_OLD) == I_OLD) & cld).sum()
                tot_cld_CROSSED = tot_cld_CROSSED + (((fls_mon & I_CROSSED) == I_CROSSED) & cld).sum()
                tot_cld_DBORNE =  tot_cld_DBORNE + (((fls_mon & I_DBORNE) == I_DBORNE) & cld).sum()
                tot_cld_DEAD = tot_cld_DEAD + (((fls_mon & I_DEAD) == I_DEAD) & cld).sum()

    #: ---- Plot ---
    #--- 1 ---
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # rects = ax.bar(range(1, 7), [tot_TOT, tot_HIT, tot_OLD, tot_CROSSED, tot_DBORNE, tot_DEAD])
    rects = ax.bar(range(1, 6), [tot_TOT, tot_HIT, tot_OLD, tot_CROSSED, tot_DBORNE])
    #: Add counts above the bars
    for rect in rects:
        hs = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%d' %hs, ha='center', va='bottom')
    ax.set_xticks(range(1, 6))
    # ax.set_xticklabels(['TOT', 'HIT', 'OLD', 'CROSSED', 'DBORNE', 'DEAD'], fontsize='large')
    ax.set_xticklabels(['TOT', 'HIT', 'OLD', 'CROSSED', 'DBORNE'], fontsize='large')
    ax.set_title('All pixels')  
    plt.tight_layout()
    figname_1 = '%s/parcel-dest_all_tot' %plotDir
    print(figname_1)
    fig.savefig(figname_1 + '.png')
    if 'oem-Latitude-5400' in socket.gethostname():
        fig.show()
    #--- 2 ---
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # rects = ax.bar(range(1, 6), [tot_HIT / tot_TOT, tot_OLD / tot_TOT, tot_CROSSED / tot_TOT, tot_DBORNE / tot_TOT, tot_DEAD / tot_TOT])
    rects = ax.bar(range(1, 5), [tot_HIT / tot_TOT, tot_OLD / tot_TOT, tot_CROSSED / tot_TOT, tot_DBORNE / tot_TOT])
    #: Add counts above the bars
    for rect in rects:
        hs = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%.4f' %hs, ha='center', va='bottom')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(['HIT', 'OLD', 'CROSSED', 'DBORNE'])
    ax.set_title('All pixels')
    plt.tight_layout()
    figname_2 = '%s/parcel-dest_all_div' %plotDir
    print(figname_2)
    fig.savefig(figname_2 + '.png')
    if 'oem-Latitude-5400' in socket.gethostname():
        fig.show()
    if dc:
        #--- 3 ---    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects = ax.bar(range(1, 7), [tot_cld_TOT, tot_cld_HIT, tot_cld_OLD, tot_cld_CROSSED, tot_cld_DBORNE, tot_cld_DEAD])
        #: Add counts above the bars
        for rect in rects:
            hs = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%d' %hs, ha='center', va='bottom')
        ax.set_xticks(range(1, 7))
        ax.set_xticklabels(['TOT', 'HIT', 'OLD', 'CROSSED', 'DBORNE', 'DEAD'], fontsize='large')
        ax.set_title('Cloudy pixels')
        plt.tight_layout()
        figname_3 = '%s/parcel-dest_cld_tot' %plotDir
        print(figname_3)
        fig.savefig(figname_3 + '.png')
        if 'oem-Latitude-5400' in socket.gethostname():
            fig.show()
        #--- 4 ---
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects = ax.bar(range(1, 6), [tot_cld_HIT / tot_cld_TOT, tot_cld_OLD / tot_cld_TOT, tot_cld_CROSSED / tot_cld_TOT, tot_cld_DBORNE / tot_cld_TOT, tot_cld_DEAD / tot_cld_TOT])
        #: Add counts above the bars
        for rect in rects:
            hs = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%.4f' %hs, ha='center', va='bottom')
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xticks(range(1, 6))
        ax.set_xticklabels(['HIT', 'OLD', 'CROSSED', 'DBORNE', 'DEAD'])
        ax.set_title('Cloudy pixels')
        plt.tight_layout()
        figname_4 = '%s/parcel-dest_cld_div' %plotDir
        print(figname_4)
        fig.savefig(figname_4 + '.png')
        if 'oem-Latitude-5400' in socket.gethostname():
            fig.show()
    
    
    pdb.set_trace()
    
    
def plotAreas(pD):

    # colours = {'asia': 'y', 'asian monsoon': 'b', 'anticyclone (AMA)': 'r', 'pacific': 'g', 'central america': 'k'}#, 'atlantic', 'africa', 'warm pool', 'nino3', 'nino4', 'nino34'}
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.coastlines()
    i = -1
    for area, latlon in areas.items():
        figa = plt.figure()
        axa = figa.add_subplot(111,projection=ccrs.PlateCarree())
        axa.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
        axa.coastlines()
        if area in ['nino34']:
            continue
        # if area not in ['asia']:
        #     continue
        # [*mcolors.TABLEAU_COLORS.keys()][i]
        if i < len([*mcolors.TABLEAU_COLORS.keys()]):
            col = mcolors.TABLEAU_COLORS[[*mcolors.TABLEAU_COLORS.keys()][i]]
        else:
            col = mcolors.BASE_COLORS[[*mcolors.BASE_COLORS.keys()][i-10]]
        # if area not in colours.keys():
        #     continue
        i = i + 1
        if 'minLon1' in latlon.keys():
            #: --- TOT Map ---
            #: Left lower
            ax.plot([latlon['minLon1'], latlon['maxLon1']], [latlon['minLat1'], latlon['minLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Left upper
            ax.plot([latlon['minLon1'], latlon['maxLon1']], [latlon['maxLat1'], latlon['maxLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Right lower            
            ax.plot([latlon['minLon2'], latlon['maxLon2']], [latlon['maxLat2'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col)
            #: Right upper
            ax.plot([latlon['minLon2'], latlon['maxLon2']], [latlon['minLat2'], latlon['minLat2']], transform=ccrs.PlateCarree(), color=col)
            
            #: Left boundary
            ax.plot([latlon['minLon1'], latlon['minLon1']], [latlon['minLat1'], latlon['maxLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Middle boundary
            ax.plot([latlon['maxLon1'], latlon['maxLon1']], [latlon['maxLat1'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col)
            ax.plot([latlon['maxLon1'], latlon['maxLon1']], [latlon['minLat1'], latlon['minLat2']], transform=ccrs.PlateCarree(), color=col)
            #: Right boundary
            ax.plot([latlon['maxLon2'], latlon['maxLon2']], [latlon['minLat2'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col, label=area)
            #: --- Single Map ---
            #: Left lower
            axa.plot([latlon['minLon1'], latlon['maxLon1']], [latlon['minLat1'], latlon['minLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Left upper
            axa.plot([latlon['minLon1'], latlon['maxLon1']], [latlon['maxLat1'], latlon['maxLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Right lower            
            axa.plot([latlon['minLon2'], latlon['maxLon2']], [latlon['maxLat2'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col)
            #: Right upper
            axa.plot([latlon['minLon2'], latlon['maxLon2']], [latlon['minLat2'], latlon['minLat2']], transform=ccrs.PlateCarree(), color=col)
            
            #: Left boundary
            axa.plot([latlon['minLon1'], latlon['minLon1']], [latlon['minLat1'], latlon['maxLat1']], transform=ccrs.PlateCarree(), color=col)
            #: Middle boundary
            axa.plot([latlon['maxLon1'], latlon['maxLon1']], [latlon['maxLat1'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col)
            axa.plot([latlon['maxLon1'], latlon['maxLon1']], [latlon['minLat1'], latlon['minLat2']], transform=ccrs.PlateCarree(), color=col)
            #: Right boundary
            axa.plot([latlon['maxLon2'], latlon['maxLon2']], [latlon['minLat2'], latlon['maxLat2']], transform=ccrs.PlateCarree(), color=col, label=area)
        else:
            #: --- TOT Map ---
            #: Lower
            ax.plot([latlon['minLon'], latlon['maxLon']], [latlon['minLat'], latlon['minLat']], transform=ccrs.PlateCarree(), color=col)
            #: Upper
            ax.plot([latlon['minLon'], latlon['maxLon']], [latlon['maxLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col)
            #: Left boundary
            ax.plot([latlon['minLon'], latlon['minLon']], [latlon['minLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col)
            #: Right boundary
            ax.plot([latlon['maxLon'], latlon['maxLon']], [latlon['minLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col, label=area)
            #: --- Single Map ---
            axa.plot([latlon['minLon'], latlon['maxLon']], [latlon['minLat'], latlon['minLat']], transform=ccrs.PlateCarree(), color=col)
            #: Upper
            axa.plot([latlon['minLon'], latlon['maxLon']], [latlon['maxLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col)
            #: Left boundary
            axa.plot([latlon['minLon'], latlon['minLon']], [latlon['minLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col)
            #: Right boundary
            axa.plot([latlon['maxLon'], latlon['maxLon']], [latlon['minLat'], latlon['maxLat']], transform=ccrs.PlateCarree(), color=col, label=area)

        # axa.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # axa.legend(loc='lower left', bbox_to_anchor=(0, 1.01))
        # axa.legend(loc='lower right', bbox_to_anchor=(1, 1.01))
        axa.legend(loc='upper left', bbox_to_anchor=(0, -0.02))
        figa.tight_layout()
        if 'oem-Latitude-5400' in socket.gethostname():
            figa.show()
        fignamea = '%s/area-map_%s' %(pD, area.replace(' ', '-').replace('(', '').replace(')', ''))
        figa.savefig(fignamea + '.png')
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.legend(loc='upper left', bbox_to_anchor=(0, -0.02), fancybox=True, ncol=2)
    fig.tight_layout()
    if 'oem-Latitude-5400' in socket.gethostname():
        fig.show()
    figname = '%s/area-map' %(pD)
    print(figname)
    fig.savefig(figname + '.png')
    
def plotGridsat(pd):
    from convsrcErikFullGridSatTropo import read_sat, read_ERA5
    import datetime
    sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
    from ECMWF_N import ECMWF  # @UnresolvedImport
    import geosat  # @UnresolvedImport
    date_end = datetime.datetime(year=2008, month=1, day=1, hour=0)
    current_time = date_end
    dstep = datetime.timedelta(hours=0)
    dtRange = datetime.timedelta(hours=0)
    dat0 = geosat.GridSat(current_time)
    dat0._get_IR0()
    dat0.var['IR0'][dat0.var['IR0']<0] = 9999
    pdb.set_trace()
    
    dat0.close()
    # remove dat and make it a view of dat0, try to avoid errors on first attempt
    print('read GridSat for ',current_time)
    get_sat = read_sat(current_time,dtRange,pre=True,vshift=0)
    get_ERA5 = read_ERA5(current_time,dtRange,pre=True)
    pdb.set_trace()
if __name__ == '__main__':
    if 'ciclad' in socket.gethostname():
        datPath = os.environ['HOME'].replace('/home/', '/data/')
        ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson'
        mainDir = '%s/flexout/STC/Calipso' %ekjDir
        outDir = '%s/flexout/STC/Calipso-OUT' %datPath
        plotDir = '%s/LMD-Traczilla/Calipso/Plots/Extra' %ekjDir
        tempDir = os.path.join(mainDir, 'TempFiles/Part')
        
    elif 'oem-Latitude-5400' in socket.gethostname():
        mainDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso'
        outDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso/Calipso-OUT'
        tempDir = os.path.join(mainDir, 'TempFiles')
        plotDir = os.path.join(mainDir, 'Plots', 'Extra')
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)
    
    

    
    plotAreas(plotDir)
    pdb.set_trace()
    plotFlexpVar(outDir, plotDir, mainDir, dc=False)
    pdb.set_trace()
    plotGridsat(plotDir)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            