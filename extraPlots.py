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
from matplotlib import pyplot as plt  # @UnresolvedImport
from analyseTraczilla import getConvfiles, missing_months, readCatalogFile,\
    getCatalogFile
from convsrcErikFullGridSatTropo import I_HIT, I_OLD, I_DEAD, I_CROSSED, I_DBORNE
import pdb



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
            rvs_mon, fls_mon, age_mon, lons_monF, lats_monF, temp_monF, pres_monF = getConvfiles(outDir, [[''], [[outname]]])  # @UnusedVariable
            
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
    rects = ax.bar(range(1, 7), [tot_TOT, tot_HIT, tot_OLD, tot_CROSSED, tot_DBORNE, tot_DEAD])
    #: Add counts above the bars
    for rect in rects:
        hs = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%d' %hs, ha='center', va='bottom')
    ax.set_xticks(range(1, 7))
    ax.set_xticklabels(['TOT', 'HIT', 'OLD', 'CROSSED', 'DBORNE', 'DEAD'], fontsize='large')
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
    rects = ax.bar(range(1, 6), [tot_HIT / tot_TOT, tot_OLD / tot_TOT, tot_CROSSED / tot_TOT, tot_DBORNE / tot_TOT, tot_DEAD / tot_TOT])
    #: Add counts above the bars
    for rect in rects:
        hs = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%.4f' %hs, ha='center', va='bottom')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticks(range(1, 6))
    ax.set_xticklabels(['HIT', 'OLD', 'CROSSED', 'DBORNE', 'DEAD'])
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
    
    
    
    
    plotFlexpVar(outDir, plotDir, mainDir, dc=False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            