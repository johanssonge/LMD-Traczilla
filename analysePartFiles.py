#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2022-02-24

Copyright (c) 2022 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""

import numpy as np
import os
import sys
import pdb
import datetime
import time
import psutil
from matplotlib import pyplot as plt  # @UnresolvedImport

from analyseTraczilla import getConvfiles, I_HIT, I_OLD, readCatalogFile, getCatalogFile
from convsrcErikFullGridSatTropo import satratio
import glob
import h5py



sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readpart107, readidx107

def convIWC(wc, p, T, qv):
    #: Convert WC from kg m-3 to kg /kg and to ppmv
    Md = 28.97
    Rv = 461
    Rd = 287.04
    #: https://www.gastec.co.jp/en/technology/knowledge/concentration/
    ppm1 = (wc) * (22.4 / Md) * (T / 273) * (101300 / p)
    
    #: https://www.cdc.gov/niosh/docs/2004-101/calc.html
    ro2 = p / (Rd * T)
    wc_kgkg2 = wc * 1/ro2
    # wc_mgkg = wc_kgkg * 1000000
    ppm2 = wc_kgkg2 * 22.41 / Md
    
    Tv = T * (1 + (Rv/Rd - 1) * qv)
    ro = p / (Rd * Tv)
    wc_kgkg = wc * 1/ro
    ppm = wc_kgkg * 1 / 0.622
    return ppm
    
def esati_murphy(p, T):
    #ei in Pa saturation vapor pressure with respect to hexagonal (most stable) ice
    #Murphy and Koop 2005, QJRMS
    lnP=9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T
    
    esati_murphy=np.exp(lnP) / p
    #ppmv
    #kg
    
    return esati_murphy


def findMin(retv, ppVF, p0IF, ppIF, use_inds, ft):
    pinds = ppIF
    temp_full = np.zeros(p0IF.shape[0]) + np.nan
    temp_full[pinds] = ppVF
    p_min = np.nanmin((retv[0,:], temp_full[use_inds]), axis=0)
    # n_nan_inds = (~np.isnan(p_min))
    # qs_min[1,:][n_nan_inds] = np.where(np.nanargmin((qs_min[0,:][n_nan_inds], temp_full[hits][n_nan_inds]), axis=0) == 1, phour, qs_min[1,:][n_nan_inds])
    retv[1,:] = np.where(np.nanargmin((retv[0,:], temp_full[use_inds]), axis=0) == 1, ft, retv[1,:])
    retv[0,:] = p_min
    return retv

def findQschange(qo, qf, qsr, qsnr, p0IF, ppIF, use_inds, ft):
    ratio = 0.5
    # ppIF_full = np.zeros(p0IF.shape[0]).astype(bool)
    # ppIF_full[ppIF] = True
    qf_full = np.zeros(p0IF.shape[0]) + np.nan
    qf_full[ppIF] = qf
    
    dp = qf_full[use_inds] / qo
    nanind = ~np.isnan(dp)
    
    abovind = nanind & (dp > ratio) & (qsr != 1)
    beloind = nanind & (dp <= ratio) & (qs != -1)
    qsr[abovind] = 1
    qsr[beloind] = -1
    qsnr[abovind] = qsnr[abovind] + 1
    qsnr[beloind] = qsnr[beloind] + 1
    return qsr, qsnr

def findFirstQschange(qo, qf, q06, q10, q12, q14, q16, p0IF, ppIF, tss, ageb, use_inds, ft):
    qf_full = np.zeros(p0IF.shape[0]) + np.nan
    qf_full[ppIF] = qf
    
    tss_full = np.zeros(p0IF.shape[0]) + np.nan
    tss_full[ppIF] = tss
    tss_use = tss_full[use_inds]

    dq = qo / qf_full[use_inds]
    nanind = ~np.isnan(dq)

    # age_ind = tss_use <= ageb
    old = (ageb < tss_use) & nanind
    
    #: Old directly
    q16[0, old & np.isnan(q16[0, :])] = 3
    q14[0, old & np.isnan(q14[0, :])] = 3
    q12[0, old & np.isnan(q12[0, :])] = 3
    q10[0, old & np.isnan(q10[0, :])] = 3
    q06[0, old & np.isnan(q06[0, :])] = 3
    #: Old Without going below
    q16[0, old & (q16[0, :] == 1)] = 2
    q14[0, old & (q14[0, :] == 1)] = 2
    q12[0, old & (q12[0, :] == 1)] = 2
    q10[0, old & (q10[0, :] == 1)] = 2
    q06[0, old & (q06[0, :] == 1)] = 2
    #: Old in any way
    # old16 = q16[0, :] > 1
    # old14 = q14[0, :] > 1
    # old12 = q12[0, :] > 1
    # old10 = q10[0, :] > 1
    # old06 = q06[0, :] > 1
    
    #: nanind & greater than th & nan :not larger than 1 i.e. old
    abov16 = nanind & (dq >= 1.6) & (np.isnan(q16[0, :]))# | (q16[0, :] == 1))
    abov14 = nanind & (dq >= 1.4) & (np.isnan(q14[0, :]))# | (q14[0, :] == 1))
    abov12 = nanind & (dq >= 1.2) & (np.isnan(q12[0, :]))# | (q12[0, :] == 1))
    abov10 = nanind & (dq >= 1.0) & (np.isnan(q10[0, :]))# | (q10[0, :] == 1))
    abov06 = nanind & (dq >= 0.6) & (np.isnan(q06[0, :]))# | (q06[0, :] == 1))
    
    #: Not above but still non nan
    #: nanind & less than th & not larger than 1 i.e. old
    below16 = nanind & (dq < 1.6) & (~(q16[0, :] > 1))
    below14 = nanind & (dq < 1.4) & (~(q14[0, :] > 1))
    below12 = nanind & (dq < 1.2) & (~(q12[0, :] > 1))
    below10 = nanind & (dq < 1.0) & (~(q10[0, :] > 1))
    below06 = nanind & (dq < 0.6) & (~(q06[0, :] > 1))
    
    q16[0, abov16] = 1 #dq[abov16]
    q14[0, abov14] = 1 #dq[abov14]
    q12[0, abov12] = 1 #dq[abov12]
    q10[0, abov10] = 1 #dq[abov10]
    q06[0, abov06] = 1 #dq[abov06]
    
    q16[0, below16] = -1 #dq[abov16]
    q14[0, below14] = -1 #dq[abov14]
    q12[0, below12] = -1 #dq[abov12]
    q10[0, below10] = -1 #dq[abov10]
    q06[0, below06] = -1 #dq[abov06]
    
    
    #: Save time for first time below   
    q16[1, (below16 & np.isnan(q16[1, :]))] = ft
    q14[1, (below14 & np.isnan(q14[1, :]))] = ft
    q12[1, (below12 & np.isnan(q12[1, :]))] = ft
    q10[1, (below10 & np.isnan(q10[1, :]))] = ft
    q06[1, (below06 & np.isnan(q06[1, :]))] = ft
    
    return q06, q10, q12, q14, q16
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    
    maindDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso'
    outDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso/Calipso-OUT'
    tempDir = os.path.join(maindDir, 'Tempfiles')
    plotDir = os.path.join(maindDir, 'Plots')
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)
    
    # outDir = '/data/ejohansson/flexout/STC/Calipso-OUT'
    forced_age_bound = None#5 * 86400 #: Seconds
    stop_day = 2 * 86400 #: Seconds
    
    
    
    years = [2008]
    months = [1]
    ticT = time.time()
    for year in years:
        for mon in months:
            tempname = 'part--%d%02d.h5' %(year, mon)
            outname = 'CALIOP-EAD-%d%02d-n-DD' %(year, mon)
            # partDir = os.path.join('/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/STC/Calipso', outname)
            partDir = os.path.join(maindDir, outname)
            initDir = os.path.join(partDir, 'Initfiles')
            catalogFile = os.path.join(initDir, 'selDardar_Calalog-%d%02d-n.pkl' %(year, mon))
            # te = readCatalogFile(partDir + '/Initfiles/selDardar_Calalog-200701-n.pkl')
            #: Read catalogfile
            cf = getCatalogFile(catalogFile, checkForNewFile=False)
            #: Read out file
            rvs, fls, age, lons, lats, temp, pres = getConvfiles(outDir, [[''], [[outname]]])
            hits = (fls & I_HIT) == I_HIT
            if forced_age_bound is not None:
                hits = np.where(age > forced_age_bound, False, hits)
            
            all_inds = np.ones(hits.shape).astype(bool)
            use_inds = hits
            #: Read part_000 file
            f0 = os.path.join(partDir, 'part_000')
            part0 = readidx107(f0, quiet=True)
            
            if stop_day is not None:
                use_inds = np.where((np.abs(part0['ir_start']) > stop_day), False, use_inds)
            
            use_inds = use_inds & (cf['iwc'] > 0) 
            #: Calc saturation
            # rvs0O = satratio(part0['p'][use_inds], part0['t'][use_inds])
            rvs0 = esati_murphy(part0['p'][use_inds], part0['t'][use_inds])
            iwc0 = convIWC(cf['iwc'][use_inds], part0['p'][use_inds], part0['t'][use_inds], cf['sh'][use_inds])
            
            #: get the desired variables and put in a 2D array
            pMin = np.zeros([2, hits.sum()])
            pMin[0, :] = part0['p'][hits]
            rvsMin = np.zeros([2, use_inds.sum()])
            rvsMin[0, :] = rvs0
            # q0 = rvs0 + iwc0
            q0 = rvs0 * 1.6
            q100 = np.zeros(use_inds.sum())
            qs = np.zeros(use_inds.sum())
            qsn = np.zeros(use_inds.sum())
            
            qs06 = np.zeros([2, use_inds.sum()]) + np.nan
            qs08 = np.zeros([2, use_inds.sum()]) + np.nan
            qs10 = np.zeros([2, use_inds.sum()]) + np.nan
            qs12 = np.zeros([2, use_inds.sum()]) + np.nan
            qs14 = np.zeros([2, use_inds.sum()]) + np.nan
            qs16 = np.zeros([2, use_inds.sum()]) + np.nan
            
            
            
            # qs = [[]] * use_inds.sum()
            
            #hmax = 18
            
            pfiles = glob.glob(partDir + '/part_*')
            pfi = []
            for pf in pfiles:
                pfi.append(int(pf.split('part_')[-1].replace('.gz', '')))
            pind = np.argsort(pfi)
            pfiles = np.asarray(pfiles)[pind][1:]
            step = np.asarray(pfi)[pind][1:][0].item()
            hmax = np.asarray(pfi)[pind][1:][-1].item()
            dstep = datetime.timedelta(hours=step)
            print('hmax=%d' %hmax)
            for phour in range(step, hmax + 1, step):
                # if phour < 1749:
                #     continue
                tic = time.time()
                #: Read part_xxx file
                partfile = f0.replace('part_000', 'part_%03d' %phour)
                if not os.path.isfile(partfile + '.gz'):
                    print("No such file")
                    print(partfile)
                    print('')
                    continue
                
                partf = readidx107(partfile, quiet=True)
                print("Read time = %.02f s" %(time.time() - tic))
                if len(partf['idx_back']) == 0:
                    continue
                #: Calc saturation
                # rvsf = satratio(partf['p'], partf['t'])
                rvsf = esati_murphy(partf['p'], partf['t'])
                #: Check memory
                pid = os.getpid()
                py = psutil.Process(pid)
                memoryUse = py.memory_info()[0]/2**30
                print('memory use: {:4.2f} gb'.format(memoryUse))
                print('')
                # temp_full = np.zeros(part0['idx_back'].shape[0]) + np.nan
                # pinds = partf['idx_back']-1
                # temp_full[pinds] = partf['p']
                # p_min = np.nanmin((qs_min[0,:], temp_full[hits]), axis=0)
                # # n_nan_inds = (~np.isnan(p_min))
                # # qs_min[1,:][n_nan_inds] = np.where(np.nanargmin((qs_min[0,:][n_nan_inds], temp_full[hits][n_nan_inds]), axis=0) == 1, phour, qs_min[1,:][n_nan_inds])
                # qs_min[1,:] = np.where(np.nanargmin((qs_min[0,:], temp_full[hits]), axis=0) == 1, phour, qs_min[1,:])
                # qs_min[0,:] = p_min
                #: Find Min
                # (q0, rvsf, part0['idx_back']-1, partf['idx_back']-1, use_inds, phour, qs, q100) 
                # rvsMin = findMin(rvsMin, rvsf, part0['idx_back']-1, partf['idx_back']-1, use_inds, phour)
                # pMin = findMin(pMin, part0['idx_back']-1, partf['idx_back']-1, hits, partf['p'], phour)
                #: Find Change
                # qs, qsn = findQschange(q0, rvsf, qs, qsn, (part0['idx_back'] - part0['idx_orgn']), (partf['idx_back'] - part0['idx_orgn']), use_inds, phour)
                if phour == 3:
                    use_inds_3 = np.zeros(use_inds.shape[0]).astype(bool)
                    use_inds_3[(partf['idx_back'] - part0['idx_orgn'])] = True
                
                qs06, qs10, qs12, qs14, qs16 = findFirstQschange(q0, rvsf, qs06, qs10, qs12, qs14, qs16, (part0['idx_back'] - part0['idx_orgn']), (partf['idx_back'] - part0['idx_orgn']), np.abs(partf['itime']-partf['ir_start']), age[use_inds], use_inds, phour)
                print("Calc time = %.02f s" %(time.time() - tic))
                #: Check memory
                pid = os.getpid()
                py = psutil.Process(pid)
                memoryUse = py.memory_info()[0]/2**30
                print('memory use: {:4.2f} gb'.format(memoryUse))
                # if (np.isnan(qs10[0, :]).sum() + np.isnan(qs12[0, :]).sum() + np.isnan(qs14[0, :]).sum() + np.isnan(qs16[0, :]).sum()) == 0:
                #     print('No more values to find')
                #     pdb.set_trace()
                #     break
                
                #: Match Indices
                # mi = np.asarray(list((set(partf['idx_back']) & set(part0['idx_back'][hits]))))
                # inds = np.searchsorted(part0['idx_back'][hits], mi)
                # pdb.set_trace()
                
                # if np.abs(partf['itime']) > (forced_age_bound + stop_day):
                #     break
                
                # idx_full = np.zeros(part0['idx_back'].shape[0]).astype(bool)
                #
                # idx_full[(partf['idx_back'] - part0['idx_orgn'])] = True
                #
                #
                # if idx_full[use_inds].sum() == 0:
                #     print('0')
                #     pdb.set_trace()
                if phour >150:
                    
                    break
            # f2r = f2.replace('.gz', '')
            h5file = h5py.File(tempname, 'w')
            h5file.create_dataset('qs06', data = qs06)
            h5file.create_dataset('qs10', data = qs10)
            h5file.create_dataset('qs12', data = qs12)
            h5file.create_dataset('qs14', data = qs14)
            h5file.create_dataset('qs16', data = qs16)
            h5file.create_dataset('age', data = age[use_inds])
            h5file.create_dataset('ir_start', data = part0['ir_start'][use_inds])
            h5file.create_dataset('use_inds', data = use_inds)
            h5file.close()
    
    heightBoundaries = [[None, None], [14,15], [15,16], [14,16], [16, 18], [18, 20], [18,19], [19,20]]
    
    title_org = 'Histogram'
    inds_org = (~np.isnan(qs16[0,:])).astype(bool)
    figname_1_org = '%s/qsn-div_total' %plotDir
    figname_2_org = '%s/qsn-total' %plotDir
    for h1, h2 in heightBoundaries:
        if not ((h1 is None) and (h2 is None)):
            title = title_org + ', %d - %d km' %(h1, h2)
            figname_1 = figname_1_org + '_%d-%d' %(h1, h2)
            figname_2 = figname_2_org + '_%d-%d' %(h1, h2)
            if h2 == heightBoundaries[-1][1]:
                inds_h = (cf['height'][use_inds] >= h1) & (cf['height'][use_inds] <= h2)
            else:
                inds_h = (cf['height'][use_inds] >= h1) & (cf['height'][use_inds] < h2)
            inds = inds_org & inds_h
        else:
            title = title_org
            figname_1 = figname_1_org
            figname_2 = figname_2_org
            inds = inds_org
    
        if inds.sum() == 0:
            continue
        abowe_until_hit_16 = (qs16[0, inds]==2).sum()
        abowe_until_hit_14 = (qs14[0, inds]==2).sum()
        abowe_until_hit_12 = (qs12[0, inds]==2).sum()
        abowe_until_hit_10 = (qs10[0, inds]==2).sum()
        abowe_until_hit_06 = (qs06[0, inds]==2).sum()
        
        
        abowe_no_hit_16 = (qs16[0, inds]==1).sum()
        abowe_no_hit_14 = (qs14[0, inds]==1).sum()
        abowe_no_hit_12 = (qs12[0, inds]==1).sum()
        abowe_no_hit_10 = (qs10[0, inds]==1).sum()
        abowe_no_hit_06 = (qs06[0, inds]==1).sum()
        
        below_before_hit_16 = (qs16[0, inds]==-1).sum()
        below_before_hit_14 = (qs14[0, inds]==-1).sum()
        below_before_hit_12 = (qs12[0, inds]==-1).sum()
        below_before_hit_10 = (qs10[0, inds]==-1).sum()
        below_before_hit_06 = (qs06[0, inds]==-1).sum()
        
        mean_time_16 = np.mean(qs16[1, inds][qs16[0, inds]==-1])
        mean_time_14 = np.mean(qs14[1, inds][qs14[0, inds]==-1])
        mean_time_12 = np.mean(qs12[1, inds][qs12[0, inds]==-1])
        mean_time_10 = np.mean(qs10[1, inds][qs10[0, inds]==-1])
        mean_time_06 = np.mean(qs06[1, inds][qs06[0, inds]==-1])
        
        total_16 = abowe_until_hit_16 + abowe_no_hit_16 + below_before_hit_16
        total_14 = abowe_until_hit_14 + abowe_no_hit_14 + below_before_hit_14
        total_12 = abowe_until_hit_12 + abowe_no_hit_12 + below_before_hit_12
        total_10 = abowe_until_hit_10 + abowe_no_hit_10 + below_before_hit_10
        total_06 = abowe_until_hit_06 + abowe_no_hit_06 + below_before_hit_06
        if inds.sum() != total_16:
            print('Wrong tot')
            pdb.set_trace()
        
        # ind06 = ~(np.isnan(qs06[0,:]) | (qs06[0,:]==-1))
        # ind10 = ~(np.isnan(qs10[0,:]) | (qs10[0,:]==-1))
        # ind12 = ~(np.isnan(qs12[0,:]) | (qs12[0,:]==-1))
        # ind14 = ~(np.isnan(qs14[0,:]) | (qs14[0,:]==-1))
        # ind16 = ~(np.isnan(qs16[0,:]) | (qs16[0,:]==-1))
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        labels = ['0.6', '1.0', '1.2', '1.4', '1.6']
        x = np.arange(len(labels))  # the label locations
        width = 0.3  # the width of the bars
        # rects1 = ax.bar(x - width/2, [(qs06[0,:]==-1).sum() / (~(np.isnan(qs06[0,:]))).sum(), (qs10[0,:]==-1).sum() / (~(np.isnan(qs10[0,:]))).sum(), (qs12[0,:]==-1).sum() / (~(np.isnan(qs12[0,:]))).sum(), (qs14[0,:]==-1).sum() / (~(np.isnan(qs14[0,:]))).sum(), (qs16[0,:]==-1).sum() / (~(np.isnan(qs16[0,:]))).sum()], width, label='Below')
        # rects2 = ax.bar(x + width/2, [ind06.sum() / (~(np.isnan(qs06[0,:]))).sum(), ind10.sum() / (~(np.isnan(qs10[0,:]))).sum(), ind12.sum() / (~(np.isnan(qs12[0,:]))).sum(), ind14.sum() / (~(np.isnan(qs14[0,:]))).sum(), ind16.sum() / (~(np.isnan(qs16[0,:]))).sum()], width, label='Convective')
        # rects1 = ax.bar(x - width/2, [(~ind06).sum() / (qs06[0,:]).shape[0], (~ind10).sum() / (qs10[0,:]).shape[0], (~ind12).sum() / (qs12[0,:]).shape[0], (~ind14).sum() / (qs14[0,:]).shape[0], (~ind16).sum() / (qs16[0,:]).shape[0]], width, label='Incitu')
        # rects2 = ax.bar(x + width/2, [ind06.sum() / (qs06[0,:]).shape[0], ind10.sum() / (qs10[0,:]).shape[0], ind12.sum() / (qs12[0,:]).shape[0], ind14.sum() / (qs14[0,:]).shape[0], ind16.sum() / (qs16[0,:]).shape[0]], width, label='Convective')
        # rects3 = ax.bar(x + width/2, [np.mean(qs10[1, ind10]), np.mean(qs12[1, ind12]), np.mean(qs14[1, ind14]), np.mean(qs16[1, ind16])], width, label='Mean Time [h]')
        
        rects1 = ax.bar(x - width/2, [abowe_until_hit_06 / total_06, abowe_until_hit_10 / total_10, abowe_until_hit_12 / total_12, abowe_until_hit_14 / total_14, abowe_until_hit_16 / total_16], width/2, label='Convective, with hit')
        rects2 = ax.bar(x, [abowe_no_hit_06 / total_06, abowe_no_hit_10 / total_10, abowe_no_hit_12 / total_12, abowe_no_hit_14 / total_14, abowe_no_hit_16 / total_16], width/2, label='Convective, no hit')
        rects3 = ax.bar(x + width/2, [below_before_hit_06 / total_06, below_before_hit_10 / total_10, below_before_hit_12 / total_12, below_before_hit_14 / total_14, below_before_hit_16 / total_16], width/2, label='In situ')
        #-1.5 0 1.5
        
        
        # Add counts above the two bar graphs
        # for rect in rects2:
        #     height = rect.get_height()
        #     plt.text(rect.get_x() + rect.get_width() / 2.0, height, '', ha='center', va='bottom')#f'{height:.0f}', ha='center', va='bottom')
    
        
        ax.text
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel('Thresholds')
        ax.set_ylabel('Num / Total')
        ax.set_title(title)
        ax.legend()
        # ax.set_title('Num of change')
        
        plt.tight_layout()
        print(figname_1)
        fig.savefig(figname_1 + '.png')
        fig.show()
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        labels = ['0.6', '1.0', '1.2', '1.4', '1.6']
        x = np.arange(len(labels))  # the label locations
        width = 0.3  # the width of the bars
        
        rects1 = ax.bar(x - width/2, [abowe_until_hit_06, abowe_until_hit_10, abowe_until_hit_12, abowe_until_hit_14, abowe_until_hit_16], width/2, label='Convective, with hit')
        rects2 = ax.bar(x, [abowe_no_hit_06, abowe_no_hit_10, abowe_no_hit_12, abowe_no_hit_14, abowe_no_hit_16], width/2, label='Convective, no hit')
        rects3 = ax.bar(x + width/2, [below_before_hit_06, below_before_hit_10, below_before_hit_12, below_before_hit_14, below_before_hit_16], width/2, label='In situ')
        
        ax.text
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel('Thresholds')
        ax.set_ylabel('Num')
        ax.set_title(title)
        ax.legend()
        # ax.set_title('Num of change')
        
        plt.tight_layout()
        print(figname_2)
        fig.savefig(figname_2 + '.png')
        fig.show()
        
    
    
    
    
    
    # plt.close(fig)
    pdb.set_trace()
    # h = ax.hist(qsn)
    # counts = np.bincount(qsn.astype(int))
    # ax.bar(range(len(counts)), counts, width=1, align='center')
    # ax.set_xticks(range(len(counts)))
    # ax.set_title('Num of change')
    # plt.tight_layout()
    # figname = '%s/qsn' %plotDir
    # fig.savefig(figname + '.png')
    # fig.show()
    # plt.close(fig)
    
    
    
    
    
    print("total time = %.02f s" %(time.time() - ticT))
    pdb.set_trace()
