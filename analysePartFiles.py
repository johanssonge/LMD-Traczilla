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
import psutil  # @UnresolvedImport
from matplotlib import pyplot as plt  # @UnresolvedImport
from matplotlib import colors  # @UnresolvedImport
import cartopy.crs as ccrs  # @UnresolvedImport
# from convsrcErikFullGridSatTropo import satratio
import glob
import h5py  # @UnresolvedImport
import socket



# sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readidx107
from analyseTraczilla import getConvfiles, getCatalogFile, missing_months, seasons, areas,\
    esati_murphy, convIWC, getAreaInds, getAreaName, getYmax
from convsrcErikFullGridSatTropo import I_HIT, I_OLD, I_DEAD, I_CROSSED, I_DBORNE

# Calculation of the saturation mixing ratio from actual temperature and pressure
def satratio2(p,T):
    """ Calculate the mass saturation ratio from pressure (in Pa) and temperature 
    (in K). Output in kg/kg """
    estar = 1.0008*np.exp(23.33086-(6111.72784/T)+0.15215*np.log(T))
    satr = 0.622 * estar/(0.01*p-estar)
    return satr


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
    beloind = nanind & (dp <= ratio) & (qsr != -1)
    qsr[abovind] = 1
    qsr[beloind] = -1
    qsnr[abovind] = qsnr[abovind] + 1
    qsnr[beloind] = qsnr[beloind] + 1
    return qsr, qsnr

def FirstQschange(q,o, d, ni, h, th):
    q[0, o & np.isnan(q[0, :])] = 3
    q[0, o & (q[0, :] == 1)] = 2
    above = ni & (d >= th) & (np.isnan(q[0, :]))# | (q16[0, :] == 1))
    below = ni & (d < th) & (~(q[0, :] > 1))
    q[0, above] = 1
    q[0, below] = -1
    q[1, (below & np.isnan(q[1, :]))] = h
    return q
def findFirstQschange(qo, qf, q06, q10, q12, q13, q14, q16, p0IF, ppIF, tss, ageb, use_inds, ft):
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
    # q16[0, old & np.isnan(q16[0, :])] = 3
    # q14[0, old & np.isnan(q14[0, :])] = 3
    # q12[0, old & np.isnan(q12[0, :])] = 3
    # q10[0, old & np.isnan(q10[0, :])] = 3
    # q06[0, old & np.isnan(q06[0, :])] = 3
    #: Old Without going below
    # q16[0, old & (q16[0, :] == 1)] = 2
    # q14[0, old & (q14[0, :] == 1)] = 2
    # q12[0, old & (q12[0, :] == 1)] = 2
    # q10[0, old & (q10[0, :] == 1)] = 2
    # q06[0, old & (q06[0, :] == 1)] = 2
    #: Old in any way
    # # old16 = q16[0, :] > 1
    # # old14 = q14[0, :] > 1
    # # old12 = q12[0, :] > 1
    # # old10 = q10[0, :] > 1
    # # old06 = q06[0, :] > 1
    
    #: nanind & greater than th & nan :not larger than 1 i.e. old
    # abov16 = nanind & (dq >= 1.6) & (np.isnan(q16[0, :]))# | (q16[0, :] == 1))
    # abov14 = nanind & (dq >= 1.4) & (np.isnan(q14[0, :]))# | (q14[0, :] == 1))
    # abov12 = nanind & (dq >= 1.2) & (np.isnan(q12[0, :]))# | (q12[0, :] == 1))
    # abov10 = nanind & (dq >= 1.0) & (np.isnan(q10[0, :]))# | (q10[0, :] == 1))
    # abov06 = nanind & (dq >= 0.6) & (np.isnan(q06[0, :]))# | (q06[0, :] == 1))
    
    #: Not above but still non nan
    #: nanind & less than th & not larger than 1 i.e. old
    # below16 = nanind & (dq < 1.6) & (~(q16[0, :] > 1))
    # below14 = nanind & (dq < 1.4) & (~(q14[0, :] > 1))
    # below12 = nanind & (dq < 1.2) & (~(q12[0, :] > 1))
    # below10 = nanind & (dq < 1.0) & (~(q10[0, :] > 1))
    # below06 = nanind & (dq < 0.6) & (~(q06[0, :] > 1))
    
    # q16[0, abov16] = 1 #dq[abov16]
    # q14[0, abov14] = 1 #dq[abov14]
    # q12[0, abov12] = 1 #dq[abov12]
    # q10[0, abov10] = 1 #dq[abov10]
    # q06[0, abov06] = 1 #dq[abov06]
    
    # q16[0, below16] = -1 #dq[abov16]
    # q14[0, below14] = -1 #dq[abov14]
    # q12[0, below12] = -1 #dq[abov12]
    # q10[0, below10] = -1 #dq[abov10]
    # q06[0, below06] = -1 #dq[abov06]
    
    
    #: Save time for first time below   
    # q16[1, (below16 & np.isnan(q16[1, :]))] = ft
    # q14[1, (below14 & np.isnan(q14[1, :]))] = ft
    # q12[1, (below12 & np.isnan(q12[1, :]))] = ft
    # q10[1, (below10 & np.isnan(q10[1, :]))] = ft
    # q06[1, (below06 & np.isnan(q06[1, :]))] = ft
    
    q16 = FirstQschange(q16,old, dq, nanind, ft, 1.6)
    q14 = FirstQschange(q14,old, dq, nanind, ft, 1.4)
    q13 = FirstQschange(q13,old, dq, nanind, ft, 1.3)
    q12 = FirstQschange(q12,old, dq, nanind, ft, 1.2)
    q10 = FirstQschange(q10,old, dq, nanind, ft, 1.0)
    q06 = FirstQschange(q06,old, dq, nanind, ft, 0.6)
    # if not (q14[~np.isnan(q14)] == q10[~np.isnan(q14)]).all():
    #     print('dame')
    #     pdb.set_trace()
    # if not (q16[~np.isnan(q16)] == q12[~np.isnan(q16)]).all():
    #     print('dame')
    #     pdb.set_trace()
    
    return q06, q10, q12, q13, q14, q16
    

def calcProcHeight(q06, q10, q12, q13, q14, q16, ohClo_h, ohClr_h, height, rh=False):
    all_heights = [*range(14000, 20121, 180)]
    inds_org = (~np.isnan(q16[0,:])).astype(bool)
    conv = np.zeros([6, len(all_heights)-1])
    conv_tot = np.zeros([6, len(all_heights)-1])
    cf = np.zeros([1, len(all_heights)-1])
    for h in range(len(all_heights)-1):
        total_clo_t = (((ohClo_h*1000) >= all_heights[h]) & ((ohClo_h*1000) < all_heights[h+1])).sum()
        total_clr_t = (((ohClr_h*1000) >= all_heights[h]) & ((ohClr_h*1000) < all_heights[h+1])).sum()
        inds_h = ((height*1000) >= all_heights[h]) & ((height*1000) < all_heights[h+1])
        inds = inds_org & inds_h
        abowe_until_hit_06 = (q06[0, inds]==2).sum()
        abowe_until_hit_10 = (q10[0, inds]==2).sum()
        abowe_until_hit_12 = (q12[0, inds]==2).sum()
        abowe_until_hit_13 = (q13[0, inds]==2).sum()
        abowe_until_hit_14 = (q14[0, inds]==2).sum()
        abowe_until_hit_16 = (q16[0, inds]==2).sum()
    
        # abowe_no_hit = (qs13[0, inds]==1).sum()
        # below_before_hit = (qs13[0, inds]==-1).sum()
        # total = abowe_until_hit_13 + abowe_no_hit + below_before_hit
        # if total != inds.sum():
        #     print('strange sum')
        #     pdb.set_trace()
        
        total_t = total_clo_t + total_clr_t
        if total_clo_t == 0:
            total_clo = np.nan
        else:
            total_clo = total_clo_t
        if total_t == 0:
            total = np.nan
        else:
            total = total_t
        
        conv[0, h] = abowe_until_hit_06 / total_clo
        conv[1, h] = abowe_until_hit_10 / total_clo
        conv[2, h] = abowe_until_hit_12 / total_clo
        conv[3, h] = abowe_until_hit_13 / total_clo
        conv[4, h] = abowe_until_hit_14 / total_clo
        conv[5, h] = abowe_until_hit_16 / total_clo
        
        conv_tot[0, h] = abowe_until_hit_06
        conv_tot[1, h] = abowe_until_hit_10
        conv_tot[2, h] = abowe_until_hit_12
        conv_tot[3, h] = abowe_until_hit_13
        conv_tot[4, h] = abowe_until_hit_14
        conv_tot[5, h] = abowe_until_hit_16
        
        cf[0, h] = total_clo / total
        
    if rh:
        return conv, conv_tot, cf, all_heights
    else:
        return conv, conv_tot, cf

def addExtra(tn, f0n):
    h5f = h5py.File(tn, 'a')
    if 'lons0' in h5f.keys():
        h5f.close()
    else:
        pf = readidx107(f0n, quiet=True)
        pf['x'] = np.where(pf['x'] > 180, pf['x'] - 360, pf['x'])
        pf['x'] = np.where(pf['x'] < -180, pf['x'] + 360, pf['x'])
        lons0_mon = pf['x']
        lats0_mon = pf['y']
        h5f.create_dataset('lons0', data = lons0_mon[h5f['use_inds'][:]])
        h5f.create_dataset('lats0', data = lats0_mon[h5f['use_inds'][:]])
        h5f.create_dataset('ohClo_lons0', data = lons0_mon[h5f['ohClo_use_inds'][:]])
        h5f.create_dataset('ohClo_lats0', data = lats0_mon[h5f['ohClo_use_inds'][:]])
        h5f.create_dataset('ohClr_lons0', data = lons0_mon[h5f['ohClr_use_inds'][:]])
        h5f.create_dataset('ohClr_lats0', data = lats0_mon[h5f['ohClr_use_inds'][:]])
        h5f.close()
                        
def readTempFile(fn):
    h5f = h5py.File(fn, 'r')   
    q10 = h5f['qs10'][:]
    q12 = h5f['qs12'][:]
    q13 = h5f['qs13'][:]
    q14 = h5f['qs14'][:]
    q16 = h5f['qs16'][:]
    q06 = h5f['qs06'][:]
    h = h5f['height'][:]
    lat = h5f['lats'][:]
    lon = h5f['lons'][:]
    # 'ohClr_height'
    extras = {'use_inds': h5f['use_inds'][:], 'ohClo_use_inds': h5f['ohClo_use_inds'][:], 'ohClr_use_inds': h5f['ohClr_use_inds'][:], \
              'ohClo_height': h5f['ohClo_height'][:], 'ohClr_height': h5f['ohClr_height'][:], \
              'ohClo_lats': h5f['ohClo_lats'][:], 'ohClo_lons': h5f['ohClo_lons'][:], \
              'ohClo_lats0': h5f['ohClo_lats0'][:],'ohClo_lons0': h5f['ohClo_lons0'][:], \
              'ohClr_lats': h5f['ohClr_lats'][:], 'ohClr_lons': h5f['ohClr_lons'][:], \
              'ohClr_lats0': h5f['ohClr_lats0'][:],'ohClr_lons0': h5f['ohClr_lons0'][:], \
              'iwc0': h5f['iwc0'][:], 'age': h5f['age'][:], \
              'part0_p': h5f['part0_p'][:], 'part0_t': h5f['part0_t'][:], \
              'lons0': h5f['lons0'][:], 'lats0': h5f['lats0'][:], \
              'pressure': h5f['pressure'][:], 'temperature': h5f['temperature'][:]}

    h5f.close()
    
    return q06, q10, q12, q13, q14, q16, h, lat, lon, extras
    
    
    
    
    
if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Program to calculate backvards through part files to find supersaturation', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-c", "--ct", action='store_false', default = True, 
                        help = "Create Temp Files, only valid if no temp file is loaded. Default = True")
    parser.add_argument("-l", "--lt", action='store_false', default = True, 
                        help = "Load Temp Files. Default = True")
    parser.add_argument("-e", "--extra", action='store_true', default = False, 
                        help = "Add extra fields to tempfile. Default = False")
    parser.add_argument("-v","--vshift",type=int,choices=[0,10], default=10, 
                        help='vertical shift. Default = 10')
    parser.add_argument("-a", "--area", type=int, choices=([*range(20)]), default=0,  
                        help = "Areas. Default = 0\n")
    parser.add_argument("-y", "--year", type=int, choices=([0]+[*range(2007, 2020)]), default=0,  
                        help = "year. Default=2018")
    parser.add_argument("-m", "--month", type=int, choices=(np.arange(-6, 13)), default=0, 
                        help = "Month. Default=6")
    
    args = parser.parse_args()
    if ('ciclad' in socket.gethostname()) or ('spirit' in socket.gethostname()):
        datPath = os.environ['HOME'].replace('/home/', '/data/')
        ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson'
        mainDir = '%s/flexout/STC/Calipso' %ekjDir
        outDir = '%s/flexout/STC/Calipso-OUT' %datPath
        plotDir = '%s/LMD-Traczilla/Calipso/Plots' %ekjDir
        tempDir = os.path.join(mainDir, 'TempFiles/Part')
        
    elif 'oem-Latitude-5400' in socket.gethostname():
        mainDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso'
        outDir = '/home/ejohansson/Projects/LMD-Traczilla/Calipso/Calipso-OUT'
        tempDir = os.path.join(mainDir, 'TempFiles')
        plotDir = os.path.join(mainDir, 'Plots')

    elif ('climserv' in socket.gethostname()) | ('polytechnique' in socket.gethostname()):
        print('climserv or polytechnique in host name')
        sys.exit()
    else:
        print ('CANNOT RECOGNIZE HOST - DO NOT RUN ON NON DEFINED HOSTS')
        sys.exit()
    
    vshift = args.vshift
    if vshift != 10:
        outDir = outDir + '/Vshift_%d' %(vshift)
        tempDir = tempDir + '/Vshift_%d' %(vshift)
        plotDir = plotDir + '/Vshift_%d' %(vshift)
    
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)
    
    forced_age_bound = 7776000 #: 90 days #864000#5 * 86400 #: Seconds
    stop_day = None#2 * 86400 #: Seconds
    # break_point = 1000
    if forced_age_bound is None:
        break_point = None
    else:
        break_point = (31 + 1 + int(forced_age_bound/86400) * 24)
    

    lt = args.lt
    ct = args.ct
   
    area = getAreaName(args.area)
    print(area.title())
    # years = [2008]#, 2008, 2009, 2010]
    # months = [2]
    if args.month == 0:
        ses = 'year'
        ses_tit = ses.title()
        months = seasons[ses]
    elif args.month == -1:
        ses = 'djf'
        ses_tit = ses
        months = seasons[ses]
    elif args.month == -2:
        ses = 'mam'
        ses_tit = ses
        months = seasons[ses]
    elif args.month == -3:
        ses = 'jja'
        ses_tit = ses
        months = seasons[ses]
    elif args.month == -4:
        ses = 'son'
        ses_tit = ses
        months = seasons[ses]
    elif args.month == -5:
        ses = 'sum'
        ses_tit = 'Summer'
        months = seasons[ses]
    elif args.month == -6:
        ses = 'win'
        ses_tit = 'Winter'
        months = seasons[ses]
    else:
        months = [args.month]
        ses = 'm%02d' %args.month
        ses_tit = 'Month %02d' %args.month
    if args.year == 0:
        # years = [*range(2007,2020)]
        years = [*range(2007,2020)]
        ses = ses
        ses_tit = ses_tit
    else:
        years = [args.year]
        ses = 'y%d-%s' %(args.year, ses)
        ses_tit = 'Year %d, %s' %(args.year, ses_tit)
    # months = seasons[ses]
        
    ticT = time.time()
    y = -1
    for year in years:
        for mon in months:
            #: Check for missing months
            if year in missing_months.keys():
                if mon in missing_months[year]:
                    continue
            print('Year = %d, Month = %02d' %(year, mon))      
            y  = y + 1
            
            outname = 'CALIOP-EAD-%d%02d-n-DD' %(year, mon)
            tempname = '%s/part-%d%02d-%s.h5' %(tempDir, year, mon, 'global')
            if forced_age_bound is not None:
                tempname = tempname.replace('.h5', '_ab-%d.h5' %forced_age_bound)
            # partDir = os.path.join('/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/STC/Calipso', outname)
            partDir = os.path.join(mainDir, outname)
            initDir = os.path.join(partDir, 'Initfiles')
            if os.path.isfile(tempname) and lt:
                qs06_mon, qs10_mon, qs12_mon, qs13_mon, qs14_mon, qs16_mon, \
                        height_mon, lats_mon, lons_mon, extras = readTempFile(tempname)
                if args.extra:
                   
                    
                    f0 = os.path.join(partDir, 'part_000')
                    addExtra(tempname, f0)
                    continue
                    
                    
            else:
                #: Read part_000 file
                f0 = os.path.join(partDir, 'part_000')
                part0 = readidx107(f0, quiet=True)
                #: Read catalogfile
                catalogFile = os.path.join(initDir, 'selDardar_Calalog-%d%02d-n.pkl' %(year, mon))
                cf = getCatalogFile(catalogFile, checkForNewFile=False)
                #: Read out file
                rvs_mon, fls_mon, age_mon, lons_monF, lats_monF, temp_monF, pres_monF = getConvfiles(outDir, ['', outname])
                lons_monF = np.where(lons_monF > 180, lons_monF - 360, lons_monF)
                lons_monF = np.where(lons_monF < -180, lons_monF + 360, lons_monF)
                hits = (fls_mon & I_HIT) == I_HIT
                olds = (fls_mon & I_OLD) == I_OLD
                olds_hits = olds | hits
                
                #: What inds should we use?
                all_inds = np.ones(hits.shape).astype(bool)
                use_inds = hits
                #: Remove inds that has age after forced_age_bound
                if forced_age_bound is not None:
                    use_inds = np.where(age_mon > forced_age_bound, False, use_inds)
                #: Remove inds that starts after stop_day
                if stop_day is not None:
                    use_inds = np.where((np.abs(part0['ir_start']) > stop_day), False, use_inds)
                
                #: IWC is required
                use_inds = use_inds & (cf['iwc'] > 0) 
                olds_hits_cloudy = olds_hits & (cf['iwc'] > 0)
                olds_hits_clear = olds_hits & (~(cf['iwc'] > 0)) 
                #: Calc saturation
                # rvs0O = satratio(part0['p'][use_inds], part0['t'][use_inds])
                rvs0 = esati_murphy(part0['p'][use_inds], part0['t'][use_inds])
                iwc0 = convIWC(cf['iwc'], part0['p'], part0['t'], cf['sh'])
                iwc0_mon = iwc0[use_inds]
                #: Calc q0
                # vodinds = (cf['vod'][use_inds] < 0.03)
                # iwcinds = (cf['iwc'][use_inds] > 0)
                # pdb.set_trace()
                q0 = rvs0 + iwc0_mon
                q0_12 = rvs0 * 1.2
                q0_16 = rvs0 * 1.6
                #: get the desired variables and put in a 2D array
                # pMin = np.zeros([2, use_inds.sum()])
                # pMin[0, :] = part0['p'][use_inds]
                # rvsMin = np.zeros([2, use_inds.sum()])
                # rvsMin[0, :] = rvs0
                # q100 = np.zeros(use_inds.sum())
                # qs = np.zeros(use_inds.sum())
                # qsn = np.zeros(use_inds.sum())
                
                # if (break_point is not None) and (y==0):
                #     lon_lat_mon = np.zeros([10, int(break_point/3) + 1])
                #     lonlatind = np.zeros(10)
                qs06_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs08_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs10_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs12_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs13_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs14_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs16_16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                
                qs06_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs08_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs10_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs12_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs13_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs14_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs16_12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                
                qs06_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs08_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs10_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs12_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs13_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs14_mon = np.zeros([2, use_inds.sum()]) + np.nan
                qs16_mon = np.zeros([2, use_inds.sum()]) + np.nan
                
                
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
                    
                    qs06_mon, qs10_mon, qs12_mon, qs13_mon, qs14_mon, qs16_mon = \
                            findFirstQschange(q0, rvsf, qs06_mon, qs10_mon, qs12_mon, qs13_mon, qs14_mon, qs16_mon, \
                                              (part0['idx_back'] - part0['idx_orgn']), \
                                              (partf['idx_back'] - part0['idx_orgn']), \
                                              np.abs(partf['itime']-partf['ir_start']), \
                                              age_mon[use_inds], use_inds, phour)
                    qs06_12_mon, qs10_12_mon, qs12_12_mon, qs13_12_mon, qs14_12_mon, qs16_12_mon = \
                            findFirstQschange(q0_12, rvsf, qs06_12_mon, qs10_12_mon, qs12_12_mon, qs13_12_mon, qs14_12_mon, qs16_12_mon, \
                                              (part0['idx_back'] - part0['idx_orgn']), \
                                              (partf['idx_back'] - part0['idx_orgn']), \
                                              np.abs(partf['itime']-partf['ir_start']), \
                                              age_mon[use_inds], use_inds, phour)
                    qs06_16_mon, qs10_16_mon, qs12_16_mon, qs13_16_mon, qs14_16_mon, qs16_16_mon = \
                            findFirstQschange(q0_16, rvsf, qs06_16_mon, qs10_16_mon, qs12_16_mon, qs13_16_mon, qs14_16_mon, qs16_16_mon, \
                                              (part0['idx_back'] - part0['idx_orgn']), \
                                              (partf['idx_back'] - part0['idx_orgn']), \
                                              np.abs(partf['itime']-partf['ir_start']), \
                                              age_mon[use_inds], use_inds, phour)
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
                    if break_point is not None:
                        if phour > break_point:
                            break
                # f2r = f2.replace('.gz', '')
                partp0_p_mon = part0['p'][use_inds]
                partp0_t_mon = part0['t'][use_inds]
                height_mon = cf['height'][use_inds]
                lons_mon = lons_monF[use_inds]
                lats_mon = lats_monF[use_inds]
                pres_mon = pres_monF[use_inds]
                temp_mon = temp_monF[use_inds]
                if ct:
                    print('Create Temp File')
                    try:
                        h5file = h5py.File(tempname, 'w')
                        h5file.create_dataset('qs06', data = qs06_mon)
                        h5file.create_dataset('qs10', data = qs10_mon)
                        h5file.create_dataset('qs12', data = qs12_mon)
                        h5file.create_dataset('qs13', data = qs13_mon)
                        h5file.create_dataset('qs14', data = qs14_mon)
                        h5file.create_dataset('qs16', data = qs16_mon)
                        h5file.create_dataset('q0', data = q0)
                        
                        h5file.create_dataset('qs06_12', data = qs06_12_mon)
                        h5file.create_dataset('qs10_12', data = qs10_12_mon)
                        h5file.create_dataset('qs12_12', data = qs12_12_mon)
                        h5file.create_dataset('qs13_12', data = qs13_12_mon)
                        h5file.create_dataset('qs14_12', data = qs14_12_mon)
                        h5file.create_dataset('qs16_12', data = qs16_12_mon)
                        h5file.create_dataset('q0_12', data = q0_12)
                        
                        h5file.create_dataset('qs06_16', data = qs06_16_mon)
                        h5file.create_dataset('qs10_16', data = qs10_16_mon)
                        h5file.create_dataset('qs12_16', data = qs12_16_mon)
                        h5file.create_dataset('qs13_16', data = qs13_16_mon)
                        h5file.create_dataset('qs14_16', data = qs14_16_mon)
                        h5file.create_dataset('qs16_16', data = qs16_16_mon)
                        h5file.create_dataset('q0_16', data = q0_16)
                        
                        h5file.create_dataset('iwc0', data = iwc0_mon)
                        h5file.create_dataset('height', data = height_mon)
                        h5file.create_dataset('pressure', data = pres_mon)
                        h5file.create_dataset('temperature', data = temp_mon)
                        h5file.create_dataset('part0_p', data = partp0_p_mon)
                        h5file.create_dataset('part0_t', data = partp0_t_mon)
                        h5file.create_dataset('lons', data = lons_mon)
                        h5file.create_dataset('lats', data = lats_mon)
                        h5file.create_dataset('age', data = age_mon[use_inds])
                        h5file.create_dataset('ir_start', data = part0['ir_start'][use_inds])
                        h5file.create_dataset('use_inds', data = use_inds)
                        
                        h5file.create_dataset('ohClo_iwc0', data = iwc0[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_height', data = cf['height'][olds_hits_cloudy])
                        h5file.create_dataset('ohClo_pressure', data = pres_monF[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_temperature', data = temp_monF[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_part0_p', data = part0['p'][olds_hits_cloudy])
                        h5file.create_dataset('ohClo_part0_t', data = part0['t'][olds_hits_cloudy])
                        h5file.create_dataset('ohClo_lons', data = lons_monF[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_lats', data = lats_monF[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_age', data = age_mon[olds_hits_cloudy])
                        h5file.create_dataset('ohClo_ir_start', data = part0['ir_start'][olds_hits_cloudy])
                        h5file.create_dataset('ohClo_use_inds', data = olds_hits_cloudy)
                        
                        h5file.create_dataset('ohClr_iwc0', data = iwc0[olds_hits_clear])
                        h5file.create_dataset('ohClr_height', data = cf['height'][olds_hits_clear])
                        h5file.create_dataset('ohClr_pressure', data = pres_monF[olds_hits_clear])
                        h5file.create_dataset('ohClr_temperature', data = temp_monF[olds_hits_clear])
                        h5file.create_dataset('ohClr_part0_p', data = part0['p'][olds_hits_clear])
                        h5file.create_dataset('ohClr_part0_t', data = part0['t'][olds_hits_clear])
                        h5file.create_dataset('ohClr_lons', data = lons_monF[olds_hits_clear])
                        h5file.create_dataset('ohClr_lats', data = lats_monF[olds_hits_clear])
                        h5file.create_dataset('ohClr_age', data = age_mon[olds_hits_clear])
                        h5file.create_dataset('ohClr_ir_start', data = part0['ir_start'][olds_hits_clear])
                        h5file.create_dataset('ohClr_use_inds', data = olds_hits_clear)
                        
                        h5file.close(tempname)
                        print()
                    except:
                        print('Something wrong in H5')
                        sys.exit()
                        # pdb.set_trace()
            
            useSource = False
            if useSource == True:
                inds_ll_m = getAreaInds(area, lats_mon, lons_mon)
                inds_ll_clo = getAreaInds(area, extras['ohClo_lats'], extras['ohClo_lons'])
                inds_ll_clr = getAreaInds(area, extras['ohClr_lats'], extras['ohClr_lons'])
            else:
                inds_ll_m = getAreaInds(area, extras['lats0'], extras['lons0'])
                inds_ll_clo = getAreaInds(area, extras['ohClo_lats0'], extras['ohClo_lons0'])
                inds_ll_clr = getAreaInds(area, extras['ohClr_lats0'], extras['ohClr_lons0'])
                
                
            if y == 0:
                qs06 = qs06_mon[:, inds_ll_m]
                qs10 = qs10_mon[:, inds_ll_m]
                qs12 = qs12_mon[:, inds_ll_m]
                qs13 = qs13_mon[:, inds_ll_m]
                qs14 = qs14_mon[:, inds_ll_m]
                qs16 = qs16_mon[:, inds_ll_m]
                height = height_mon[inds_ll_m]
                lats = lats_mon[inds_ll_m]
                lons = lons_mon[inds_ll_m]
                age = extras['age'][inds_ll_m]
                ohClo_height = extras['ohClo_height'][inds_ll_clo]
                ohClr_height = extras['ohClr_height'][inds_ll_clr]
                pres_conv = extras['pressure'][inds_ll_m]
                temp_conv = extras['temperature'][inds_ll_m]
                pres_0 = extras['part0_p'][inds_ll_m]
                temp_0 = extras['part0_t'][inds_ll_m]
                lons_0 = extras['lons0'][inds_ll_m]
                lats_0 = extras['lats0'][inds_ll_m]
                ohClo_lons_0 = extras['ohClo_lons0'][inds_ll_clo]
                ohClo_lats_0 = extras['ohClo_lats0'][inds_ll_clo]
                # height_ohClo = height_ohClo
                conv_m, conv_tot_m, cf_m = calcProcHeight(qs06_mon[:, inds_ll_m], qs10_mon[:, inds_ll_m], qs12_mon[:, inds_ll_m], qs13_mon[:, inds_ll_m], qs14_mon[:, inds_ll_m], qs16_mon[:, inds_ll_m], extras['ohClo_height'][inds_ll_clo], extras['ohClr_height'][inds_ll_clr], height_mon[inds_ll_m])
                conv_m = conv_m.reshape(1, conv_m.shape[0], conv_m.shape[1])
                conv_tot_m = conv_tot_m.reshape(1, conv_tot_m.shape[0], conv_tot_m.shape[1])
            else:
                qs06 = np.hstack((qs06, qs06_mon[:, inds_ll_m]))
                qs10 = np.hstack((qs10, qs10_mon[:, inds_ll_m]))
                qs12 = np.hstack((qs12, qs12_mon[:, inds_ll_m]))
                qs13 = np.hstack((qs13, qs13_mon[:, inds_ll_m]))
                qs14 = np.hstack((qs14, qs14_mon[:, inds_ll_m]))
                qs16 = np.hstack((qs16, qs16_mon[:, inds_ll_m]))
                height = np.hstack((height, height_mon[inds_ll_m]))
                lats = np.hstack((lats, lats_mon[inds_ll_m]))
                lons = np.hstack((lons, lons_mon[inds_ll_m]))
                age = np.hstack((age, extras['age'][inds_ll_m]))
                ohClo_height = np.hstack((ohClo_height, extras['ohClo_height'][inds_ll_clo]))
                ohClr_height = np.hstack((ohClr_height, extras['ohClr_height'][inds_ll_clr]))
                pres_conv = np.hstack((pres_conv, extras['pressure'][inds_ll_m]))
                temp_conv = np.hstack((temp_conv, extras['temperature'][inds_ll_m]))
                pres_0 = np.hstack((pres_0, extras['part0_p'][inds_ll_m]))
                temp_0 = np.hstack((temp_0, extras['part0_t'][inds_ll_m]))
                lons_0 = np.hstack((lons_0, extras['lons0'][inds_ll_m]))
                lats_0 = np.hstack((lats_0, extras['lats0'][inds_ll_m]))
                ohClo_lons_0 = np.hstack((ohClo_lons_0, extras['ohClo_lons0'][inds_ll_clo]))
                ohClo_lats_0 = np.hstack((ohClo_lats_0, extras['ohClo_lats0'][inds_ll_clo]))
                
                conv_me, conv_tot_me, cf_me = calcProcHeight(qs06_mon[:, inds_ll_m], qs10_mon[:, inds_ll_m], qs12_mon[:, inds_ll_m], qs13_mon[:, inds_ll_m], qs14_mon[:, inds_ll_m], qs16_mon[:, inds_ll_m], extras['ohClo_height'][inds_ll_clo], extras['ohClr_height'][inds_ll_clr], height_mon[inds_ll_m])
                conv_me = conv_me.reshape(1, conv_me.shape[0], conv_me.shape[1])
                conv_m = np.concatenate((conv_m, conv_me), axis=0)
                conv_tot_me = conv_tot_me.reshape(1, conv_tot_me.shape[0], conv_tot_me.shape[1])
                conv_tot_m = np.concatenate((conv_tot_m, conv_tot_me), axis=0)
                cf_m = np.concatenate((cf_m, cf_me), axis=0)
                # conv_m[0] = np.vstack((conv_m[0], conv_me[0]))
                # conv_m[1] = np.vstack((conv_m[1], conv_me[1]))
                # conv_m[2] = np.vstack((conv_m[2], conv_me[2]))
                # conv_m[3] = np.vstack((conv_m[3], conv_me[3]))
                # conv_m[4] = np.vstack((conv_m[4], conv_me[4]))
                # conv_m[5] = np.vstack((conv_m[5], conv_me[5]))
                
    print('Tot calc time = %d' %(time.time() - ticT))
    if lt == False:
        sys.exit()
    if args.extra:
        print('Extra added')
        sys.exit()
    
    #: Calc
    if np.isnan(qs16[0,:]).any():
        print('Hmmm nan')
        pdb.set_trace()
        
    def calcPotT(p, t):
        #: Convert p (in Pa) to hPa
        hpa = (p/100.)
        #: Calc potential temperature using 1000hPa as base
        #: R / cp = 0.286
        #: t in K
        pt = t * (hpa/1000) **(0.286)
        return pt
    
    pt_conv = calcPotT(pres_conv, temp_conv)
    pt_0 = calcPotT(pres_0, temp_0)
    conv, conv_tot, cf, all_heights = calcProcHeight(qs06, qs10, qs12, qs13, qs14, qs16, ohClo_height, ohClr_height, height, rh=True)
    if conv_m.shape[0] > 1:
        conv_std = np.nanstd(conv_m, axis=0)
        conv_tot_std = np.nanstd(conv_tot_m, axis=0)
        cf_std = np.nanstd(cf_m, axis=0)
    else:
        conv_std = np.zeros(conv_m[0,:,:].shape)
        conv_tot_std = np.zeros(conv_tot_m[0,:,:].shape)
        cf_std = np.zeros(cf_m[0,:].shape)
    
    #: -- Plot --
    title_end = '%s, %s' %(ses_tit, area.title())
    figname_end = '%s_%s' %(ses, area.replace(' ', '_').replace('(', '').replace(')', ''))
    if vshift != 10:
        title_end = title_end + ', Vshift=%d' %vshift                  
    #: Plot height
    fig = plt.figure()
    fig.suptitle('%s' %(title_end))
    ax = fig.add_subplot(1,2,1)
    ax.plot(conv[0, :]*100, all_heights[0:-1], 'b', label='0.6')
    
    
    ax.fill_betweenx(all_heights[0:-1], conv[0, :]*100, conv[0, :]*100 - conv_std[0, :] * 100, facecolor='b', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[0, :]*100, conv[0, :]*100 + conv_std[0, :] * 100, facecolor='b', alpha=0.3)
    
    
    ax.plot(conv[1, :]*100, all_heights[0:-1], 'm', label='1.0')
    ax.fill_betweenx(all_heights[0:-1], conv[1, :]*100, conv[1, :]*100 - conv_std[1, :] * 100, facecolor='m', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[1, :]*100, conv[1, :]*100 + conv_std[1, :] * 100, facecolor='m', alpha=0.3)
    ax.plot(conv[2, :]*100, all_heights[0:-1], 'r', label='1.2')
    ax.fill_betweenx(all_heights[0:-1], conv[2, :]*100, conv[2, :]*100 - conv_std[2, :] * 100, facecolor='r', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[2, :]*100, conv[2, :]*100 + conv_std[2, :] * 100, facecolor='r', alpha=0.3)
    ax.plot(conv[3, :]*100, all_heights[0:-1], 'c', label='1.3')
    ax.fill_betweenx(all_heights[0:-1], conv[3, :]*100, conv[3, :]*100 - conv_std[3, :] * 100, facecolor='c', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[3, :]*100, conv[3, :]*100 + conv_std[3, :] * 100, facecolor='c', alpha=0.3)
    ax.plot(conv[4, :]*100, all_heights[0:-1], 'g', label='1.4')
    ax.fill_betweenx(all_heights[0:-1], conv[4, :]*100, conv[4, :]*100 - conv_std[4, :] * 100, facecolor='g', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[4, :]*100, conv[4, :]*100 + conv_std[4, :] * 100, facecolor='g', alpha=0.3)
    ax.plot(conv[5, :]*100, all_heights[0:-1], 'y', label='1.6')
    ax.fill_betweenx(all_heights[0:-1], conv[5, :]*100, conv[5, :]*100 - conv_std[5, :] * 100, facecolor='y', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], conv[5, :]*100, conv[5, :]*100 + conv_std[5, :] * 100, facecolor='y', alpha=0.3)
    
    ax.set_ylim(14000, 20000)
    ax.legend(fontsize='large')
    ax.set_xlabel('Convective [%]')
    ax.set_ylabel('Height [m]')
    yt = ax.get_yticks()
    # ax.set_title('Convective origion [%%], %s' %(title_end))
    
    ax = fig.add_subplot(1,2,2)
    ax.plot(cf[0, :]*100, all_heights[0:-1], 'b', label='CF')
    ax.fill_betweenx(all_heights[0:-1], cf[0, :]*100, cf[0, :]*100 - cf_std * 100, facecolor='b', alpha=0.3)
    ax.fill_betweenx(all_heights[0:-1], cf[0, :]*100, cf[0, :]*100 + cf_std * 100, facecolor='b', alpha=0.3)
    ax.set_ylim(14000, 20000)
    
    ytl = [''] * len(yt)
    ax.legend(fontsize='large')
    ax.set_xlabel('CF [%]')
    ax.set_yticklabels(ytl)
    
    figname = '%s/conv-div_height_%s' %(plotDir, figname_end)
    print(figname)
    # fig.savefig('test' + '.png')
    fig.savefig(figname + '.png')
    
    
    
    
    qss = [qs06, qs10, qs12, qs13, qs14, qs16]
    qssname = ['All pixel with hit', 'Threshold 0.6', 'Threshold 1.0', 'Threshold 1.2', 'Threshold 1.3', 'Threshold 1.4', 'Threshold 1.6']
    lenqss = len(qssname)
    fig = plt.figure(figsize=(8, 24))
    fig.suptitle('Histograms of time since convection, %s' %(title_end))
    for i in range(lenqss):
        f = i + 1
        qs_i = i - 1
        ax = fig.add_subplot(lenqss, 1, f)
        if i == 0:
            plot_age = age / (24*3600)
        else:
            qs = qss[qs_i]
            abowe_until_hit = (qs[0, :] == 2)
            plot_age = (age / (24*3600))[abowe_until_hit]
        #: Calculate Histogram
        #: 15 = range with 15 days steps
        #: 1 = range with 1 days steps
        h15 = np.histogram(plot_age, bins=range(0,91,15), density=True)
        #: Plot histogram
        #: * with np.diff(h15[1]) i.e. the width get the bars sum to be one
        # h = ax.hist(extras['age'] / (24*3600), bins=range(91), density=True)
        rects15 = ax.bar((np.diff(h15[1]) / 2) + h15[1][0:-1], h15[0] * np.diff(h15[1]), width=(14.8), color='r', alpha=0.5, label='15 days')
        h1 = np.histogram(plot_age, bins=range(91), density=True)
        rects1 = ax.bar((np.diff(h1[1]) / 2) + h1[1][0:-1], h1[0], color='b', label='1 day')
        ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlabel('Days from convection')
        ax.set_title(qssname[i])
        ax.legend()
        fig.tight_layout(rect=[0, 0, 1, 0.99])
    figname = '%s/day-since-conv_div_hist_%s' %(plotDir, figname_end)
    print(figname)
    fig.savefig(figname + '.png')
    if 'oem-Latitude-5400' in socket.gethostname():
        fig.show()
    
    fig = plt.figure(figsize=(24, 8))
    fig.suptitle('Potential Temperature, %s' %(title_end))
    for i in range(lenqss):
        f = i + 1
        qs_i = i - 1
        ax = fig.add_subplot(1, lenqss, f)
        if i == 0:
            abowe_until_hit = np.ones(qs16.shape[1]).astype(bool)
        else:
            qs = qss[qs_i]
            abowe_until_hit = (qs[0, :] == 2)
        hh, xedges, yedges = np.histogram2d(pt_0[abowe_until_hit], pt_conv[abowe_until_hit], bins=[[*range(85, 134)],  [*range(83, 194)]], density=True)
        # print(hh.max())
        # xedges 85.85816108, 132.81934959
        # yedges 83.11706566, 192.31496401
        im = ax.imshow(hh.T, origin ='lower', vmin=0, vmax=0.003)#, norm=colors.Normalize())#, vmin=0, vmax=6000)
        ax.set_xticks([5, 15, 25, 35, 45])
        ax.set_xticklabels([xedges[5], xedges[15], xedges[25], xedges[35], xedges[45]])
        ax.set_yticks([2, 17, 32, 47, 62, 77, 92, 107])
        ax.set_yticklabels([yedges[2], yedges[17], yedges[32], yedges[47], yedges[62], yedges[77], yedges[92], yedges[107]])
        ax.set_title(qssname[i])
        ax.set_xlabel('Satellite')

        # fig.subplots_adjust(right=0.89)
        # cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
        # cbar = fig.colorbar(im1, cbar_ax)
        
        if f == 1:
            ax.set_ylabel('Convective Source')
        if f == lenqss:
            axpos = ax.get_position()
            pos_x = axpos.x0 + axpos.width + 0.01
            # cbar_ax = fig.add_axes([pos_x, axpos.y0, 0.01, axpos.height])
            cbar_ax = fig.add_axes([pos_x, 0.12, 0.01, 0.71])
            cbar = fig.colorbar(im, orientation='vertical', cax=cbar_ax)#, ticks=barticks)  # @UnusedVariable
    
    # ax1 = fig.add_subplot(2,1,2)
    # ax1.hist2d(pt_0, pt_conv, bins=[[*range(85, 134)],  [*range(83, 194)]])
    fig.tight_layout(rect=[0, 0, 0.91, 0.99])
    figname = '%s/pt_div_2D-hist_%s' %(plotDir, figname_end)
    print(figname)
    fig.savefig(figname + '.png')
    
    figname_org = '%s/maps_2D-hist_%s' %(plotDir, figname_end)
    title_org = 'Convective, %s' %(title_end)
    #: Height
    heightBoundaries = [[None, None]]#, [14,15], [15,16], [14,16], [16, 18], [18, 20], [18,19], [19,20]] 
    for h1, h2 in heightBoundaries:
        if not ((h1 is None) and (h2 is None)):
            title = title_org + ', %d - %d km' %(h1, h2)
            figname = figname_org + '_%d-%d' %(h1, h2)
            if h2 == heightBoundaries[-1][1]:
                inds_h = (height >= h1) & (height <= h2)
            else:
                inds_h = (height >= h1) & (height < h2)
            inds = inds_h
        else:
            title = title_org
            figname = figname_org
            inds = np.ones(height.shape[0]).astype(bool)
        if inds.sum() == 0:
            continue
        #: 2D Histogram
        hh1, xedges1, yedges1 = np.histogram2d(lons[inds], lats[inds], bins=[180, 90])
        hh2, xedges2, yedges2 = np.histogram2d(lons_0[inds], lats_0[inds], bins=[180, 90])
        #: Latitude boundary
        ymax = getYmax(yedges1, yedges2)
        aspect = float('%.2f' %((1/3) / (ymax / 180.)))
        fig = plt.figure()
        fig.suptitle(title)
        #: Subplot 1
        ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
        ax.coastlines()
        im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]])
        ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
        ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.set_title('Convective source')
        fig.subplots_adjust(right=0.89)
        cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
        cbar = fig.colorbar(im1, cbar_ax)
        #: Subplot 2
        ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
        ax.coastlines()
        im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]])
        ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
        ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.set_title('Cloud')
        fig.subplots_adjust(right=0.89)
        cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
        cbar = fig.colorbar(im2, cax=cbar_ax)
        print(figname)
        fig.savefig(figname + '.png')
        plt.close(fig)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    sys.exit()
    
    #: Add counts above the bars
    
    heightBoundaries = [[None, None], [14,15], [15,16], [14,16], [16, 18], [18, 20], [18,19], [19,20]]
    
    title_org = 'Histogram'
    inds_org = (~np.isnan(qs16[0,:])).astype(bool)
    figname_1_org = '%s/qsn-div_total' %plotDir
    figname_2_org = '%s/qsn-total' %plotDir
    for h1, h2 in heightBoundaries:
        if not ((h1 is None) and (h2 is None)):
            # title = title_org + ', %d - %d km' %(h1, h2)
            title = '%d - %d km' %(h1, h2)
            figname_1 = figname_1_org + '_%d-%d' %(h1, h2)
            figname_2 = figname_2_org + '_%d-%d' %(h1, h2)
            if h2 == heightBoundaries[-1][1]:
                inds_h = (height >= h1) & (height <= h2)
            else:
                inds_h = (height >= h1) & (height < h2)
            inds = inds_org & inds_h
        else:
            title = ' Above 14 km'
            # title = title_org
            figname_1 = figname_1_org
            figname_2 = figname_2_org
            inds = inds_org
        
        # title = '%s, %s, %s' %(title, ses.title(), area.title())
        figname_1 = '%s_%s_%s' %(figname_1, ses, area.replace(' ', '_'))
        figname_2 = '%s_%s_%s' %(figname_2, ses, area.replace(' ', '_'))
        
        #: Area
        if area != 'global':
            latmin = areas[area]['minLat']
            latmax = areas[area]['maxLat']
            lonmin = areas[area]['minLon']
            lonmax = areas[area]['maxLon']
            inds_lat = (lats > latmin) & (lats <= latmax)
            inds_lon = (lons > lonmin) & (lons <= lonmax)
            inds = inds & inds_lat & inds_lon
        
        if inds.sum() == 0:
            continue
        abowe_until_hit_16 = (qs16[0, inds]==2).sum()
        abowe_until_hit_14 = (qs14[0, inds]==2).sum()
        abowe_until_hit_13 = (qs13[0, inds]==2).sum()
        abowe_until_hit_12 = (qs12[0, inds]==2).sum()
        abowe_until_hit_10 = (qs10[0, inds]==2).sum()
        abowe_until_hit_06 = (qs06[0, inds]==2).sum()
        
        
        abowe_no_hit_16 = (qs16[0, inds]==1).sum()
        abowe_no_hit_14 = (qs14[0, inds]==1).sum()
        abowe_no_hit_13 = (qs13[0, inds]==1).sum()
        abowe_no_hit_12 = (qs12[0, inds]==1).sum()
        abowe_no_hit_10 = (qs10[0, inds]==1).sum()
        abowe_no_hit_06 = (qs06[0, inds]==1).sum()
        
        below_before_hit_16 = (qs16[0, inds]==-1).sum()
        below_before_hit_14 = (qs14[0, inds]==-1).sum()
        below_before_hit_13 = (qs13[0, inds]==-1).sum()
        below_before_hit_12 = (qs12[0, inds]==-1).sum()
        below_before_hit_10 = (qs10[0, inds]==-1).sum()
        below_before_hit_06 = (qs06[0, inds]==-1).sum()
        
        mean_time_16 = np.mean(qs16[1, inds][qs16[0, inds]==-1])
        mean_time_14 = np.mean(qs14[1, inds][qs14[0, inds]==-1])
        mean_time_13 = np.mean(qs13[1, inds][qs13[0, inds]==-1])
        mean_time_12 = np.mean(qs12[1, inds][qs12[0, inds]==-1])
        mean_time_10 = np.mean(qs10[1, inds][qs10[0, inds]==-1])
        mean_time_06 = np.mean(qs06[1, inds][qs06[0, inds]==-1])
        
        total_16 = abowe_until_hit_16 + abowe_no_hit_16 + below_before_hit_16
        total_14 = abowe_until_hit_14 + abowe_no_hit_14 + below_before_hit_14
        total_13 = abowe_until_hit_13 + abowe_no_hit_13 + below_before_hit_13
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
    
        labels = ['0.6', '1.0', '1.2', '1.3', '1.4', '1.6']
        x = np.arange(len(labels))  # the label locations
        width = 0.3  # the width of the bars

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # rects1 = ax.bar(x - width/2, [(qs06[0,:]==-1).sum() / (~(np.isnan(qs06[0,:]))).sum(), (qs10[0,:]==-1).sum() / (~(np.isnan(qs10[0,:]))).sum(), (qs12[0,:]==-1).sum() / (~(np.isnan(qs12[0,:]))).sum(), (qs14[0,:]==-1).sum() / (~(np.isnan(qs14[0,:]))).sum(), (qs16[0,:]==-1).sum() / (~(np.isnan(qs16[0,:]))).sum()], width, label='Below')
        # rects2 = ax.bar(x + width/2, [ind06.sum() / (~(np.isnan(qs06[0,:]))).sum(), ind10.sum() / (~(np.isnan(qs10[0,:]))).sum(), ind12.sum() / (~(np.isnan(qs12[0,:]))).sum(), ind14.sum() / (~(np.isnan(qs14[0,:]))).sum(), ind16.sum() / (~(np.isnan(qs16[0,:]))).sum()], width, label='Convective')
        # rects1 = ax.bar(x - width/2, [(~ind06).sum() / (qs06[0,:]).shape[0], (~ind10).sum() / (qs10[0,:]).shape[0], (~ind12).sum() / (qs12[0,:]).shape[0], (~ind14).sum() / (qs14[0,:]).shape[0], (~ind16).sum() / (qs16[0,:]).shape[0]], width, label='Incitu')
        # rects2 = ax.bar(x + width/2, [ind06.sum() / (qs06[0,:]).shape[0], ind10.sum() / (qs10[0,:]).shape[0], ind12.sum() / (qs12[0,:]).shape[0], ind14.sum() / (qs14[0,:]).shape[0], ind16.sum() / (qs16[0,:]).shape[0]], width, label='Convective')
        # rects3 = ax.bar(x + width/2, [np.mean(qs10[1, ind10]), np.mean(qs12[1, ind12]), np.mean(qs14[1, ind14]), np.mean(qs16[1, ind16])], width, label='Mean Time [h]')
        
        if abowe_no_hit_06 != 0:
            rects1 = ax.bar(x - width/2, [abowe_until_hit_06 / total_06, abowe_until_hit_10 / total_10, abowe_until_hit_12 / total_12, abowe_until_hit_13 / total_13, abowe_until_hit_14 / total_14, abowe_until_hit_16 / total_16], width/2, color="b", label='Convective')
            rects2 = ax.bar(x, [abowe_no_hit_06 / total_06, abowe_no_hit_10 / total_10, abowe_no_hit_12 / total_12, abowe_no_hit_13 / total_13, abowe_no_hit_14 / total_14, abowe_no_hit_16 / total_16], width/2, color="g", label='Convective, no hit')
            rects3 = ax.bar(x + width/2, [below_before_hit_06 / total_06, below_before_hit_10 / total_10, below_before_hit_12 / total_12, below_before_hit_13 / total_13, below_before_hit_14 / total_14, below_before_hit_16 / total_16], width/2, color="r", label='In situ')
        else:
            rects1 = ax.bar(x - width/2, [abowe_until_hit_06 / total_06, abowe_until_hit_10 / total_10, abowe_until_hit_12 / total_12, abowe_until_hit_13 / total_13, abowe_until_hit_14 / total_14, abowe_until_hit_16 / total_16], width, color="b", label='Convective')
            rects3 = ax.bar(x + width/2, [below_before_hit_06 / total_06, below_before_hit_10 / total_10, below_before_hit_12 / total_12, below_before_hit_13 / total_13, below_before_hit_14 / total_14, below_before_hit_16 / total_16], width, color="r", label='In situ')
        #-1.5 0 1.5
        
        
        #: Add counts above the two bar graphs
        for rect in rects1:
            hs = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2.0, hs, '%.2f' %hs, ha='center', va='bottom')#f'{height:.0f}', ha='center', va='bottom')
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize='large')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize='large')
        ax.set_xlabel('Thresholds', fontsize='x-large')
        # ax.set_ylabel('%', fontsize='x-large')
        ax.set_title(title, fontsize='x-large')
        ax.legend(fontsize='large')
        # ax.set_title('Num of change')
        
        plt.tight_layout()
        print(figname_1)
        fig.savefig(figname_1 + '.png')
        if 'oem-Latitude-5400' in socket.gethostname():
            fig.show()
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if abowe_no_hit_06 != 0:
            rects1 = ax.bar(x - width/2, [abowe_until_hit_06, abowe_until_hit_10, abowe_until_hit_12, abowe_until_hit_13, abowe_until_hit_14, abowe_until_hit_16], width/2, color="b", label='Convective')
            rects2 = ax.bar(x, [abowe_no_hit_06, abowe_no_hit_10, abowe_no_hit_12, abowe_no_hit_13, abowe_no_hit_14, abowe_no_hit_16], width/2, color="g", label='Convective, no hit')
            rects3 = ax.bar(x + width/2, [below_before_hit_06, below_before_hit_10, below_before_hit_12, below_before_hit_13, below_before_hit_14, below_before_hit_16], width/2, color="r", label='In situ')
        else:
            rects1 = ax.bar(x - width/2, [abowe_until_hit_06, abowe_until_hit_10, abowe_until_hit_12, abowe_until_hit_13, abowe_until_hit_14, abowe_until_hit_16], width, color="b", label='Convective')
            rects3 = ax.bar(x + width/2, [below_before_hit_06, below_before_hit_10, below_before_hit_12, below_before_hit_13, below_before_hit_14, below_before_hit_16], width, color="r", label='In situ')
        
        ax.text
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize='large')
        ax.set_xlabel('Thresholds', fontsize='large')
        ax.set_ylabel('Total number', fontsize='large')
        ax.set_title(title, fontsize='large')
        ax.legend()
        # ax.set_title('Num of change')
        
        plt.tight_layout()
        print(figname_2)
        fig.savefig(figname_2 + '.png')
        if 'oem-Latitude-5400' in socket.gethostname():
            fig.show()
        
    
    
    
    
    
    # plt.close(fig)
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
