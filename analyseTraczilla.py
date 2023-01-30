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
import flammkuchen as fk  # @UnresolvedImport

# sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readidx107
from convsrcErikFullGridSatTropo import I_HIT, I_OLD, I_DEAD, I_CROSSED, I_DBORNE
import calendar
# I_DEAD = 0x200000 #: Dead
# I_HIT = 0x400000 #: Hit a cloud
# I_CROSSED = 0x2000000 #: outside domain
# I_DBORNE =  0x1000000 #: Lost after first step
# I_OLD = 0x800000 #: Reached end of time without encounter a cloud
I_STOP = I_HIT + I_DEAD

# areas = {'asia': {'minLat': -30, 'maxLat': 20, 'minLon': 80, 'maxLon': 150}, \
#          'asian monsoon': {'minLat': -30, 'maxLat': 30, 'minLon': 150, 'maxLon': 180}, \
#          'anticyclone (AMA)': {'minLat': 20, 'maxLat': 30, 'minLon': -50, 'maxLon': 150}, \
#          'pacific': {'minLat': -30, 'maxLat': 30, 'minLon': -180, 'maxLon': -115}, \
#          'central america': {'minLat': -30, 'maxLat': 30, 'minLon': -115, 'maxLon': -40}, \
#          'atlantic': {'minLat': -30, 'maxLat': 30, 'minLon': 150, 'maxLon': 180}, \
#          'africa': {'minLat': -30, 'maxLat': 20, 'minLon': 150, 'maxLon': 180}}
areas = {'asia': {'num': 1,'minLat1': -30, 'maxLat1': 20, 'minLon1': 80, 'maxLon1': 150, \
                           'minLat2': -30, 'maxLat2': 30, 'minLon2': 150, 'maxLon2': 180}, \
         'asian monsoon anticyclone (AMA)': {'num': 2,'minLat': 20, 'maxLat': 30, 'minLon': 50, 'maxLon': 150}, \
         'pacific': {'num': 3,'minLat': -30, 'maxLat': 30, 'minLon': -180, 'maxLon': -115}, \
         'central america': {'num': 4,'minLat': -30, 'maxLat': 30, 'minLon': -115, 'maxLon': -40}, \
         'africa': {'num': 5, 'minLat1': -30, 'maxLat1': 30, 'minLon1': -40, 'maxLon1': 50, \
                              'minLat2': -30, 'maxLat2': 20, 'minLon2': 50, 'maxLon2': 80}, \
         'asian monsoon': {'num': 6, 'minLat': -20, 'maxLat': 40, 'minLon': 65, 'maxLon': 95}, \
         'pacific 2': {'num': 7, 'minLat': -30, 'maxLat': 30, 'minLon': -180, 'maxLon': -115}, \
         'central america 2': {'num': 8, 'minLat': -30, 'maxLat': 30, 'minLon': -115, 'maxLon': -40}, \
         'atlantic': {'num': 9, 'minLat': -30, 'maxLat': 30, 'minLon': -40, 'maxLon': 0}, \
         'africa 2': {'num': 10, 'minLat': -30, 'maxLat': 30, 'minLon': 0, 'maxLon': 45}, \
         'warm pool': {'num': 11, 'minLat': -30, 'maxLat': 30, 'minLon': 70, 'maxLon': 160}, \
         'nino3':{'num': 12, 'minLat': -5, 'maxLat': 5, 'minLon': -150, 'maxLon': -90}, \
         'nino4':{'num': 13, 'minLat1': -5, 'maxLat1': 5, 'minLon1': 150, 'maxLon1': 180, \
                             'minLat2': -5, 'maxLat2': 5, 'minLon2': -180, 'maxLon2': -150}, \
         'nino34':{'num': 14, 'minLat1': -5, 'maxLat1': 5, 'minLon1': 150, 'maxLon1': 180, \
                              'minLat2': -5, 'maxLat2': 5, 'minLon2': -180, 'maxLon2': -90}}

        
seasons = {'year': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'djf': [12, 1, 2], 'mam': [3, 4, 5], \
           'jja': [6, 7, 8], 'son': [9, 10, 11], 'sum': [3, 4, 5, 6, 7, 8], 'win': [9, 10, 11, 12, 1, 2]}

missing_months = {2007: [1,2], 2011: [1, 5, 6, 7, 8, 9, 10, 11, 12], 2012: [1, 2, 3, 4], 2016: [2], 2017: [7], 2018: [1, 2, 3, 4, 5], 2019: [7, 8, 9, 10, 11, 12]}
magic = 0.84


def getAreaInds(arn, lat, lon):
    if arn.lower() in ['global']:
        llInd = np.ones(lat.shape).astype(bool)
    elif arn in areas.keys():
        if 'minLat1' in areas[arn].keys():
            latmin1 = areas[arn]['minLat1']
            latmax1 = areas[arn]['maxLat1']
            lonmin1 = areas[arn]['minLon1']
            lonmax1 = areas[arn]['maxLon1']
            latmin2 = areas[arn]['minLat2']
            latmax2 = areas[arn]['maxLat2']
            lonmin2 = areas[arn]['minLon2']
            lonmax2 = areas[arn]['maxLon2']
            indsLat1 = (lat > latmin1) & (lat <= latmax1)
            indsLon1 = (lon > lonmin1) & (lon <= lonmax1)
            llInd1 = indsLat1 & indsLon1
            indsLat2 = (lat > latmin2) & (lat <= latmax2)
            indsLon2 = (lon > lonmin2) & (lon <= lonmax2)
            llInd2 = indsLat2 & indsLon2
            llInd = llInd1 | llInd2
        else:
            latmin = areas[arn]['minLat']
            latmax = areas[arn]['maxLat']
            lonmin = areas[arn]['minLon']
            lonmax = areas[arn]['maxLon']
            indsLat = (lat > latmin) & (lat <= latmax)
            indsLon = (lon > lonmin) & (lon <= lonmax)
            llInd = indsLat & indsLon
    else:
        print('Wrong Area')
        sys.exit()
    return llInd


def getAreaName(an):
    ans = []
    if an == 0:
        area = 'global'
    else:
        for arrname, arrdata in areas.items():
            ans.append(arrdata['num'])
            if an == arrdata['num']:
                area = arrname
                break
    try:
        return area
    except:
        print('wrong area num')
        print('Area num avail')
        print(ans)
        print('Area num used = %d' %an)
        sys.exit()


def convIWC(wc, p, T, qv):
    #: Convert WC from kg m-3 to kg /kg and to ppmv
    Md = 28.97
    Rv = 461
    Rd = 287.04
    # #: https://www.gastec.co.jp/en/technology/knowledge/concentration/
    # ppm1 = (wc) * (22.4 / Md) * (T / 273) * (101300 / p)
    # #: https://www.cdc.gov/niosh/docs/2004-101/calc.html
    # ro2 = p / (Rd * T)
    # wc_kgkg2 = wc * 1/ro2
    # # wc_mgkg = wc_kgkg * 1000000
    # ppm2 = wc_kgkg2 * 22.41 / Md
    
    Tv = T * (1 + (Rv/Rd - 1) * qv)
    ro = p / (Rd * Tv)
    wc_kgkg = wc * 1/ro
    ppm = wc_kgkg# * 1 / 0.622
    return ppm
    
def esati_murphy(p, T):
    #ei in Pa saturation vapor pressure with respect to hexagonal (most stable) ice
    #Murphy and Koop 2005, QJRMS
    lnP=9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T
    
    esati_murphy=np.exp(lnP) / p
    #: Convert to kg/kg
    esati_murphy = esati_murphy * 0.622
    #ppmv
    #kg
    
    return esati_murphy






#:----------------------------------------------------


def getYmax(y1, y2):
    ym = np.max((y1, y2))
    if ym > 60:
        ym = 90
    else:
        ym = 60
    return ym




#:----------------------------------------------------
def compareFiles(f1, f1n, f2, f2n, ft=0):
    #: 0 => Both files Catalog
    #: 1 => f1 = Catalog, f2 = parad
    f1l = 0
    f2l = 0
    f12l = 0
    for d in f1.keys():
        for o in f1[d].keys():
            for an in f1[d][o].keys():
                
                # try:
                #     isinstance(f1[d][o][an][0], (datetime.datetime))
                # except:
                #     pdb.set_trace()
                if isinstance(f1[d][o][an], (int)):
                    if not (f1[d][o][an] == f2[d][o][an]):
                        print('not same int content')
                        pdb.set_trace()
                    continue
                if isinstance(f1[d][o][an][0], (datetime.datetime)):
                    f1isnan = np.zeros(np.asarray(f1[d][o][an]).shape).astype(bool)
                    f2isnan = np.zeros(np.asarray(f2[d][o][an]).shape).astype(bool)
                else:
                    f1isnan = np.isnan(f1[d][o][an])
                    f2isnan = np.isnan(f2[d][o][an])
                if not f1isnan.shape == f2isnan.shape:
                    print('wrong nan shape')
                    pdb.set_trace()
                if not (np.asarray(f1[d][o][an])[~f1isnan] == np.asarray(f2[d][o][an])[~f2isnan]).all():
                    print('not same content')
                    pdb.set_trace()
    
            if len(f1[d][o].keys()) == len(f2[d][o].keys()):
                f12l = f12l + 1
            elif len(f1[d][o].keys()) > len(f2[d][o].keys()):
                f1l = f1l + 1
            elif len(f1[d][o].keys()) < len(f2[d][o].keys()):
                f2l = f2l + 1
    
                
                
    #: TODO: Kolla
    # datetime.datetime(2008, 10, 1, 0, 0)
    print("")
    print("All data in f1 exist (and are the same) in f2 ")
    if (f12l > 0) and (f1l == 0) and (f2l == 0):
        print('f1 and f2 contains same amount of data')
    elif (f12l == 0) and (f1l > 0) and (f2l == 0):
        print('f1 contains more data')
    elif (f12l == 0) and (f1l == 0) and (f2l > 0):    
        print('f2 contains more data')
    else:
        print('f1 and f2 contains defferent length of data')
        print('f12l = %d' %f12l)
        print('f1l = %d' %f1l)
        print('f2l = %d' %f2l)
        pdb.set_trace()
    
    print("f1 was created: %s" %time.ctime(os.path.getctime(f1n)))
    print("f2 was created: %s" %time.ctime(os.path.getctime(f2n)))
    if os.path.getctime(f1n) == os.path.getctime(f2n):
        print("Files created at the same time")
    if os.path.getctime(f1n) > os.path.getctime(f2n):
        print("f1 is newer")
    if os.path.getctime(f1n) < os.path.getctime(f2n):
        print("f2 is newer")
    print("Delete one file")
    print("")
    print("f1 = %s" %f1n)
    print("f2 = %s" %f2n)
    sys.exit()


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


def readCatalogFile(fn):
    retv = []
    with (gzip.open(fn, "rb")) as openfile:
        while True:
            try:
                retv.append(pickle.load(openfile))
            except EOFError:
                break
    return retv[0]


def readConvFile(fn, usefk=False):
    """ 
    Read the resultfile from convsrcErikFullGridSatTropo
    Deafault is to use h5py but flamkuchen is an option, set usefk = True
    """
    if usefk:
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


def getCatalogFile(fn, p0=None, checkForNewFile=True):
    
    objects = readCatalogFile(fn)
    if checkForNewFile:
        tf = '/data/ejohansson/flexout/STC/Calipso/CALIOP-EAD/Catalog/%s' %os.path.basename(fn)
        if os.path.isfile(tf):
            to = readCatalogFile(tf)
            compareFiles(objects, fn, to, tf, 0)
    
    nameDic = {'CM': 'Cloud_Mask', 'TpH': 'Tropopause_Height', \
                'SZA': 'SZA', 'SC': 'Simplified_Categorization', \
                'lensel': 'lensel', 'npart': 'npart', 'vod': 'vis_optical_depth', \
                'height': 'Track_Height', 'iwc': 'iwc', 'sh': 'Specific_humidity', 'utc': 'utc'}
    retv = {}
    for rn in nameDic.keys():
        retv.update({rn: np.array([])})
    sp = 0
    for da in objects.keys():
        for orbit in objects[da].keys():
            ep = sp + objects[da][orbit]['npart']
            if p0 is not None:
                #: Test to make sure the param korrespond to the part0
                if not objects[da][orbit]['longitudes'][0] == p0['x'][sp]:
                    print('Lon 0 test')
                    pdb.set_trace()
                #: Not sure why -1
                if not objects[da][orbit]['longitudes'][1] == p0['x'][ep-1]:
                    print('Lon 1 test')
                    pdb.set_trace()
            
            # (np.reshape(pres, npart).reshape(npart//nlev, nlev) == pres).all()
            # (np.reshape(pres, npart).reshape(len(sel), nlev) == pres).all()
            # npart = nlev * len(sel)
            # nlev = objects[da][orbit]['npart'] // objects[da][orbit]['lensel']
            # npart = objects[da][orbit]['npart']
            lensel = objects[da][orbit]['lensel']
            for retvn in retv.keys():
                if retvn in ['lensel', 'npart']:
                    retv[retvn] = np.concatenate((retv[retvn], np.repeat(objects[da][orbit][nameDic[retvn]], (ep-sp))))
                elif retvn in ['height']:
                    #: Create 2D
                    # np.repeat(altx.reshape((1, -1)), len(sel), axis=0)
                    #: reshape to 1D
                    # np.reshape(np.repeat(altx.reshape((1, -1)), len(sel), axis=0), npart)
                    #: Both in 1 line
                    # np.reshape(np.repeat(objects[da][orbit][nameDic[retvn]].reshape((1, -1)), lensel, axis=0), npart)
                    #: This is the same
                    htemp = np.tile(objects[da][orbit][nameDic[retvn]], lensel)
                    retv[retvn] = np.concatenate((retv[retvn], htemp))
                else:
                    retv[retvn] = np.concatenate((retv[retvn], objects[da][orbit][nameDic[retvn]]))
            sp = ep
    
    
    # if p0 is not None:
    #     retv = {'CM': np.zeros(p0['numpart']), 'TpH': np.zeros(p0['numpart']), \
    #             'SZA': np.zeros(p0['numpart']), 'SC': np.zeros(p0['numpart']), \
    #             'lensel': np.zeros(p0['numpart']), 'npart': np.zeros(p0['numpart'])}
    # objects
    # sp = 0
    # for da in objects.keys():
    #     for orbit in objects[da].keys():
    #         ep = sp + objects[da][orbit]['npart']
    #         retv['CM'][sp:ep] = objects[da][orbit]['Cloud_Mask']
    #         retv['TpH'][sp:ep] = objects[da][orbit]['Tropopause_Height']
    #         retv['SZA'][sp:ep] = objects[da][orbit]['SZA']
    #         retv['SC'][sp:ep] = objects[da][orbit]['Simplified_Categorization']
    #         retv['lensel'][sp:ep] = np.repeat(objects[da][orbit]['lensel'], (ep-sp))
    #         retv['npart'][sp:ep] = np.repeat(objects[da][orbit]['npart'], (ep-sp))
    #         #: Test to make sure the param korrespond to the part0
    #         if not objects[da][orbit]['longitudes'][0] == p0['x'][sp]:
    #             print('Lon 0 test')
    #             pdb.set_trace()
    #         #: Not sure why -1
    #         if not objects[da][orbit]['longitudes'][1] == p0['x'][ep-1]:
    #             print('Lon 1 test')
    #             pdb.set_trace()
    #         # pdb.set_trace()
    #         sp = ep
    # else:
    #     retv = objects
    
    return retv


# def getInitfiles(mD, years, se, dn, uD, lt=True):
#     if (not isinstance(years, list)):# or (not isinstance(months, list)):
#         print('Wrong type')
#         pdb.set_trace()
#     j = -1
#     if isinstance(se, int):
#         months = [se]
#     else:
#         months = seasons[se]
#
#     for y in years:
#         j = j + 1
#         #: Decide the name for tempfile
#         if not isinstance(se, int): #(len(months) == 12):
#             tempname = os.path.join(mD,'TempFiles', 'CALIOP-EAD-%d-%s' %(y, dn))
#             if uD:
#                 tempname = tempname + '-DD'
#             # if not ((se=='year') and (year == 2009)):
#             tempname = tempname + '-%s' %se
#             tempname = tempname + '-init'
#             loadname = tempname + '-param.h5'
#             # loadname = tempname + '.npy'
#         else:
#             loadname = ''
#             tempname = ''
#         #: If tempfile exist and lt = True, load the tempfile
#         if (lt and os.path.isfile(loadname)):
#             # retvPL, retvCL, retvML, outnamesTL = np.load(loadname, allow_pickle=True)
#             retvPL = fk.load(tempname + '-part.h5')
#             retvCL = fk.load(tempname + '-catalog.h5')
#             retvML = fk.load(tempname + '-param.h5')
#             outnamesTL = retvPL.pop('on')
#
#
#         #: if tempfile do not exist or lt = False. Do calculations
#         else:
#             i = -1
#             outnamesTL = []
#             for m in months:
#                 if y in missing_months.keys():
#                     if m in missing_months[y]:
#                         continue
#                 i = i + 1
#                 #: Directories of the backward trajectories and name of the output file 77193056,
#                 outname = 'CALIOP-EAD-%d%02d-%s' %(y, m, dn)
#                 if uD:
#                     outname = outname + '-DD'
#                 outnamesTL.append(outname)
#                 trajDir = os.path.join(mD,outname) 
#                 initDir = os.path.join(trajDir, 'Initfiles')
#
#                 paramname = 'selDardar_Params-%s.pkl' %'-'.join(outname.split('-')[2:4])
#                 catalogname = paramname.replace('_Params-', '_Calalog-')
#                 paramFile = os.path.join(initDir, paramname)
#                 catalFile = os.path.join(initDir, catalogname)
#                 pf = readidx107(os.path.join(trajDir,'part_000'),quiet=False)
#                 cf = getCatalogFile(catalFile, pf, checkForNewFile=False)
#                 mf = readParamFile(paramFile)
#                 pf['x'] = np.where(pf['x'] > 180, pf['x'] - 360, pf['x'])
#                 if i == 0:
#                     retvPL = pf.copy()
#                     retvCL = cf.copy()
#                     retvML = mf.copy()
#                 else:
#                     for arname in retvPL.keys():
#                         if isinstance(retvPL[arname], int):
#                             retvPL[arname] = [retvPL[arname]]
#                         if isinstance(pf[arname], int):
#                             pf[arname] = [pf[arname]]
#                         retvPL[arname] = np.concatenate((np.asarray(retvPL[arname]), np.asarray(pf[arname])))
#
#                     for arname in retvCL.keys():
#                         if isinstance(retvCL[arname], int):
#                             retvCL[arname] = [retvCL[arname]]
#                         if isinstance(cf[arname], int):
#                             cf[arname] = [cf[arname]]
#                         retvCL[arname] = np.concatenate((np.asarray(retvCL[arname]), np.asarray(cf[arname])))
#
#                     for arname in retvML.keys():
#                         if isinstance(retvML[arname], (int, np.int64, str, datetime.datetime)):
#                             retvML[arname] = [retvML[arname]]
#                         if isinstance(mf[arname], (int, np.int64, str, datetime.datetime)):
#                             mf[arname] = [mf[arname]]
#                         retvML[arname] = np.concatenate((np.asarray(retvML[arname]), np.asarray(mf[arname])))
#             #: If there are no tempfile, create one
#             if (loadname != '') and (not os.path.isfile(loadname)):
#                 retvPL.update({'on': outnamesTL})
#                 fk.save(tempname + '-part.h5', retvPL)
#                 fk.save(tempname + '-catalog.h5', retvCL)
#                 fk.save(tempname + '-param.h5', retvML)
#                 dummy = retvPL.pop('on')
#                 dummy = None
#                 # np.save(tempname, [retvPL, retvCL, retvML, outnamesTL])
#
#
#         if j == 0:
#             retvP = retvPL.copy()
#             retvC = retvCL.copy()
#             retvM = retvML.copy()
#             outnames = [[tempname], [outnamesTL]]
#         else:
#             for arname in retvP.keys():
#                 if isinstance(retvP[arname], int) or isinstance(retvP[arname], np.int64):
#                     retvP[arname] = [retvP[arname]]
#                 if isinstance(retvPL[arname], int) or isinstance(retvPL[arname], np.int64):
#                     retvPL[arname] = [retvPL[arname]]
#                 try:
#                     retvP[arname] = np.concatenate((np.asarray(retvP[arname]), np.asarray(retvPL[arname])))
#                 except:
#                     pdb.set_trace()
#             for arname in retvC.keys():
#                 if isinstance(retvC[arname], int):
#                     retvC[arname] = [retvC[arname]]
#                 if isinstance(retvCL[arname], int):
#                     retvCL[arname] = [retvCL[arname]]
#                 retvC[arname] = np.concatenate((np.asarray(retvC[arname]), np.asarray(retvCL[arname])))
#
#             for arname in retvM.keys():
#                 if isinstance(retvM[arname], (int, np.int64, str, datetime.datetime)):
#                     retvM[arname] = [retvM[arname]]
#                 if isinstance(retvML[arname], (int, np.int64, str, datetime.datetime)):
#                     retvML[arname] = [retvML[arname]]
#                 retvM[arname] = np.concatenate((np.asarray(retvM[arname]), np.asarray(retvML[arname])))
#
#             outnames[0].append(tempname)
#             outnames[1].append(outnamesTL)
#
#     return retvP, retvC, retvM, outnames



def convertIrStartToUTC(partf, reshape_col = None):
    #: if the time is repeted in steps it might be 
    #: unessecary to calculate same time several times 
    if reshape_col is None:
        ir_start = partf['ir_start']
    else:
        ir_start = partf['ir_start'].reshape(-1, reshape_col)[:, 0]
    utc = [datetime.datetime(int(str(partf['stamp_date'])[0:4]), int(str(partf['stamp_date'])[4:6]), \
                             int(str(partf['stamp_date'])[6:8]), int(str(partf['stamp_date'])[8:10]), \
                             int(str(partf['stamp_date'])[10:12]), int(str(partf['stamp_date'])[12:14])) + \
                             datetime.timedelta(seconds=int(ist)) for ist in ir_start]
    return utc
def readRedmaskDatetime(dt):
    fn = '/data/legras/SVC/%d/redmask-%d_%02d_%02d.pkl'%(dt.year, dt.year, dt.month, dt.day)
    if not os.path.isfile(fn):
        return {'Nofile': np.nan}
    with gzip.open(fn ,'rb') as f:
        redmask = pickle.load(f)
    return redmask
def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)
def addVarToTemp(fn, var_data, var_name):
    h5file = h5py.File(fn, 'a')
    if var_name in h5file.keys():
        print('%s already exist in tempfile' %var_name)
        pdb.set_trace()
    else:
        try:
            h5file.create_dataset(var_name, data = var_data)
        except:
            print('H5 kan ej spara')
            pdb.set_trace()
    h5file.close()


def findSVC(y, m, utc, pf, mon0_lon, mon0_lat):
    svc_mask_hmask = np.zeros([len(utc), 100])
    svc_mask_lat = np.zeros([len(utc)])
    svc_mask_lon = np.zeros([len(utc)])
    svc_mask_time = np.zeros([len(utc)]).astype(datetime.datetime)
    svc_mask_time_diff = np.zeros([len(utc)])
    used_utc = np.zeros([len(utc)])

    #: Redmask is updated in days while data from flexpart is monthly
    weekDay, nrofdays = calendar.monthrange(y, m)  # @UnusedVariable
    for dom in range(1, nrofdays + 1)[::-1]:
        print('%d of %d' %(dom, nrofdays))
        redmask = {}
        dt_dom = datetime.datetime(y, m, dom)
        #: Add a day before and after to make sure we get the right match at day boundaries
        redmask1 = readRedmaskDatetime(dt_dom - datetime.timedelta(1))
        redmask2 = readRedmaskDatetime(dt_dom)
        redmask3 = readRedmaskDatetime(dt_dom + datetime.timedelta(1))
        for arrname in redmask2.keys():
            #: redmask2, important mask must exist. 
            #: It is the one for the specific day
            if arrname == 'Nofile':
                redmask = {'time': [np.nan], 'hmask': [np.nan], 'time': [np.nan], 'longitude': [np.nan], 'latitude': [np.nan]}
                break
            elif arrname == 'time':
                #: Time in UTC                            
                mask2 = [redmask2['day'] + datetime.timedelta(seconds=int(ist)) for ist in redmask2['time']]
            else:
                mask2 = redmask2[arrname]
            #: This one is 1D and is therfore not important at the moment
            #: Needs to be after mask2 otherwise redmask is not uppdated for Nofile
            #: Do not mather if it is before mask1/3 but saves "some time" if it is before
            if isinstance(redmask2[arrname], datetime.datetime):
                continue
            #: redmask1 and readmask 3 day before and after. 
            #: Not as important, so continu without if not exist
            if 'Nofile' in redmask1.keys():
                #: 2D
                if arrname in ['hmask']:
                    mask1 = np.array([], dtype=np.int64).reshape(0,100)
                else:
                    mask1 = []
            elif arrname == 'time':
                #: Time in UTC
                mask1 = [redmask1['day'] + datetime.timedelta(seconds=int(ist)) for ist in redmask1['time']]
            else:
                mask1 = redmask1[arrname]
            if 'Nofile' in redmask3.keys():
                #: 2D
                if arrname in ['hmask']:
                    mask3 = np.array([], dtype=np.int64).reshape(0,100)
                else:
                    mask3 = []
            elif arrname == 'time':
                #: Time in UTC
                mask3 = [redmask3['day'] + datetime.timedelta(seconds=int(ist)) for ist in redmask3['time']]
            else:
                mask3 = redmask3[arrname]
            redmask.update({arrname: np.concatenate((mask1, mask2, mask3))})
        
        day_ind = ((np.asarray(utc) >= datetime.datetime(y, m, dom, 0, 0, 0)) & \
                   (np.asarray(utc) < (datetime.datetime(y, m, dom, 0, 0, 0) + datetime.timedelta(1))))
        used_utc[day_ind] = used_utc[day_ind] + 1
        if day_ind.sum() == 0:
            continue
        tic = time.time()
        for di in np.where(day_ind)[0]:
            utctid = utc[di]
            if (len(redmask['time']) == 1) and np.isnan(redmask['time'][0]):
                close_ind = 0
                #: This one needs to be for each hight
                #: 1D here but 2D under else:
                svc_mask_hmask[di, : ] = np.repeat(redmask['hmask'][close_ind], 100)
                #: can not take .total_seconds() on a nan value
                svc_mask_time_diff[di] = redmask['time'][close_ind]
            else:
                close_ind = nearest_ind(redmask['time'], utctid)
                svc_mask_hmask[di, :] = redmask['hmask'][close_ind, :]
                svc_mask_time_diff[di] = (utctid - redmask['time'][close_ind]).total_seconds()
            svc_mask_time[di] = redmask['time'][close_ind]
            svc_mask_lon[di] = redmask['longitude'][close_ind]
            svc_mask_lat[di] = redmask['latitude'][close_ind]
        print(time.time() - tic)
    svc_mask_hmask[(used_utc  == 0), :] = -9999
    svc_mask_lon[(used_utc  == 0)] = -9999
    svc_mask_lat[(used_utc  == 0)] = -9999
    svc_mask_time_diff[(used_utc  == 0)] = -9999
    svc_mask_time[(used_utc  == 0)] = -9999
    print('done with SVC')
    return svc_mask_hmask, svc_mask_lon, svc_mask_lat, svc_mask_time_diff, svc_mask_time


def getInitFiles(mD, y, m, dn, uD, lt=True, ct=True):

    outname = 'CALIOP-EAD-%d%02d-%s' %(y, m, dn)
    if uD:
        outname = outname + '-DD'
    # outnamesTL.append(outname)
    trajDir = os.path.join(mD,outname)
    initDir = os.path.join(trajDir, 'Initfiles')
    
    tempname = os.path.join(mD,'TempFiles', '%s-init.h5' %(outname))
    if (lt and os.path.isfile(tempname)):
        lons0_mon, lats0_mon, p0_mon, t0_mon, sh_mon, vod_mon, height_mon, cm_mon, sc_mon, rvs0_mon, iwc0_mon, irs_mon, svc_mon = readTempFileInit(tempname)
        if (irs_mon is None) and (y in [2007, 2008, 2009, 2010]):
            print('Add irs to tempfile')
            pf = readidx107(os.path.join(trajDir,'part_000'),quiet=False)
            utc = convertIrStartToUTC(pf, reshape_col=34)
            svc_hmask, svc_lon, svc_lat, svc_tdiff, svc = findSVC(y, m, utc, pf, lons0_mon, lats0_mon)
            addVarToTemp(tempname, svc_hmask, 'svc_hmask')
            addVarToTemp(tempname, svc_lon, 'svc_lon')
            addVarToTemp(tempname, svc_lat, 'svc_lat')
            addVarToTemp(tempname, svc_tdiff, 'svc_tdiff')
#             addVarToTemp(tempname, svc_utc, 'svc_utc')
            addVarToTemp(tempname, pf['ir_start'], 'ir_start')
            addVarToTemp(tempname, [pf['stamp_date']], 'stamp_date')
            print('SVC is added to the tempfile')
            print('Closing program')
            sys.exit()
            
    else:
        paramname = 'selDardar_Params-%s.pkl' %'-'.join(outname.split('-')[2:4])
        catalogname = paramname.replace('_Params-', '_Calalog-')
        # paramFile = os.path.join(initDir, paramname)
        catalFile = os.path.join(initDir, catalogname)
        pf = readidx107(os.path.join(trajDir,'part_000'),quiet=False)
        cf = getCatalogFile(catalFile, pf, checkForNewFile=False)
        # mf = readParamFile(paramFile)
        pf['x'] = np.where(pf['x'] > 180, pf['x'] - 360, pf['x'])
        pf['x'] = np.where(pf['x'] < -180, pf['x'] + 360, pf['x'])
        lons0_mon = pf['x']
        lats0_mon = pf['y']
        p0_mon = pf['p']
        t0_mon = pf['t']
        sh_mon = cf['sh']
        vod_mon = cf['vod']
        height_mon = cf['height']
        cm_mon = cf['CM']
        sc_mon = cf['SC']
        rvs0_mon = esati_murphy(p0_mon, t0_mon)
        iwc0_mon = convIWC(cf['iwc'], p0_mon, t0_mon, sh_mon)
        
        utc = convertIrStartToUTC(pf, reshape_col=34)
        svc_hmask, svc_lon, svc_lat, svc_tdiff, svc = findSVC(y, m, utc, pf, lons0_mon, lats0_mon)
        
#         def nearest_ind(items, pivot):
#             time_diff = np.abs([date - pivot for date in items])
#             return time_diff.argmin(0)
        
#         fig = plt.figure()
#         fig.suptitle('2008-07-03, Time 01:02 - 01:12')
#         ax = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
#         ax.set_extent([0, 30, -30, 30], crs=ccrs.PlateCarree())
#         ax.coastlines()
#         ax.plot(lons0_mon[org_cut_inds_file], lats0_mon[org_cut_inds_file], 'r+', label='DARDAR')
#         ax.plot(np.asarray(redmask['longitude'])[file_cut_inds], np.asarray(redmask['latitude'])[file_cut_inds], 'b+', label='SVC-MASK')
#         ax.set_xticks([0, 15, 30], crs=ccrs.PlateCarree())
#         ax.set_xticklabels(['0', '15', '30'], fontsize='large')
#         ax.set_yticks([-30, 0, 30], crs=ccrs.PlateCarree())
#         ax.set_yticklabels(['-30', '0', '30'], fontsize='large')
#         ax.set_xlabel('Longitude')
#         ax.set_ylabel('Latitude')
#         ax.legend()
#         
#         ax = fig.add_subplot(1,2,2,projection=ccrs.PlateCarree())
#         ax.set_extent([5, 15, -5, 20], crs=ccrs.PlateCarree())
#         ax.coastlines()
#         ax.plot(lons0_mon[org_cut_inds_file][1:][::18], lats0_mon[org_cut_inds_file][1:][::18], 'r+', label='DARDAR')
#         ax.plot(np.asarray(redmask['longitude'])[file_cut_inds][1:][::18], np.asarray(redmask['latitude'])[file_cut_inds][1:][::18], 'b+', label='SVC-MASK')
#         ax.set_xticks([5, 10, 15], crs=ccrs.PlateCarree())
#         ax.set_xticklabels(['5', '10', '15'], fontsize='large')
#         ax.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
#         ax.set_yticklabels(['0', '10', '20'], fontsize='large')
#         ax.set_xlabel('Longitude')
#         ax.set_ylabel('Latitude')
# #         ax.legend()
#         fig.savefig('svc_mask' + '.png')
        
        if ct:
            print('Create Temp File')
            try:
                h5file = h5py.File(tempname, 'w')
                h5file.create_dataset('lons', data = lons0_mon)
                h5file.create_dataset('lats', data = lats0_mon)
                h5file.create_dataset('p0', data = p0_mon)
                h5file.create_dataset('t0', data = t0_mon)
                
                h5file.create_dataset('sh', data = sh_mon)
                h5file.create_dataset('vod', data = vod_mon)
                h5file.create_dataset('height', data = height_mon)
                h5file.create_dataset('cm', data = cm_mon)
                h5file.create_dataset('sc', data = sc_mon)
                
                h5file.create_dataset('rvs0', data = rvs0_mon)
                h5file.create_dataset('iwc0', data = iwc0_mon)
                
                h5file.create_dataset('svc_hmask', data = svc_hmask)
                h5file.create_dataset('svc_lon', data = svc_lon)
                h5file.create_dataset('svc_lat', data = svc_lat)
                h5file.create_dataset('svc_tdiff', data = svc_tdiff)
                h5file.create_dataset('ir_start', data = pf['ir_start'])
                h5file.create_dataset('stamp_date', data = [pf['stamp_date']])
                
                h5file.close()
            except:
                print('Something wrong in H5')
                sys.exit()
        
    return [tempname, outname], lons0_mon, lats0_mon, p0_mon, t0_mon, sh_mon, vod_mon, height_mon, cm_mon, sc_mon, rvs0_mon, iwc0_mon



def getConvfiles(oD, inames, lt=True, ct=True):
    tempname = inames[0]
    outname = inames[1]
    if tempname != '':
        tempname = tempname.replace('-init', '-conv')
    if lt and os.path.isfile(tempname):
        retvRvsL, retvFlsL, retvAgeL, retvLonsL, retvLatsL, retvTempL, retvPresL = readTempFileConv(tempname)
    else:
        i = -1
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
            retvRvs = np.concatenate((retvRvsL, rvs))
            retvFls = np.concatenate((retvFlsL, fls))
            retvAge = np.concatenate((retvAgeL, age))
            retvLons = np.concatenate((retvLonsL, lons))
            retvLats = np.concatenate((retvLatsL, lats))
            retvTemp = np.concatenate((retvTempL, temp))
            retvPres = np.concatenate((retvPresL, pres))
        
        if ct and (tempname != ''):
            np.save(tempname, [retvRvsL, retvFlsL, retvAgeL, retvLonsL, retvLatsL, retvTempL, retvPresL])
    # if ln == 0:
    #     retvRvs = retvRvsL.copy()
    #     retvFls = retvFlsL.copy()
    #     retvAge = retvAgeL.copy()
    #     retvLons = retvLonsL.copy()
    #     retvLats = retvLatsL.copy()
    #     retvTemp = retvTempL.copy()
    #     retvPres = retvPresL.copy()
    # else:
    #     retvRvs = np.concatenate((retvRvs, retvRvsL))
    #     retvFls = np.concatenate((retvFls, retvFlsL))
    #     retvAge = np.concatenate((retvAge, retvAgeL))
    #     retvLons = np.concatenate((retvLons, retvLonsL))
    #     retvLats = np.concatenate((retvLats, retvLatsL))
    #     retvTemp = np.concatenate((retvTemp, retvTempL))
    #     retvPres = np.concatenate((retvPres, retvPresL))
        
    return retvRvs, retvFls.astype(int), retvAge, retvLons, retvLats, retvTemp, retvPres


def normalizeLons(lon, lat):
    print('Normalize Longitudes')
    retv = lon * np.cos(np.deg2rad(lat))
    return retv


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


def conv2DHistToImage(h, x, y, yb=60):
    xstep = x[3] - x[2]  # @UnusedVariable
    ystep = y[3] - y[2]
    yim = np.arange(-1*yb, yb, ystep)
    yargmin = np.argmin(np.abs(yim - y[0]))
    yargmax = np.argmin(np.abs(yim - y[-1]))
    if (y[0] < -179) and (yim.shape[0] > h.shape[1]):
        use_yshape = h.shape[1]
        # yargmax = use_yshape
    else:
        use_yshape = yim.shape[0]
    try:
        him = np.zeros([h.shape[0], use_yshape])
        him[:, yargmin: yargmax] = h
    except:
        pdb.set_trace()
    return him

def readTempFileInit(fn):
    h5f = h5py.File(fn, 'r')   
    lons = h5f['lons'][:]
    lats = h5f['lats'][:]
    p0 = h5f['p0'][:]
    t0 = h5f['t0'][:]
    
    sh = h5f['sh'][:]
    vod = h5f['vod'][:]
    height = h5f['height'][:]
    cm = h5f['cm'][:]
    sc = h5f['sc'][:]
    
    rvs0 = h5f['rvs0'][:]
    iwc0 = h5f['iwc0'][:]
    
    #: updating the tempfiles with new variables
    if 'ir_start' in h5f.keys():
        irs = h5f['ir_start'][:]
        svc = {'svc_hmask': h5f['svc_hmask'][:], 'svc_lon': h5f['svc_lon'][:], \
               'svc_lat': h5f['svc_lat'][:], 'svc_tdiff': h5f['svc_tdiff'][:], \
               'ir_start': h5f['ir_start'][:], 'stamp_date': h5f['stamp_date'][:]}
    else:
        irs = None
        svc = None
    
    
    h5f.close()
    
    
    # extras = {'use_inds': h5f['use_inds'][:], 'ohClo_use_inds': h5f['ohClo_use_inds'][:], 'ohClr_use_inds': h5f['ohClr_use_inds'][:], \
    #           'ohClo_height': h5f['ohClo_height'][:], 'ohClr_height': h5f['ohClr_height'][:], \
    #           'ohClo_lats': h5f['ohClo_lats'][:], 'ohClo_lons': h5f['ohClo_lons'][:], \
    #           'iwc0': h5f['iwc0'][:], 'age': h5f['age'][:], \
    #           'part0_p': h5f['part0_p'][:], 'part0_t': h5f['part0_t'][:], \
    #           'pressure': h5f['pressure'][:], 'temperature': h5f['temperature'][:]}

    return lons, lats, p0, t0, sh, vod, height, cm, sc, rvs0, iwc0, irs, svc

def readTempFileConv(fn):
    h5f = h5py.File(fn, 'r')   
    lons = h5f['lons'][:]
    lats = h5f['lats'][:]
    p0 = h5f['p0'][:]
    t0 = h5f['t0'][:]
    
    sh = h5f['sh'][:]
    vod = h5f['vod'][:]
    height = h5f['height'][:]
    cm = h5f['cm'][:]
    sc = h5f['sc'][:]
    
    rvs0 = h5f['rvs0'][:]
    iwc0 = h5f['iwc0'][:]
    
    h5f.close()
    
    
    # extras = {'use_inds': h5f['use_inds'][:], 'ohClo_use_inds': h5f['ohClo_use_inds'][:], 'ohClr_use_inds': h5f['ohClr_use_inds'][:], \
    #           'ohClo_height': h5f['ohClo_height'][:], 'ohClr_height': h5f['ohClr_height'][:], \
    #           'ohClo_lats': h5f['ohClo_lats'][:], 'ohClo_lons': h5f['ohClo_lons'][:], \
    #           'iwc0': h5f['iwc0'][:], 'age': h5f['age'][:], \
    #           'part0_p': h5f['part0_p'][:], 'part0_t': h5f['part0_t'][:], \
    #           'pressure': h5f['pressure'][:], 'temperature': h5f['temperature'][:]}

    return lons, lats, p0, t0, sh, vod, height, cm, sc, rvs0, iwc0



if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Program to calculate backvards through part files to find supersaturation', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-c", "--ct", action='store_false', default = True, 
                        help = "Create Temp Files, only valid if no temp file is loaded. Default = True")
    parser.add_argument("-l", "--lt", action='store_false', default = True, 
                        help = "Load Temp Files. Default = True")
    parser.add_argument("-a", "--area", type=int, choices=([*range(20)]), default=0,  
                        help = "Areas. Default = 0\n")
    parser.add_argument("-t", "--temp", action='store_true', default = False, 
                        help = "Create Temp Files. Default = False")
    parser.add_argument("-y", "--year", type=int, choices=([0]+[*range(2007, 2020)]), default=0,  
                        help = "year. Default=2018")
    parser.add_argument("-m", "--month", type=int, choices=(np.arange(-6, 13)), default=0, 
                        help = "Month. Default=6")
    
    args = parser.parse_args()
    
    # if 'ciclad' in socket.gethostname():
        #main_sat_dir = '/data/legras/flexpart_in/SAFNWC'
            #SVC_Dir = '/bdd/CFMIP/SEL2'
    datPath = os.environ['HOME'].replace('/home/', '/data/')
    scrPath = os.environ['HOME'].replace('/home/', '/scratchu/')
    ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson'
    mainDir = '%s/flexout/STC/Calipso' %ekjDir
    outDir = '%s/flexout/STC/Calipso-OUT' %datPath
    plotDir = '%s/LMD-Traczilla/Plots' %ekjDir
    
    
    compare = False
    dn = 'n'
    useDardar = True
    clType = 'clr'
    area = getAreaName(args.area)
    
    lt = args.lt
    
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
        if clType == 'clr':
            years = [*range(2007,2011)]
        else:
            years = [*range(2007,2020)]
        ses = ses
        ses_tit = ses_tit
    else:
        years = [args.year]
        ses = 'y%d-%s' %(args.year, ses)
        ses_tit = 'Year %d, %s' %(args.year, ses_tit)
    #: Create Temp Files
    if args.temp:
        years = [*range(2007,2020)]
        months = [*range(1,13)]
        ct = True
        
    y = -1
    for year in years:
        for mon in months:
            #: Check for missing months
            if year in missing_months.keys():
                if mon in missing_months[year]:
                    continue
            print('Year = %d, Month = %02d' %(year, mon))      
            y  = y + 1
            
            tic=time.time()
            outnames, lons0_mon, lats0_mon, p0_mon, t0_mon, sh0_mon, vod0_mon, height0_mon, cm0_mon, sc0_mon, rvs0_mon, iwc0_mon = getInitFiles(mainDir, year, mon, dn, useDardar, lt=lt)
            if args.temp:
                continue
            
            fname = os.path.join(outDir, outnames[1] + '.h5')
            rvs_mon, fls_mon, age_mon, lons_mon, lats_mon, temp_mon, pres_mon = readConvFile(fname)
            print(time.time() - tic)
            if False:
                lons_mon = normalizeLons(lons_mon, lats_mon)
                lons0_mon = normalizeLons(lons0_mon, lats0_mon)
            useType = 'Hit'
            if useType == 'Hit':
                typeInd = (fls_mon & I_HIT) == I_HIT
                typeName = 'hits'
            elif useType == 'Old':
                typeInd = (fls_mon & I_OLD) == I_OLD
                typeName = 'olds'
            else:
                typeInd = np.ones(fls_mon.shape).astype(bool)
                typeName = 'all'
            lons_mon = checkLons(lons_mon, typeInd)
            
            if clType == 'cld':
                clInd = iwc0_mon > 0
            elif clType == 'cld_thin':
                clInd = (iwc0_mon > 0) & (vod0_mon < 0.03)
            elif clType == 'clr':
                clInd = ~(iwc0_mon > 0)
            else:
                clInd = np.ones(fls_mon.shape).astype(bool)
                clType = 'all'
                
                
            llInd = getAreaInds(area, lats0_mon, lons0_mon)
                
            useInds = typeInd & clInd & llInd
            
            if y == 0:
                lons0 = lons0_mon[useInds]
                lats0 = lats0_mon[useInds]
                vod0 = vod0_mon[useInds]
                height0 = height0_mon[useInds]
            
                lons = lons_mon[useInds]
                lats = lats_mon[useInds]
                fls = fls_mon[useInds]
                age = age_mon[useInds]
                hh1f, xedges, yedges = np.histogram2d(age / 86400, height0*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)
                hh3f = hh1f * np.diff(yedges)
                hh2f, xedges, yedges = np.histogram2d(age / 86400, height0*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]])
            else:
                lons0 = np.concatenate((lons0, lons0_mon[useInds]))
                lats0 = np.concatenate((lats0, lats0_mon[useInds]))
                vod0 = np.concatenate((vod0, vod0_mon[useInds]))
                height0 = np.concatenate((height0, height0_mon[useInds]))
            
                fls = np.concatenate((fls, fls_mon[useInds]))
                lons = np.concatenate((lons, lons_mon[useInds]))
                lats = np.concatenate((lats, lats_mon[useInds]))
                hh1a, xedges, yedges = np.histogram2d(age_mon[useInds] / 86400, height0_mon[useInds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)
                hh3a = hh1a * np.diff(yedges)
                hh2a, xedges, yedges = np.histogram2d(age_mon[useInds] / 86400, height0_mon[useInds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]])
                age = np.concatenate((age, age_mon[useInds]))
                
                
                hh1t, xedges, yedges = np.histogram2d(age / 86400, height0*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)
                hh3t = hh1t * np.diff(yedges)
                hh2t, xedges, yedges = np.histogram2d(age / 86400, height0*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]])
                pdb.set_trace()
                # = np.concatenate(())
                # = np.concatenate(())
                # = np.concatenate(())
            
            
    if args.temp:
        sys.exit()       
    
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

    #: --- Plot ---
    inds_org = np.ones(age.shape[0]).astype(bool)
    if inds_org.sum() == 0:
        print("No match")
        pdb.set_trace()
    #
    # 2D-Histograms, Height vs Age 
    print('Plot 2D Histograms - height')
    # for ctyp in ['all', 'cld', 'cld_thin', 'clr']:
    figname_end = '%s_%s_%s' %(ses, area.replace(' ', '_').replace('(', '').replace(')', ''),  useType.lower())
    if clType == 'all':
        title_org = '%s pixels' %useType.title()
        figname_end = '%s_all' %(figname_end)
        vmax = 100000
    elif clType == 'cld':
        title_org = '%s and Cloudy pixels' %useType.title()
        figname_end = '%s_cld' %(figname_end)
        vmax1 = 0.003
        vmax2 = 10000
    elif clType == 'cld_thin':
        title_org = '%s and Thin Cloudy pixels' %useType.title()
        figname_end = '%s_cld_thin' %(figname_end)
        vmax = 10000
    elif clType == 'clr':
        title_org = '%s and Clear pixels' %useType.title()
        figname_end = '%s_clr' %(figname_end)
        vmax1 = 0.003
        vmax2 = 10000
        vmax = 50000
    figname_org = '%s/2dhist_height_%s' %(plotDir, figname_end)
    inds = np.ones(age.shape[0]).astype(bool)
    figname = figname_org
    title = title_org
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(2,1,1)
    # hh = ax.hist2d(age[inds] / 86400, height0[inds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)#, vmin=0, vmax=vmax)#, bins=400)
    hh1e, xedges, yedges = np.histogram2d(age[inds] / 86400, height0[inds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)
    hh1 = (hh1e * np.diff(yedges))         
    aspect = float('%.2f' %((1.) / ((yedges[0] - yedges[-1]) / (xedges[0] - xedges[-1]))))
    im = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]], vmin=0, vmax=vmax1)
    # ax.set_xlabel('Age [days]', fontsize='x-large')
    ax.set_ylabel('Height [km]', fontsize='x-large')
    ax.set_yticklabels((ax.get_yticks()/1000).astype(int))
    # fig.subplots_adjust(right=0.89)
    cbar = fig.colorbar(im, label='[1 / day]')#, cax=cbar_ax)
    ax = fig.add_subplot(2,1,2)
    # hh = ax.hist2d(age[inds] / 86400, height0[inds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]], density=True)#, vmin=0, vmax=vmax)#, bins=400)
    hh2, xedges, yedges = np.histogram2d(age[inds] / 86400, height0[inds]*1000, bins=[[*range(0,201)], [*range(14000,20001, 500)]])
    aspect = float('%.2f' %((1.) / ((yedges[0] - yedges[-1]) / (xedges[0] - xedges[-1]))))
    # aspect = float('%.2f' %((1.) / ((yedges.shape[0]) / (xedges.shape[0]))))
    im = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]], vmin=0, vmax=vmax2)
    ax.plot((np.diff(xedges)/2) + xedges[0:-1], np.argmax(hh1, axis=1)*500+14000, 'r')
    ax.set_xlabel('Age [days]', fontsize='x-large')
    ax.set_ylabel('Height [km]', fontsize='x-large')
    ax.set_yticklabels((ax.get_yticks()/1000).astype(int))
    fig.subplots_adjust(right=0.89)
    cbar = fig.colorbar(im)#, label='[1 / day]')#, cax=cbar_ax)


    
    # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    # cbar = fig.colorbar(hh[-1])#, label='[1 / day]')#, cax=cbar_ax)
    
    print(figname)
    fig.savefig(figname + '.png')
    pdb.set_trace()

    heightBoundaries = [[None, None], [14,16], [16, 18], [18, 20], [18,19], [19,20]]
    #: https://stackoverflow.com/questions/50611018/cartopy-heatmap-over-openstreetmap-background/50638691
    print('Plot Age Histograms')

    figname_org = '%s/hist_age_%s' %(plotDir, figname_end)
    for h1, h2 in heightBoundaries:
        if not ((h1 is None) and (h2 is None)):
            title = title_org + ', %d - %d km' %(h1, h2)
            figname = figname_org + '_%d-%d' %(h1, h2)
            if h2 == heightBoundaries[-1][1]:
                inds_h = (height0 >= h1) & (height0 <= h2)
            else:
                inds_h = (height0 >= h1) & (height0 < h2)
            inds = inds_h
        else:
            title = title_org
            figname = figname_org
            inds = np.ones(age.shape[0]).astype(bool)

        if inds.sum() == 0:
            continue
           
        #: --- Plot ---
        fig = plt.figure()
        ax = fig.add_subplot(2,1,1)
        # fig.suptitle('Age histogram')
        h = ax.hist(age[inds] / 86400, bins=range(0,201), density=True)
        # h = ax.hist(age[hits] / 86400, bins=np.logspace(np.log10(0.001),np.log10(42.0), 400))#bins=400)#, density=True)
        # ax.set_xscale('log')
        # fig.gca().set_xscale("log")
        # ax.set_xlabel('Age [days]')
        ax.set_title(title)
        ax.text(0.7, 0.9,'total = %d' %inds.sum(),
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        ax.set_ylabel('[1 / day]')
        ax = fig.add_subplot(2,1,2)
        h = ax.hist(age[inds] / 86400, bins=range(0,201), density=True)
        ax.set_yscale('log')
        ax.set_xlabel('Age [days]')
        # ax.set_title('Hits pixels')
        plt.tight_layout()
        print(figname)
        fig.savefig(figname + '.png')
        plt.close(fig)

    
    #: --- Plot ---
    #: trajectory
    if False:
        tic = time.time()
        fig = plt.figure()
        ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 190, -90, 90], crs=ccrs.PlateCarree())
        ax.coastlines()
        # [0::500]
        step = int(lons.shape[0]/250)
        plt.plot([lons[0::step], lons0[0::step]], [lats[0::step], lats0[0::step]],
             linewidth=1, transform=ccrs.Geodetic())
        plt.tight_layout()
        fig.savefig('%s/traj.png' %(plotDir))
        toc = time.time()
        print(toc-tic)
    
    #: TODO: Height, cloudmask, 
    
    #: --- Plot ---
    # from matplotlib import ticker
    #: --- All pixels ---
    if False:
        hh1, xedges1, yedges1 = np.histogram2d(lons, lats, bins=[180, 90])
        hh2, xedges2, yedges2 = np.histogram2d(lons0, lats0, bins=[180, 90])#, bins=400)#, density=True)
        ymax = getYmax(yedges1, yedges2)
        aspect = float('%.2f' %((1/3) / (ymax / 180.)))
        # hhim1 = conv2DHistToImage(hh1, xedges1, yedges1, ymax)
        # hhim2 = conv2DHistToImage(hh2, xedges2, yedges2, ymax)
        
        fig = plt.figure()
        fig.suptitle('2D Hist - All pixels')
        ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
        ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax], fontsize='large')
        ax.tick_params(axis=u'both', which=u'both',length=0)
        im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
        ax.set_title('Trazilla', fontsize='x-large')
        # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False)#, linewidth=0.001, color='k',alpha=0)
        # gl.right_labels = gl.top_labels = gl.bottom_labels = gl.ylines = gl.xlines = False
        # gl.ylocator = ticker.FixedLocator([-1*ymax+1, -1*ymax//2, 0, ymax//2, ymax-1])
        # glLS = '%d%s' %(ymax, gl.ylabel_artists[0].get_text()[2:])
        # gl.ylabel_artists[0].set_text(glLS)
        # gl.ylabel_artists[-2].set_text('%d%s' %(ymax, gl.ylabel_artists[-2].get_text()[2:]))
        fig.subplots_adjust(right=0.89)
        cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
        cbar = fig.colorbar(im1, cbar_ax)
        ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
        ax.set_extent([-180, 190, -1*ymax, ymax], crs=ccrs.PlateCarree())
        ax.coastlines()
        im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
        ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
        ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax], fontsize='large')
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.set_title('Dardar', fontsize='x-large')
        fig.subplots_adjust(right=0.89)
        # cbar_ax = fig.add_axes([0.90, 0.06, 0.01, 0.40])
        cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
        cbar = fig.colorbar(im2, cax=cbar_ax)
        # fig.tight_layout(rect=[0, 0.03, 0.97, 0.97])
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        # plt.tight_layout()
        figname = '%s/2dhist_%s_%s_all' %(plotDir, ses, area.replace(' ', '_'))
        fig.savefig(figname + '.png')
    
    #: --- Plot ---
    #: --- Maps ---
    print('Maps - 2D-Histograms')
    figname_org = '%s/maps_%s' %(plotDir, figname_end)
    #: Height
    for h1, h2 in heightBoundaries:
        if not ((h1 is None) and (h2 is None)):
            title = title_org + ', %d - %d km' %(h1, h2)
            figname = figname_org + '_%d-%d' %(h1, h2)
            if h2 == heightBoundaries[-1][1]:
                inds_h = (height0 >= h1) & (height0 <= h2)
            else:
                inds_h = (height0 >= h1) & (height0 < h2)
            inds = inds_org & inds_h
        else:
            title = title_org
            figname = figname_org
            inds = inds_org
        if inds.sum() == 0:
            continue
        #: 2D Histogram
        hh1, xedges1, yedges1 = np.histogram2d(lons[inds], lats[inds], bins=[180, 90])
        hh2, xedges2, yedges2 = np.histogram2d(lons0[inds], lats0[inds], bins=[180, 90])
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
        ax.set_title('Trazilla')
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
        ax.set_title('Dardar')
        fig.subplots_adjust(right=0.89)
        cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
        cbar = fig.colorbar(im2, cax=cbar_ax)
        print(figname)
        fig.savefig(figname + '.png')
        plt.close(fig)
    
    
    
    
    
    