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


def getCatalogFile(fn, p0=None):
    tf = '/data/ejohansson/flexout/STC/Calipso/CALIOP-EAD/Catalog/%s' %os.path.basename(fn)
    
    objects = readCatalogFile(fn)
    if os.path.isfile(tf):
        to = readCatalogFile(tf)
        compareFiles(objects, fn, to, tf, 0)
    
    nameDic = {'CM': 'Cloud_Mask', 'TpH': 'Tropopause_Height', \
                'SZA': 'SZA', 'SC': 'Simplified_Categorization', \
                'lensel': 'lensel', 'npart': 'npart', 'vod': 'vis_optical_depth', \
                'height': 'Track_Height'}
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


def getInitfiles(mD, years, months, dn, uD, lt=True):
    if (not isinstance(years, list)) or (not isinstance(months, list)):
        print('Wrong type')
        pdb.set_trace()
    
    #: Decide the name for tempfile
    if (len(years) == 1) and (len(months) == 12):
        tempname = os.path.join(mD,'TempFiles', 'CALIOP-EAD-%d-%s' %(years[0], dn))
        if uD:
            tempname = tempname + '-DD'
        tempname = tempname + '-init'
        loadname = tempname + '.npy'
    else:
        loadname = ''
    #: If tempfile exist and lt = True, load the tempfile
    if lt and os.path.isfile(loadname):
        retvP, retvC, retvM, outnamesT = np.load(loadname, allow_pickle=True)
        outnames = [tempname, outnamesT]
    #: if tempfile do not exist or lt = False. Do calculations
    else:
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
                cf = getCatalogFile(catalFile, pf)
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
        #: If there are no tempfile, create one
        if (loadname != '') and (not os.path.isfile(loadname)):
            np.save(tempname, [retvP, retvC, retvM, outnames])
    return retvP, retvC, retvM, outnames


def getConvfiles(oD, inames, lt=True):
    i = -1
    if (len(inames) == 2) and (isinstance(inames[1], list)):
        tempname = inames[0].replace('-init', '-conv')
        loadname = tempname + '.npy'
        onames = inames[1]
    else:
        onames = inames
        tempname = ''
        loadname = ''
    if lt and os.path.isfile(loadname):
        retvRvs, retvFls, retvAge, retvLons, retvLats, retvTemp, retvPres = np.load(loadname, allow_pickle=True)
    else:
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
        if not os.path.isfile(loadname) and (loadname != ''):
            np.save(tempname, [retvRvs, retvFls, retvAge, retvLons, retvLats, retvTemp, retvPres])
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
    year = [2007]
    month = [1]
    #month = [*range(1, 13)]
    dn = 'n'
    useDardar = True
    lt = False
    # outDir = '/data/ejohansson/flexout/STC/Calipso-OUT/Coldpoint'
    #: filename of outfiles
    # outnames = 'CALIOP-'+advect+'-'+date_end.strftime('%b%Y')+diffus+'-%s' %args.night 
    # if args.use_dardar:
    #     outnames = outnames + '-DD'
    print('Read Init-files')
    tic = time.time()
    part0, catalog, params, outnames = getInitfiles(mainDir, year, month, dn, useDardar, lt=lt)
    print('It took %d sec to read Init-files' %(time.time() - tic))
    #: Read the index file that contains the initial positions
    # print('numpart',part0['numpart'])
    #:
    
    lons_0 = part0['x']#np.where(part0['x'] > 180, part0['x'] - 360, part0['x'])
    lats_0 = part0['y']
    
    print('Read Conv-files')
    tic = time.time()
    rvs, fls, age, lons, lats, temp, pres = getConvfiles(outDir, outnames, lt=lt)
    print('It took %d sec to read Conv-files' %(time.time() - tic))
    
    # print(time.time() - tic)
    olds = (fls & I_OLD) == I_OLD
    hits = (fls & I_HIT) == I_HIT
    
    lons = checkLons(lons, hits)
    
    if False:
        lons = normalizeLons(lons, lats)
        lons_0 = normalizeLons(lons_0, lats_0)
    # originDate, idxDates = getDates(part0, hits)
    # cloudy = (catalog['CM'] > 0)
    # cloudy = ((catalog['SC']>0) & (catalog['SC']<5))
    cloudy = (catalog['vod'] > 0)
    cloudy_thin = cloudy & (catalog['vod'] < 0.03)
    hits_cld = hits & cloudy
    hits_clr = hits & ~cloudy
    hits_cldt = hits & cloudy_thin
    
    
    olds_cld = olds & cloudy
    olds_clr = olds & ~cloudy
    olds_cldt = olds & cloudy_thin
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
    #: --- Plot ---
    #
    print('Plot 2D Histograms - height')
    for ctyp in ['all', 'cld', 'cld_thin', 'clr']:
        if ctyp == 'all':
            title_org = '2D Age - Hits pixels'
            figname_org = '%s/2dhist_height_hits_all' %plotDir
            inds_org = hits
            vmax = 100000
        elif ctyp == 'cld':
            title_org = '2D Age - Hits and Cloudy pixels'
            figname_org = '%s/2dhist_height_hits_cld' %plotDir
            inds_org = hits_cld
            vmax = 100000
        elif ctyp == 'cld_thin':
            title_org = '2D Age - Hits and Thin Cloudy pixels'
            figname_org = '%s/2dhist_height_hits_cld_thin' %plotDir
            inds_org = hits_cldt
            vmax = 10000
        elif ctyp == 'clr':
            title_org = '2D Hist - Hits and Clear pixels'
            figname_org = '%s/2dhist_height_hits_clr' %plotDir
            inds_org = hits_clr
            vmax = 50000
        inds = inds_org
        figname = figname_org
        title = title_org
        fig = plt.figure()
        ax = fig.add_subplot(111)
        hh = ax.hist2d(age[inds] / 86400, catalog['height'][inds], bins=[100, 10], vmin=0, vmax=vmax)#, bins=400)
        ax.set_ylabel('Height')
        ax.set_xlabel('Age [days]')
        ax.set_title(title)
        fig.subplots_adjust(right=0.89)
        # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
        cbar = fig.colorbar(hh[3])#, cax=cbar_ax)
        fig.savefig(figname + '.png')
    pdb.set_trace()
    
    #: https://stackoverflow.com/questions/50611018/cartopy-heatmap-over-openstreetmap-background/50638691
    heightBoundaries = [[None, None], [14,16], [16, 18], [18,20]]
    print('Plot Age Histograms')
    for ctyp in ['all', 'cld', 'cld_thin', 'clr']:
        if ctyp == 'all':
            title_org = 'Hits pixels'
            figname_org = '%s/hist_age_all' %plotDir
            inds_org = hits
        elif ctyp == 'cld':
            title_org = 'Hits and Cloudy pixels'
            figname_org = '%s/hist_age_cld' %plotDir
            inds_org = hits_cld
        elif ctyp == 'cld_thin':
            title_org = 'Hits and Thin Cloudy pixels'
            figname_org = '%s/hist_age_cld_thin' %plotDir
            inds_org = hits_cldt
        elif ctyp == 'clr':
            title_org = 'Hits pixels Clear pixels'
            figname_org = '%s/hist_age_clr' %plotDir
            inds_org = hits_clr
        
        for h1, h2 in heightBoundaries:
            if not ((h1 is None) and (h2 is None)):
                title = title_org + ', %d - %d km' %(h1, h2)
                figname = figname_org + '_%d-%d' %(h1, h2)
                if h2 == heightBoundaries[-1][1]:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] <= h2)
                else:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] < h2)
                inds = inds_org & inds_h
            else:
                title = title_org
                figname = figname_org
                inds = inds_org
                
            #: --- Plot ---
            fig = plt.figure()
            ax = fig.add_subplot(2,1,1)
            # fig.suptitle('Age histogram')
            h = ax.hist(age[inds] / 86400, bins=400)#, density=True)
            # h = ax.hist(age[hits] / 86400, bins=np.logspace(np.log10(0.001),np.log10(42.0), 400))#bins=400)#, density=True)
            # ax.set_xscale('log')
            # fig.gca().set_xscale("log")
            # ax.set_xlabel('Age [days]')
            ax.set_title(title)
            ax.text(0.7, 0.9,'total = %d' %inds.sum(),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = ax.transAxes)
            ax = fig.add_subplot(2,1,2)
            h = ax.hist(age[inds] / 86400, bins=400)#, density=True)
            ax.set_yscale('log')
            ax.set_xlabel('Age [days]')
            # ax.set_title('Hits pixels')
            plt.tight_layout()
            fig.savefig(figname + '.png')
            plt.close(fig)
    # pdb.set_trace()

    # #: --- Plot ---
    # fig = plt.figure()
    # ax = fig.add_subplot(2,1,1)
    # ax.hist(age[hits_cld] / 86400, bins=400)#, density=True)
    # ax.set_title('Hits and Cloudy pixels')
    # ax.text(0.7, 0.9,'total = %d' %hits_cld.sum(),
    #         horizontalalignment='center',
    #         verticalalignment='center',
    #         transform = ax.transAxes)
    # ax = fig.add_subplot(2,1,2)
    # ax.hist(age[hits_cld] / 86400, bins=400)#, density=True)
    # ax.set_yscale('log')
    # ax.set_xlabel('Age [days]')
    # plt.tight_layout()
    # figname = '%s/hist_age_cld' %plotDir
    # fig.savefig(figname + '.png')
    #
    # #: --- Plot ---
    # fig = plt.figure()
    # # fig.suptitle('Age histogram')
    # ax = fig.add_subplot(2,1,1)
    # ax.hist(age[hits_clr] / 86400, bins=400)#, density=True)
    # ax.set_title('Hits and Clear pixels')
    # ax.text(0.7, 0.9,'total = %d' %hits_clr.sum(),
    #         horizontalalignment='center',
    #         verticalalignment='center',
    #         transform = ax.transAxes)
    # ax = fig.add_subplot(2,1,2)
    # ax.hist(age[hits_clr] / 86400, bins=400)#, density=True)
    # ax.set_yscale('log')
    # ax.set_xlabel('Age [days]')
    # plt.tight_layout()
    # figname = '%s/hist_age_clr' %plotDir
    # fig.savefig(figname + '.png')
    #
    # #: --- Plot ---
    # fig = plt.figure()
    # ax = fig.add_subplot(2,1,1)
    # ax.hist(age[hits_cldt] / 86400, bins=400)#, density=True)
    # ax.set_title('Hits and Thin Cloudy pixels')
    # ax.text(0.7, 0.9,'total = %d' %hits_cldt.sum(),
    #         horizontalalignment='center',
    #         verticalalignment='center',
    #         transform = ax.transAxes)
    # ax = fig.add_subplot(2,1,2)
    # ax.hist(age[hits_cldt] / 86400, bins=400)#, density=True)
    # ax.set_yscale('log')
    # ax.set_xlabel('Age [days]')
    # plt.tight_layout()
    # figname = '%s/hist_age_cld_thin' %plotDir
    # fig.savefig(figname + '.png')

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
    
    #: TODO: Height, cloudmask, 
    
    #: --- Plot ---
    # from matplotlib import ticker
    #: --- All pixels ---

    hh1, xedges1, yedges1 = np.histogram2d(lons, lats, bins=[180, 90])
    hh2, xedges2, yedges2 = np.histogram2d(lons_0, lats_0, bins=[180, 90])#, bins=400)#, density=True)
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
    ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
    ax.set_title('Trazilla')
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
    ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.set_title('Dardar')
    fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.06, 0.01, 0.40])
    cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    cbar = fig.colorbar(im2, cax=cbar_ax)
    # fig.tight_layout(rect=[0, 0.03, 0.97, 0.97])
    # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.tight_layout()
    figname = '%s/2dhist_all' %plotDir
    fig.savefig(figname + '.png')
    
    #: --- Plot ---
    #: --- Hits pixels ---
    print('Plot 2D Histograms - Hits')
    for ctyp in ['all', 'cld', 'cld_thin', 'clr']:
        if ctyp == 'all':
            title_org = '2D Hist - Hits pixels'
            figname_org = '%s/2dhist_hits_all' %plotDir
            inds_org = hits
        elif ctyp == 'cld':
            title_org = '2D Hist - Hits and Cloudy pixels'
            figname_org = '%s/2dhist_hits_cld' %plotDir
            inds_org = hits_cld
        elif ctyp == 'cld_thin':
            title_org = '2D Hist - Hits and Thin Cloudy pixels'
            figname_org = '%s/2dhist_hits_cld_thin' %plotDir
            inds_org = hits_cldt
        elif ctyp == 'clr':
            title_org = '2D Hist - Hits and Clear pixels'
            figname_org = '%s/2dhist_hits_clr' %plotDir
            inds_org = hits_clr
        #: Height
        for h1, h2 in heightBoundaries:
            if not ((h1 is None) and (h2 is None)):
                title = title_org + ', %d - %d km' %(h1, h2)
                figname = figname_org + '_%d-%d' %(h1, h2)
                if h2 == heightBoundaries[-1][1]:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] <= h2)
                else:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] < h2)
                inds = inds_org & inds_h
            else:
                title = title_org
                figname = figname_org
                inds = inds_org
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
            im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])
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
            im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])
            ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
            ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_title('Dardar')
            fig.subplots_adjust(right=0.89)
            cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
            cbar = fig.colorbar(im2, cax=cbar_ax)
            fig.savefig(figname + '.png')
            plt.close(fig)
    # pdb.set_trace()
        
    #: --- Plot ---
    #: --- OLD Pixels ---
    print('Plot 2D Histograms - Old')
    for ctyp in ['all', 'cld', 'cld_thin', 'clr']:
        if ctyp == 'all':
            title_org = '2D Hist - Old pixels'
            figname_org = '%s/2dhist_olds_all' %plotDir
            inds_org = olds
        elif ctyp == 'cld':
            title_org = '2D Hist - Old and Cloudy pixels'
            figname_org = '%s/2dhist_olds_cld' %plotDir
            inds_org = olds_cld
        elif ctyp == 'cld_thin':
            title_org = '2D Hist - Old and Thin Cloudy pixels'
            figname_org = '%s/2dhist_olds_cld_thin' %plotDir
            inds_org = olds_cldt
        elif ctyp == 'clr':
            title_org = '2D Hist - Old and Clear pixels'
            figname_org = '%s/2dhist_olds_clr' %plotDir
            inds_org = olds_clr
        #: Height
        for h1, h2 in heightBoundaries:
            if not ((h1 is None) and (h2 is None)):
                title = title_org + ', %d - %d km' %(h1, h2)
                figname = figname_org + '_%d-%d' %(h1, h2)
                if h2 == heightBoundaries[-1][1]:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] <= h2)
                else:
                    inds_h = (catalog['height'] >= h1) & (catalog['height'] < h2)
                inds = inds_org & inds_h
            else:
                title = title_org
                figname = figname_org
                inds = inds_org
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
            im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])
            ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
            ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_title('Trazilla')
            fig.subplots_adjust(right=0.89)
            cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
            cbar = fig.colorbar(im1, cbar_ax)
            ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
            ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
            ax.coastlines()
            im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])
            ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
            ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_title('Dardar')
            fig.subplots_adjust(right=0.89)
            cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
            cbar = fig.colorbar(im2, cax=cbar_ax)
            # figname = '%s/2dhist_olds_all' %plotDir
            fig.savefig(figname + '.png') 
            plt.close(fig)
    
    
    
#     hh1, xedges1, yedges1 = np.histogram2d(lons[hits], lats[hits], bins=[180, 90])
#     hh2, xedges2, yedges2 = np.histogram2d(lons_0[hits], lats_0[hits], bins=[180, 90])#, bins=400)#, density=True)
#     ymax = getYmax(yedges1, yedges2)
#     aspect = float('%.2f' %((1/3) / (ymax / 180.)))
#     # hhim1 = conv2DHistToImage(hh1, xedges1, yedges1, ymax)
#     # hhim2 = conv2DHistToImage(hh2, xedges2, yedges2, ymax)
#     fig = plt.figure()
#     fig.suptitle('2D Hist - Hits pixels')
#     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Trazilla')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
#     cbar = fig.colorbar(im1, cbar_ax)
#
#     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Dardar')
#     fig.subplots_adjust(right=0.89)
#     # cbar_ax = fig.add_axes([0.90, 0.06, 0.01, 0.40])
#     cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
#     cbar = fig.colorbar(im2, cax=cbar_ax)
#     # fig.tight_layout(rect=[0, 0.03, 0.97, 0.97])
#     # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#     # plt.tight_layout()
#     figname = '%s/2dhist_hits_all' %plotDir
#     fig.savefig(figname + '.png')
# #     fig = plt.figure()
# #     fig.suptitle('2D Hist - Hits pixels')
# #     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
# #     ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
# #     ax.coastlines()
# #     # ax.gridlines(draw_labels=True)
# #     # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, linewidth=2.001, color='k',alpha=0)
# #     # gl.right_labels = gl.top_labels = False
# #     # gl.ylocator = ticker.FixedLocator([60, -30, 0, 30, 60])
# #     # gl.xlocator = ticker.FixedLocator([-180, -90, 0, 90, 180])
# #     # do coordinate conversion of (x,y)
# #     # xynps = ax.projection.transform_points(ccrs.Geodetic(), x, y)
# #     # im1 = ax.imshow(hh1, vmin=vmin, vmax=vmax)
# #     hh = ax.hist2d(lons[hits], lats[hits], bins=[180, 90])#, bins=400)#, density=True)
# #     fig.subplots_adjust(right=0.9)
# #     # cbar_ax = fig.add_axes([0.91, 0.54, 0.008, 0.24])
# #     # fig.colorbar(hh[3], cax=cbar_ax)#, orientation='horizontal')
# #     # import matplotlib as mpl
# #     # cb1 = mpl.colorbar.ColorbarBase(ax, orientation='horizontal')
# #     ax.set_xlabel('Longitude')
# #     ax.set_ylabel('Latitude')
# #     ax.set_title('Trazilla')
# #     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
# #     ax.set_extent([-180, 190, -60, 60], crs=ccrs.PlateCarree())
# #     ax.coastlines()
# #     # im2 = ax.imshow(hh2, vmin=vmin, vmax=vmax, extent = [-180, 180, -60, 60])
# #     hh = ax.hist2d(lons_0[hits], lats_0[hits], bins=[180, 90])#, bins=400)#, density=True)
# #     # cbar_ax = fig.add_axes([pos_x, axpos.y0, 0.01, axpos.height])
# #     # cbar = fig.colorbar(im, orientation='vertical', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
# # #   cbar=fig.colorbar(im)#, cax=pos_cax)
# #     # cbar.ax.tick_params(labelsize=fs)
# #     # divider = make_axes_locatable(ax)
# #     # cax = divider.append_axes("right", size="5%", pad=0.05)
# #     fig.subplots_adjust(right=0.9)
# #     # cbar_ax = fig.add_axes([0.91, 0.04, 0.008, 0.24])
# #     # cbar = fig.colorbar(hh[3], cax=cbar_ax)#, orientation='horizontal')
# #     ax.set_xlabel('Longitude')
# #     ax.set_ylabel('Latitude')
# #     ax.set_title('Dardar')
# #
# #     # plt.tight_layout()
# #     fig.savefig('%s/2dhist_hits_all.png' %plotDir)
# #
#
#     hh1, xedges1, yedges1 = np.histogram2d(lons[hits_cld], lats[hits_cld], bins=[180, 90])
#     hh2, xedges2, yedges2 = np.histogram2d(lons_0[hits_cld], lats_0[hits_cld], bins=[180, 90])#, bins=400)#, density=True)
#     ymax = getYmax(yedges1, yedges2)
#     aspect = float('%.2f' %((1/3) / (ymax / 180.)))
#     fig = plt.figure()
#     fig.suptitle('2D Hist - Hits and Cloudy pixels')
#     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Trazilla')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
#     cbar = fig.colorbar(im1, cbar_ax)
#     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Dardar')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
#     cbar = fig.colorbar(im2, cax=cbar_ax)
#     figname = '%s/2dhist_hits_cld' %plotDir
#     fig.savefig(figname + '.png')
#
#     hh1, xedges1, yedges1 = np.histogram2d(lons[hits_clr], lats[hits_clr], bins=[180, 90])
#     hh2, xedges2, yedges2 = np.histogram2d(lons_0[hits_clr], lats_0[hits_clr], bins=[180, 90])#, bins=400)#, density=True)
#     ymax = getYmax(yedges1, yedges2)
#     aspect = float('%.2f' %((1/3) / (ymax / 180.)))
#     fig = plt.figure()
#     fig.suptitle('2D Hist - Hits and Clear pixels')
#     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Trazilla')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
#     cbar = fig.colorbar(im1, cbar_ax)
#     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Dardar')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
#     cbar = fig.colorbar(im2, cax=cbar_ax)
#     figname = '%s/2dhist_hits_clr' %plotDir
#     fig.savefig(figname + '.png')
#
#     hh1, xedges1, yedges1 = np.histogram2d(lons[hits_cldt], lats[hits_cldt], bins=[180, 90])
#     hh2, xedges2, yedges2 = np.histogram2d(lons_0[hits_cldt], lats_0[hits_cldt], bins=[180, 90])#, bins=400)#, density=True)
#     ymax = getYmax(yedges1, yedges2)
#     aspect = float('%.2f' %((1/3) / (ymax / 180.)))
#     fig = plt.figure()
#     fig.suptitle('2D Hist - Hits and Thin Cloudy pixels')
#     ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Trazilla')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
#     cbar = fig.colorbar(im1, cbar_ax)
#     ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
#     ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
#     ax.coastlines()
#     im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
#     ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
#     ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
#     ax.tick_params(axis=u'both', which=u'both',length=0)
#     ax.set_title('Dardar')
#     fig.subplots_adjust(right=0.89)
#     cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
#     cbar = fig.colorbar(im2, cax=cbar_ax)
#     figname = '%s/2dhist_hits_cld_thin' %plotDir
#     fig.savefig(figname + '.png')
#
#
#     pdb.set_trace()
#
    # #: --- Plot ---
    # #: --- OLD Pixels ---
    # hh1, xedges1, yedges1 = np.histogram2d(lons[olds], lats[olds], bins=[180, 90])
    # hh2, xedges2, yedges2 = np.histogram2d(lons_0[olds], lats_0[olds], bins=[180, 90])#, bins=400)#, density=True)
    # ymax = getYmax(yedges1, yedges2)
    # aspect = float('%.2f' %((1/3) / (ymax / 180.)))
    # fig = plt.figure()
    # fig.suptitle('2D Hist - Old pixels')
    # ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Trazilla')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
    # cbar = fig.colorbar(im1, cbar_ax)
    # ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Dardar')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    # cbar = fig.colorbar(im2, cax=cbar_ax)
    # figname = '%s/2dhist_olds_all' %plotDir
    # fig.savefig(figname + '.png')
    #
    # hh1, xedges1, yedges1 = np.histogram2d(lons[olds_cld], lats[olds_cld], bins=[180, 90])
    # hh2, xedges2, yedges2 = np.histogram2d(lons_0[olds_cld], lats_0[olds_cld], bins=[180, 90])#, bins=400)#, density=True)
    # ymax = getYmax(yedges1, yedges2)
    # aspect = float('%.2f' %((1/3) / (ymax / 180.)))
    # fig = plt.figure()
    # fig.suptitle('2D Hist - Old and Cloudy pixels')
    # ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Trazilla')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
    # cbar = fig.colorbar(im1, cbar_ax)
    # ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Dardar')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    # cbar = fig.colorbar(im2, cax=cbar_ax)
    # figname = '%s/2dhist_olds_cld' %plotDir
    # fig.savefig(figname + '.png')
    #
    # hh1, xedges1, yedges1 = np.histogram2d(lons[olds_clr], lats[olds_clr], bins=[180, 90])
    # hh2, xedges2, yedges2 = np.histogram2d(lons_0[olds_clr], lats_0[olds_clr], bins=[180, 90])#, bins=400)#, density=True)
    # ymax = getYmax(yedges1, yedges2)
    # aspect = float('%.2f' %((1/3) / (ymax / 180.)))
    # fig = plt.figure()
    # fig.suptitle('2D Hist - Old and Clear pixels')
    # ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Trazilla')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
    # cbar = fig.colorbar(im1, cbar_ax)
    # ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Dardar')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    # cbar = fig.colorbar(im2, cax=cbar_ax)
    # figname = '%s/2dhist_olds_clr' %plotDir
    # fig.savefig(figname + '.png')
    #
    # hh1, xedges1, yedges1 = np.histogram2d(lons[olds_cldt], lats[olds_cldt], bins=[180, 90])
    # hh2, xedges2, yedges2 = np.histogram2d(lons_0[olds_cldt], lats_0[olds_cldt], bins=[180, 90])#, bins=400)#, density=True)
    # ymax = getYmax(yedges1, yedges2)
    # aspect = float('%.2f' %((1/3) / (ymax / 180.)))
    # fig = plt.figure()
    # fig.suptitle('2D Hist - Old and Thin Cloudy pixels')
    # ax = fig.add_subplot(2,1,1,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im1 = ax.imshow(hh1.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges1[0], yedges1[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Trazilla')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.55, 0.01, 0.30])
    # cbar = fig.colorbar(im1, cbar_ax)
    # ax = fig.add_subplot(2,1,2,projection=ccrs.PlateCarree())
    # ax.set_extent([-180, 180, -1*ymax, ymax], crs=ccrs.PlateCarree())
    # ax.coastlines()
    # im2 = ax.imshow(hh2.T, origin ='lower', aspect=aspect, extent = [-180, 180, yedges2[0], yedges2[-1]])#, vmin=vmin, vmax=vmax)
    # ax.set_yticks([-1*ymax, -1*ymax//2, 0, ymax//2, ymax], crs=ccrs.PlateCarree())
    # ax.set_yticklabels([r'$%d\degree S$' %ymax, r'$%d\degree S$' %(ymax//2), r'$0\degree$', r'$%d\degree N$' %(ymax//2), r'$%d\degree N$' %ymax])
    # ax.tick_params(axis=u'both', which=u'both',length=0)
    # ax.set_title('Dardar')
    # fig.subplots_adjust(right=0.89)
    # cbar_ax = fig.add_axes([0.90, 0.13, 0.01, 0.30])
    # cbar = fig.colorbar(im2, cax=cbar_ax)
    # figname = '%s/2dhist_olds_cld_thin' %plotDir
    # fig.savefig(figname + '.png')
    
    
    
    
    
    
    
    
    
    
    
    
    sys.exit()    
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
    