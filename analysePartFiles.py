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


from analyseTraczilla import getConvfiles, I_HIT, I_OLD, readCatalogFile
from convsrcErikFullGridSatTropo import satratio

sys.path.append(os.environ['HOME'] + '/Projects/STC/pylib')
from io107 import readpart107, readidx107


def findMin(retv, p0IF, ppIF, use_inds, ppVF, ft):
    pinds = ppIF
    temp_full = np.zeros(p0IF.shape[0]) + np.nan
    temp_full[pinds] = ppVF
    p_min = np.nanmin((retv[0,:], temp_full[use_inds]), axis=0)
    # n_nan_inds = (~np.isnan(p_min))
    # qs_min[1,:][n_nan_inds] = np.where(np.nanargmin((qs_min[0,:][n_nan_inds], temp_full[hits][n_nan_inds]), axis=0) == 1, phour, qs_min[1,:][n_nan_inds])
    retv[1,:] = np.where(np.nanargmin((retv[0,:], temp_full[use_inds]), axis=0) == 1, ft, retv[1,:])
    retv[0,:] = p_min
    return retv
if __name__ == '__main__':
    
    outDir = '/data/ejohansson/flexout/STC/Calipso-OUT'
    years = [2008]
    months = [1]
    ticT = time.time()
    for year in years:
        for mon in months:
            outname = 'CALIOP-EAD-%d%02d-n-DD' %(year, mon)
            partDir = os.path.join('/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/STC/Calipso', outname)
            # te = readCatalogFile(partDir + '/Initfiles/selDardar_Calalog-200701-n.pkl')
            #: Read out file
            rvs, fls, age, lons, lats, temp, pres = getConvfiles(outDir, [outname])
            hits = (fls & I_HIT) == I_HIT
            all_inds = np.ones(hits.shape).astype(bool)
            #: Read part_000 file
            f0 = os.path.join(partDir, 'part_000')
            part0 = readidx107(f0, quiet=True)
            #: Calc saturation
            rvs0 = satratio(part0['p'], part0['t'])
            #: get the desired variables and put in a 2D array
            pMin = np.zeros([2, hits.sum()])
            pMin[0, :] = part0['p'][hits]
            rvsMin = np.zeros([2, hits.shape[0]])
            rvsMin[0, :] = rvs0
            #hmax = 18
            step = 3
            hmax = 1800#3240 #1800 #: 75 days (75 * 24 = 1800)
            dstep = datetime.timedelta(hours=step)
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
                rvsf = satratio(partf['p'], partf['t'])
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
                rvsMin = findMin(rvsMin, part0['idx_back']-1, partf['idx_back']-1, all_inds, rvsf, phour)
                # pMin = findMin(pMin, part0['idx_back']-1, partf['idx_back']-1, hits, partf['p'], phour)
                print("Calc time = %.02f s" %(time.time() - tic))
                #: Check memory
                pid = os.getpid()
                py = psutil.Process(pid)
                memoryUse = py.memory_info()[0]/2**30
                print('memory use: {:4.2f} gb'.format(memoryUse))
                #: Match Indices
                # mi = np.asarray(list((set(partf['idx_back']) & set(part0['idx_back'][hits]))))
                # inds = np.searchsorted(part0['idx_back'][hits], mi)
                # pdb.set_trace()
            # f2r = f2.replace('.gz', '')
    print("total time = %.02f s" %(time.time() - ticT))
    pdb.set_trace()
