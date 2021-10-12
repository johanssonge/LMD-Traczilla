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
import sys


I_DEAD = 0x200000 #: Dead
I_HIT = 0x400000 #: Hit a cloud
I_CROSSED = 0x2000000 #: outside domain
I_DBORNE =  0x1000000 #: Lost after first step
I_OLD = 0x800000 #: Reached end of time without encounter a cloud
I_STOP = I_HIT + I_DEAD


# import flammkuchen as fk  # @UnresolvedImport

if __name__ == '__main__':
    
    outDir = '/data/ejohansson/flexout/STC/Calipso-OUT/Coldpoint'
    
    fname = '%s/CALIOP-EAD-May2019-n.h5' %(outDir)
    # tic=time.time()
    # fkf = fk.load(fname)#,prod0,compression='zlib')
    # age = fkf['src']['age']
    # rvs = fkf['rvs']
    #
    # lons = fkf['src']['x'][:]
    # lats = fkf['src']['y'][:]
    # temp = fkf['src']['t'][:]
    # pres = fkf['src']['p'][:]
    # print(time.time() - tic)
    # tic=time.time()
    h5f = h5py.File(fname, 'r')
    rvs = h5f['data/rvs'][:]
    fls = h5f['data/flag_source'][:]
    age = h5f['data/src']['age'][:]
    lons = h5f['data/src']['x'][:]
    lats = h5f['data/src']['y'][:]
    temp = h5f['data/src']['t'][:]
    pres = h5f['data/src']['p'][:]
    h5f.close()
    # print(time.time() - tic)
    old = (fls & I_OLD) == I_OLD
    hits = (fls & I_HIT) == I_HIT
    
    if lons[hits].max() > 180:
        print('Check lons')
        print('lons max = %f' %lons.max())
        sys.exit()
    else:
        lons = np.where(lons > 180, lons-360, lons)
    
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(age[hits] // 86400)
    # ax.plot(age[hits], pres[hits])
    fig.savefig('test.png')
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    # rotated_pole = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
    # ax.plot(0, 0, 'o', transform=rotated_pole)
    # ax.coastlines()
    # # ax.plot(age[argsort], rvs[argsort])
    #
    # fig.savefig('test.png')
    pdb.set_trace()
    