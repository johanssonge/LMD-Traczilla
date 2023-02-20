#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2023-02-10

Copyright (c) 2023 Erik Johansson

@author:     Erik Johansson
@contact:    <erik.johansson@smhi.se>
 
'''

import numpy as np
import h5py
import pdb

from matplotlib import pyplot as plt  # @UnresolvedImport
import cartopy.crs as ccrs  # @UnresolvedImport

from analyseTraczilla import readRedmaskDatetime, readTempFileInit
import datetime
from selCaliop import readDARDAR
import os
import glob
if __name__ == '__main__':
    dirDardar = '/home/ejohansson/Data/Satellite/Cloudsat/DARDAR-MASK.v2.23'
    dirSat = dirDardar
    
    # Latitude band
    latmin = -30
    latmax = 45
    
    ddf = ''
    tempf = '/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/STC/Calipso/TempFiles/CALIOP-EAD-200808-n-DD-init.h5'
    svcf = ''
    
    y = 2008
    m = 8
    d = 31
    dt = datetime.datetime(y,m,d)

    lons0_mon, lats0_mon, p0_mon, t0_mon, sh_mon, vod_mon, height_mon, cm_mon, sc_mon, rvs0_mon, iwc0_mon, irs_mon, svc_mon = readTempFileInit(tempf)
    
    d_lon = lons0_mon.reshape(-1, 34)[:, 0]
    d_lat = lats0_mon.reshape(-1, 34)[:, 0]
    

    svc_tdiff = svc_mon['svc_tdiff']
    svc_lon = svc_mon['svc_lon']
    svc_lat = svc_mon['svc_lat']
    ti = np.where(svc_tdiff == 0)[0]

    redmask = readRedmaskDatetime(dt)

    dirday = os.path.join(dirSat, dt.strftime('%Y/%Y_%m_%d'))
    fic = sorted(glob.glob(dirday+'/DARDAR-MASK*.nc'))
    for filename in fic:
        altx, lats, lons, pres, temp, utc, extras = readDARDAR(filename, 'n')
        sel = np.where((lats>latmin) & (lats<latmax))[0][0:-1:10]
        lats = lats[sel]
        lons180 = np.where(lons > 180, lons-180, lons)[sel]
    
        np.where((lats == d_lat[0]) & (lons180 == d_lon[0]))
        lats[0:100] == d_lat[0:100]
        lons180[0:100] == d_lon[0:100]
        pdb.set_trace()
    
    fig = plt.figure()
    fig.suptitle('2008-07-03, Time 01:02 - 01:12')
    ax = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
    ax.set_extent([0, 30, -30, 30], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.plot(lons0_mon[org_cut_inds_file], lats0_mon[org_cut_inds_file], 'r+', label='DARDAR')
    ax.plot(np.asarray(redmask['longitude'])[file_cut_inds], np.asarray(redmask['latitude'])[file_cut_inds], 'b+', label='SVC-MASK')
    ax.set_xticks([0, 15, 30], crs=ccrs.PlateCarree())
    ax.set_xticklabels(['0', '15', '30'], fontsize='large')
    ax.set_yticks([-30, 0, 30], crs=ccrs.PlateCarree())
    ax.set_yticklabels(['-30', '0', '30'], fontsize='large')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.legend()
    
    ax = fig.add_subplot(1,2,2,projection=ccrs.PlateCarree())
    ax.set_extent([5, 15, -5, 20], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.plot(lons0_mon[org_cut_inds_file][1:][::18], lats0_mon[org_cut_inds_file][1:][::18], 'r+', label='DARDAR')
    ax.plot(np.asarray(redmask['longitude'])[file_cut_inds][1:][::18], np.asarray(redmask['latitude'])[file_cut_inds][1:][::18], 'b+', label='SVC-MASK')
    ax.set_xticks([5, 10, 15], crs=ccrs.PlateCarree())
    ax.set_xticklabels(['5', '10', '15'], fontsize='large')
    ax.set_yticks([0, 10, 20], crs=ccrs.PlateCarree())
    ax.set_yticklabels(['0', '10', '20'], fontsize='large')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.legend()
    fig.savefig('svc_mask' + '.png')
    
    
    
    