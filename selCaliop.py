#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exploring the CALIOP database to select orbits and initialize parcels for
backward trajectory calculations on a regular grid.

This is meant to run on ICARE.

Generates part_000 initialization file for a year and a month and two auxiliary files
containing the list of parameters and the catalog of retained orbits.

In the initial version,
- night orbits are retained
- orbits are retained every 3 days in the month (interdate=3)
- levels from 20 to 14 kms are retained with a sampling of ~ 200 m (inlev=3)
- On the horizontal the resolution is about 10 km (ns=2) between 30S (latmin) and
45N (latmax)
 
@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta  # @UnresolvedImport
import pickle,gzip
#import matplotlib.pyplot as plt
import os
import netCDF4  # @UnresolvedImport
import h5py  # @UnresolvedImport
from pyhdf.SD import SD  # @UnresolvedImport
from pyhdf import HDF, VS, V  # @UnresolvedImport @UnusedImport
import glob
from astropy.time import Time, TimeDelta  # @UnresolvedImport
import argparse
import io107
import pdb


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--use_dardar",action='store_true',default=False,help="Use DARDAR data")
    parser.add_argument("-y","--year",type=int,default=2017,help="year")
    parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),default=8,help="month")
    #parser.add_argument("-d","--day",type=int,choices=1+np.arange(31),help="day0")
    
    #: Default parameters
    #: End day of the selection
    # year = 2017
    # month = 8
    day = 1
    
    args = parser.parse_args()
    year = args.year
    month = args.month
    useDardar = args.use_dardar
    # if args.year is not None: year = args.year
    # if args.month is not None: month = args.month
    # Define dates and date interval
    endDate = datetime(year, month, 1, 0)
    originDate = endDate + relativedelta(months=1)
    interdate = 3
    
    # Main irectories for aerosol profiles and L1 data
    if (year > 2020) or ((year == 2020) and (month >= 7)):
        dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v4.21'
    else:
        dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v4.20'
    # dirDardar = '/DATA/LIENS/CLOUDSAT/DARDAR-CLOUD.v3.10'
    dirDardar = '/DATA/LIENS/CLOUDSAT/DARDAR-MASK.v2.23/'
    #dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
    
    
    # Latitude band
    latmin = -30
    latmax = 45
    
    # Horizontal spacing in 5 km units
    ns = 2 # that is 10 km
    
    # # first element below 20 km
    # l20 = 57
    # # last element above 14 km
    # l14 = 156
    # # The interval is 59.8755 m
    # intlev = 3
    # # range of altitudes every 3 points (about 180 m)
    # altidx = np.arange(l20,l14+1,intlev)
    # nlev = len(altidx)
    
    altx_ref = pickle.load(open('alx_ref.pkl','rb'))
    altx_ref = np.array(altx_ref)
    
    # Create catalog
    catalog = {}
    catalog_file = 'selCaliop_Calalog'+endDate.strftime('-%b%Y.pkl')
    params = {}
    params = {'enddate':endDate,'originDate':originDate,'latmin':latmin,'latmax':latmax,
              'ns':ns,
              'altx':altx_ref,'type':'night','interdate':interdate}
    params_file =  'selCaliop_Params'+endDate.strftime('-%b%Y.pkl')
    
    # Generate the dictionary to be used to write part_000
    part0 = {}
    part0_file = endDate.strftime('part_000-%b%Y')
    # Heading data
    part0['lhead'] = 3
    part0['outnfmt'] = 107
    part0['mode'] = 3   # modify that
    part0['stamp_date'] = originDate.year*10**10 + originDate.month*10**8 + \
        originDate.day*10**6 + originDate.hour*10**4 + originDate.minute*100
    part0['itime'] = 0
    part0['step'] = 450
    part0['idx_orgn'] = 1
    part0['nact_lastO'] = 0
    part0['nact_lastNM'] = 0
    part0['nact_lastNH'] = 0
    part0['flag'] = np.empty(0,dtype=int)
    part0['ir_start'] = np.empty(0,dtype=int)
    part0['x'] = np.empty(0,dtype=float)
    part0['y'] = np.empty(0,dtype=float)
    part0['t'] = np.empty(0,dtype=float)
    part0['p'] = np.empty(0,dtype=float)
    part0['idx_back'] = np.empty(0,dtype=int)
    numpart = 0
    
    # Browse dates
    # Starts from the last day of the month
    date = originDate - timedelta(days=1)
    
    while date >= endDate:
        # Generate names of daily directories
        dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
        dirdayD = os.path.join(dirDardar,date.strftime('%Y/%Y_%m_%d'))
        # List the content of the daily aeorosol directory
        # fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro-*.hdf'))
        fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro-*.h*'))
        # ficD = sorted(glob.glob(dirdayD+'/DARDAR-CLOUD*.nc'))
        ficD = sorted(glob.glob(dirdayD+'/DARDAR-MASK*.nc'))
        print(dirday, len(fic))
        print(dirdayD, len(ficD))
        # process all the half-orbits in the aerosol directory
        catalog[date] = {}
        for i in range(len(fic)):
            # pop file
            file = fic.pop()
            #print(file)
            # skip day files
            if ('ZD' in file): continue
            # open file
            try:
                if file.split('.')[-1] == 'hdf':
                    hdf = SD(file)
                    hh = HDF.HDF(file,HDF.HC.READ)
                else:
                    h5 = h5py.File(file, 'r')
            except :
                print('%s Error -> giving up' %file.split('.')[-1].upper())
                continue
            # Check altitude
            ncf = netCDF4.Dataset(ficD[i],'r')
            meta = hh.vstart().attach('metadata')
            altx = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
            if np.max((altx - altx_ref)**2)>1.e-5:
                print('ACHTUNG ALARM! NON CONFORM ALTITUDE')
                continue
            # altxD = ncf.variables['height'][:].data
            altxD = ncf.variables['CS_TRACK_Height'][:].data
            
            #: calculate height levels
            #: first element below 20 km
            l20 = np.where(altx<20)[0][0]
            l20D = np.where(altxD<20)[0][0]
            #: last element above 14 km
            l14 = np.where(altx>14)[0][-1]
            l14D = np.where(altxD>14)[0][-1]
            #: The interval is 59.8755 m
            #: Use 180 m height points
            intlev = round(0.180 / (altx[0] - altx[1]))
            intlevD = round(0.180 / (altxD[0] - altxD[1]))
            #: range of altitudes
            altidx = np.arange(l20,l14+1,intlev)
            altidxD = np.arange(l20D,l14D+1,intlevD)
            nlev = len(altidx)
            nlevD = len(altidxD)
            
            
            
            # Reads latitudes and longitudes
            lats = hdf.select('Latitude').get()[:,1]
            lons = hdf.select('Longitude').get()[:,1] % 360
            
            # latsD = ncf.variables['latitude'][:].data
            # lonsD = ncf.variables['longitude'][:].data % 360
            latsD = ncf.variables['CLOUDSAT_Latitude'][:].data
            lonsD = ncf.variables['CLOUDSAT_Longitude'][:].data % 360
            
            # Read pressure (hPa), temperature (degree C)
            # Conversion to Pa and K
            pres = hdf.select('Pressure').get()[:]
            temp = hdf.select('Temperature').get()[:]
    
            pres *= 100
            temp += 273.15
            
            # tempD = ncf.variables['temperature'][:].data
            presD = ncf.variables['Pressure'][:].data
            tempD = ncf.variables['Temperature'][:].data
            # Read time
            tai = hdf.select('Profile_Time').get()[:,1]
            tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
            utc = tt.utc.datetime
    
            # sy, sm, sd, st = ncf.variables['time'].units.split(' ')[2:-1]
            # sh, smin, ss = st.split(':')
            # datetime(int(sy), int(sm), int(sd), int(sh), int(smin), int(ss))
            # np.datetime64('%s-%s-%s %s' %(sy, sm, sd, st))
            # utcD = np.datetime64('%s-%s-%s %s' %(sy, sm, sd, st)) + ncf.variables['time'][:].data.astype('timedelta64[s]')
            # utcD = np.asarray(utcD.tolist())
            taiD = ncf.variables['CLOUDSAT_TAI_Time'][:].data
            ttD = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(taiD, format='sec')
            utcD = ttD.utc.datetime
            # Selection of the latitude range and sampling every ns
            sel = np.where((lats>latmin) & (lats<latmax))[0][0:-1:ns]
            if len(sel) == 0: continue
            selD = np.where((latsD>latmin) & (latsD<latmax))[0][0:-1:ns]
            if len(selD) == 0: continue
            # Selection of the 1D data
            lats = lats[sel]
            lons = lons[sel]
            utc = utc[sel]
            latsD = latsD[selD]
            lonsD = lonsD[selD]
            utcD = utcD[selD]

            ir_start = np.array([int((utc[i] - originDate).total_seconds()) for i in range(len(utc))])
            ir_startD = np.array([int((utcD[i] - originDate).total_seconds()) for i in range(len(utcD))])
            # Selection of the 2D data
            pres = pres[sel,:][:,altidx]
            temp = temp[sel,:][:,altidx]
            presD = presD[selD,:][:,altidxD]
            tempD = tempD[selD,:][:,altidxD]
            # Expand the 1D fields
            ir_start = np.repeat(ir_start,nlev).astype(int)
            lats = np.repeat(lats,nlev)
            lons = np.repeat(lons,nlev)
            ir_startD = np.repeat(ir_startD,nlevD).astype(int)
            latsD = np.repeat(latsD,nlevD)
            lonsD = np.repeat(lonsD,nlevD)
            # Unidimensionalize the 2D fields
            npart = nlev * len(sel)
            pres = np.reshape(pres,npart)
            temp = np.reshape(temp,npart)
            npartD = nlevD * len(selD)
            presD = np.reshape(presD,npartD)
            tempD = np.reshape(tempD,npartD)
    
            pdb.set_trace()
            # Enrich the catalog
            fname = os.path.basename(file)
            # extract orbit
            orbit = fname[35:-4]
            catalog[date][orbit] = {'type':'night','longitudes':[lons[0],lons[-1]],
                                    'utc':[utc[0],utc[-1]],'selection':sel,'lensel':len(sel),
                                    'npart':npart}
            print(date,orbit,len(sel),npart)
            # fill part0
            idx1 = numpart
            numpart += npart
            part0['x'] = np.append(part0['x'],lons)
            part0['y'] = np.append(part0['y'],lats)
            part0['t'] = np.append(part0['t'],temp)
            part0['p'] = np.append(part0['p'],pres)
            part0['ir_start'] = np.append(part0['ir_start'],ir_start)
            part0['idx_back'] = np.append(part0['idx_back'],np.arange(idx1+1,numpart+1,dtype=int))
            part0['flag'] = np.append(part0['flag'],np.full(npart,127,dtype=int))
    
        date -= timedelta(days=interdate)
        
    # store the dictionary of traces
    # final size information
    print('End of generation, particles:',numpart)
    part0['numpart'] = numpart
    part0['nact'] = numpart
    params['numpart'] = numpart
    params['lensel'] = len(sel)
    params['npart'] = npart
    params.update({'toplev':l20})
    params.update({'botlev':l14})
    params.update({'intlev':intlev})
    params.update({'nlev':nlev})
    with gzip.open(catalog_file,'wb') as f:
        pickle.dump(catalog,f)
    with open(params_file,'wb') as f:
        pickle.dump(params,f)
    io107.writeidx107(part0_file,part0)