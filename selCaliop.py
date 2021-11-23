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
import sys

def readDARDAR(fname, dn):
    # open file
    ncf = netCDF4.Dataset(fname,'r')
    
    if 'DARDAR-MASK' in fname:
        altx = ncf.variables['CS_TRACK_Height'][:].data
        # Reads latitudes and longitudes
        lats = ncf.variables['CLOUDSAT_Latitude'][:].data
        lons = ncf.variables['CLOUDSAT_Longitude'][:].data % 360
        # Read pressure (Pa), temperature (degree K)
        pres = ncf.variables['Pressure'][:].data
        temp = ncf.variables['Temperature'][:].data
        # Read time
        tai = ncf.variables['CLOUDSAT_TAI_Time'][:].data
        tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
        utc = tt.utc.datetime

        #: Extras
        #: Cloud Mask
        clm = ncf.variables['CLOUDSAT_2B_GEOPROF_CPR_Cloud_Mask'][:].data
        clv = ncf.variables['CALIOP_Land_Water_Mask'][:].data
        cst = ncf.variables['CALIOP_IGBP_Surface_Type'][:].data
        th = ncf.variables['Tropopause_Height'][:].data
        sp = ncf.variables['Surface_pressure'][:].data
        skt = ncf.variables['Skin_temperature'][:].data
        t2 = ncf.variables['Temperature_2m'][:].data
        sst = ncf.variables['Sea_surface_temperature'][:].data
        sh = ncf.variables['Specific_humidity'][:].data
        msz = ncf.variables['MODIS_Solar_zenith'][:].data
        dsc = ncf.variables['DARMASK_Simplified_Categorization'][:].data
        
        extras = {'Cloud_Mask': clm, 'Land_Water_Mask': clv, 'Surface_Type': cst, \
                  'Tropopause_Height': th, 'Surface_pressure': sp, 'Skin_temperature': skt, \
                  'Temperature_2m': t2, 'SST': sst, 'Specific_humidity': sh, 'SZA': msz, 'Simplified_Categorization': dsc}
        
        #: Day/night
        dn_var = ncf.variables['CALIOP_Day_Night_Flag'][:].data
        #: 1 = night
        #: 0 = day
        if dn == 'n':
            dnf = dn_var.astype(bool)
        elif dn == 'd':
            dnf = ~(dn_var.astype(bool))
        else:
            dnf = np.ones(dn_var.shape).astype(bool)
    else:
        altx = ncf.variables['height'][:].data
        print('check for unit. Needs to be km')
        pdb.set_trace()
        lats = ncf.variables['latitude'][:].data
        lons = ncf.variables['longitude'][:].data % 360
        temp = ncf.variables['temperature'][:].data
        print('No pressure')
        pdb.set_trace()
        sy, sm, sd, st = ncf.variables['time'].units.split(' ')[2:-1]
        # sh, smin, ss = st.split(':')
        # datetime(int(sy), int(sm), int(sd), int(sh), int(smin), int(ss))
        # np.datetime64('%s-%s-%s %s' %(sy, sm, sd, st))
        utc = np.datetime64('%s-%s-%s %s' %(sy, sm, sd, st)) + ncf.variables['time'][:].data.astype('timedelta64[s]')
        utc = np.asarray(utc.tolist())

    #: Cut day/nigh/all
    #: 1D
    lats = lats[dnf]
    lons = lons[dnf]
    utc = utc[dnf]
    #: 2D
    pres = pres[dnf, :]
    temp = temp[dnf, :]
    
    #: Extras
    for arname, val in extras.items():
        if val.ndim == 1:
            val = val[dnf]
        elif val.ndim == 2:
            val = val[dnf, :]
        else:
            print('Wrong ndim')
            sys.exit()
        extras[arname] = val


    ncf.close()
    return altx, lats, lons, pres, temp, utc, extras
    
    
            
def readCalipso(fname):
    # open file
    try:
        if fname.split('.')[-1] == 'hdf':
            hdf = SD(fname)
            hh = HDF.HDF(fname,HDF.HC.READ)
        else:
            h5 = h5py.File(fname, 'r')  # @UnusedVariable
            print('not prepared for h5')
            pdb.set_trace()
    except :
        print('%s Error -> giving up' %fname.split('.')[-1].upper())
        return -1, -1, -1, -1, -1, -1

    # Check altitude
    meta = hh.vstart().attach('metadata')
    altx = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    if np.max((altx - altx_ref)**2)>1.e-5:
        print('ACHTUNG ALARM! NON CONFORM ALTITUDE')
        return -1, -1, -1, -1, -1, -1
        
    # Reads latitudes and longitudes
    lats = hdf.select('Latitude').get()[:,1]
    lons = hdf.select('Longitude').get()[:,1] % 360
    
    # Read pressure (hPa), temperature (degree C)
    # Conversion to Pa and K
    pres = hdf.select('Pressure').get()[:]
    temp = hdf.select('Temperature').get()[:]

    pres *= 100
    temp += 273.15
        
    # Read time
    tai = hdf.select('Profile_Time').get()[:,1]
    tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
    utc = tt.utc.datetime
        
    return altx, lats, lons, pres, temp, utc
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--use_dardar", action='store_true', default=False, 
                        help = "Use DARDAR data")
    parser.add_argument("-y", "--year", type = int, default = 2018,  
                        help = "year. Default = 2018")
    parser.add_argument("-m", "--month", type = int, choices=np.arange(1, 13), default = 6, 
                        help = "Month. Default = 6")
    parser.add_argument("-n", "--night", type=str, choices=["d", "n", "a"], default = 'a', 
                        help = "pixlar with day (d), night (n) or all (a). Default = a")
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
    #: day (d) / Night (n) / All (a)
    dn = args.night
    if (dn == 'n'):
        dntype = 'night'
    elif (dn == 'd'):
        dntype = 'day'
    else:
        dntype = '24h'
    
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
    if useDardar:
        dirSat = dirDardar
    else:
        dirSat = dirAProf
    
    # Latitude band
    latmin = -30
    latmax = 45
    
    #: Use every second datapoint. 
    #: This mean 10km Horizontal spacing in 5 km units
    #: and 2km Horizontal spacing in 1 km units
    if useDardar:
        ns = 10
    else:
        ns = 2    
    
    catalogdir = './Catalog'
    if not os.path.isdir(catalogdir):
        os.makedirs(catalogdir)
    paramdir = './Param'
    if not os.path.isdir(paramdir):
        os.makedirs(paramdir)
    partdir = './Part'
    if not os.path.isdir(partdir):
        os.makedirs(partdir)
    
    altx_ref = pickle.load(open('alx_ref.pkl','rb'))
    altx_ref = np.array(altx_ref)
    
    # Create catalog
    catalog = {}
    params = {}
    params = {'enddate':endDate,'originDate':originDate,'latmin':latmin,'latmax':latmax,
              'ns':ns,
              'altx':altx_ref,'type': dntype,'interdate':interdate}
    catalog_file = '%s/selCaliop_Calalog-%s-%s.pkl' %(catalogdir, endDate.strftime('%b%Y'), dn)
    params_file =  '%s/selCaliop_Params-%s-%s.pkl' %(paramdir, endDate.strftime('%b%Y'), dn)
    part0_file = '%s/part_000-%s-%s' %(partdir, endDate.strftime('%b%Y'), dn)
    if useDardar:
        catalog_file = catalog_file.replace('/selCaliop_', '/selDardar_')
        params_file = params_file.replace('/selCaliop_', '/selDardar_')
        part0_file = part0_file + '-DD'
    # Generate the dictionary to be used to write part_000
    part0 = {}
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
        dirday = os.path.join(dirSat,date.strftime('%Y/%Y_%m_%d'))
        # List the content of the daily aeorosol directory
        if useDardar:
            # fic = sorted(glob.glob(dirday+'/DARDAR-CLOUD*.nc'))
            fic = sorted(glob.glob(dirday+'/DARDAR-MASK*.nc'))
        else:
            # fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro-*.hdf'))
            fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro-*.h*'))
        print(dirday, len(fic))
        # process all the half-orbits in the aerosol directory
        catalog[date] = {}
        for filename in fic:
            #print(file)
            if useDardar:
                altx, lats, lons, pres, temp, utc, extras = readDARDAR(filename, dn)
            else:
                # skip day/night files
                if (dn == 'n') and ('ZD' in filename):
                    continue
                if (dn == 'd') and ('ZN' in filename):
                    continue
                altx, lats, lons, pres, temp, utc = readCalipso(filename)
            
            if isinstance(altx, int) and (altx == -1):
                continue
            
            #: calculate height levels
            #: first element below 20 km
            l20 = np.where(altx<20)[0][0]
            #: last element above 14 km
            l14 = np.where(altx>14)[0][-1]
            #: The interval is 59.8755 m
            #: Use 180 m height points
            intlev = round(0.180 / (altx[0] - altx[1]))
            #: range of altitudes
            altidx = np.arange(l20,l14+1,intlev)
            nlev = len(altidx)
            
            # Selection of the latitude range and sampling every ns
            sel = np.where((lats>latmin) & (lats<latmax))[0][0:-1:ns]
            if len(sel) == 0:
                continue
            # Selection of the 1D data
            lats = lats[sel]
            lons = lons[sel]
            utc = utc[sel]

            ir_start = np.array([int((utc[i] - originDate).total_seconds()) for i in range(len(utc))])
            # Selection of the 2D data
            pres = pres[sel,:][:,altidx]
            temp = temp[sel,:][:,altidx]
            # Expand the 1D fields
            ir_start = np.repeat(ir_start,nlev).astype(int)
            lats = np.repeat(lats,nlev)
            lons = np.repeat(lons,nlev)
            # Unidimensionalize the 2D fields
            npart = nlev * len(sel)
            pres = np.reshape(pres, npart)
            temp = np.reshape(temp, npart)
            #: Extras
            for arname, val in extras.items():
                if val.ndim == 1:
                    #: Selection of the 1D data
                    val = val[sel]
                    #: Expand the 1D fields
                    val = np.repeat(val, nlev)
                elif val.ndim == 2:
                    #: Selection of the 2D data
                    val = val[sel,:][:,altidx]
                    #: Unidimensionalize the 2D fields
                    val = np.reshape(val, npart)
                extras[arname] = val
            
            # Enrich the catalog
            fname = os.path.basename(filename)
            # extract orbit
            if useDardar:
                #: Really granule
                orbit = fname.split('_')[2]
            else:
                #: Really date and time
                orbit = fname[35:-4]
                
            catalog[date][orbit] = {'longitudes':[lons[0],lons[-1]],
                                    'utc':[utc[0],utc[-1]],'selection':sel,'lensel':len(sel),
                                    'npart':npart, 
                                    'Cloud_Mask': extras['Cloud_Mask'], 'Tropopause_Height': extras['Tropopause_Height'], 
                                    'SZA': extras['SZA'], 'Simplified_Categorization': extras['Simplified_Categorization']}
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
    if numpart == 0:
        print('No points to save ()')
        sys.exit()
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
