# -*- coding: utf-8 -*-
"""

This script makes the ERA5 AVAILABLE file for one year or for 
a selection of months during that year

Created on Thu May 21 23:02:19 2020

@author: Bernard Legras
"""
import argparse
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta  # @UnresolvedImport
import numpy as np
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("-y","--year",type=int,help="year")
parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="first month")
parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="last month")
parser.add_argument("-t","--type",type=str,choices=('uvwt','hr'),help="uvwt (default) or hr")

extended = True
year = 2017
m1 = 1
m2 = 12
d1 = 1
typ = 'uvwt'

args = parser.parse_args()
if args.year is not None: year = args.year
if args.month1 is not None: m1 = args.month1
if args.month2 is not None: m2 = args.month2
if args.type is not None: typ = args.type


date = datetime(year,m1,d1,0)
date_stop = date + relativedelta(months=m2-m1+1)
if extended:
    date_stop = date_stop + relativedelta(days=1)
    date_start = date - relativedelta(months=2)
    saveDir = '/data/ejohansson/ERA5/indexes'
else:
    date_start = date
    saveDir = './'

filename = '%s/AVAILABLE-%s-%s' %(saveDir, str(year), typ)
fid = open(filename, 'w')

if typ == 'uvwt':
    # Heading
    fid.write('AVAILABLE '+typ+' file for '+str(year)+' '+str(m1)+'-'+str(m2)+'\n')
    fid.write('ldat  ltime  fname spec (hh in fname)\n')
    fid.write('i8,1x,i6,3x,a16,3s,a10 (hh decoded in read_era5.f90)\n')
    while date_start < date_stop:
        fid.write(date_start.strftime('%Y%m%d %H0000      %Y%m%d '))
        fid.write('{0:4}   ON_DISK\n'.format(100*date_start.hour))
        date_start += timedelta(hours=3)
    fid.close()
    
elif typ == 'hr':
    # Heading
    fid.write('AVAILABLE '+typ+' file for '+str(year)+' '+str(m1)+'-'+str(m2)+'\n')
    fid.write('ldat  ltime  fname  spec (time and hh in fname)\n')
    fid.write('i8,1x,i6,3x,a16,3x,a10 (time and hh decoded in read_era5_diab)\n')
    date2 = date_start - timedelta(hours=6)
    date_start += timedelta(minutes=30)
    time = 1800
    step = 7
    while date_start < date_stop:
        fid.write(date_start.strftime('%Y%m%d %H%M00   '))
        fid.write(date2.strftime('%Y%m%d '))
        fid.write('{0:4} '.format(time))
        fid.write('{0:2}   ON_DISK\n'.format(step))
        date_start += timedelta(hours=1)
        date2 += timedelta(hours=1)
        step += 1
        if step ==13:
            step = 1
            time = (time +1200) %2400
    fid.close()
        
print(filename)


