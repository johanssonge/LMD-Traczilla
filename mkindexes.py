#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generates the index for a full month of ECMWF ERA5 data, either for the EN (uvwt) or DI (hr) files 

Created on Thu Jun  4 00:00:09 2020

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import os
from subprocess import call
from os.path import join
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-y","--year",type=int,help="year")
parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="first month")
parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="last month")
parser.add_argument("-t","--typ",type=str,choices=('uvwt','hr','all'),help="all (default), uvwt or hr")

year = 2010
m1 = 1
m2 = 4
typ = "all"

args = parser.parse_args()
if args.year is not None: year = args.year
if args.month1 is not None: m1 = args.month1
if args.month2 is not None: m2 = args.month2
if args.typ is not None: typ = args.typ

os.chdir(os.path.join("/home/legras/ERA5/indexes",str(year)))
for mm in range(m1,m2+1):
    date = datetime(year,mm,1)
    # call(["grib_index_build","-o",date.strftime("uvwt-%Y-%m.gribidx"),
    #       date.strftime("/data/legras/flexpart_in/ERA5/EN-true/%Y/ERA5EN%Y%m*")])
    # call(["grib_index_build","-o",date.strftime("hr-%Y-%m.gribidx"),
    #       date.strftime("/data/legras/flexpart_in/ERA5/DI-true/%Y/ERA5DI%Y%m*")])
    if typ in ["all","uvwt"]:
        call("grib_index_build -o "+date.strftime("uvwt-%Y-%m.gribidx ")+
             date.strftime("/data/legras/flexpart_in/ERA5/EN-true/%Y/ERA5EN%Y%m*"),shell=True)
    if typ in ["all","hr"]:
        call("grib_index_build -o "+date.strftime("hr-%Y-%m.gribidx ")+
             date.strftime("/data/legras/flexpart_in/ERA5/DI-true/%Y/ERA5DI%Y%m*"),shell=True)

