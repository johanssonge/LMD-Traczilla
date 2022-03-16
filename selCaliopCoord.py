#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2022-03-16

Copyright (c) 2022 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""

import numpy as np
import pdb
import datetime



def getTempPressure(H):
    #: the universal gas constant, N m kmol⁻¹ K⁻¹
    R = 8.31432*(10**3)
    #: the gravitational acceleration, m/s²
    g0 = 9.80665
    #: the molar mass of Earth’s air, kg/mol
    M = 0.0289644
    if H < 0:
        Lb = 0
        Hb = 0
        Tb = 288.15
        Pb = 101325
    elif H < 11:
        Lb = -6.5
        Hb = 0
        Tb = 288.15
        Pb = 101325
    elif H < 20:
        Lb = 0
        Hb = 11
        Tb = 216.65
        Pb = 22632.06
    elif H < 32:
        Lb = 1
        Hb = 20
        Tb = 216.65
        Pb = 5474.889
    else:
        print('Add more height layers')
        print('https://www.translatorscafe.com/unit-converter/en-US/calculator/altitude/#international-standard-atmosphere-ISA')
    Tm = Tb + Lb * (H - Hb)
    if Lb == 0:
        P = Pb * np.exp((-1 * g0 * M * (H - Hb)) / (R * Tb))
    else:
        P = Pb * ((Tb / Tm)**((g0 * M) / (R * Lb) ))
    return(Tm, P)
    

def readTextFile(fn):
    cf = open(fn, 'r')
    cfl = cf.readlines()
    cf.close()
    datearr = []
    latarr = []
    lonarr = []
    heightarr = []
    otarr = []
    temparr = []
    presurearr = []
    for cdl in cfl:
        if (cdl[0] == '#') or (cdl == '\n'):
            continue
        datestr, latstr, lonstr, heightstr, otstr = cdl.split(' ')
        datearr.append(datetime.datetime.strptime(datestr, "%Y%m%d_%H%M%S"))
        latarr.append(float(latstr))
        lonf = float(lonstr)
        if lonf > 180:
            lonf = lonf - 360
        lonarr.append(lonf)
        heightf = float(heightstr)
        heightarr.append(heightf)
        otarr.append(otstr)
        Tm, P = getTempPressure(heightf)
        temparr.append(Tm)
        presurearr.append(P)
    return datearr, latarr, lonarr, heightarr, otarr, temparr, presurearr

if __name__ == '__main__':
    coordfile = 'coord_high_thin_layers_17.5.txt'
    datum, lats, lons, heights, ots, temps, pressures = readTextFile(coordfile)
    
    print('hej')
    pdb.set_trace()
    
    
