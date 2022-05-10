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
import sys
import io107
import os
# def create2D(da, lat, lon,  ):
dateWithNans = {'20200111': '20200113_003758', \
                '20220102': '20220101_080814', \
                '20220103': '20220101_080814', \
                '20220104': '20220107_092916', \
                '20220116': '20220115_123549', \
                '20220117': '20220115_123549'}
def getTempPressure(Z):
    #: the universal gas constant, N m kmol⁻¹ K⁻¹
    R = 8.31432#*(10**3)
    #: the gravitational acceleration, m/s²
    g0 = 9.80665
    #: the molar mass of Earth’s air, kg/mol
    M = 0.0289644
    #:Radius of Earth [km]
    #Re = 6356.766
    Re = 6371.0088
    H = (Re * Z) / (Re + Z)
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
        P = Pb * np.exp((-1 * g0 * M * ((H - Hb)*1000)) / (R * Tb))
    else:
        P = Pb * ((Tb / Tm)**((g0 * M) / (R * (Lb / 1000) ) ))
    if Z < 100:
        T = Tm
    else:
        print("Function is not valid for Altitudes > 100 km")
        print("Altitude = %d km" %Z)
        sys.exit()
    return(T, P)
    

def readTextFile(fn):
    cf = open(fn, 'r')
    cfl = cf.readlines()
    cf.close()
    
    datedic = {}
    latdic = {}
    londic = {}
    heightdic = {}
    otdic = {}
    tempdic = {}
    presuredic = {}
    
    for cdl in cfl:
        if (cdl[0] == '#') or (cdl == '\n'):
            continue
        datestr, latstr, lonstr, heightstr, otstr = cdl.split(' ')
        day = datestr.split('_')[0]
        datum = datetime.datetime.strptime(datestr, "%Y%m%d_%H%M%S")
        lonf = float(lonstr)
        latf = float(latstr)
        heightf = float(heightstr)
        T, P = getTempPressure(heightf)
        if day not in datedic.keys():
            datedic.update({day: [datum]})
            latdic.update({day: [latf]})
            londic.update({day: [lonf]})
            heightdic.update({day: [heightf]})
            otdic.update({day: [float(otstr.replace('\n', ''))]})
            tempdic.update({day: [T]})
            presuredic.update({day: [P]})
        else:
            datedic[day].append(datum)
            latdic[day].append(latf)
            londic[day].append(lonf)
            heightdic[day].append(heightf)
            otdic[day].append(float(otstr.replace('\n', '')))
            tempdic[day].append(T)
            presuredic[day].append(P)
        
        # if lonf > 180:
        #     lonf = lonf - 360
    return datedic, latdic, londic, heightdic, otdic, tempdic, presuredic


def readTextFileExtend(fn):
    cf = open(fn, 'r')
    cfl = cf.readlines()
    cf.close()
    
    datedic = {}
    latdic = {}
    londic = {}
    heightdic = {}
    otdic = {}
    tempdic = {}
    presuredic = {}
    
    for cdl in cfl:
        if (cdl[0] == 't') or (cdl[0] == 'f') or (cdl[0] == '#') or (cdl == '\n'):
            continue
        if len(cdl.split(' ')) == 11:
            flight, datestr, latstr, lonstr, heightTopstr, T_topstr, P_topstr, heightBotstr, T_botstr, P_botstr, otstr = cdl.split(' ')  # @UnusedVariable
        elif len(cdl.split(' ')) == 1:
            datestr, latstr, lonstr, heightTopstr, T_topstr, P_topstr, heightBotstr, T_botstr, P_botstr, otstr = cdl.split(' ')
        day = datestr.split('_')[0]
        datum = datetime.datetime.strptime(datestr, "%Y%m%d_%H%M%S")
        lonf = float(lonstr)
        latf = float(latstr)
        h_top = float(heightTopstr)
        h_bot = float(heightBotstr)
        T_top = float(T_topstr)
        T_bot = float(T_botstr)
        P_top = float(P_topstr) * 100
        P_bot = float(P_botstr) * 100
        if day not in datedic.keys():
            datedic.update({day: [datum]})
            latdic.update({day: [latf]})
            londic.update({day: [lonf]})
            heightdic.update({day: [[h_top, h_bot]]})
            tempdic.update({day: [[T_top, T_bot]]})
            presuredic.update({day: [[P_top, P_bot]]})
            otdic.update({day: [float(otstr.replace('\n', ''))]})
        else:
            datedic[day].append(datum)
            latdic[day].append(latf)
            londic[day].append(lonf)
            heightdic[day].append([h_top, h_bot])
            tempdic[day].append([T_top, T_bot])
            presuredic[day].append([P_top, P_bot])
            otdic[day].append(float(otstr.replace('\n', '')))
        
        # if lonf > 180:
        #     lonf = lonf - 360
    return datedic, latdic, londic, heightdic, otdic, tempdic, presuredic


def getTempertureExtended(pres):
    fn = '220425_T_profiles_high_thin_layers.txt'
    tf = open(fn, 'r')
    tfl = tf.readlines()
    tf.close()
    retv = {}
    for tdl in tfl:
        tdls = tdl.split(' ')
        if (tdls[0] == 'time'):# or (tdl[0] == 'f') or (tdl[0] == '#') or (tdl == '\n'):
            overhead = []
            i = -1
            for p in tdls:
                if (p[0] == 't'):
                    continue
                i = i + 1
                pi = int(p.replace('hPa', ''))
                if pi != (pres[i] / 100):
                    print('Wrong pressure')
                    pdb.set_trace()
                else:
                    overhead.append(pi)
        # if len(tdl.split(' ')) == 11:
        else:
            j = -1
            temp = []
            for t in tdls:
                j = j + 1
                if j == 0:
                    tid = datetime.datetime.strptime(t, "%Y%m%d_%H%M%S")
                else:
                    temp.append(float(t))
            retv.update({tid: temp})
    
    #: Check for nans
    for time in retv.keys():
        if time.strftime("%Y%m%d") in dateWithNans.keys():
            retv[time] = retv[datetime.datetime.strptime(dateWithNans[time.strftime("%Y%m%d")], "%Y%m%d_%H%M%S")]
        if np.isnan(retv[time][0]):
            print(time)
            print(time.strftime("%Y%m%d_%H%M%S"))
            pdb.set_trace()
    return retv



def create2D(dats, lats, lons, hs, os, ts, ps):
    pr = [*range(5000, 9500+1, 100)]
    lpr = len(pr)
    if isinstance(hs[0], list) and (len(hs[0]) == 2):
        lpr = lpr + 2
        isExtend = True
        #: Same length as pr
        
        tr = getTempertureExtended(pr)
        # tr = np.repeat(216.65, lpr - 2)
    elif isinstance(hs[0], float):
        lpr = lpr + 1
        isExtend = False
    
    rdate = []
    rlat = []
    rlon = []
    rh = []
    ro = []
    rT = []
    rP = []
    for i in range(len(ps)):
        rdate = np.hstack((rdate, np.repeat(dats[i], lpr)))
        rlat = np.hstack((rlat, np.repeat(lats[i], lpr)))
        rlon = np.hstack((rlon, np.repeat(lons[i], lpr)))
        ro = np.hstack((ro, np.repeat(os[i], lpr)))
        if isExtend:
            if i == 0:
                rh = np.repeat(np.asarray(hs[i]).reshape(-1, 1), lpr, axis=1)
            else:
                rh = np.hstack((rh, np.repeat(np.asarray(hs[i]).reshape(-1, 1), lpr, axis=1)))
            prs = pr.copy()
            trs = tr[dats[i]].copy()
            for j in range(2):
                ind = np.searchsorted(pr, ps[i][j])
                prs = np.insert(prs, ind, ps[i][j])
                trs = np.insert(trs, ind, ts[i][j])
        
        else:
            rh = np.hstack((rh, np.repeat(hs[i], lpr)))
            prs = np.insert(pr, np.searchsorted(pr, ps[i]), ps[i])
            tr = np.repeat(ts[i], lpr)
            
        rT = np.hstack((rT, trs))
        rP = np.hstack((rP, prs))
        
        # np.hstack((rlat, np.repeat(lats[i], lpr + 1)))
    return rdate, rlat, rlon, rh, ro, rT, rP    


if __name__ == '__main__':
    #coordfile = 'coord_high_thin_layers_17.5.txt'
    coordfile = 'coord_high_thin_layers_PT_17.5.txt'
    if '_PT_' in coordfile:
        datums_dict, lats_dict, lons_dict, heights_dict, ots_dict, temps_dict, pressures_dict = readTextFileExtend(coordfile)
    else:
        datums_dict, lats_dict, lons_dict, heights_dict, ots_dict, temps_dict, pressures_dict = readTextFile(coordfile)
    
    partdir = './Part_Coord'
    if not os.path.isdir(partdir):
        os.makedirs(partdir)
    catalog = {}
    for day in datums_dict.keys():
        print(day)
        part0_file = '%s/part_000-Coord_%s' %(partdir, day)
        
        datums = datums_dict[day]
        lats = lats_dict[day]
        lons = lons_dict[day]
        heights = heights_dict[day]
        ots = ots_dict[day]
        temps = temps_dict[day]
        pressures = pressures_dict[day]
        datums2D, lats2D, lons2D, heights2D, ots2D, temps2D, pressures2D = \
                        create2D(datums, lats, lons, heights, ots, temps, pressures)
    
    
        maxdatum = np.max(datums)
        originDate = datetime.datetime(maxdatum.year, maxdatum.month, maxdatum.day) + datetime.timedelta(days=1)
        ir_start = np.array([int((datums[i] - originDate).total_seconds()) for i in range(len(datums))])
        ir_start2D = np.array([int((datums2D[i] - originDate).total_seconds()) for i in range(len(datums2D))])
    
        # catalog[day] = {'longitudes':[lons[0],lons[-1]],
        #                             'utc':[utc[0],utc[-1]],'selection':sel,'lensel':len(sel),
        #                             'npart':npart, 'Track_Height': altx}
        # Generate the dictionary to be used to write part_000
        part0 = {}
        numpart = len(lons2D)
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
        part0['ir_start'] = ir_start2D
        part0['x'] = np.asarray(lons2D)
        part0['y'] = np.asarray(lats2D)
        part0['t'] = np.asarray(temps2D)
        part0['p'] = np.asarray(pressures2D)
        part0['idx_back'] = np.arange(1,numpart+1,dtype=int)
        part0['flag'] = np.full(numpart,127,dtype=int)
        part0['numpart'] = numpart
        part0['nact'] = numpart
        
        io107.writeidx107(part0_file,part0)
    print('Done')
    pdb.set_trace()
    
    
