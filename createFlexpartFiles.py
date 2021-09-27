#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

Created on 2021-09-20

Copyright (c) 2021 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

'''

import numpy as np  # @UnusedImport
import pdb
import calendar
import os
import sys
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--use_dardar",action='store_true',default=False,help="Use DARDAR data")
    parser.add_argument("-y","--year", type=int, default=2018,  
                        help="year. Default=2018")
    parser.add_argument("-m","--month", type=int, default=6, 
                        help="Month. Default=6")
    parser.add_argument("-n","--night",type=str, default='a', 
                        help="pixlar with day (d), night (n) or all (a). Default = a")
   
    args = parser.parse_args()
    year = args.year
    mon = args.month
    useDardar = args.use_dardar
    #: Day / night / 24h
    dnf = args.night
    
    if mon == 12:
        year_p1 = year + 1
        mon_p1 = 1
    else:
        year_p1 = year
        mon_p1 = mon + 1

    if mon in [1,2]:
        year_m2 = year - 1
        if mon == 1:
            mon_m2 = 11
        elif min == 2:
            mon_m2 = 12
    else:
        year_m2 = year
        mon_m2 = mon - 2

    abr = calendar.month_abbr[mon]
    
    runDir = '/home/ejohansson/Projects/flexpart/work/STC/Calipso'
    optDir = '/home/ejohansson/Projects/flexpart/traczilla/options/STC/Calipso'
    outDir = '/data/ejohansson/flexout/STC/Calipso'
    
    #: Do file to start tracilla with qsub
    dofn = '%s/do-CALIOP-EAD-%d-%02d-%s' %(runDir, year, mon, dnf)
    #: File with paths used by tracilla
    pathfn = '%s/path-CALIOP-EAD-%d-%02d-%s'  %(runDir, year, mon, dnf)

    #: Result dir
    outpath = '%s/CALIOP-EAD-%s%d-%s' %(outDir, abr, year, dnf)
    #: part link from selCaliop. 
    #: The program will create a link (outlink) to this one in the outpath
    srclink = '%s/CALIOP-EAD/part_000-%s%d-%s' %(outDir, abr, year, dnf)
    
    #: Path for RELEASES and COMMAND files
    optPath = '%s/CALIOP-EAD-%s%d-%s' %(optDir, abr, year, dnf)
    #: Yhe RELEASES file. This one is same for everyone
    src_Relise = '%s/RELEASES' %optDir
    
    #: File with the printstatment traxilla is creating during run
    trazilla_outfn = '%s/traczilla-CALIOP-EAD-%d-%02d-%s' %(runDir, year, mon, dnf)
    #: Jobanme for qsub
    jobbname = '%s/O-${PBS_JOBNAME}-%s' %(runDir, dnf)
    
    #: If dardar is used as input, add DD to keep seperated
    if useDardar:
        dofn = dofn + '-DD'
        pathfn = pathfn + '-DD'
        jobbname = jobbname + 'DD'
        trazilla_outfn = trazilla_outfn + '-DD'
        srclink = srclink + '-DD'
        outpath = outpath + '-DD'
        optPath = optPath + '-DD'
        
    #: Links to files used by tracilla
    #: startfile i.e. file from selCaliop
    outlink = '%s/part_000' %outpath
    #: Link to RELEASES file.
    link_Relise = '%s/RELEASES' %optPath
    #: No link. File with commands for tracilla. Is created here
    commandfn = '%s/COMMAND' %optPath
    pdb.set_trace()
    print(dofn)
    dof = open(dofn, 'w')
    dof.writelines('#!/bin/bash \n')
    dof.writelines('#PBS -q day \n')
    dof.writelines('#PBS -l nodes=1:ppn=16 -l mem=5gb -l vmem=5gb \n')
    dof.writelines('#PBS -N %s%d \n' %(abr, year))
    dof.writelines('#PBS -o %s \n' %jobbname)
    dof.writelines('#PBS -j oe \n')
    dof.writelines('\n')
    dof.writelines('module load pgi/2016 \n')
    dof.writelines('module load netcdf4/4.3.3.1-pgf2016 \n')
    dof.writelines('\n')
    dof.writelines('export OMP_NUM_THREADS=16 \n')
    dof.writelines('export OMP_SCHEDULE="GUIDED,10000" \n')
    dof.writelines('/home/legras/flexpart/new6-devel/TRACZILLA-TT-par-ng %s > %s \n' %(pathfn, trazilla_outfn))
    dof.writelines('\n')
    dof.close()
    
    
    print(pathfn)
    pathf = open(pathfn, 'w')
#     pathf.writelines('/home/ejohansson/Projects/flexpart/traczilla/options/STC/Calipso/CALIOP-EAD-Jul2008/ \n')
    pathf.writelines('%s/ \n' %optPath)
    pathf.writelines('%s/ \n' %outpath)
    pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
    pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/AVAILABLE-%d-uvwt \n' %year)
    pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
    pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/AVAILABLE-%d-hr \n' %year)
    pathf.close()
    
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    if not os.path.exists(outlink):
        os.symlink(srclink, outlink)
    if not os.path.exists(outlink):
        print('broken link')
        print('src = %s' %srclink)
        print('link = %s' %outlink)
        sys.exit()
        
    if not os.path.isdir(optPath):
        os.makedirs(optPath)
    if not os.path.exists(link_Relise):
        os.symlink(src_Relise, link_Relise)
    if not os.path.exists(link_Relise):
        print('broken link')
        print('src = %s' %src_Relise)
        print('link = %s' %link_Relise)
        sys.exit()
    
    print(commandfn)
    commandf = open(commandfn, 'w')
    commandf.writelines("# COMMAND\n")
    commandf.writelines("#*******************************************************************************\n")
    commandf.writelines("#                                                                              *\n")
    commandf.writelines("#      Input file for the Lagrangian particle dispersion model TRACZILLA       *\n")
    commandf.writelines("#                           Please select your options       *\n")
    commandf.writelines("#                                           *\n")
    commandf.writelines("#*******************************************************************************\n")
    commandf.writelines(" &COMMAND\n")
    commandf.writelines(" ibdate=%d%02d15 ibtime=000000\n" %(year_m2, mon_m2))
    commandf.writelines(" iedate=%d%02d01 ietime=000000\n" %(year_p1, mon_p1))
    commandf.writelines(" diffus=0.\n")
    commandf.writelines(" release_plan='StratoClim'\n")
    commandf.writelines(" diftype=1\n")
    commandf.writelines(" loutstep=10800\n")
    commandf.writelines(" loutprint=10800\n")
    commandf.writelines(" lsynctime=450\n")
    commandf.writelines(" loutsav=86400\n")
    commandf.writelines(" loffset=0\n")
    commandf.writelines(" loffset2=-0\n")
    commandf.writelines(" restart=.true.\n")
    commandf.writelines(" hrstart=-1\n")
    commandf.writelines(" ldirect=-1\n")
    commandf.writelines(" isentropic_motion=.false.\n")
    commandf.writelines(" z_motion=.false.\n")
    commandf.writelines(" theta_bounds=.false.\n")
    commandf.writelines(" vert_interpol='log'\n")
    commandf.writelines(" diabatic_w=.true.\n")
    commandf.writelines(" ecmwf_diabatic_w=.false.\n")
    commandf.writelines(" era5_data=.true.\n")
    commandf.writelines(" era5_diab=.true.\n")
    commandf.writelines(" cloud_sky=.true.\n")
    commandf.writelines(" debug_out=.false.\n")
    commandf.writelines(" hour_accu=1\n")
    commandf.writelines(" lower_theta_level=18\n")
    commandf.writelines(" upper_theta_level=118\n")
    commandf.writelines(" &END\n")
    commandf.writelines("\n")
    commandf.close()

    
    
    
    pdb.set_trace()