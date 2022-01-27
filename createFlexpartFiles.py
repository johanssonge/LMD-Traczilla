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
import os
import sys
import datetime
from dateutil.relativedelta import relativedelta  # @UnresolvedImport
import shutil

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cont", action='store_true', default = False, 
                        help = "Continue a run on an already started but broken run. Default = False")
    parser.add_argument("-d", "--use_dardar", action='store_true', default=False, 
                        help = "Use DARDAR data")
    parser.add_argument("-y", "--year", type=int, default=2018,  
                        help = "year. Default=2018")
    parser.add_argument("-m", "--month", type=int, choices=np.arange(1, 13), default=6, 
                        help = "Month. Default=6")
    parser.add_argument("-n", "--night", type=str, choices=["d", "n", "a"], default='a', 
                        help = "pixlar with day (d), night (n) or all (a). Default = a")
    parser.add_argument("-q", "--qsub", action='store_true', default = False, 
                        help = "Send direktly to qsub. Default = False")
    
    args = parser.parse_args()
    year = args.year
    mon = args.month
    useDardar = args.use_dardar
    #: Day / night / 24h
    dnf = args.night
    
    startDate = datetime.datetime(year, mon, 1, 0, 0) + relativedelta(months=1)
    endDate = startDate - relativedelta(days=75)

    # year_p1 = startDate.year
    # mon_p1 = startDate.month
    # year_m2 = endDate.year
    # mon_m2 = endDate.month
    
    #: Dirs
    runDir = '/home/ejohansson/Projects/flexpart/work/STC/Calipso/%d' %year
    optDir = '/home/ejohansson/Projects/flexpart/traczilla/options/STC/Calipso'
    dataDir = '/data/ejohansson/flexout/STC/Calipso'
    ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/STC/Calipso'
    if not os.path.isdir(runDir):
        os.makedirs(runDir)
    
    #: Do file to start tracilla with qsub
    dofn = '%s/do-CALIOP-EAD-%d%02d-%s' %(runDir, year, mon, dnf)
    #: File with paths used by tracilla
    pathfn = '%s/path-CALIOP-EAD-%d%02d-%s'  %(runDir, year, mon, dnf)
    #: Result dir
    resultDir = '%s/CALIOP-EAD-%d%02d-%s' %(ekjDir, year, mon, dnf)
    #: Path for RELEASES and COMMAND files
    optPath = '%s/CALIOP-EAD-%d%02d-%s' %(optDir, year, mon, dnf)
    #: The RELEASES file. This one is same for everyone
    src_Relise = '%s/RELEASES' %optDir
    #: File with the printstatment traxilla is creating during run
    trazilla_outfn = '%s/traczilla-CALIOP-EAD-%d%02d-%s' %(runDir, year, mon, dnf)
    #: Part file
    partName = 'part_000-%d%02d-%s' %(year, mon, dnf)
    #: Param adn Catalog start of filename
    iniStartName = 'selCaliop'
    
    #: Qsub
    #:Outfile for qsub
    # qsub_outfile = '%s/O-${PBS_JOBNAME}-%s' %(runDir, dnf)
    qsub_outfile = '%s/O-%d%02d-%s' %(runDir, year, mon, dnf)
    #: Jobname
    jobname = '%d%02d%s' %(year, mon, dnf)
    
    #: If dardar is used as input, add DD to keep seperated
    if useDardar:
        dofn = dofn + '-DD'
        pathfn = pathfn + '-DD'
        qsub_outfile = qsub_outfile + 'DD'
        trazilla_outfn = trazilla_outfn + '-DD'
        partName = partName + '-DD'
        resultDir = resultDir + '-DD'
        optPath = optPath + '-DD'
        jobname = jobname + 'DD'
        iniStartName = 'selDardar'
    
    #: Dir with initiation files
    iniDir = '%s/Initfiles' %resultDir
    #: Initiation files
    paramName = partName.replace('part_000-', '%s_Params-' %iniStartName).replace('-DD', '') + '.pkl'
    catalogName = partName.replace('part_000-', '%s_Calalog-' %iniStartName).replace('-DD', '') + '.pkl'
    partFile = '%s/CALIOP-EAD/Part/%s' %(dataDir, partName)
    paramFile = partFile.replace('/Part/', '/Param/').replace(partName, paramName)
    catalogFile = partFile.replace('/Part/', '/Catalog/').replace(partName, catalogName) 

    #: Links to files used by tracilla
    #: startfile i.e. file from selCaliop
    outlink = '%s/part_000' %resultDir
    #: Link to RELEASES file.
    link_Relise = '%s/RELEASES' %optPath
    #: No link. File with commands for tracilla. Is created here
    commandfn = '%s/COMMAND' %optPath

    #: Restartung an old run
    if args.cont:
        #: Make sure the old run exist
        if not os.path.isfile(resultDir + '/savpos'):
            print('No savpos file from tracilla')
            print('Argument (-c --cont) cant be used')
            sys.exit()
        #: 0 means restart
        restart = 0
        #: Move old output file just in case
        newname = trazilla_outfn
        i = 0
        check_exist = True
        while check_exist:
            i = i + 1
            if os.path.isfile(newname):
                newname = '%s-%d' %(trazilla_outfn, i)
            else:
                os.rename(trazilla_outfn, newname)
                check_exist = False
    else:
        #: -1 means new run
        restart = -1
        
        #: Create and move some Dir/Files/Linkes
        if not os.path.isdir(resultDir):
            os.makedirs(resultDir)
        if not os.path.isdir(iniDir):
            os.makedirs(iniDir)
        if not os.path.isfile(iniDir + '/' + partName):
            shutil.move(partFile, iniDir + '/' + partName)
            shutil.move(paramFile, iniDir + '/' + paramName)
            shutil.move(catalogFile, iniDir + '/' + catalogName)
        
        if not os.path.islink(outlink):
            #: Only controls if link exist not if broken
            os.symlink('Initfiles/' + partName, outlink)
        if not os.path.exists(outlink):
            #: Both controls if link exist and if broken
            print('broken link')
            print('src = %s' %(iniDir + '/' + partName))
            print('link = %s' %outlink)
            sys.exit()
            
        if not os.path.isdir(optPath):
            os.makedirs(optPath)
        if not os.path.islink(link_Relise):
            os.symlink(src_Relise, link_Relise)
        if not os.path.exists(link_Relise):
            print('broken link')
            print('src = %s' %src_Relise)
            print('link = %s' %link_Relise)
            sys.exit()
        
        #: Create do-file
        print(dofn)
        dof = open(dofn, 'w')
        dof.writelines('#!/bin/bash \n')
        dof.writelines('#PBS -q day \n')
        dof.writelines('#PBS -l nodes=1:ppn=16 -l mem=5gb -l vmem=5gb \n')
        dof.writelines('#PBS -N %s \n' %(jobname))
        dof.writelines('#PBS -o %s.out \n' %qsub_outfile)
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
        availDir = '/data/ejohansson/ERA5/indexes'
        if not (os.path.isfile('%s/AVAILABLE-%d-uvwt' %(availDir, year)) and os.path.isfile('%s/AVAILABLE-%d-hr' %(availDir, year))):
            print('ERA5-indexes are missing. Use mkAVAILABLE.py to create')
            print('%s/AVAILABLE-%d-uvwt \n' %(availDir, year))
            print('%s/AVAILABLE-%d-hr \n' %(availDir, year))
            sys.exit()
        print(pathfn)
        pathf = open(pathfn, 'w')
    #     pathf.writelines('/home/ejohansson/Projects/flexpart/traczilla/options/STC/Calipso/CALIOP-EAD-Jul2008/ \n')
        pathf.writelines('%s/ \n' %optPath)
        pathf.writelines('%s/ \n' %resultDir)
        pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
        pathf.writelines('%s/AVAILABLE-%d-uvwt \n' %(availDir, year))
        pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
        pathf.writelines('%s/AVAILABLE-%d-hr \n' %(availDir, year))
        pathf.close()
        

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
    commandf.writelines(" ibdate=%d%02d%02d ibtime=000000\n" %(endDate.year, endDate.month, endDate.day))
    commandf.writelines(" iedate=%d%02d%02d ietime=000000\n" %(startDate.year, startDate.month, startDate.day))
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
    commandf.writelines(" hrstart=%d\n" %restart)
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

    print(resultDir)
    
    if args.qsub:
        cmd = "qsub %s" %dofn
        print("Submitt qsub")
        print(cmd)
        c1 = os.system(cmd)
    
    pdb.set_trace()