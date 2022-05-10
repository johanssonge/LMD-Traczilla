#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on 2022-04-07

Copyright (c) 2022 Erik Johansson


@author:     Erik Johansson
@contact: <erik.johansson@lmd.ipsl.fr>

"""

import numpy as np  # @UnusedImport
import pdb
import os
import sys
import datetime
from dateutil.relativedelta import relativedelta  # @UnresolvedImport
import shutil
import glob

#: ---------------------------------------------------------------
max_life_time = 15

#: ---------------------------------------------------------------

def rmfile(fn):
    if os.path.islink(fn):
        os.unlink(fn)
    elif os.path.isfile(fn):
        os.remove(fn)

#: ---------------------------------------------------------------

if __name__ == '__main__':
    # import argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-c", "--cont", action='store_true', default = False, 
    #                     help = "Continue a run on an already started but broken run. Default = False")
    # parser.add_argument("-d", "--use_dardar", action='store_true', default=False, 
    #                     help = "Use DARDAR data")
    # parser.add_argument("-y", "--year", type=int, default=2018,  
    #                     help = "year. Default=2018")
    # parser.add_argument("-m", "--month", type=int, choices=np.arange(1, 13), default=6, 
    #                     help = "Month. Default=6")
    # parser.add_argument("-n", "--night", type=str, choices=["d", "n", "a"], default='a', 
    #                     help = "pixlar with day (d), night (n) or all (a). Default = a")
    # parser.add_argument("-q", "--qsub", action='store_true', default = False, 
    #                     help = "Send direktly to qsub. Default = False")
    #
    # args = parser.parse_args()
    # year = args.year
    # mon = args.month
    # useDardar = args.use_dardar
    # #: Day / night / 24h
    # dnf = args.night
    
    
    
    
    #: Dirs
    #: Work dir whith do-, path-, traczilla- and O- files
    runDir = '/home/ejohansson/Projects/flexpart/work/Coord'
    #: Dir for option files i.e. COMMAND  RELEASES
    optDir = '/home/ejohansson/Projects/flexpart/traczilla/options/Coord'
    #: Location for Part files
    dataDir = '/data/ejohansson/flexout/Coord'
    #: Flexout results
    ekjDir = '/proju/flexpart/flexpart_in/EKJ/ejohansson/flexout/Coord'
    
    outputHour = 3
    if outputHour != 3:
        runDir = '%s/%d' %(runDir, outputHour * 3600)
        optDir = '%s/%d' %(optDir, outputHour * 3600)
        ekjDir = '%s/%d' %(ekjDir, outputHour * 3600)
    if not os.path.isdir(runDir):
        os.makedirs(runDir)
    
    
    
    # coordfile = 'coord_high_thin_layers_17.5.txt'
    coordfile = 'coord_high_thin_layers_PT_17.5.txt'
    cf = open(coordfile, 'r')
    cfl = cf.readlines()
    cf.close()
    days = []
    for cdl in cfl:
        if (cdl[0] == 't') or (cdl[0] == 'f') or (cdl[0] == '#') or (cdl == '\n'):
            continue
        if '_PT_' in coordfile:
            datestr = cdl.split(' ')[1].split('_')[0]
        else:
            datestr = cdl.split(' ')[0].split('_')[0]
        if datestr in days:
            continue
        else:
            days.append(datestr)
    # pfs = glob.glob('%s/Part/part_000-Coord_*' %dataDir)
    # pfs.sort()
    # pf = pfs[0]
    for day in days:
        # if day in ['20200111']:
        #     continue
        if int(day[0:4]) in [2020]:
            continue
        
        print(day)
        
        # day = pf.split('_')[-1]
        year = int(day[0:4])
        mon = int(day[4:6])
        d = int(day[6:8])
    
    
        startDate = datetime.datetime(year, mon, d, 0, 0) + datetime.timedelta(days=1)
        endDate = startDate - relativedelta(days=max_life_time + 1)

    # year_p1 = startDate.year
    # mon_p1 = startDate.month
    # year_m2 = endDate.year
    # mon_m2 = endDate.month
    
    
        #: Do file to start tracilla with qsub
        dofn = '%s/do-Coord-%s' %(runDir, day)
        #: File with paths used by tracilla
        pathfn = '%s/path-Coord-%s' %(runDir, day)
    
        #: Result dir
        resultDir = '%s/Coord-%s' %(ekjDir, day)
        #: Path for RELEASES and COMMAND files
        optPath = '%s/Coord-%s' %(optDir, day)
    
        #: File with the printstatment traxilla is creating during run
        trazilla_outfn = '%s/traczilla-Coord-%s' %(runDir, day)
        #: Part file
        partName = 'part_000-Coord_%s' %day #os.path.basename(pf)
        
        #: Qsub
        #:Outfile for qsub
        # qsub_outfile = '%s/O-${PBS_JOBNAME}-%s' %(runDir, dnf)
        qsub_outfile = '%s/O-Coord-%s' %(runDir, day)
        #: Jobname
        jobname = 'Co-%s' %(day)
    
    
        #: Dir with initiation files
        iniDir = '%s/Initfiles' %resultDir
        #: Initiation files
        # paramName = partName.replace('part_000-', '%s_Params-' %iniStartName).replace('-DD', '') + '.pkl'
        # catalogName = partName.replace('part_000-', '%s_Calalog-' %iniStartName).replace('-DD', '') + '.pkl'
        partFile = '%s/Part/%s' %(dataDir, partName) #pf
        # paramFile = partFile.replace('/Part/', '/Param/').replace(partName, paramName)
        # catalogFile = partFile.replace('/Part/', '/Catalog/').replace(partName, catalogName) 

        #: Links to files used by tracilla
        #: startfile i.e. file from selCaliop
        outlink = '%s/part_000' %resultDir
        #: Link to RELEASES file.
        # link_Relise = '%s/RELEASES' %optPath
        releasfn = '%s/RELEASES' %optPath
        #: No link. File with commands for tracilla. Is created here
        commandfn = '%s/COMMAND' %optPath
    
        #: Restartung an old run
        if False:#args.cont:
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
                # shutil.move(paramFile, iniDir + '/' + paramName)
                # shutil.move(catalogFile, iniDir + '/' + catalogName)
            
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
        # if not os.path.islink(link_Relise):
        #     os.symlink(src_Relise, link_Relise)
        # if not os.path.exists(link_Relise):
        #     print('broken link')
        #     print('src = %s' %src_Relise)
        #     print('link = %s' %link_Relise)
        #     sys.exit()
        
            #: Create do-file
            print(dofn)
            rmfile(dofn)
            dof = open(dofn, 'w')
            dof.writelines('#!/bin/bash \n')
            dof.writelines('#PBS -q short \n')
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
            #: Create path-file
            print(pathfn)
            rmfile(pathfn)
            pathf = open(pathfn, 'w')
        #     pathf.writelines('/home/ejohansson/Projects/flexpart/traczilla/options/STC/Calipso/CALIOP-EAD-Jul2008/ \n')
            pathf.writelines('%s/ \n' %optPath)
            pathf.writelines('%s/ \n' %resultDir)
            pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
            pathf.writelines('%s/AVAILABLE-%d-uvwt \n' %(availDir, year))
            pathf.writelines('/data/legras/flexpart_in/ERA5/indexes/ \n')
            pathf.writelines('%s/AVAILABLE-%d-hr \n' %(availDir, year))
            pathf.close()
        
            #: Create releas-file
            print(releasfn)
            rmfile(releasfn)
            releasf = open(releasfn, 'w')
            releasf.writelines("#RELEASES\n")
            releasf.writelines("#************************************************************************\n")
            releasf.writelines("#                                                                       *\n")
            releasf.writelines("#   Input file for the Lagrangian particle dispersion model TRACZILLA   *\n")
            releasf.writelines("#                        Please select your option    *\n")
            releasf.writelines("#                                                                       *\n")
            releasf.writelines("#************************************************************************\n")
            releasf.writelines("&StratoClim\n")
            releasf.writelines("pcut=50000\n")
            releasf.writelines("plowcut=50\n")
            releasf.writelines("NearRealTime=.false.\n")
            releasf.writelines("startfrom0=.true.\n")
            releasf.writelines("setxylim=.false.\n")
            releasf.writelines("max_life_time=%d\n" %max_life_time)
            releasf.writelines("&END\n")
            # releasf.writelines("\n")
            releasf.close()
        
        #: Create command-file
        print(commandfn)
        rmfile(commandfn)
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
        commandf.writelines(" loutstep=%d\n" %(outputHour * 3600))
        commandf.writelines(" loutprint=%d\n" %(outputHour * 3600))
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
    
        if True:#args.qsub:
            cmd = "qsub %s" %dofn
            print("Submitt qsub")
            print(cmd)
            c1 = os.system(cmd)
            # pdb.set_trace()
        else:
            pdb.set_trace()
    
    