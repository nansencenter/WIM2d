import numpy as np
import os
import sys
import shutil
import struct
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
w2d = os.getenv("WIM2D_PATH")
dd  = os.path.abspath(".")
sys.path.append(dd+"/bin")
sys.path.append(w2d+"/fortran/py_funs")

import run_WIM2d     as Rwim
import WIM2d_f2py    as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

testing  = 1
grid_dir = 'grid'

if testing is 1:
    # run and plot figs later (no I/O)
    RUN_OPT           = 0
    out_fields,outdir = Rwim.do_run(RUN_OPT)
    sys.exit()

elif testing is 2:

    RUN_OPT  = 2
    # check passing in of 'in_fields'
    # - read in inputs from saved files:
    # (need to run without I/O first)
    if 1:
        in_dir = 'out_py/binaries'
    else:
        in_dir = 'out_py_io/binaries'

    gf                     = Fdat.fn_check_grid(in_dir)
    grid_prams             = gf
    ice_fields,wave_fields = Fdat.fn_check_init(in_dir)

    # merge ice and wave fields:
    ice_fields.update(wave_fields)
    in_fields = ice_fields

    # other inputs:
    int_prams   = None # default parameters
    real_prams  = None # default parameters

    if 1:
        # check passing in of integer parameters:
        SCATMOD     = 0
        ADV_DIM     = 2
        CHECK_FINAL = 1
        CHECK_PROG  = 1
        CHECK_INIT  = 1
        int_prams   = np.array([SCATMOD,ADV_DIM,
                                        CHECK_FINAL,CHECK_PROG,CHECK_INIT])

    if 1:
        # check passing in of real parameters:
        young          = 5.0e9
        visc_rp        = 0.0
        duration_hours = 12.0
        duration       = duration_hours*60*60
        real_prams     = np.array([young,visc_rp,duration])

    
    out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                                    int_prams=int_prams,
                                                    real_prams=real_prams)

elif testing is 3:
    # compare runs with and without input/output:
    if 1:
        # need to rerun:
        out_fields,outdir  = Rwim.do_run(RUN_OPT=0)
    else:
        # use saved results:
        out_fields,outdir  = Rwim.do_run(RUN_OPT=1)

    if 1:
        # need to rerun:
        out_fields,outdir3 = Rwim.do_run(RUN_OPT=2)
    else:
        # use saved results:
        out_fields2,outdir3  = Rwim.do_run(RUN_OPT=3)


if (testing is 1) or (testing is 2):

    # plot results
    results  = Fdat.wim_results(outdir)
    results.plot_initial() # look at initial  fields:
    results.plot_final()   # look at final     fields:
    results.plot_prog()    # look at progress fields:

elif testing is 3:
    # compare binaries from different runs (wim2d_run & wim2d_io)
    # NB need same grid & initial conditions
    bindir1 = outdir+'/binaries'
    bindir2 = outdir3+'/binaries'

    gp1 = Fdat.fn_check_grid(bindir1)
    gp2 = Fdat.fn_check_grid(bindir2)
    if 0:
        # check grids are the same
        print('Checking grids are the same...')
        ##
        keys = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
        Key  = '            '
        for key in keys:
            diff     = np.abs(gp1[key]-gp2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            #
            ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
            print(' '+key+Key[len(key):]+ss)
    elif 0:
        # check initial fields are the same
        print('Checking initial fields are the same...')
        ##
        if1,wf1 = Fdat.fn_check_init(bindir1)
        if2,wf2 = Fdat.fn_check_init(bindir2)

        Key = '            '
        for key in if1.keys():
            diff     = np.abs(if1[key]-if2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            #
            ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
            print(' '+key+Key[len(key):]+ss)

        for key in wf1.keys():
            diff     = np.abs(wf1[key]-wf2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            #
            ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
            print(' '+key+Key[len(key):]+ss)

    elif 1:
        # check out fields are the same
        print('Checking final fields are the same...')
        ##
        of1 = Fdat.fn_check_out_bin(bindir1)
        of2 = Fdat.fn_check_out_bin(bindir2)

        Key = '            '
        for key in of1.keys():
            diff      = np.abs(of1[key]-of2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            #
            ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
            print(' '+key+Key[len(key):]+ss)
