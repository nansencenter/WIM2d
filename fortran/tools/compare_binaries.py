import numpy as np
import os,sys
import shutil
import struct

##
## NB run from 'run' directory !!
##
w2d    = os.getenv('WIM2D_PATH')
dd     = w2d+'/fortran'
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

import fns_get_data  as Fdat
import fns_plot_data as Fplt

mdir  = w2d+'/matlab/main/'
fdir  = w2d+'/fortran/run/'

# choose first set of results to compare
# odir     = [fdir+'out'    ] # fortran, non-IO
odir     = [fdir+'out_io'] # fortran, IO
# odir     = [fdir+'out_2' ] # fortran, IO, full spec
# odir     = [mdir+'m_out' ] # pure matlab
# odir     = [mdir+'out_io'] # mex, IO
# odir     = [mdir+'out_2' ] # mex, full spec

# choose second set of results to compare
# odir.append(fdir+'out'    ) # fortran, non-IO
# odir.append(fdir+'out_io') # fortran, IO
# odir.append(fdir+'out_2' ) # fortran, IO, full spec
odir.append(mdir+'m_out' ) # pure matlab
# odir.append(mdir+'out_io') # mex, IO
# odir.append(mdir+'out_2' ) # mex, full spec

print('Comparing directories:')
print(odir)
print('\n')

OPT    = 1    # 1: initial conditions; 2: final results; 3: prog results
if OPT==1:
    print("\n************************************************")
    print("Checking initial conditions...")
    print("************************************************\n")
elif OPT==2:
    print("\n************************************************")
    print("Checking final results...")
    print("************************************************\n")
elif OPT==3:
    print("\n************************************************")
    print("Checking progress files...")
    print("************************************************\n")
    n_prog    = 1

arrays    = 2*[0]
for j in range(2):
    outdir    = odir[j]
    bindir    = outdir+'/binaries'

    if OPT==1:

        # Look at initial fields:
        ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
        arrays[j]                    = ice_fields
        arrays[j].update(wave_fields)

    elif OPT==2:

        # Look at end results:
        out_fields  = Fdat.fn_check_out_bin(bindir)
        arrays[j]    = out_fields

    elif OPT==3:

        # Plot progress files (if they exist)
        # - as colormaps
        prog_files  = os.listdir(bindir+'/prog')
        steps         = []
        for pf in prog_files:
            if '.a'==pf[-2:]:
                stepno    = pf.strip('wim_prog').strip('.a')
                steps.append(stepno)

        print('Available time steps:')
        print(steps)

        stepno    = steps[n_prog]
        print("Checking results at time step "+stepno+" ...")
        prog_fields = Fdat.fn_check_prog(outdir,stepno)
        arrays[j]    = prog_fields
        #############################################################

keys  = arrays[0].keys()
for key in keys:
    print('Comparing field: '+key)
    diff  = abs(arrays[0][key]-arrays[1][key])
    print('Maximum difference = '+str(np.max(diff))+'\n')
