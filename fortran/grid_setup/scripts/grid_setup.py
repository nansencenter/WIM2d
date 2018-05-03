import os,sys
import matplotlib.pyplot as plt
from getopt import getopt

sys.path.append('../py_funs')
import fns_grid_setup as gs

CHANGE_GRID     = 1
CHANGE_WAVES    = 1

# check arguments given to script in command line
GRID_OPT     = None
opts,args    = getopt(sys.argv[1:],"",["GRID_OPT="])
for opt,arg in opts:
    if opt=='--GRID_OPT':
        GRID_OPT = int(arg)
        # GRID_OPT = 0 # standard 1d grid configuration
        # GRID_OPT = 1 # standard 2d grid configuration
        # GRID_OPT = 2 # Philipp's small-square


if GRID_OPT is None:
    # no input given in command line
    infile    = 'infile_grid.txt'
    if os.path.exists(infile):
        # read in GRID_OPT from file:
        GRID_OPT = infile
    else:
        ss = "either give 'GRID_OPT' argument to script\n"\
             +"eg python --GRID_OPT=1\n"\
             +"or get infile_grid.txt from \n"\
             +os.getenv('WIM2D_PATH')+"/grid_setup/infiles\n"\
             +"(edit if necessary)"
        raise ValueError(ss)

LAND_OPT = 0
grid_fields,grid_arrays = gs.grid_setup(GRID_OPT=GRID_OPT,LAND_OPT=LAND_OPT)

if CHANGE_GRID:
    print('**************************************************')
    print('Testing values of grid_fields\n')
    keys  = grid_fields.keys()

    for key in keys:
        fld    = grid_fields[key]
        if hasattr(fld,'shape'):
            # array:
            print('field: '+key)
            print('field min: '+str(fld.min()))
            print('field max: '+str(fld.max())+'\n')
        else:
            # scalar:
            print(key+' = '+str(fld)+'\n')

    print('**************************************************')

if CHANGE_WAVES:
    # change number of wave frequencies and directions

    infile    = 'infile_waves.txt'
    hfil      = 'wave_info.h'

    if os.path.exists(infile):
        # get info from infile
        fid        = open(infile)
        lines     = fid.readlines()
        fid.close()

        vers_no  = int  (lines[0].split()[0])
        nfreq     = int  (lines[1].split()[0])
        ndir      = int  (lines[2].split()[0])
        Tmin      = float(lines[3].split()[0])
        Tmax      = float(lines[4].split()[0])

    else:
        # set here
        Tmin  = 2.5 # min period
        Tmax  = 25  # max period
        nfreq = 1    # no of freq's
        ndir  = 16  # no of dir's

    print(' ')
    print('**************************************************')
    print('Changing number of wave frequencies and directions:')
    print('nfreq : '+str(nfreq))
    print('ndir  : '+str(ndir))
    print('editing '+hfil)
    print('**************************************************')
    print(' ')

    """
    wave_info.h should look like this:
    integer,parameter :: n_wavdir     = 16
    integer,parameter :: n_wave_freq = 1!should be odd (Simpson's rule)

    real,parameter :: Tmin  = 2.5
    real,parameter :: Tmax  = 25.0
    """

    spc    = 6*' ' # F77 indent
    str1  = 'integer,parameter :: n_wavdir     = '
    str2  = 'integer,parameter :: n_wave_freq = '
    str3  = "!should be odd (Simpson's rule)"
    str4  = 'real,parameter :: Tmin  = '
    str5  = 'real,parameter :: Tmax  = '


    # write file
    hf = open(hfil,'w')
    hf.write(spc+str1+str(ndir)+'\n')
    hf.write(spc+str2+str(nfreq)+str3+'\n')
    hf.write('\n')
    hf.write(spc+str4+str(Tmin)+'\n')
    hf.write(spc+str5+str(Tmax)+'\n')
    hf.close()
