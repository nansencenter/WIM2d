import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt

gsd = os.path.join(
        os.getenv("WIM2D_PATH"),
        'fortran',
        'grid_setup',
        )
sys.path.append(gsd+"/bin")
sys.path.append(gsd+"/py_funs")
import save_grid_f2py as gs
import fns_get_data   as Fdat
import fns_plot_data  as Fplt

def _get_grid_arrays_SmallSquare(diag_length,resolution):

    dx  = resolution
    nx  = int(np.floor(diag_length/resolution))
    x0  = -(nx*dx)/2.
    out = get_grid_arrays(x0=x0,y0=x0,nx=nx,ny=nx,dx=dx,dy=dx,LAND_OPT=0)

    # fix landmask
    gf   = out[1]
    xmin = np.min(gf['X'])
    xmax = np.max(gf['X'])
    ymin = np.min(gf['Y'])
    ymax = np.max(gf['Y'])
    #
    X        = gf['X']
    Y        = gf['Y']
    LANDMASK = np.zeros((nx,nx))
    
    # land outside the diamond of Philipp:
    LANDMASK[Y>(X-xmin)]  = 1.
    LANDMASK[Y<-(X-xmin)] = 1.
    LANDMASK[Y<(X-xmax)]  = 1.
    LANDMASK[Y>-(X-xmax)] = 1.

    out[1]['LANDMASK'] = LANDMASK
    out[0][:,:,6]      = LANDMASK

    return out # = [grid_arrays,grid_fields]


def get_grid_arrays(x0=0.,y0=0.,nx=100,ny=4,dx=4.e3,dy=4.e3,LAND_OPT=0):

    # get centres of grid
    vx = np.arange(.5,nx)*dx
    vy = np.arange(.5,ny)*dy
    oo = np.ones((nx,ny))

    grid_arrays = np.zeros((nx,ny,7))
    grid_fields = {}
    keys        = ['nx','ny','dx','dy']
    vals        = [nx,ny,dx,dy]

    # make placeholders for arrays
    keys.extend(['X','Y','scuy','scvx','scp2','scp2i','LANDMASK'])
    vals.extend(7*[0])
    for n,key in enumerate(keys):
        val = vals[n]
        grid_fields.update({key:val})

    gf             = grid_fields  # pointer (with short name) to grid_fields
    gf['X']        = np.zeros((nx,ny))
    gf['Y']        = np.zeros((nx,ny))
    gf['scuy']     = dy*oo
    gf['scvx']     = dx*oo
    gf['scp2']     = (dx*dy)*oo
    gf['scp2i']    = (1./dx/dy)*oo
    gf['LANDMASK'] = np.zeros((nx,ny))

    # X,Y:
    for j in range(ny):
        gf['X'][:,j] = vx
    for i in range(nx):
        gf['Y'][i,:] = vy

    # LANDMASK:
    if LAND_OPT is 0:
        print('\nNo land present\n')
    elif LAND_OPT is 1:
        # standard 2d setup (works with the standard 1d/2d grids):
        xl                                = x0+.9*(vx.max()-x0)
        gf['LANDMASK'][vx>xl,:]    = 1.0
        print('\nLand in last 10% of domain\n')
    elif LAND_OPT is 2:
        # standard 2d, but with an island
        xc = vx.mean()
        yc = vy.mean()
        Rc = .1*(vy.max()-y0)

        for i in range(nx):
            for j in range(ny):
                x  = gf['X'][i,j]
                y  = gf['Y'][i,j]
                R  = np.sqrt(pow(x-xc,2)+pow(y-yc,2))
                if R<Rc:
                    gf['LANDMASK'][i,j]  = 1.
        print('\nMaking island\n')
    elif LAND_OPT is 3:
        print('\nLand on 3 out of 4 boundaries (not left-most one)\n')
        gf['LANDMASK'][-1,: ] = 1. # far right cells
        gf['LANDMASK'][: ,0 ] = 1.# bottom cells
        gf['LANDMASK'][: ,-1] = 1. # top cells
    else:
        raise ValueError('Invalid value of LAND_OPT ('+str(LAND_OPT)+' - should be 0/1/2/3)')


    # assign integers to sort the fields when
    # passing stuff to grid_arrays
    gf0 = {'X':0,'Y':1,'scuy':2,'scvx':3,'scp2':4,
                'scp2i':5,'LANDMASK':6}

    for key in gf0:
        grid_arrays[:,:,gf0[key]] = gf[key]

    return grid_arrays,grid_fields


def grid_setup(GRID_OPT=1,LAND_OPT=0):

    outdir  = '.' # wim_grid.[a,b] saved here
    outdir2 = '.' # wave_info.h,grid_info.h saved here
    dd2     = os.path.abspath(outdir)
    print('\n')
    print("Saving grid files to:")
    print(dd2)
    print('\n')
    if not os.path.exists(dd2):
        os.mkdir(dd2)

    # set grid configurations
    if type(GRID_OPT)==type('string'):
        # read in parameters from infile
        infile = GRID_OPT
        print('using infile '+infile+' for parameters\n')

        fid   = open(infile)
        lines = fid.readlines()
        fid.close()

        # check version number
        Vno = 1
        vno = int(lines[0].split()[0])
        if vno!=Vno:
            s1 = '\nWrong version of '+infile
            s2 = '\nCurrent version number is '+str(vno)
            s3 = '\nVersion number should be  '+str(Vno)+'\n'
            raise ValueError(s1+s2+s3)

        # get values as floats 1st, then convert some to integers
        params = []
        for lin in lines[1:]:
            lin2 = lin.split()
            params.append(float(lin2[0]))
        nx,ny,dx,dy,x0,y0,LAND_OPT = params

        # convert some parameters to integers
        nx       = int(nx)
        ny       = int(ny)
        LAND_OPT = int(LAND_OPT)
        #
        grid_arrays,grid_fields = get_grid_arrays(x0=x0,y0=y0,nx=nx,ny=ny,\
                dx=dx,dy=dy,LAND_OPT=LAND_OPT)

    elif GRID_OPT is 0:
        # standard 1d setup:
        nx = 150
        ny = 4
        dx = 4.0e3
        dy = 10*dx # to make plots clearer
        x0 = 0.
        y0 = 0.
        #
        grid_arrays,grid_fields = get_grid_arrays(x0=x0,y0=y0,nx=nx,ny=ny,\
                dx=dx,dy=dy,LAND_OPT=LAND_OPT)

    elif GRID_OPT is 1:
        # standard 2d setup:
        nx = 150
        ny = 50
        dx = 4.0e3
        dy = 4.0e3
        x0 = 0.
        y0 = 0.
        #
        grid_arrays,grid_fields = get_grid_arrays(x0=x0,y0=y0,nx=nx,ny=ny,\
                dx=dx,dy=dy,LAND_OPT=LAND_OPT)

    elif GRID_OPT is 2:
        # to use with Philipp's "small-square" grid
        diag_length = 96e3    # m
        resolution  = 0.75e3  # m

        # this gives a rotated (x',y') grid to align with the sides of the square
        grid_arrays,grid_fields = _get_grid_arrays_SmallSquare(diag_length,resolution)
        #
        nx = grid_fields['nx']
        ny = grid_fields['ny']
        dx = grid_fields['dx']
        dy = grid_fields['dy']
        print('nx = '+str(nx))
        print('ny = '+str(ny))
        print('dx = '+str(dx/1.e3)+'km')
        print('dy = '+str(dy/1.e3)+'km')

    print(' ')
    print(60*'*')
    gs.save_grid_info_hdr_f2py(outdir2,nx,ny,dx,dy)
    gs.save_grid_f2py(outdir,grid_arrays)
    print(60*'*')
    print(' ')

    if 0:
        # check difference between binaries
        # saved by pure fortran and f2py:
        outdir1 = 'test/out'    # fortran binaries here
                                        # run grid_setup.sh with
                                        # testing=1 in p_save_grid.F (recompile)
        gf1  = Fdat.fn_check_grid(outdir1)
        gf2  = Fdat.fn_check_grid(outdir)
        keys = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']

        print('Comparing fortran to python:\n')
        for key in keys:
            diff     = np.abs(gf1[key]-gf2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            print('max difference in : '+key+' = '+str(diff_max))
            print('min difference in : '+key+' = '+str(diff_min)+'\n')
    elif 1:
        # check difference between binaries
        # saved by f2py and the input fields:
        gf   = grid_fields
        gf2  = Fdat.fn_check_grid(outdir)
        keys = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']

        print('Comparing python in to python out:\n')
        for key in keys:
            diff     = np.abs(gf[key]-gf2[key])
            diff_max = diff.max()
            diff_min = diff.min()
            print('max difference in : '+key+' = '+str(diff_max))
            print('min difference in : '+key+' = '+str(diff_min)+'\n')

    SV_FIG = 1
    if SV_FIG:
        # save test figure:
        if not os.path.exists('out_py'):
            os.mkdir('out_py')
        fig = 'out_py/land.png'
        #
        gf = grid_fields
        nx = gf['nx']
        ny = gf['ny']
        dx = gf['dx']
        dy = gf['dy']
        x  = gf['X'][:,0]+dx/2.
        x  = np.concatenate([[x[0]-dx],x])/1.e3
        y  = gf['Y'][0,:]+dy/2.
        y  = np.concatenate([[y[0]-dy],y])/1.e3

        if gf['ny']>1:
            # 2d grid => make colormap of LANDMASK
            z    = gf['LANDMASK'].transpose()
            f,ax = Fplt.cmap_3d(x,y,z,['$x$, km','$y$, km','LANDMASK'])
        else:
            # 1d grid => make 1d plot of LANDMASK
            f,ax,line = Fplt.plot_1d(gf['X'][:,0]/1.0e3, gf['LANDMASK'][:,0],\
                                        labs=['$x$, km', 'LANDMASK'])

        print('Saving test figure : '+fig)
        f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
        plt.close(f)

    # print(' ')
    # print(60*'*')
    # print("Now compile WIM code in .Build")
    # print("Run in                        ..")
    # print(60*'*')
    # print(' ')

    return grid_fields,grid_arrays
