import numpy as np
import os
import sys
import matplotlib.pyplot as plt
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/misc_py")

import save_grid_f2py   as gs
import fns_get_data     as Fdat
import fns_plot_data    as Fplt

if 0:
   # test
   outdir   = 'test/out_py'
   outdir2  = 'test/out_py'
   nc       = len(outdir)
   nc2      = len(outdir2)
else:
   # proper places
   outdir   = '../run/inputs'
   outdir2  = '../header_files'
   nc       = len(outdir)
   nc2      = len(outdir2)


dd    = os.path.abspath(".")
dd2   = dd+'/'+outdir
print(dd2)
if not (os.path.exists(dd2)):
   os.makedirs(dd2)

###########################################################
# grid size
GRID_OPT = 1
if GRID_OPT is -1:
   # test:
   nx = 10
   ny = 4
   dx = 3.0e3
   dy = 3.0e3

elif GRID_OPT is 0:
   # standard 1d setup:
   nx = 150
   ny = 1
   dx = 4.0e3
   dy = 4.0e3

elif GRID_OPT is 1:
   # standard 2d setup:
   nx = 150
   ny = 50
   dx = 4.0e3
   dy = 4.0e3
###########################################################

###########################################################
def get_grid_arrays(nx,ny,dx,dy):

   vx = np.array(range(0,nx))*dx
   vx = vx-vx.mean()
   vy = np.array(range(0,ny))*dy
   vy = vy-vy.mean()
   oo = np.ones((nx,ny))

   grid_arrays = np.zeros((nx,ny,7))
   grid_fields = Fdat.Grid_Prams(nx=nx,ny=ny,dx=dx,dy=dy)

   gf             = grid_fields  # pointer (with short name) to grid_fields
   gf['X']        = np.zeros((nx,ny))
   gf['Y']        = np.zeros((nx,ny))
   gf['scuy']     = oo*dy
   gf['scvx']     = oo*dx
   gf['scp2']     = oo*dx*dy
   gf['scp2i']    = oo/dx/dy
   gf['LANDMASK'] = np.zeros((nx,ny))

   # X,Y:
   for j in range(ny):
      gf['X'][:,j]   = vx
   for i in range(nx):
      gf['Y'][i,:]   = vy

   # LANDMASK:
   LAND_OPT = 1
   if LAND_OPT is 1:
      # standard 2d setup (works with the standard 2d grid):
      gf['LANDMASK'][113:,:]  = 1.0

   ###########################################################

   ###########################################################
   # assign integers to sort the fields when
   # passing stuff to grid_arrays
   gf0   = {'X':0,'Y':1,'scuy':2,'scvx':3,'scp2':4,
            'scp2i':5,'LANDMASK':6}

   for key in gf0.keys():
      grid_arrays[:,:,gf0[key]]  = gf[key]

   return grid_arrays,grid_fields
###########################################################

###########################################################
grid_arrays,grid_fields = get_grid_arrays(nx,ny,dx,dy)
print(' ')
print(60*'*')
gs.save_grid_info_hdr(outdir2,nx,ny,dx,dy,nc2)
gs.save_grid(outdir,grid_arrays,nc)
print(60*'*')
print(' ')
###########################################################

###########################################################
if 0:
   # check difference between binaries
   # saved by pure fortran and python:
   outdir1  = 'test/out'   # fortran binaries here
                           # run grid_setup.sh with
                           # testing=1 in p_save_grid.F (recompile)
   gf1      = Fdat.fn_check_grid(outdir1)
   gf2      = Fdat.fn_check_grid(outdir)
   keys     = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']

   print('Comparing fortran to python:\n')
   for key in keys:
      diff     = np.abs(gf1[key]-gf2[key])
      diff_max = diff.max()
      diff_min = diff.min()
      print('max difference in : '+key+' = '+str(diff_max))
      print('min difference in : '+key+' = '+str(diff_min)+'\n')
elif 1:
   # check difference between binaries
   # saved by python and the input fields:
   gf    = grid_fields
   gf2   = Fdat.fn_check_grid(outdir)
   keys  = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']

   print('Comparing python in to python out:\n')
   for key in keys:
      diff     = np.abs(gf[key]-gf2[key])
      diff_max = diff.max()
      diff_min = diff.min()
      print('max difference in : '+key+' = '+str(diff_max))
      print('min difference in : '+key+' = '+str(diff_min)+'\n')
###########################################################

###########################################################
if 1:
   # save test figure:
   fig   = 'test/out_py/land.png'
   print('Saving test figure : '+fig)
   #
   gf = grid_fields
   f  = Fplt.cmap_3d(gf['X']/1.0e3,gf['Y']/1.0e3,
                     gf['LANDMASK'],
                     ['$x$, km','$y$, km','LANDMASK'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()
###########################################################
