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

if 1:
   # test
   outdir   = 'test/out'
   outdir2  = 'test/out'
   nc       = len(outdir)
   nc2      = len(outdir2)
else:
   # proper places
   outdir   = '../run/inputs'
   outdir2  = '../header_files'
   nc       = len(outdir)
   nc2      = len(outdir2)

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
# make arrays
vx = np.array(range(0,nx))*dx
vx = vx-vx.mean()
vy = np.array(range(0,ny))*dy
vy = vy-vy.mean()

grid_arrays = np.zeros((nx,ny,7))
grid_fields = {'X'         :0.0,
               'Y'         :0.0,
               'scuy'      :0.0,
               'scvx'      :0.0,
               'scp2'      :0.0,
               'scp2i'     :0.0,
               'LANDMASK'  :0.0}

gf = grid_fields  # pointer to grid_fields
oo = np.ones((nx,ny))

gf['X']        = np.zeros((nx,ny))
gf['Y']        = np.zeros((nx,ny))
gf['scuy']     = oo*dy
gf['scvx']     = oo*dx
gf['scp2']     = oo*dx*dy
gf['scp2i']    = oo/dx/dy
gf['LANDMASK'] = np.zeros((nx,ny))

for j in range(ny):
   gf['X'][:,j]   = vx
for i in range(nx):
   gf['Y'][i,:]   = vy

LAND_OPT = 1
if LAND_OPT is 1:
   # standard 2d setup (works with the standard 2d grid):
   gf['LANDMASK'][113:,:]  = 1.0
   print('max LANDMASK:')
   print(gf['LANDMASK'].max())

###########################################################

###########################################################
keys  = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
n     = 0
for key in keys:
   grid_arrays[:,:,n]   = gf[key]
   n                    = n+1
###########################################################

###########################################################
if 1:
   # save test figure:
   fig   = 'test/out/land.png'
   print('Saving test figure : '+fig)
   #
   f     = Fplt.cmap_3d(gf['X']/1.0e3,gf['Y']/1.0e3,
                        gf['LANDMASK'],
                        ['$x$, km','$y$, km','LANDMASK'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()
###########################################################

###########################################################
if 1:
   # finish by saving files
   gs.save_grid_info_hdr(outdir2,nx,ny,dx,dy,nc2)
   gs.save_grid(outdir,grid_arrays,nc)
###########################################################


# if 1:
#    RUN_OPT     = 2
#    out_fields  = Rwim.do_run(RUN_OPT)
# elif 1:
#    # check passing in of 'in_fields'
#    # - read in inputs from saved files:
#    in_dir                  = 'out'
#    grid_prams              = Fdat.fn_check_grid(in_dir)
#    ice_fields,wave_fields  = Fdat.fn_check_init(in_dir)
# 
#    # merge ice and wave fields:
#    ice_fields.update(wave_fields)
#    in_fields   = ice_fields
# 
#    out_fields  = Rwim.do_run(RUN_OPT=2,in_fields=in_fields)
# elif 0:
#    out_fields  = Rwim.do_run(0)
#    out_fields2 = Rwim.do_run(2)
# elif 0:
#    out_fields  = Rwim.do_run(1)
#    out_fields2 = Rwim.do_run(3)
# 
# if 0:
#    # plot results
#    ## look at initial fields:
#    print("Plotting initial conditions...")
#    grid_prams              = Fdat.fn_check_grid(outdir) # load grid from binaries
#    ice_fields,wave_fields  = Fdat.fn_check_init(outdir) # load initial conditions from binaries
#    ##
#    Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
#    print("Plots in "+figdir+"/init")
#    print(" ")
# 
#    ## look at results:
#    print("Plotting results...")
#    Fplt.fn_plot_final(grid_prams,out_fields,figdir)
#    print("Plots in "+figdir+"/final")
# elif 0:
#    # compare binaries from different runs (wim2d_run & wim2d_io)
#    # NB need same grid & initial conditions
#    outdir1  = 'out'
#    outdir2  = 'out_io'
# 
#    gp1   = Fdat.fn_check_grid(outdir1)
#    gp2   = Fdat.fn_check_grid(outdir2)
#    if 0:
#       # check grids are the same
#       print('Checking grids are the same...')
#       ##
#       keys  = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
#       Key   = '         '
#       for key in keys:
#          diff     = np.abs(gp1[key]-gp2[key])
#          diff_max = diff.max()
#          diff_min = diff.min()
#          #
#          ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
#          print(' '+key+Key[len(key):]+ss)
#    elif 0:
#       # check initial fields are the same
#       print('Checking initial fields are the same...')
#       ##
#       if1,wf1  = Fdat.fn_check_init(outdir1)
#       if2,wf2  = Fdat.fn_check_init(outdir2)
# 
#       Key   = '         '
#       for key in if1.keys():
#          diff     = np.abs(if1[key]-if2[key])
#          diff_max = diff.max()
#          diff_min = diff.min()
#          #
#          ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
#          print(' '+key+Key[len(key):]+ss)
# 
#       for key in wf1.keys():
#          diff     = np.abs(wf1[key]-wf2[key])
#          diff_max = diff.max()
#          diff_min = diff.min()
#          #
#          ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
#          print(' '+key+Key[len(key):]+ss)
# 
#    elif 1:
#       # check out fields are the same
#       print('Checking final fields are the same...')
#       ##
#       of1   = Fdat.fn_check_out_bin(outdir1)
#       of2   = Fdat.fn_check_out_bin(outdir2)
# 
#       Key   = '         '
#       for key in of1.keys():
#          diff     = np.abs(of1[key]-of2[key])
#          diff_max = diff.max()
#          diff_min = diff.min()
#          #
#          ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
#          print(' '+key+Key[len(key):]+ss)
# 
# # elif 1:
# #    # compare exponential decay of Hs for SOLVER = 1,0
# # 
