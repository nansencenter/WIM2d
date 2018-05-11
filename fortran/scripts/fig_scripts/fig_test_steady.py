import numpy as np
import os
import sys
import shutil
import struct
from matplotlib import pyplot as plt
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
fd = os.path.join(os.getenv('WIM2D_PATH'),
          'fortran')
sys.path.append(os.path.join(fd,"bin"))
sys.path.append(os.path.join(fd,"py_funs"))

# interface to fortran code
# - compile with single frequency
import run_WIM2d as Rwim

# other python modules
import WIM2d_f2py    as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# steady state results to compare to:
import fns_boltzmann_steady    as Fbs

RUN_OPT  = 2 # rerun then plot
# RUN_OPT  = 3 # plot saved results

gf             = Fdat.fn_check_grid('inputs')
gfl            = gf['LANDMASK']
ICEMASK      = 1.-gfl
WAVEMASK     = 1.-gfl
grid_prams  = gf

# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

    # ice band of finite width
    xe                  = -220.e3
    ICEMASK             = 1+0*gf['X']
    ICEMASK[gf['X']<xe] = 0.
    ICEMASK[gfl>0]      = 0.

    # right hand ice edge
    ice_width                       = 100.e3
    ICEMASK[gf['X']>(xe+ice_width)] = 0.

    # edge of wave mask
    xw                   = -260.e3
    WAVEMASK             = 1+0*gf['X']
    WAVEMASK[gf['X']>xw] = 0.
    WAVEMASK[gfl>0]      = 0.

    c_in      = .7   # conc
    h_in      = 2.   # thickness
    D_in      = 100. # initial floe size
    Hs_in     = 3.   # initial Hs
    Tp_in     = 12.  # initial Tp
    mwd_in    = -90  # initial mean wave direction

    in_fields    = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK,
                        'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}

int_prams  = None # default integer parameters
real_prams = None # default real parameters

if 1:
    # change integer parameters:
    SCATMOD     = 1
    ADV_DIM     = 1
    CHECK_FINAL = 1
    CHECK_PROG  = 1
    CHECK_INIT  = 1
    DO_BREAKING = 0 # no breaking - testing convergence of Hs to steady-state
    int_prams   = np.array([SCATMOD,ADV_DIM,
                                    CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                                    DO_BREAKING])

if 1:
    # change real parameters:
    young          = 5.49e9
    visc_rp        = 0.0
    duration_hours = 48.0
    duration       = duration_hours*60*60
    real_prams     = np.array([young,visc_rp,duration])


# out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
#                                                 int_prams=int_prams,
#                                                 real_prams=real_prams)

cols   = ['k','b','r','g','m','c']
lstil  = ['-','--','-.',':']
Nc     = len(cols)
loop_c = -1
loop_s = 0

ndirs = [16,32,64,128,256]
fig   = None
labs  = ['$x$, km','$H_s$, m']

for N in ndirs:
    # get steady state solution
    om = 2*np.pi/Tp_in
    # gravity = 9.81
    atten_in  = np.array([h_in,om,young,visc_rp])
    atten_out = Mwim.atten_youngs(atten_in)
    alp       = c_in/D_in*atten_out[4]    # scattering "attenuation" [m^{-1}]
    alp_dis   = 2*c_in*atten_out[0]        # damping [m^{-1}]
    kwtr      = atten_out[2]
    cp        = om/kwtr # phase vel (open water) [m/s]
    cg        = cp/2.    # group vel (open water, inf depth relation) [m/s]
    #
    #N      = pow(2,4)
    out    = Fbs.solve_boltzmann_ft(width=ice_width,
                alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs_in)

    # alp2  = out['inputs'][0]
    # print(alp,alp2)
    # sys.exit('fig_test_convergence2steady')

    nx0       = 1.e3
    x_ice     = np.linspace(0.,ice_width,nx0)
    E_n       = Fbs.calc_energy(out,x_ice,L=ice_width,n_test=[0])
    Hs_steady = 4*np.sqrt(E_n[:,0].real)

    # print(x_ice)
    # print(Hs_steady)
    # fig        = Fplt.plot_1d(1.e-3*(xe+x_ice),Hs_steady,labs=labs,color='y',linewidth=3)
    #out_fields  = Fdat.fn_check_prog(outdir,nstep) # load ice/wave conditions from binaries
    #Hs_n          = out_fields['Hs'][:,0]
    #
    if loop_c==Nc-1:
        loop_c = 0
        loop_s = loop_s+1
    else:
        loop_c = loop_c+1

    fig = Fplt.plot_1d(1.e-3*(xe+x_ice),Hs_steady,labs=labs,f=fig,color=cols[loop_c],linestyle=lstil[loop_s])
#
figdir  = 'fig_scripts/figs'
figname = figdir+'/test_steady.png'
print('saving to '+figname+'...')
plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
plt.close()
fig.clf()
# 
# if 1:
#     # print to dat-file
#     dfil  = figdir+'/test_steady1.dat'
#     fid    = open(dfil,'w')
#     blk    = 4*' '
#     for loop_x in range(len(xx)):
#         lin    =              ('%f' % (1.e3*xx[loop_x])  )
#         lin    = lin+blk + ('%f' % (Hs_n[loop_x])      )
#         fid.write(lin+'\n')
#     fid.close()
#     print('saving to '+dfil+'...')
#     ##
#     dfil  = figdir+'/test_steady2.dat'
#     fid    = open(dfil,'w')
#     blk    = 4*' '
#     for loop_x in range(len(x_ice)):
#         lin    =              ('%f' % (xe+x_ice[loop_x])  )
#         lin    = lin+blk + ('%f' % (Hs_steady[loop_x])      )
#         fid.write(lin+'\n')
#     fid.close()
#     print('saving to '+dfil+'...')
# 
# # sys.exit('fig_test_convergence2steady.py, L94')
# 
# # figdir1                      = figdir+'/init/'
# # ice_fields,wave_fields  = Fdat.fn_check_init(bindir)
# # Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir1) # plot initial conditions
# # print("Plots in "+figdir+"/init")
# # print(" ")
# # 
# # # Look at end results:
# # print("Plotting results...")
# # figdir2  = figdir+'/final/'
# # Fplt.fn_plot_final(grid_prams,out_fields,figdir2)
# # print("Plots in "+figdir2+'\n')
# # print(" ")
# # 
# # if 1:
# #     # Plot progress files (if they exist)
# #     figdir3      = figdir+'/prog'
# #     prog_files  = os.listdir(bindir+'/prog')
# #     steps         = []
# #     for pf in prog_files:
# #         if '.a'==pf[-2:]:
# #             stepno    = pf[-5:-2]
# #             steps.append(stepno)
# # 
# #     # make dir for progress plots
# #     if (not os.path.exists(figdir3)) and len(steps)>0:
# #         os.mkdir(figdir3)
# # 
# #     # clear old progress plots
# #     old_dirs = os.listdir(figdir3)
# #     for od in old_dirs:
# #         # need this for macs
# #         if od!='.DS_Store':
# #             # os.rmdir(figdir3+'/'+od)
# #             shutil.rmtree(figdir3+'/'+od)
# # 
# #     for stepno in steps:
# #         print("Plotting results at time step "+stepno+" ...")
# #         prog_fields = Fdat.fn_check_prog(outdir,int(stepno))
# #         figdir3_0    = figdir3+'/'+stepno
# #         Fplt.fn_plot_final(grid_prams,prog_fields,figdir3_0)
# #         print("Plots in "+figdir3_0+'\n')
