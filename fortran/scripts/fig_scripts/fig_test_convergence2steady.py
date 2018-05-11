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
dd    = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

# interface to fortran code
# - compile with single frequency
import run_WIM2d      as Rwim

# other python modules
import WIM2d_f2py     as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# steady state results to compare to:
import fns_boltzmann_steady as Fbs

# RUN_OPT  = 1 # plot from saved dir (not in default location)
# RUN_OPT  = 2 # rerun then plot
RUN_OPT  = 3 # plot saved results (from default location)

do_legend         = 1 # show times as legends
cmp_steady        = 1 # compare to steady state solution
cmp_steady_num  = 0 # compare to steady state Galerkin solution from matlab

gf             = Fdat.fn_check_grid('inputs')
gfl            = gf['LANDMASK']
ICEMASK      = 1.-gfl
WAVEMASK     = 1.-gfl
grid_prams  = gf

# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

    # ice band of finite width
    xe                         = 50.e3
    ICEMASK                  = 1+0*gf['X']
    ICEMASK[gf['X']<xe]  = 0.
    ICEMASK[gfl>0]         = 0.

    # right hand ice edge
    ice_width                                = 150.e3
    ICEMASK[gf['X']>(xe+ice_width)]  = 0.

    # edge of wave mask
    xw                         = 30.e3
    WAVEMASK                 = 1+0*gf['X']
    WAVEMASK[gf['X']>xw] = 0.
    WAVEMASK[gfl>0]        = 0.

    c_in      = .7      # conc
    h_in      = 2.      # thickness
    D_in      = 100.    # initial floe size
    Hs_in     = 3.      # initial Hs
    Tp_in     = 12.     # initial Tp
    mwd_in    = -90     # initial mean wave direction

    in_fields    = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK,
                        'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}

int_prams    = None # default integer parameters
real_prams  = None # default real parameters

if 1:
    # change integer parameters:
    SCATMOD      = 1
    ADV_DIM      = 1
    ADV_OPT      = 2
    CHECK_FINAL = 1
    CHECK_PROG  = 1
    CHECK_INIT  = 1
    STEADY        = 1
    DO_BREAKING = 0 # no breaking - testing convergence of Hs to steady-state
    DO_ATTEN     = 1 # no breaking - testing convergence of Hs to steady-state
    int_prams    = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                                    CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                                    STEADY,DO_BREAKING,DO_ATTEN])

if 1:
    # change real parameters:
    young             = 5.49e9
    visc_rp          = 0.0
    duration_hours = 72.0
    duration         = duration_hours*60*60
    CFL                = 0.7
    real_prams      = np.array([young,visc_rp,duration,CFL])

if RUN_OPT==1:
    # specify manually where results are
    # outdir    = '/Volumes/Tim_Ext_HD2/WORK/Model-Results/Boltzmann/convergence/16dirs/500m'
    # outdir    = '/Volumes/Tim_Ext_HD2/WORK/Model-Results/Boltzmann/convergence/16dirs/1km'
    # outdir    = '/Volumes/Tim_Ext_HD2/WORK/Model-Results/Boltzmann/convergence/16dirs/2km'
    outdir    = '/Volumes/Tim_Ext_HD2/WORK/Model-Results/Boltzmann/convergence/16dirs/4km'
    figdir    = outdir
else:
    figdir0  = 'fig_scripts/figs/'
    figdir    = figdir0+'TC2S'
    if not os.path.exists(figdir0):
        os.mkdir(figdir0)
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    if RUN_OPT==3:
        outdir    = 'out_io'                # usual place
        # outdir    = '../../matlab/main/m_out' # matlab results
    else:
        out_fields,outdir = \
            Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,\
                            int_prams=int_prams,\
                            real_prams=real_prams)

# Make plots
bindir    = outdir+'/binaries'

grid_prams  = Fdat.fn_check_grid(bindir) # load grid from binaries
xx             = 1.e-3*grid_prams['X'][:,0]
# labs1         = ['x, km','$H_s$, m']
labs1         = ['','$H_s$, m']
labs2         = ['$x$, km',r'$\tau_x$, Pa']

fig    = plt.figure()
ax1    = fig.add_subplot(2,1,1)
ax2    = fig.add_subplot(2,1,2)
# 
if do_legend:
    lines_h      = []
    lines_t      = []
    text_leg    = []

if cmp_steady:
    # get steady state solution
    om             = 2*np.pi/Tp_in
    atten_in     = np.array([h_in,om,young,visc_rp])
    atten_out    = Mwim.atten_youngs(atten_in)
    alp            = c_in/D_in*atten_out[4]    # scattering "attenuation" [m^{-1}]
    alp_dis      = 2*c_in*atten_out[0]        # damping [m^{-1}]
    kwtr          = atten_out[2]
    cp             = om/kwtr # phase vel (open water) [m/s]
    cg             = cp/2.    # group vel (open water, inf depth relation) [m/s]
    print('alp_scat = '+str(alp)+'; alp_dis = '+str(alp_dis))
    #
    N              = pow(2,7);
    nx0            = 500
    x_ice         = np.linspace(0.,ice_width,nx0)

    # solve in position space
    out    = Fbs.solve_boltzmann(width=ice_width,
                alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs_in)
    #
    E_th          = Fbs.calc_expansion(out,x_ice,L=ice_width)
    dtheta        = 2*np.pi/float(N)
    E0             = dtheta*E_th.sum(1)                             # integral over directions (freq spec)
    Hs_steady    = 4*np.sqrt(E0.real)
    #
    S_th          = E_th.dot(out['solution']['Rmat'].transpose())
    S_cos         = S_th.dot(dtheta*np.cos(out['angles']))  # integral with cos over directions
    rhow          = 1025 # kg/m^3
    gravity      = 9.81 # m/s^2
    tx_steady    = -(rhow*gravity/cp)*S_cos.real

    # plot Hs
    pobj  = Fplt.plot_1d(1.e-3*(xe+x_ice),Hs_steady,pobj=[fig,ax1],plot_steps=False,\
                labs=labs1,color='y',linewidth=3)
    lines_h.append(pobj[-1])
    
    # plot taux
    pobj  = Fplt.plot_1d(1.e-3*(xe+x_ice),tx_steady,pobj=[fig,ax2],plot_steps=False,\
                labs=labs1,color='y',linewidth=3)
    lines_t.append(pobj[-1])
    text_leg.append('Steady')

    if cmp_steady_num:
        # compare to numerical scheme from matlab
        # NB needs to be run for the same parameters
        w2d    = os.getenv('WIM2D_PATH')
        tdir  = w2d+'/matlab/boltzmann/out/'
        tfil  = tdir+'boltzmann_steady_numeric.dat'
        #
        t_out,t_prams  = Fdat.read_datfile(tfil)
        xm                 = t_out['x'] .data
        Hm                 = t_out['Hs'].data
        ax1.plot((xe+xm)/1.e3,Hm,'.k')

pdir      = bindir+'/prog/'
pfiles    = os.listdir(pdir)
afiles    = []
psteps    = []
psteps_i = []
tsteps    = []

for pf in pfiles:
    if pf[-2:]=='.a':
        afiles.append(pf)
    elif pf[-2:]=='.b':
        binfo,vlist = Fdat.fn_bfile_info(pdir+pf)
        tsteps.append(binfo['t_out'])
        stepno    = pf.strip('.b').strip('wim_prog')
        psteps.append(stepno)
        psteps_i.append(int(stepno))

# # get approx time step from log file - now get from .b file
# logfil    = outdir+'/log/wim2d.log'
# fid        = open(logfil)
# lines     = fid.readlines()
# fid.close()
# dt         = float(lines[25].split()[-1]) # time step (s)
# psteps_i = np.array(psteps_i)
# tsteps    = dt*psteps_i

def find_nearest(nparr,val):
    R              = list(abs(nparr-val))
    Rs             = sorted([(r,i) for i,r in enumerate(R)])
    Rval,ival    = Rs[0]
    return int(ival)

################################################ 
Nprog = np.max(np.array(psteps_i))
if 0:
    if 1:
        # compare every stp steps
        stp            = 80
        steps2plot  = range(0,Nprog,stp)
        steps2plot.extend([Nprog])
        times2plot  = []
        for step in steps2plot:
            ip = steps2plot.index(step)
            times2plot.append(tsteps[ip])
    else:
        # just compare final prog file
        steps2plot  = [Nprog]
        times2plot  = [tsteps[-1]]
else:
    # specify the times manually
    times = [0.,1.,1.5,2.,2.5,3.] # times in h
    times.extend([6.,12.,18.,24.,48,72.])
    steps2plot  = []
    times2plot  = []
    tsteps        = np.array(tsteps)
    for tval in times:
        i_t    = find_nearest(tsteps,tval*3600)
        steps2plot.append(psteps_i[i_t])
        times2plot.append(tsteps[i_t])
################################################ 


cols      = ['k','b','r','g','m','c']
lstil     = ['-','--','-.',':']
Nc         = len(cols)
Ns         = len(lstil)
loop_c    = -1
loop_s    = -1


print('\n')
print('prog files in '+outdir+'\n')
for i,nstep in enumerate(steps2plot):
    out_fields  = Fdat.fn_check_prog(outdir,int(nstep)) # load ice/wave conditions from binaries
    ny = grid_prams['ny']
    if ny>1:
        Hs_n  = out_fields['Hs']    [:,int(np.ceil(ny/2.))]
        tx_n  = out_fields['taux'] [:,int(np.ceil(ny/2.))]
    else:
        Hs_n  = out_fields['Hs'][:,0]
        tx_n  = out_fields['taux'][:,0]

    #
    loop_c    = np.mod(loop_c+1,Nc)
    if loop_c==0:
        loop_s    = np.mod(loop_s+1,Ns)
    
    # plot and set up for legend
    pobj  = Fplt.plot_1d(xx,Hs_n,labs=labs1,pobj=(fig,ax1),\
                                color=cols[loop_c],linestyle=lstil[loop_s],linewidth=2)
    if do_legend:
        lines_h.append(pobj[-1])

    pobj  = Fplt.plot_1d(xx,tx_n,labs=labs2,pobj=(fig,ax2),\
                                color=cols[loop_c],linestyle=lstil[loop_s],linewidth=2)
    tt     = times2plot[i]/3600.
    print('n, t(h): '+str(nstep)+', '+str(tt))
    if do_legend:
        lines_t.append(pobj[-1])
        ss = '%5.1fh' % (tt)
        text_leg.append(ss)

ax1.set_xlim([0,300])
ax2.set_ylim([0,0.3])
ax2.set_xlim([0,300])

# legends
if do_legend:
    if len(lines_h)>0:
        ax1.legend(lines_h,text_leg,loc='upper right', bbox_to_anchor=(1.3, .78))
    # if len(lines_t)>0:
    #     ax2.legend(lines_t,text_leg)

# figname  = figdir+'/convergence2steady.png'
figname  = figdir+'/convergence2steady.eps'
print('\nSaving to '+figname+'...')
fig.savefig(figname,bbox_inches='tight',pad_inches=0.05)
plt.close(fig)

if 1:
    # print  time-dep results to dat-file
    dfil  = figdir+'/test_steady1.dat'
    blk    = 4*' '
    fid    = open(dfil,'w')

    # header:
    lin    = []
    lin.append('# Time-dependant results (end of simulation)\n')
    lin.append('# Time (h): '+str(duration_hours)+'\n')
    lin.append('# x (m), Hs (m), taux, Pa'+'\n')
    lin.append('#####################################################\n')
    lin.append('\n')

    for loop_h in range(len(lin)):
        fid.write(lin[loop_h])
    
    # results:
    for loop_x in range(len(xx)):
        lin    =              ('%f' % (1.e3*xx[loop_x])  )
        lin    = lin+blk + ('%f' % (Hs_n[loop_x])      )
        lin    = lin+blk + ('%f' % (tx_n[loop_x])      )
        fid.write(lin+'\n')
    fid.close()
    print('saving to '+dfil+'...\n')

if cmp_steady:
    # print  steady-state results to dat-file
    dfil  = figdir+'/test_steady2.dat'
    fid    = open(dfil,'w')

    # header:
    lin    = []
    lin.append('# Steady-state results\n')
    lin.append('# x (m), Hs (m)'+'\n')
    lin.append('#####################################################\n')
    lin.append('\n')

    for loop_h in range(len(lin)):
        fid.write(lin[loop_h])

    # results:
    for loop_x in range(len(x_ice)):
        lin    =              ('%f' % (xe+x_ice[loop_x])  )
        lin    = lin+blk + ('%f' % (Hs_steady[loop_x])      )
        fid.write(lin+'\n')
    fid.close()
    print('saving to '+dfil+'...')
