import os,sys
import fns_get_data  as Fdat

# run from root results directory eg out, out_io
outdir   = os.getcwd()
results  = Fdat.wim_results(outdir=outdir)

PLOT_INIT   = 1
PLOT_FINAL  = 1
PLOT_PROG   = 1


##########################################################################
if PLOT_INIT:
   results.plot_initial()
##########################################################################


################################################################
if PLOT_FINAL:
   results.plot_final()
################################################################


################################################################
if PLOT_PROG:
   results.plot_prog()


   ################################################################
   wim2d_path  = os.getenv('WIM2D_PATH')

   print('\n**********************************************************************')
   print('to make movie, type')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs '+outdir+'/figs/prog')
   print('or')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Dmax '+outdir+'/figs/prog')
   print('**********************************************************************\n')
   ################################################################


################################################################
