import os,sys
import fns_get_data  as Fdat
from getopt import getopt

# run from root results directory eg out, out_io
outdir      = os.getcwd()
PLOT_INIT   = 1
PLOT_FINAL  = 1
PLOT_PROG   = 1

# =======================================================================
# command line inputs
opts,args   = getopt(sys.argv[1:],"",["outdir=","init=","final=","prog="])
for opt,arg in opts:
   if opt=='--outdir':
      outdir   = arg 
   if opt=='--init':
      PLOT_INIT  = int(arg)
   if opt=='--final':
      PLOT_FINAL  = int(arg)
   if opt=='--prog':
      PLOT_PROG  = int(arg)
# =======================================================================

results  = Fdat.wim_results(outdir=outdir)



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
