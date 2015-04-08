import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
# from matplotlib.font_manager import FontProperties

import fns_get_data as Fdat

#######################################################################
def plot_1d(x,y,labs=None,f=None,**kwargs):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)

   # fontname = 'cursive'
   fontname = 'serif'
   # fontname = 'sans-serif'

   if f is None:
      # no figure open so open a new one
      f  = plt.figure()

   ax = plt.plot(x,y,**kwargs)
   if not(labs is None):
      xl = plt.xlabel(labs[0], fontsize=16)
      xl.set_fontname(fontname)
      yl = plt.ylabel(labs[1], fontsize=16)
      yl.set_fontname(fontname)

#  # fonts of axes:
#  for tick in ax.axes.xaxis.get_major_ticks():
#     tick.label.set_fontsize(14)
#     tick.label.set_fontname(fontname)
#  for tick in ax.axes.yaxis.get_major_ticks():
#     tick.label.set_fontsize(14)
#     tick.label.set_fontname(fontname)

   return f
#######################################################################

#######################################################################
def cmap_3d_V1d(x,y,z,labs,ADD_CONTS=1,fmt='%4.1f'):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)

   # fontname = 'cursive'
   fontname = 'serif'
   # fontname = 'sans-serif'

   f  = plt.figure()
   ax = plt.pcolor(x,y,z)
   xl = plt.xlabel(labs[0], fontsize=16)
   xl.set_fontname(fontname)
   yl = plt.ylabel(labs[1], fontsize=16)
   yl.set_fontname(fontname)

   #colorbar:
   #cbar = plt.colorbar(ax, extend='neither', spacing='proportional',
                   #orientation='vertical', format="%4.2f")
   
   cbar  = plt.colorbar(ax)#,ticks=np.arange(0,1+dc,dc))
   cbar.set_label(labs[2], size=14)
   cbar.ax.tick_params(labelsize=14)
   #plt.locator_params(nbins=4)

   ####################################################################
   if ADD_CONTS is 1:
      # add contours:
      levels0  = cbar.ax.get_yticklabels()
      levels   = []
      for j in range(1,len(levels0)-2,2):
      # for j in range(1,len(levels0),1):
         lev   = levels0[j]
         levels.append(float(lev.get_text()))

      CS = plt.contour(x,y,z,levels=levels,
                        colors='k',linestyles='--',linewidths=1)
      plt.clabel(CS,inline=1,fontsize=10,fmt=fmt)
   ####################################################################

   tick_locator = ticker.MaxNLocator(nbins=5)
   cbar.locator = tick_locator
   cbar.update_ticks()

   # fonts of axes:
   for tick in ax.axes.xaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)
   for tick in ax.axes.yaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)

   return f
#######################################################################

#######################################################################
def cmap_3d(x,y,z,labs):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)

   # fontname = 'cursive'
   fontname = 'serif'
   # fontname = 'sans-serif'

   f  = plt.figure()
   ax = plt.pcolor(x,y,z)
   ax.axes.set_aspect('equal')
   xl = plt.xlabel(labs[0], fontsize=16)
   xl.set_fontname(fontname)
   yl = plt.ylabel(labs[1], fontsize=16)
   yl.set_fontname(fontname)

   ############################################################
   #colorbar:
   if z.max()!=z.min():
      # only have colorbar if not constant
      # - doesn't work on linux otherwise
      cbar  = plt.colorbar(ax)#,ticks=np.arange(0,1+dc,dc))
      cbar.set_label(labs[2], size=14)
      cbar.ax.tick_params(labelsize=14)
      #plt.locator_params(nbins=4)
      #
      cpos     = cbar.ax.get_position().extents
      cpos[2]  = cpos[2]+.15  # cbar width
      cpos[1]  = cpos[1]+.21  # lower height
      cpos[3]  = cpos[3]*.38  # colorbar height
      cbar.ax.set_position(cpos)         

      tick_locator = ticker.MaxNLocator(nbins=5)
      cbar.locator = tick_locator
      cbar.update_ticks()
   ############################################################

   # fonts of axes:
   for tick in ax.axes.xaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)
   for tick in ax.axes.yaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)

   return f
#######################################################################

#######################################################################
def fn_plot_init(grid_prams,ice_fields,wave_fields,figdir):

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   x  = grid_prams['X']/1.0e3
   y  = grid_prams['Y']/1.0e3
   nx = grid_prams['nx']
   ny = grid_prams['ny']

   # ice fields
   # dictionary of figure names
   figs   = {'icec':'icec.png','iceh':'iceh.png','dfloe':'Dmax0.png'}
   # dictionary of labels for colorbars
   labs   = {'icec':'$c$'     ,'iceh':'$h$, m'  ,'dfloe':'$D_{max}$, m'}
   keys   = labs.keys()

   for key in keys:
      fig   = figdir+figs[key]
      if ny>1:
         f  = cmap_3d(x,y,ice_fields[key],
                  ['$x$, km','$y$, km',labs[key]])
      else:
         f  = plot_1d(x,ice_fields[key],
                  ['$x$, km',labs[key]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()

   # wave fields
   # dictionary of figure names
   figs     = {'Hs':'Hs0.png'    ,'Tp':'Tp0.png'   ,'mwd':'mwd0.png'}
   # dictionary of labels for colorbars
   labs  = {'Hs':'$H_{s}$, m' ,'Tp':'$T_p$, s'  ,'mwd':'mwd, degrees'}
   keys   = labs.keys()

   for key in keys:
      fig   = figdir+'/'+figs[key]
      if ny>1:
         f  = cmap_3d(x,y,wave_fields[key],
               ['$x$, km','$y$, km',labs[key]])
      else:
         f  = plot_1d(x,wave_fields[key],
               ['$x$, km',labs[key]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################

############################################################################
def fn_plot_final(grid_prams,out_fields,figdir):

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   keys  = out_fields.keys()
   x     = grid_prams['X']/1.0e3
   y     = grid_prams['Y']/1.0e3
   nx    = grid_prams['nx']
   ny    = grid_prams['ny']

   # dictionary of figure names
   figs  = {'dfloe':'Dmax.png','taux':'taux.png','tauy':'tauy.png',
            'Hs':'Hs.png','Tp':'Tp.png'}

   # dictionary of labels for colorbars
   labs  = {'dfloe':'$D_{max}$, m','taux':'Stress (x dir.), Pa','tauy':'Stress (y dir.), Pa',
            'Hs':'$H_s$, m','Tp':'$T_p$, s'}

   for key in keys:
      fig   = figdir+'/'+figs[key]
      if ny>1:
         f     = cmap_3d(x,y,out_fields[key],
                  ['$x$, km','$y$, km',labs[key]])
      else:
         f     = plot_1d(x,out_fields[key],
                  ['$x$, km',labs[key]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################

############################################################################
def fn_plot_final_V1d(grid_prams,Tp_vec,out_fields,figdir):
   # plot against T_p instead of y

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   # keys  = out_fields.keys()
   keys  = ['dfloe','taux','Hs']
   x     = grid_prams['X']/1.0e3

   nx,ny = x.shape
   y     = np.ones((nx,ny))
   for j in range(0,ny):
      y[:,j]   = Tp_vec[j]*y[:,j]

   # dictionary of figure names
   figs  = Fdat.Out_Fields('Dmax.png','taux.png','tauy.png',
                           'Hs.png','Tp.png',)
   # dictionary of labels for colorbars
   labs  = Fdat.Out_Fields('$D_{max}$, m',
                           'Stress (x dir.), Pa','Stress (y dir.), Pa',
                           '$H_s$, m','$T_p$, s')
   AC    = Fdat.Out_Fields(0,0,0,1,0)
   FMT   = Fdat.Out_Fields(0,0,0,'%2.1f',0)

   for key in keys:
      fig   = fdir+figs[key]
      f     = cmap_3d_V1d(x,y,out_fields[key],
               ['$x$, km','$T_p$, s',labs[key]],
               ADD_CONTS=AC[key],fmt=FMT[key])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################
