import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker

import fns_get_data as Mdat

#######################################################################
def cmap_3d(x,y,z,labs):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)
   f  = plt.figure()
   ax = plt.pcolor(x,y,z)
   ax.axes.set_aspect('equal')
   xl = plt.xlabel(labs[0], fontsize=16)
   xl.set_fontname('times')
   yl = plt.ylabel(labs[1], fontsize=16)
   yl.set_fontname('times')

   #colorbar:
   #cbar = plt.colorbar(ax, extend='neither', spacing='proportional',
                   #orientation='vertical', format="%4.2f")
   
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

   # fonts of axes:
   for tick in ax.axes.xaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname('times')
   for tick in ax.axes.yaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname('times')

   return f
#######################################################################

#######################################################################
def fn_plot_init(grid_prams,ice_fields,wave_fields,figdir):

   fdir  = figdir+'/init/'
   x     = grid_prams['X']/1.0e3
   y     = grid_prams['Y']/1.0e3

   # ice fields
   figs   = Mdat.Ice_Fields(icec='icec.png',iceh='iceh.png',dfloe='Dmax0.png')
   labs   = Mdat.Ice_Fields(icec='$c$',iceh='$h$, m',dfloe='$D_{max}$, m')
   keys   = ['icec','iceh','dfloe']

   for key in keys:
      fig   = fdir+figs[key]
      f     = cmap_3d(x,y,ice_fields[key],
                      ['$x$, km','$y$, km','$c$'])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()

   # wave fields
   figs  = Mdat.Wave_Fields(Hs='Hs0.png',Tp='Tp0.png',mwd='mwd0.png')
   labs  = Mdat.Wave_Fields('$H_{s}$, m','$T_p$, s','mwd, degrees')
   keys   = ['Hs','Tp','mwd']

   for key in keys:
      fig   = fdir+figs[key]
      f     = cmap_3d(x,y,wave_fields[key],
                      ['$x$, km','$y$, km',labs[key]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################

############################################################################
def fn_plot_final(grid_prams,out_fields,figdir):

   fdir  = figdir+'/final/'
   keys  = out_fields.keys()
   x     = grid_prams['X']/1.0e3
   y     = grid_prams['Y']/1.0e3

   # dictionary of figure names
   figs  = Mdat.Out_Fields('Dmax.png','taux.png','tauy.png',
                           'Hs.png','Tp.png',)
   # dictionary of labels for colorbars
   labs  = Mdat.Out_Fields('$D_{max}$, m',
                           'Stress (x dir.), Pa','Stress (y dir.), Pa',
                           '$H_s$, m','$T_p$, s')
   for key in keys:
      fig   = fdir+figs[key]
      f     = cmap_3d(x,y,out_fields[key],
               ['$x$, km','$y$, km',labs[key]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################
