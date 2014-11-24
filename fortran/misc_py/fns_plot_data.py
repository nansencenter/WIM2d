import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker

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

   fig   = fdir+'icec.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,ice_fields.icec,
            ['$x$, km','$y$, km','$c$'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   fig   = fdir+'iceh.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,ice_fields.iceh,
                   ['$x$, km','$y$, km','$h$, m'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   fig   = fdir+'Dmax0.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,ice_fields.dfloe,
                   ['$x$, km','$y$, km','$h$, m'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   fig   = fdir+'Hs0.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,wave_fields.Hs,
            ['$x$, km','$y$, km','$H_s$, m'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   fig   = fdir+'Tp0.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,wave_fields.Tp,
                   ['$x$, km','$y$, km','$T_p$, s'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   fig   = fdir+'mwd0.png'
   f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,wave_fields.mwd,
                   ['$x$, km','$y$, km','mwd$, degrees'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()
############################################################################

#######################################################################
def fn_plot_final(grid_prams,out_arrays,figdir):

   fdir  = figdir+'/final/'

   figs  = ['Dmax.png','Hs.png','Tp.png','taux.png','tauy.png']
   labs  = ['$D_{max}$','$H_s$','$T_p$','Stress (x dir.), Pa','Stress (y dir.), Pa']
   for j in range(0,5):
      fig   = fdir+figs[j]
      f     = cmap_3d(grid_prams.X/1.0e3,grid_prams.Y/1.0e3,out_arrays[:,:,j],
               ['$x$, km','$y$, km',labs[j]])
      plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      f.clf()
############################################################################
