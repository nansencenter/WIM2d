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
def cmap_3d(x,y,z,labs,pobj=None,zlims=None):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)

   # fontname = 'cursive'
   fontname = 'serif'
   # fontname = 'sans-serif'

   ################################################
   if zlims is not None:
      vmin  = zlims[0]
      vmax  = zlims[1]
   else:
      vmax  = z.max()
      vmin  = z.min()
      dv    = vmax-vmin

      if dv==0:
         # z is const
         vv    = z.mean()
         vmin  = vmin-.1*vv
         vmax  = vmax+.1*vv
      else:
         vmin  = vmin-.1*dv
         vmax  = vmax+.1*dv
   ################################################

   ################################################
   if pobj is None:
      f     = plt.figure()
      ax    = f.add_subplot(1,1,1)
      pobj  = f,ax
   else:
      f,ax  = pobj
   ################################################

   P  = ax.pcolor(x,y,z,vmin=vmin,vmax=vmax)
   ax.axes.set_aspect('equal')
   ax.set_xlabel(labs[0], fontsize=16,fontname=fontname)
   ax.set_ylabel(labs[1], fontsize=16,fontname=fontname)
   ax.set_xlim((x[0],x[-1]))
   ax.set_ylim((y[0],y[-1]))

   ############################################################
   #colorbar:
   if z.max()!=z.min():
      # only have colorbar if not constant
      # - doesn't work on linux otherwise
      cbar  = f.colorbar(P)#,ticks=np.arange(0,1+dc,dc))
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
   for tick in P.axes.xaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)
   for tick in P.axes.yaxis.get_major_ticks():
      tick.label.set_fontsize(14)
      tick.label.set_fontname(fontname)

   return pobj
#######################################################################

#######################################################################
def fn_plot_gen(grid_prams,fields,figdir):

   from matplotlib import pyplot as plt

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   nx = grid_prams['nx']
   ny = grid_prams['ny']
   dx = grid_prams['dx']
   dy = grid_prams['dy']
   x  = grid_prams['X'][:,0]+dx/2.
   x  = np.concatenate([[x[0]-dx],x])/1.e3
   y  = grid_prams['Y'][0,:]+dy/2.
   y  = np.concatenate([[y[0]-dy],y])/1.e3

   # dictionary of figure names
   figs   = {'icec':'icec.png','iceh':'iceh.png','dfloe':'Dmax.png',\
             'taux':'taux.png','tauy':'tauy.png',\
             'Hs':'Hs.png'    ,'Tp':'Tp.png'   ,'mwd':'mwd.png',\
             'ICE_MASK':'ice_mask.png','WAVE_MASK':'wave_mask.png'}

   # dictionary of labels for colorbars
   labs   = {'icec':'$c$'     ,'iceh':'$h$, m'  ,'dfloe':'$D_{max}$, m',\
             'taux':r'$\tau_x$, Pa','tauy':r'$\tau_y$, Pa',\
             'Hs':'$H_{s}$, m' ,'Tp':'$T_p$, s'  ,'mwd':'mwd, degrees',\
             'ICE_MASK':'Ice mask','WAVE_MASK':'Wave mask'}

   # allow for other variations of key names
   aliases  = {'Dmax':'dfloe',\
               'cice':'icec',\
               'hice':'iceh',\
               'tau_x':'taux',\
               'tau_y':'tauy'}
   for key in aliases.keys():
      figs.update({key:figs[aliases[key]]})
      labs.update({key:labs[aliases[key]]})

   # make plots
   for key in fields.keys():
      fig   = figdir+'/'+figs[key]
      if ny>1:
         f,ax  = cmap_3d(x,y,fields[key].transpose(),\
                        ['$x$, km','$y$, km',labs[key]])
      else:
         f,ax  = plot_1d(x,fields[key],\
                        ['$x$, km',labs[key]])
      f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close(f)
############################################################################

#######################################################################
def fn_plot_init(grid_prams,ice_fields,wave_fields,figdir):

   from matplotlib import pyplot as plt

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   nx = grid_prams['nx']
   ny = grid_prams['ny']
   dx = grid_prams['dx']
   dy = grid_prams['dy']
   x  = grid_prams['X'][:,0]+dx/2.
   x  = np.concatenate([[x[0]-dx],x])/1.e3
   y  = grid_prams['Y'][0,:]+dy/2.
   y  = np.concatenate([[y[0]-dy],y])/1.e3

   # ice fields
   # dictionary of figure names
   figs   = {'icec':'icec.png','iceh':'iceh.png','dfloe':'Dmax.png'}
   # dictionary of labels for colorbars
   labs   = {'icec':'$c$'     ,'iceh':'$h$, m'  ,'dfloe':'$D_{max}$, m'}

   for key in labs.keys():
      fig   = figdir+figs[key]
      if ny>1:
         f,ax  = cmap_3d(x,y,ice_fields[key].transpose(),\
                        ['$x$, km','$y$, km',labs[key]])
      else:
         f,ax  = plot_1d(x,ice_fields[key],\
                        ['$x$, km',labs[key]])
      f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close(f)

   # wave fields
   # dictionary of figure names
   figs     = {'Hs':'Hs.png'    ,'Tp':'Tp.png'   ,'mwd':'mwd.png'}
   # dictionary of labels for colorbars
   labs  = {'Hs':'$H_{s}$, m' ,'Tp':'$T_p$, s'  ,'mwd':'mwd, degrees'}

   for key in labs.keys():
      fig   = figdir+'/'+figs[key]
      if ny>1:
         f,ax  = cmap_3d(x,y,wave_fields[key].transpose(),\
                        ['$x$, km','$y$, km',labs[key]])
      else:
         f,ax  = plot_1d(x,wave_fields[key],\
               ['$x$, km',labs[key]])
      f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      ax.cla()
      f.clf()
############################################################################

############################################################################
def fn_plot_final(grid_prams,out_fields,figdir):

   from matplotlib import pyplot as plt

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   nx = grid_prams['nx']
   ny = grid_prams['ny']
   dx = grid_prams['dx']
   dy = grid_prams['dy']
   x  = grid_prams['X'][:,0]+dx/2.
   x  = np.concatenate([[x[0]-dx],x])/1.e3
   y  = grid_prams['Y'][0,:]+dy/2.
   y  = np.concatenate([[y[0]-dy],y])/1.e3

   # dictionary of figure names
   figs  = {'dfloe':'Dmax.png','taux':'taux.png','tauy':'tauy.png',
            'Hs':'Hs.png','Tp':'Tp.png'}

   # dictionary of labels for colorbars
   labs  = {'dfloe':'$D_{max}$, m','taux':'Stress (x dir.), Pa','tauy':'Stress (y dir.), Pa',
            'Hs':'$H_s$, m','Tp':'$T_p$, s'}

   for key in out_fields.keys():
      fig   = figdir+'/'+figs[key]
      if ny>1:
         f,ax  = cmap_3d(x,y,out_fields[key].transpose(),\
                  ['$x$, km','$y$, km',labs[key]])
      else:
         f,ax  = plot_1d(x,out_fields[key],\
                  ['$x$, km',labs[key]])
      f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close(f)
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
