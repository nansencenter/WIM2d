import os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
# from matplotlib.font_manager import FontProperties

import fns_get_data as Fdat

#######################################################################
def step_1d(x,y):
   # vectors for constant-panel plotting

   nx       = len(x)
   dx       = x[1]-x[0]
   #
   xx       = np.zeros(2*nx)
   xx[::2]  = x-.5*dx
   xx[1::2] = x+.5*dx
   #
   yy       = np.zeros(2*nx)
   yy[::2]  = y
   yy[1::2] = y

   return xx,yy
#######################################################################

#######################################################################
def plot_1d(x,y,labs=None,plot_steps=True,pobj=None,**kwargs):
   # f  = plt.figure(figsize=[6,6],dpi=50)
   # plt.pcolor(grid_prams.X,grid_prams.Y,ice_fields.icec,cmap=cm.jet,vmax=Vmax,vmin=Vmin)

   # fontname = 'cursive'
   fontname = 'serif'
   # fontname = 'sans-serif'

   if pobj is None:
      # no figure open so open a new one
      f     = plt.figure()
      ax    = f.add_subplot(1,1,1)
      pobj  = f,ax
   elif len(pobj)==2:
      f,ax  = pobj
   else:
      f,ax,line = pobj

   if plot_steps:
      x,y   = step_1d(x,y)

   line ,= ax.plot(x,y,**kwargs)
   pobj  = f,ax,line

   if not(labs is None):
      ax.set_xlabel(labs[0], fontsize=16,fontname=fontname)
      ax.set_ylabel(labs[1], fontsize=16,fontname=fontname)

#  # fonts of axes:
#  for tick in ax.axes.xaxis.get_major_ticks():
#     tick.label.set_fontsize(14)
#     tick.label.set_fontname(fontname)
#  for tick in ax.axes.yaxis.get_major_ticks():
#     tick.label.set_fontsize(14)
#     tick.label.set_fontname(fontname)

   return pobj
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
def fn_plot_gen(grid_prams,fields,figdir,zlims_in=None,text='',vlist=None):

   from matplotlib import pyplot as plt

   if not os.path.exists(figdir):
      os.mkdir(figdir)

   nx = grid_prams['nx']
   ny = grid_prams['ny']
   dx = grid_prams['dx']
   dy = grid_prams['dy']
   if ny==1:
      x  = grid_prams['X']/1.e3
      if x.ndim==2:
         x  = x[:,0]
   else:
      x  = grid_prams['X'][:,0]+dx/2.
      x  = np.concatenate([[x[0]-dx],x])/1.e3
      y  = grid_prams['Y'][0,:]+dy/2.
      y  = np.concatenate([[y[0]-dy],y])/1.e3

   # dictionary of figure names
   figs   = {'icec':'icec'+text+'.png',\
             'iceh':'iceh'+text+'.png',\
             'dfloe':'Dmax'+text+'.png',\
             'taux':'taux'+text+'.png',\
             'tauy':'tauy'+text+'.png',\
             'Hs':'Hs'+text+'.png',\
             'Tp':'Tp'+text+'.png',\
             'mwd':'mwd'+text+'.png',\
             'ICE_MASK':'ice_mask'+text+'.png',\
             'WAVE_MASK':'wave_mask'+text+'.png',\
             'LANDMASK':'land_mask'+text+'.png'}

   # dictionary of labels for colorbars
   labs   = {'icec':'$c$',\
         'iceh':'$h$, m',\
         'dfloe':'$D_{max}$, m',\
         'taux':r'$\tau_x$, Pa',\
         'tauy':r'$\tau_y$, Pa',\
         'Hs':'$H_{s}$, m',\
         'Tp':'$T_p$, s',\
         'mwd':'mwd, degrees',\
         'ICE_MASK':'Ice mask',\
         'WAVE_MASK':'Wave mask',\
         'LANDMASK':'Land mask'}

   # Default is let python choose range for variables
   zlims  = {'icec'     : None,   \
             'iceh'     : None,   \
             'dfloe'    : None,   \
             'taux'     : None,   \
             'tauy'     : None,   \
             'Hs'       : None,   \
             'Tp'       : None,   \
             'mwd'      : None,   \
             'ICE_MASK' : None,   \
             'WAVE_MASK': None,   \
             'LANDMASK' : None}

   # allow for other variations of key names
   aliases  = Fdat.key_aliases()

   for key in aliases.keys():
      key2  = aliases[key]
      figs.update ({key:figs [key2]})
      labs.update ({key:labs [key2]})
      zlims.update({key:zlims[key2]})

   if zlims_in is not None:
      for key in zlims_in.keys():
         zlims[key]  = zlims_in[key]

   # make plots
   if vlist is None:
      vlist = fields.keys()

   for key in vlist:
      fig   = figdir+'/'+figs[key]
      zlim  = zlims[key]
      if ny>1:
         f,ax  = cmap_3d(x,y,fields[key].transpose(),\
                         labs=['$x$, km','$y$, km',labs[key]],\
                         zlims=zlim)
      else:
         F  = fields[key]
         if F.ndim==2:
            F  = F[:,0]
         f,ax,line  = plot_1d(x,F,labs=['$x$, km',labs[key]])
         ax.set_ylim(zlim)
      f.savefig(fig,bbox_inches='tight',pad_inches=0.05)
      plt.close(f)
   return
############################################################################

#######################################################################
def fn_plot_init(grid_prams,ice_fields,wave_fields,figdir):

   ice_fields.update(wave_fields)
   fn_plot_gen(grid_prams,ice_fields,figdir)

   return
############################################################################

############################################################################
def fn_plot_final(grid_prams,out_fields,figdir):

   fn_plot_gen(grid_prams,out_fields,figdir)

   return
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

   return
############################################################################
