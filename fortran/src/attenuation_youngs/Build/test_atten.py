import numpy as np
import sys

sys.path.append("../../bin")

import atten_young_f2py as AT

iceh     = 2.0
drag_rp  = 13.0
visc_ws  = 0.0
young    = 2.0e9
gravity  = 9.81
#
pi    = np.pi
Tmin  = 2.5
Tmax  = 25.0
nw    = 25;
#
fmin  = 1.0/Tmax
fmax  = 1.0/Tmin
df    = (fmax-fmin)/(nw-1.0)
#
for w in range(0,nw):
   freq  = fmin+w*df
   om    = 2*pi*freq
   T     = 1.0/freq

   ## initial guess for ice wavenumber
   kw_inf   = om**2/gravity
   if (w is 0):
      print(' ')
      guess = kw_inf
      print('guess='+str(guess)+'\n')

   ## get atten coeff
   inputs   = np.array([iceh,om,young,drag_rp,visc_ws])
   damping,kice,kwtr,int_adm,alp_nd,modT,argR,argT = AT.atten_young(inputs,guess)

   ## update guess for ice wavenumber
   guess = kice

   ## check attenuation coeff's, wavenumbers
   if (w is 0):
      print('***********************************************')
      print('check outputs from RTparam_outer:')
      print('atten : '+'%f,%f' %(alp_nd,damping))
      print('ki,kw,2pi/wlng_wtr/guess : '+ '%f,%f,%f,%f'%(kice,kwtr,kw_inf,guess))
      print('|T|,argRT_s : '+'%f,%f,%f,%f' %(modT,argR,argT,int_adm))
      print('***********************************************')
      print(' ')

   print('T (s), atten (per floe), damping (/m) : '+ '%f,%e,%e'%(T,alp_nd,damping))
