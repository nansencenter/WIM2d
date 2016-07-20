import numpy as np
import os
import sys

dd    = os.path.abspath("../..")
dirs  = ["bin","misc_py"] # directories relative to "fortran"
for n in range(0,len(dirs)):
   dd2   = dd+"/"+dirs[n]
   print('adding path : '+dd2)
   sys.path.append(dd2)

import RTparam_outer as RT

iceh     = 2.0
visc_rp  = 13.0
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
   inputs   = np.array([iceh,om,young,visc_rp])
   damping,kice,kwtr,int_adm,alp_nd,modT,argR,argT = RT.rtparam_outer(inputs,guess)

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
