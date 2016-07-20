from matplotlib import pyplot as plt
import numpy as np
import os

for loop_i in range(2):
   if loop_i==0:
      figname  = 'fig_egHsVsX_steady.png'
      dname    = 'test_steady2.dat'
      Nh       = 4
   elif loop_i==1:
      figname  = 'fig_egHsVsX_steady_FT.png'
      dname    = 'test_steady2_FT.dat'
      Nh       = 4
   else:
      figname  = 'fig_egHsVsX_timedep.png'
      dname    = 'test_steady1.dat'
      Nh       = 5

   if os.path.exists(dname):
      print(' ')
      print('Opening '+dname)
      fid      = open(dname,'r')

      lines    = fid.readlines()
      fid.close()

      Nlines   = len(lines)
      Nx       = Nlines-Nh
      xx       = []
      yy       = []

      x0 = -220.e3
      x1 = -120.e3

      for i in range(Nx):
         line  = lines[Nh+i].split()
         x     = float(line[0])
         y     = float(line[1])
         if (x0<=x) and (x<=x1):
            xx.append((x-x0)/1e3)
            yy.append(y)

      xx = np.array(xx)
      yy = np.array(yy)

      plt.plot(xx,yy/yy[0])
      plt.xlabel('x, km')
      plt.ylabel('$H_s(x)/H_s(0)$')
      plt.savefig(figname)
      plt.close()

      print('Saving fig to '+figname)
      print(' ')
