import os,sys
import matplotlib.pyplot as plt

sys.path.append('../py_funs')
import fns_get_data as Fdat
import fns_plot_data as Fplt

if 1:
   gdir  = '../run/inputs/'
else:
   gdir  = 'test/out/'

gf       = Fdat.fn_check_grid(gdir)
fields   = {'LANDMASK':gf['LANDMASK']}
Fplt.fn_plot_gen(gf,fields,'test/out/')

figname = 'land_mask.png'
print('Saved to '+figname)
