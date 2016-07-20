import os,sys
import matplotlib.pyplot as plt

w2d   = os.getenv('WIM2D_PATH')
sys.path.append(w2d+'/fortran/py_funs')
import fns_get_data as Fdat
import fns_plot_data as Fplt

gdir     = '.'
gf       = Fdat.fn_check_grid(gdir)
fields   = {'LANDMASK':gf['LANDMASK']}
Fplt.fn_plot_gen(gf,fields,'out/')

figname = 'land_mask.png'
print('Saved to '+figname)
