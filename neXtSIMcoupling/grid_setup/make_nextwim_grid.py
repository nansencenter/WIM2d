## make_nextwim_grid.py
## Author: Timothy Williams
## Date: 20160916
import matplotlib
# matplotlib.use('Agg')
import os,sys
import numpy as np
import pyproj
import struct
from matplotlib import pyplot as plt

# ===============================================================
# set up projection
nsd      = os.getenv('NEXTSIMDIR')+'/data'
# nsd      = '.'
mppfile  = nsd+'/NpsNextsim.mpp'
mf       = open(mppfile)
lines    = mf.readlines()
mf.close()

TEST_PLOT   = 0
USE_LL      = 1
TEST_PROJ   = 0 # test projection against mapxy and exit
LARGE_GRID  = 0

# shape of earth
ecc   = float(lines[-1].split()[0])
a     = 1e3*float(lines[-2].split()[0]) # convert to m
b     = a*np.sqrt(1-pow(ecc,2))

# stere info
lat0,lon0,lat_ts  = lines[1].split()[:3]
lon0              = float(lon0)
lat0              = float(lat0)
lat_ts            = float(lat_ts)

rotation = float(lines[2].split()[0])

print(a,b,ecc)
print(rotation,lat0,lat_ts)
mapx  = pyproj.Proj(proj='stere',a=a,b=b,\
               lon_0=rotation,lat_0=lat0,lat_ts=lat_ts)

if TEST_PROJ:
   # test projection against mapll/mapxy
   # - works for north-pole-centered projections,
   # but maybe not general centers & rotations
   # (TODO may need to add rotation manually - if we can get its definition!)
   x  = np.linspace(-20e3,20e3,41)
   y  = 12.e3
   print("x   y   lon   lat\n")
   for xi in x:
      lon,lat  = mapx(xi,y,inverse=True)
      ss       = "%f   %f   %f   %f\n" %(xi/1.e3,y/1.e3,lon,lat)
      print(ss)
   sys.exit()
# ===============================================================


# ===============================================================
# get limits of grid
# - NB outer cells only used for boundary conditions
# (fixed or free)
res      = 2e3;# resolution in m
vertices = []
if not LARGE_GRID:
   # too small?
   figname  = "test_grid_small"
   gridname = "wim_grid_full_ONR_Oct2015_%ikm_small" %(int(res/1.e3))
   vertices.append(( -153.7   ,71.5 ))
   vertices.append(( -157     ,74.1 ))
   vertices.append(( -149.4   ,74.8 ))
   vertices.append(( -147.2   ,72   ))
else:
   figname  = "test_grid_big"
   gridname = "wim_grid_full_ONR_Oct2015_%ikm_big" %(int(res/1.e3))
   vertices.append(( -160  ,72.5 ))
   vertices.append(( -160  ,77 ))
   vertices.append(( -143  ,77 ))
   vertices.append(( -143  ,72.5 ))

xc = []
yc = []
for lon,lat in vertices:
   x,y   = mapx(lon,lat)
   xc.append(x)
   yc.append(y)
#    Lon,Lat=mapx(x,y,inverse=True)
#    print(lon,lat)
#    print(Lon,Lat)
#    print('\n')
# sys.exit()

xmin  = np.min(xc)
xmax  = np.max(xc)
ymin  = np.min(yc)
ymax  = np.max(yc)
Dx    = xmax-xmin
Dy    = ymax-ymin

# reduce xmax,ymax to give correct resolution
nx    = int(np.floor(Dx/res))
ny    = int(np.floor(Dy/res))
xmax  = xmin+nx*res
ymax  = ymin+ny*res
Dx    = xmax-xmin
Dy    = ymax-ymin

print('Resolution (km): '+str(res/1.e3))
print('x range (km)   : '+str(Dx /1.e3))
print('y range (km)   : '+str(Dy /1.e3))
print('nx             : '+str(nx))
print('ny             : '+str(ny))
# ===============================================================


# ===============================================================
def sav2bin(aid,arr,fields,nbytes=8):
   key   = arr.keys()[0]
   fields.append(key)

   X     = arr[key]
   nX    = X.size

   if nbytes==8:
      ss = nX*'d'
   else:
      ss = nX*'f'

   # reshape to fortran order
   A  = X.reshape((nX),order='F')
   aid.write(struct.pack(ss,*A))
   
   return
# ===============================================================


# ===============================================================
def write_bfile(bid,fields,nx,ny,nbytes=8):

   nrecs = len(fields)

   bid.write("%2.2i       Nrecs    # Number of records\n" %(nrecs))
   bid.write("01       Norder   # Storage order [column-major (F/matlab) = 1, row-major (C) = 0]\n")
   bid.write("%4.4i     nx       # Record length in x direction (elements)\n" %(nx))
   bid.write("%4.4i     ny       # Record length in y direction (elements)\n" %(ny))
   bid.write("%2.2i       nbytes       # bytes per entry\n" %(nbytes))
   bid.write("\n")
   bid.write("Record number and name:\n")

   """
   Rest of file:
   01       qlon
   02       qlat
   03       plon
   04       plat
   05       ulon
   06       ulat
   07       vlon
   08       vlat
   09       scuy
   10       scvx
   11       LANDMASK
   """

   for recno,vname in enumerate(fields):
      ss = "%2.2i       %s\n" %(recno+1,vname)
      bid.write(ss)

   return
# ===============================================================
   
# ===============================================================
# write to binary
aid      = open(gridname+'.a','wb')
fields   = [] # list of variable names

# corners: (nx+1)*(ny+1)
# X increases down rows, Y increases along columns
qx    = np.linspace(xmin,xmax,nx+1)
qy    = np.linspace(ymin,ymax,ny+1)
qY,qX = np.meshgrid(qy,qx)
if not USE_LL:
   sav2bin(aid,{'qx':qX},fields)
   sav2bin(aid,{'qy':qY},fields)
else:
   qlon,qlat   = mapx(qX,qY,inverse=True)
   sav2bin(aid,{'qlon':qlon},fields)
   sav2bin(aid,{'qlat':qlat},fields)

# centres: nx*ny
px    = .5*(qx[1:]+qx[:-1])
py    = .5*(qy[1:]+qy[:-1])
pY,pX = np.meshgrid(py,px)
if not USE_LL:
   sav2bin(aid,{'px':pX},fields)
   sav2bin(aid,{'py':pY},fields)
else:
   plon,plat   = mapx(pX,pY,inverse=True)
   sav2bin(aid,{'plon':plon},fields)
   sav2bin(aid,{'plat':plat},fields)

# side points: (nx+1)*ny
uY,uX = np.meshgrid(py,qx)
if not USE_LL:
   sav2bin(aid,{'ux':uX},fields)
   sav2bin(aid,{'uy':uY},fields)
else:
   ulon,ulat   = mapx(uX,uY,inverse=True)
   sav2bin(aid,{'ulon':ulon},fields)
   sav2bin(aid,{'ulat':ulat},fields)

# up/down points: nx*(ny+1)
vY,vX = np.meshgrid(qy,px)
if not USE_LL:
   sav2bin(aid,{'vx':vX},fields)
   sav2bin(aid,{'vy':vY},fields)
else:
   vlon,vlat   = mapx(vX,vY,inverse=True)
   sav2bin(aid,{'vlon':vlon},fields)
   sav2bin(aid,{'vlat':vlat},fields)

# ================================================
# grid cell sizes
# TODO should these be distance on the sphere?
scuy     = np.zeros((nx,ny))+res
sav2bin(aid,{'scuy':scuy},fields)

scvx     = np.zeros((nx,ny))+res
sav2bin(aid,{'scvx':scvx},fields)
# ================================================

# land mask
landmask = np.zeros((nx,ny))
# landmask[:,-1] = 1 # uncomment for testing
sav2bin(aid,{'LANDMASK':landmask},fields)

aid.close()

bid   = open(gridname+'.b','w')
write_bfile(bid,fields,nx,ny)
bid.close()
print("\n grid-files "+gridname+'.[a,b] saved\n')
# ===============================================================

# ===============================================================

if TEST_PLOT:
   # plot grid outline on AMSR2 conc data
   from mpl_toolkits.basemap import Basemap
   import mod_reading as mr
   wmsc  = '/work/shared/nersc/msc'
   bmap  = None
   if 0:
      figname += '_AMSR2.png'
      ddir     = wmsc+'/AMSR2_3125'
      ncfil    = ddir+'/Arc_2015/Arc_20151012_res3.125_pyres.nc'
      nci      = mr.nc_getinfo(ncfil,lonlat_file=ddir+'/LongitudeLatitudeGrid_3.125km_Arctic.nc')
      po,bmap  = nci.plot_var('fice',HYCOMreg='beau',show=False)
   elif 1:
      figname += '_OSISAF.png'
      #ddir     = wmsc+'/OSI-SAF/2015_nh_polstere'
      ddir     = 'OSI-SAF'
      ncfil    = ddir+'/ice_conc_nh_polstere-100_multi_201510121200.nc'
      nci      = mr.nc_getinfo(ncfil)
      po,bmap  = nci.plot_var('fice',HYCOMreg='beau',show=False)

   # plot new selection
   xyVerts  = []
   xyVerts.append((xmin,ymin))
   xyVerts.append((xmin,ymax))
   xyVerts.append((xmax,ymax))
   xyVerts.append((xmax,ymin))

   Verts = []
   for x,y in xyVerts:
      lon,lat  = mapx(x,y,inverse=True)
      Verts.append((lon,lat))
      # print(x,y)
      # print(lon,lat)
      # X,Y   = mapx(lon,lat)
      # print(X,Y)
      # print('\n')

   Verts.append(Verts[0])
   print(Verts)
   Vlon,Vlat   = np.array(Verts).transpose()

   if bmap is not None:
      # plot original selection
      vertices.append(vertices[0])
      vlon,vlat   = np.array(vertices).transpose()
      bmap.plot(vlon,vlat,'m',linewidth=4,ax=po.ax,latlon=True)

      bmap.plot(Vlon,Vlat,'c',linewidth=2,ax=po.ax,latlon=True)
      
      # buoy coord:
      bmap.plot(-150.6,72.76,'^r',ax=po.ax,latlon=True)

      po.fig.savefig(figname)
      # plt.show(po.fig)
