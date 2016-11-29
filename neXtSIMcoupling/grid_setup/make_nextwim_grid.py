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
from getopt import getopt

# =======================================================================
def read_grid_info_file(filename):
   f     = open(filename)
   lines = f.readlines()
   f.close()
   res      = float(lines[1].split()[0])
   gridname = lines[0].split()[0]
   vertices = []
   for lin in lines[2:]:
      lon   = float(lin.split()[0])
      lat   = float(lin.split()[1])
      vertices.append((lon,lat))
   return gridname,res,vertices
# =======================================================================


# =======================================================================
def get_depth_mask(plon,plat,**kwargs):
   # ncfil = 'ETOPO/ETOPO_Arctic_10arcmin.nc' # ~20km
   ncfil = 'ETOPO/ETOPO_Arctic_5arcmin.nc'  # ~10km
   # ncfil = 'ETOPO/ETOPO_Arctic_2arcmin.nc'  # ~4km
   # ncfil = 'ETOPO/ETOPO_Arctic_1arcmin.nc'  # ~2km

   #interp depth to plon,plat (centres of wim grid cells)
   import mod_reading as mr
   nci      = mr.nc_getinfo(ncfil)
   depth    = nci.interp2points('z', (plon,plat),**kwargs) # masked array
   landmask = np.ones(depth.shape)
   good     = np.logical_not(depth.mask)
   data     = depth.data[good]
   lmgood   = landmask[good] 

   lmgood[data<0.]   = 0. # z<0 is water
   landmask[good]    = lmgood

   return depth,landmask
# =======================================================================


# =======================================================================
# command line inputs
opts,args      = getopt(sys.argv[1:],"",["TEST_PLOT=","grid_info_file="])
grid_info_file = None
TEST_PLOT      = 1 #TODO get HYCOMreg automatically from vertices
for opt,arg in opts:
   if opt=='--grid_info_file':
      grid_info_file  = arg
   if opt=='--TEST_PLOT':
      TEST_PLOT   = int(arg)

if grid_info_file is not None:
   # get gridname, resolution (m) and 4 vertices (lon/lat) from grid info file
   gridname,res,vertices   = read_grid_info_file(grid_info_file)
   print('Grid file name : '+gridname)
   print('Resolution     : %0.1fkm' %(res/1.e3))
   print('Vertices       :')
   for v in vertices:
      print(v)
   print('\n')
else:
   raise ValueError('please specify input file with --grid_info_file=')
# =======================================================================


# =======================================================================
# determine where to save output files
tmp   = grid_info_file.split('/')
Nlev  = len(tmp)-1 # no of directory levels
if len(tmp)==0:
   outdir   = '.'
else:
   outdir   = ''
   for j in range(Nlev):
      outdir   = outdir+tmp[j]+'/'
# =======================================================================


# ===============================================================
# set up projection
nsd      = os.getenv('NEXTSIMDIR')+'/data'
# nsd      = '.'
USE_MPROJ   = 0
if not USE_MPROJ:
   mppfile  = nsd+'/NpsNextsim.mpp'
else:
   mppfile  = nsd+'/NpsMproj.mpp'

mf       = open(mppfile)
lines    = mf.readlines()
mf.close()

USE_LL      = 1
TEST_PROJ   = 0 # test projection against mapxy and exit

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
geod  = pyproj.Geod(a=a,b=b)

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


#for test plots
HYCOMregions   = {}
HYCOMregions.update({'wim_grid_FS_7km':'gre'})
HYCOMregions.update({'wim_grid_FS_4km':'gre'})
HYCOMregions.update({'wim_grid_full_ONR_Oct2015_2km_big':'beau'})
HYCOMregions.update({'wim_grid_full_ONR_Oct2015_4km_big':'beau'})
HYCOMregions.update({'wim_grid_full_ONR_Oct2015_2km_small':'beau'})
HYCOMreg = HYCOMregions[gridname]


# ===============================================================
if TEST_PLOT==1:
   ddir        = 'OSI-SAF'
   ncfil       = ddir+'/ice_conc_nh_polstere-100_multi_201510121200.nc'
   vbl         = 'fice'
   lonlat_file = None
   time_index  = 0
   dlabel      = 1
   figname     = outdir+"/test_"+gridname+"_OSISAF.png"

elif TEST_PLOT==2:
      ddir        = os.getenv('NEXTSIMDIR')+'/data'
      ncfil       = ddir+'/SWARP_WW3_ARCTIC-12K_20151217.nc'
      vbl         = 'swh'
      time_index  = 6
      dlabel      = 2
      lonlat_file = None
      figname     = outdir+"/test_"+gridname+"_WW3.png"
# ===============================================================

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

   print("Range of "+key,X.min(),X.max())

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
aid      = open(outdir+'/'+gridname+'.a','wb')
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
scuy  = np.zeros((nx+1,ny))+res
scvx  = np.zeros((nx,ny+1))+res
scpy  = np.zeros((nx,ny))+res
scpx  = np.zeros((nx,ny))+res

if USE_LL:
   # CHANGE TO DISTANCES ON THE ELLIPSOID

   # uy,py:
   # y increases in j dirn, so get distance between cols
   for j in range(ny):
      # loop over cols
      scuy[:,j]   = geod.inv(qlon[:,j+1],qlat[:,j+1],qlon[:,j],qlat[:,j])[2]
      scpy[:,j]   = geod.inv(vlon[:,j+1],vlat[:,j+1],vlon[:,j],vlat[:,j])[2]

   # vx:
   # x increases in i dirn, so get distance between rows
   for i in range(nx):
      # loop over rows
      scvx[i,:]   = geod.inv(qlon[i+1,:],qlat[i+1,:],qlon[i,:],qlat[i,:])[2]
      scpx[i,:]   = geod.inv(ulon[i+1,:],ulat[i+1,:],ulon[i,:],ulat[i,:])[2]

scp2  = scpx*scpy
sav2bin(aid,{'scuy':scuy},fields)
sav2bin(aid,{'scvx':scvx},fields)
sav2bin(aid,{'scp2':scp2},fields)
# ================================================

# land mask
USE_BATHY   = 1
if USE_BATHY:
   # use bathymetry
   # where depth is "missing", it is land
   depth,landmask = get_depth_mask(plon,plat,mapping=mapx)
else:
   landmask = np.zeros((nx,ny))

# landmask[:,-1] = 1 # uncomment for testing
sav2bin(aid,{'LANDMASK':landmask},fields)

aid.close()

bid   = open(outdir+'/'+gridname+'.b','w')
write_bfile(bid,fields,nx,ny)
bid.close()
print("\n grid-files "+outdir+'/'+gridname+'.[a,b] saved\n')
# ===============================================================

# ===============================================================
if TEST_PLOT>0:
   # plot grid outline on top of data from netcdf file
   import mod_reading as mr
   import fns_plotting as Fplt
   bmap  = Fplt.start_HYCOM_map(HYCOMreg)
   print('Using '+ncfil)
   nci      = mr.nc_getinfo(ncfil,lonlat_file=lonlat_file)
   print(' - contains:')
   print(nci.variables)
   po,bmap  = nci.plot_var(vbl,bmap=bmap,show=False,time_index=time_index,date_label=dlabel,clabel=vbl)

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

   if 1:
      # plot original selection
      vertices.append(vertices[0])
      vlon,vlat   = np.array(vertices).transpose()
      bmap.plot(vlon,vlat,'m',linewidth=4,ax=po.ax,latlon=True)

      bmap.plot(Vlon,Vlat,'c',linewidth=2,ax=po.ax,latlon=True)
      
      # buoy coord:
      bmap.plot(-150.6,72.76,'^r',ax=po.ax,latlon=True)

   print('\nSaving '+figname+'\n')
   po.fig.savefig(figname)
   po.ax.cla()


   if USE_BATHY:
      # plot LANDMASK & depth
      from mpl_toolkits.axes_grid1 import make_axes_locatable
      cax   = po.cbar.ax

      for vbl in ['depth','landmask']:
         if vbl=='depth':
            PC    = bmap.pcolor(plon,plat,depth,latlon=True,ax=po.ax)
            cols  = ['b']
         else:
            PC    = bmap.pcolor(plon,plat,landmask,latlon=True,ax=po.ax)
            cols  = ['c']

         cb = po.fig.colorbar(mappable=PC,cax=cax)

         # plot original selection
         vertices.append(vertices[0])
         vlon,vlat   = np.array(vertices).transpose()
         bmap.plot(vlon,vlat,'m',linewidth=4,ax=po.ax,latlon=True)

         if 0:
            # finish with coast etc
            Fplt.finish_map(bmap,ax=po.ax)

         # depth=0 contour
         bmap.contour(plon,plat,depth,colors=cols,linewidth=4,latlon=True,ax=po.ax,levels=[0.])
         
         # save fig
         figname  = outdir+"/test_"+gridname+"_"+vbl+".png"
         print('\nSaving '+figname+'\n')
         po.fig.savefig(figname)
         po.ax.cla()

   plt.close(po.fig)
