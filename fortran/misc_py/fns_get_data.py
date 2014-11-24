import numpy as np
import os
import sys
import struct

##############################################################
def get_array(fid,nx,ny):
   #routine to get the array from the .a (binary) file
   fmt_size = 4               # real*4 so 4B per number 
   recs     = nx*ny
   rec_size = recs*fmt_size
   ##
   data  = fid.read(rec_size)
   fld   = struct.unpack(recs*'f',data)
   fld   = np.array(fld)
   fld   = fld.reshape((ny,nx)).transpose()  # need to transpose because of differences between
                                             # python/c and fortran/matlab 

   return fld
##############################################################

##############################################################
#define Grid_Prams object:
class Grid_Prams(object):
   def __init__(self,nx=0,ny=0,dx=0.0,dy=0.0,X=0.0,Y=0.0,
                  scuy=0.0,scvx=0.0,scp2=0.0,scp2i=0.0,LANDMASK=0.0):
      self.nx        = nx
      self.ny        = ny
      self.dx        = dx
      self.dy        = dy
      self.X         = X
      self.Y         = Y
      self.scuy      = scuy
      self.scvx      = scvx
      self.scp2      = scp2
      self.scp2i     = scp2i
      self.LANDMASK  = LANDMASK
##############################################################

##############################################################
#define Ice_Fields object:
class Ice_Fields(object):
   def __init__(self,icec=0,iceh=0,dfloe=0.0,ICE_MASK=0.0):
      self.icec      = icec
      self.iceh      = iceh
      self.dfloe     = dfloe
      self.ICE_MASK  = ICE_MASK
##############################################################

##############################################################
#define Wave_Fields object:
class Wave_Fields(object):
   def __init__(self,Hs=0,Tp=0,mwd=0.0,WAVE_MASK=0.0):
      self.Hs        = Hs
      self.Tp        = Tp
      self.mwd       = mwd
      self.WAVE_MASK = WAVE_MASK
##############################################################

##############################################################
#define Out_Fields object:
class Out_Fields(object):
   def __init__(self,dfloe=0.0,taux=0.0,tauy=0.0,Hs=0,Tp=0):
      self.dfloe     = dfloe
      self.taux      = taux
      self.tauy      = tauy
      self.Hs        = Hs
      self.Tp        = Tp
##############################################################

##############################################################
def fn_check_grid(outdir):
   # routine to get grid and other parameters
   # from binary files

   ###########################################################
   # get dimensions from .b file
   afile    = outdir+'/binaries/wim_grid.a'
   bfile    = outdir+'/binaries/wim_grid.b'
   #
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   nxline   = lines[1].split()
   nyline   = lines[2].split()
   nx       = int(nxline[0])
   ny       = int(nyline[0])

   grid_prams  = Grid_Prams()
   s1          = grid_prams   #pointer with shorter name (NB change s1=> change grid_prams)
   s1.nx       = nx
   s1.ny       = ny
   #
   aid         = open(afile,'rb')
   s1.X        = get_array(aid,nx,ny)
   s1.Y        = get_array(aid,nx,ny)
   s1.scuy     = get_array(aid,nx,ny)
   s1.scvx     = get_array(aid,nx,ny)
   s1.scp2     = get_array(aid,nx,ny)
   s1.scp2i    = get_array(aid,nx,ny)
   s1.LANDMASK = get_array(aid,nx,ny)
   #
   aid.close()
   s1.dx = s1.X[1,0]-s1.X[0,0]
   s1.dy = s1.Y[0,1]-s1.Y[0,0]
   ###########################################################
   
   # output
   return grid_prams
##############################################################

##############################################################
def fn_check_init(outdir):
   # routine to get initial fields from binary files:
   afile    = outdir+'/binaries/wim_init.a'
   bfile    = outdir+'/binaries/wim_init.b'

   ###########################################################
   # get dimensions from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   nxline   = lines[1].split()
   nyline   = lines[2].split()
   nx       = int(nxline[0])
   ny       = int(nyline[0])
   ###########################################################

   ###########################################################
   # can now read data from .a file
   aid   = open(afile,'rb')
   ##
   ice_fields  = Ice_Fields()
   wave_fields = Wave_Fields()
   s2          = ice_fields
   s3          = wave_fields
   ##
   s2.icec  = get_array(aid,nx,ny)
   s2.iceh  = get_array(aid,nx,ny)
   s2.dfloe = get_array(aid,nx,ny)
   ##
   s3.Hs  = get_array(aid,nx,ny)
   s3.Tp  = get_array(aid,nx,ny)
   s3.mwd = get_array(aid,nx,ny)
   aid.close()
   ##
   s2.ICE_MASK    = np.zeros((nx,ny))
   s2.ICE_MASK[s2.icec>0.05]  = 1.0
   s3.WAVE_MASK   = np.zeros((nx,ny))
   s3.WAVE_MASK[s3.Hs>0.0]  = 1.0
   ###########################################################
   
   # outputs
   return ice_fields,wave_fields
##############################################################

##############################################################
def fn_check_out_bin(outdir):
   # routine to get output fields from binary files:
   afile    = outdir+'/binaries/wim_out.a'
   bfile    = outdir+'/binaries/wim_out.b'

   ###########################################################
   # get dimensions from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   nxline   = lines[1].split()
   nyline   = lines[2].split()
   nx       = int(nxline[0])
   ny       = int(nyline[0])
   ###########################################################

   ###########################################################
   # can now read data from .a file
   aid   = open(afile,'rb')
   ##
   out_fields  = Out_Fields()
   s2          = out_fields
   ##
   s2.dfloe = get_array(aid,nx,ny)
   s2.taux  = get_array(aid,nx,ny)
   s2.tauy  = get_array(aid,nx,ny)
   s2.Hs    = get_array(aid,nx,ny)
   s2.Tp    = get_array(aid,nx,ny)
   aid.close()
   ###########################################################
   
   # outputs
   return out_fields
##############################################################

##############################################################
def fn_check_out_arr(out_arrays):
   # routine to convert out_arrays to Out_Fields object
   out_fields  = Out_Fields()
   s2          = out_fields
   ##
   s2.dfloe = out_arrays[:,:,1]
   s2.taux  = out_arrays[:,:,2]
   s2.tauy  = out_arrays[:,:,3]
   s2.Hs    = out_arrays[:,:,4]
   s2.Tp    = out_arrays[:,:,5]
   aid.close()
   ###########################################################
   
   # outputs
   return out_fields
##############################################################

##############################################################
def fn_check_prog(outdir,n):
   # routine to get progress fields from binary files:
   cts      = '%3.3d'   % (n) # time step
   afile    = outdir+'/binaries/prog/wim_prog'+cts+'.a'
   bfile    = outdir+'/binaries/prog/wim_prog'+cts+'.b'

   ###########################################################
   # get dimensions from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   nxline   = lines[1].split()
   nyline   = lines[2].split()
   nx       = int(nxline[0])
   ny       = int(nyline[0])
   ###########################################################

   ###########################################################
   # can now read data from .a file
   aid   = open(afile,'rb')
   ##
   out_fields  = Out_Fields()
   s2          = out_fields
   ##
   s2.dfloe = get_array(aid,nx,ny)
   s2.taux  = get_array(aid,nx,ny)
   s2.tauy  = get_array(aid,nx,ny)
   s2.Hs    = get_array(aid,nx,ny)
   s2.Tp    = get_array(aid,nx,ny)
   aid.close()
   ###########################################################
   
   # outputs
   return out_fields
##############################################################
