import numpy as np
import os
import sys
import struct

##############################################################
def fn_check_init(outdir):
   # routine to get initial fields from binary files:

   ###########################################################
   def get_array(fid,nx,ny):
      #routine to get the array from the .a (binary) file
      fmt_size = 4               # real*4 so 4B per number 
      recs     = nx*ny
      rec_size = recs*fmt_size
      ##
      data  = fid.read(rec_size)
      fld   = struct.unpack(recs*'f',data)
      fld   = np.array(fld)
      fld   = fld.reshape((ny,nx)).transpose()

      return fld
   ###########################################################

   ###########################################################
   #define grid_prams object:
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
   ###########################################################

   ###########################################################
   #define ice_fields object:
   class Ice_Fields(object):
      def __init__(self,icec=0,iceh=0,dfloe=0.0,ICE_MASK=0.0):
         self.icec      = icec
         self.iceh      = iceh
         self.dfloe     = dfloe
         self.ICE_MASK  = ICE_MASK
   ###########################################################

   ###########################################################
   #define wave_fields object:
   class Wave_Fields(object):
      def __init__(self,Hs=0,Tp=0,mwd=0.0,WAVE_MASK=0.0):
         self.Hs        = Hs
         self.Tp        = Tp
         self.mwd       = mwd
         self.WAVE_MASK = WAVE_MASK
   ###########################################################

   ###########################################################
   # get dimensions from .b file
   afile    = outdir+'/binaries/wim_grid.a'
   bfile    = outdir+'/binaries/wim_grid.b'
   afile2   = outdir+'/binaries/wim_init.a'
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

   ###########################################################
   # can now read data from .a file
   aid   = open(afile2,'rb')
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
   return grid_prams,ice_fields,wave_fields
##############################################################
