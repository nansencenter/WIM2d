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
#initialise grid_prams dictionary:
def Grid_Prams():
   grid_prams  = {'nx'        :0,
                  'ny'        :0,
                  'dx'        :0.0,
                  'dy'        :0.0,
                  'X'         :0.0,
                  'Y'         :0.0,
                  'scuy'      :0.0,
                  'scvx'      :0.0,
                  'scp2'      :0.0,
                  'scp2i'     :0.0,
                  'LANDMASK'  :0.0}
   return grid_prams
##############################################################

##############################################################
# initialise ice_fields dictionary:
def Ice_Fields(icec=0.0,iceh=0.0,dfloe=0.0,ICE_MASK=0.0,):

   ice_fields  = {'icec'      :icec,
                  'iceh'      :iceh,
                  'dfloe'     :dfloe,
                  'ICE_MASK'  :ICE_MASK}

   return ice_fields
##############################################################

##############################################################
# initialise wave_fields dictionary:
def Wave_Fields(Hs=0.0,Tp=0.0,mwd=0.0,WAVE_MASK = 0.0):

   wave_fields = {'mwd'       :mwd,
                  'WAVE_MASK' :WAVE_MASK,
                  'Hs'        :Hs,
                  'Tp'        :Tp}

   return wave_fields
##############################################################

##############################################################
# initialise out_fields dictionary:
def Out_Fields(dfloe=0.0,taux=0.0,tauy=0.0,Hs=0.0,Tp=0.0):

   out_fields  = {'dfloe'  :dfloe,
                  'taux'   :taux,
                  'tauy'   :tauy,
                  'Hs'     :Hs,
                  'Tp'     :Tp}

   return out_fields
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
   s1['nx']    = nx
   s1['ny']    = ny
   aid         = open(afile,'rb')
   #
   keys  = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
   for key in keys:
      s1[key]  = get_array(aid,nx,ny)
   #
   aid.close()
   s1['dx'] = s1['X'][1,0]-s1['X'][0,0]
   s1['dy'] = s1['Y'][0,1]-s1['Y'][0,0]
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

   ## ice fields
   keys  = ['icec','iceh','dfloe']
   for key in keys:
      s2[key]  = get_array(aid,nx,ny)

   ## wave fields
   keys  = ['Hs','Tp','mwd']
   for key in keys:
      s3[key]  = get_array(aid,nx,ny)
   aid.close()

   ## masks
   s2['ICE_MASK']    = np.zeros((nx,ny))
   s3['WAVE_MASK']   = np.zeros((nx,ny))
   s2['ICE_MASK'] [s2['icec']>0.05] = 1.0
   s3['WAVE_MASK'][s3['Hs']>0.0]    = 1.0
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
   keys        = ['dfloe','taux','tauy','Hs','Tp']
   ##
   for key in keys:
      s2[key]  = get_array(aid,nx,ny)

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

   keys  = ['dfloe','taux','tauy','Hs','Tp'] # can be got from s2.keys(),
   n     = 0
   for key in keys:
      s2[key]  = out_arrays[:,:,n]
      n        = n+1
   
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
   keys        = ['dfloe','taux','tauy','Hs','Tp']
   ##
   for key in keys:
      s2[key]  = get_array(aid,nx,ny)

   aid.close()
   ###########################################################
   
   # outputs
   return out_fields
##############################################################
