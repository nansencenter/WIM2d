import numpy as np
import os
import sys
import struct

##############################################################
def key_aliases(inverse=False):
   if not inverse:
      aliases  = {'Dmax'   :'dfloe' ,\
                  'cice'   :'icec'  ,\
                  'hice'   :'iceh'  ,\
                  'tau_x'  :'taux'  ,\
                  'tau_y'  :'tauy'  }
   else:
      aliases  = {'dfloe'  :'Dmax'  ,\
                  'icec'   :'cice'  ,\
                  'iceh'   :'hice'  ,\
                  'taux'   :'tau_x' ,\
                  'tauy'   :'tau_y' }
   return aliases
##############################################################

##############################################################
def get_array(fid,nx,ny,fmt_size=4,order='F'):
   # routine to get the array from the .a (binary) file
   # * fmt_size = size in bytes of each entry)
   #   > default = 4 (real*4/single precision)
   recs     = nx*ny
   rec_size = recs*fmt_size
   #
   if fmt_size==4:
      fmt_py   = 'f' # python string for single
   else:
      fmt_py   = 'd' # python string for double


   data  = fid.read(rec_size)
   fld   = struct.unpack(recs*fmt_py,data)
   fld   = np.array(fld)
   fld   = fld.reshape((nx,ny),order=order)

   return fld
##############################################################

##############################################################
#initialise grid_prams dictionary:
def Grid_Prams(nx=0,ny=0,dx=0,dy=0,
      X=0,Y=0,scuy=0,scvx=0,scp2=0,
      scp2i=0,LANDMASK=0):
   grid_prams  = {'nx'        :nx,
                  'ny'        :ny,
                  'dx'        :dx,
                  'dy'        :dy,
                  'X'         :X,
                  'Y'         :Y,
                  'scuy'      :scuy,
                  'scvx'      :scvx,
                  'scp2'      :scp2,
                  'scp2i'     :scp2i,
                  'LANDMASK'  :LANDMASK}
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
   afile    = outdir+'/wim_grid.a'
   bfile    = outdir+'/wim_grid.b'
   fields   = fn_read_general_binary(afile)
   aliases  = key_aliases(inverse=True)

   grid_prams  = {}
   keys        = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      grid_prams.update({key:fields[key2]})
   #
   nx,ny = grid_prams['X'].shape
   grid_prams.update({'nx':nx})
   grid_prams.update({'ny':ny})
   ###########################################################
   
   # output
   return grid_prams
##############################################################

##############################################################
def fn_check_init(outdir):
   # routine to get initial fields from binary files:
   afile    = outdir+'/wim_init.a'
   bfile    = outdir+'/wim_init.b'
   fields   = fn_read_general_binary(afile)
   aliases  = key_aliases(inverse=True)
   
   ###########################################################
   ## ice fields
   keys        = ['icec','iceh','dfloe']
   ice_fields  = {}
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      ice_fields.update({key:fields[key2]})

   # ice mask
   ice_fields.update({'ICE_MASK':0*ice_fields['icec']})
   ice_fields['ICE_MASK'][ice_fields['icec']>0.05] = 1.0
   ###########################################################

   ###########################################################
   ## wave fields
   keys        = ['Hs','Tp','mwd']
   wave_fields = {}
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      wave_fields.update({key:fields[key]})

   # wave mask
   wave_fields.update({'WAVE_MASK':0*wave_fields['Hs']})
   wave_fields['WAVE_MASK'][wave_fields['Hs']>0.0] = 1.0
   ###########################################################
   
   # outputs
   return ice_fields,wave_fields
##############################################################

##############################################################
def fn_read_general_binary(afile):
   # routine to get output fields from binary files:
   bfile = afile[:-2]+'.b'

   ###########################################################
   # get dimensions and variable names from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   Nrecs    = int(lines[0].split()[0])
   Nord     = int(lines[1].split()[0])
   nx       = int(lines[2].split()[0])
   ny       = int(lines[3].split()[0])
   if Nord==1:
      order = 'fortran'
   else:
      order = 'C'

   Nlines   = len(lines)
   for n in range(2,Nlines):
      lin   = lines[n]
      if 'Record number and name:' in lin:
         n0 = n+1
         break

   keys  = []
   for n in range(n0,Nlines):
      lin   = lines[n]
      vname = lin.split()[1]
      keys.append(vname)
   ###########################################################

   if len(keys)!=Nrecs:
      raise ValueError('Inconsistent number of records in file: '+bfile)

   ###########################################################
   # can now read data from .a file
   aid   = open(afile,'rb')
   ##
   out   = {}
   ##
   for key in keys:
      out.update({key:get_array(aid,nx,ny,order=order)})

   aid.close()
   ###########################################################
   
   # outputs
   return out
##############################################################

##############################################################
def fn_check_out_bin(outdir):
   # routine to get output fields from binary files:
   afile    = outdir+'/wim_out.a'
   bfile    = outdir+'/wim_out.b'
   fields   = fn_read_general_binary(afile)
   aliases  = key_aliases(inverse=True)
   
   ###########################################################
   ## out fields
   keys        = ['dfloe','taux','tauy','Hs','Tp']
   out_fields  = {}
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      out_fields.update({key:fields[key2]})
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
def fn_check_prog(outdir,cts):
   # routine to get progress fields from binary files:
   # cts is a string eg '010' or '0010' corresponding to the time step
   afile    = outdir+'/binaries/prog/wim_prog'+cts+'.a'
   bfile    = outdir+'/binaries/prog/wim_prog'+cts+'.b'
   fields   = fn_read_general_binary(afile)
   aliases  = key_aliases(inverse=True)
   
   ###########################################################
   ## out fields
   keys        = ['dfloe','taux','tauy','Hs','Tp']
   out_fields  = {}
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      out_fields.update({key:fields[key2]})
   ###########################################################
   
   # outputs
   return out_fields
##############################################################
