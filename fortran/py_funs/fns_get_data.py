import numpy as np
import os
import sys
import struct

#####################################################
class var_info:
   # simple object
   # scalars: obj.value,obj.unit
   # numpy arrays: obj.data,obj.unit,obj.length
   def __init__(self,value,unit=None):

      # set units
      if unit is not '':
         self.unit   = unit
      else:
         self.unit   = None

      if hasattr(value,'ndim'):
         # an array
         self.data   = value
         self.length = len(value)
         self.min    = np.min(value)
         self.max    = np.max(value)
         self.first  = value[0]
         self.last   = value[-1]
      else:
         # a number
         self.value  = value

      return
#####################################################

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
def fn_check_grid(outdir):
   # routine to get grid and other parameters
   # from binary files

   ###########################################################
   afile       = outdir+'/wim_grid.a'
   bfile       = outdir+'/wim_grid.b'
   fields,info = fn_read_general_binary(afile)
   aliases     = key_aliases(inverse=True)

   grid_prams  = {}
   keys        = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
   for key in keys:
      if key in fields.keys():
         key2  = key
      else:
         key2  = aliases[key]
      grid_prams.update({key:fields[key2]})
   ###########################################################
   
   ###########################################################
   # extra info
   nx,ny = grid_prams['X'].shape
   grid_prams.update({'nx':nx})
   grid_prams.update({'ny':ny})
   #
   dx = np.mean(grid_prams['scvx'])
   dy = np.mean(grid_prams['scuy'])
   grid_prams.update({'dx':dx})
   grid_prams.update({'dy':dy})
   ###########################################################
   
   # output
   return grid_prams
##############################################################

##############################################################
def fn_check_init(outdir):
   # routine to get initial fields from binary files:
   afile       = outdir+'/wim_init.a'
   bfile       = outdir+'/wim_init.b'
   fields,info = fn_read_general_binary(afile)
   aliases     = key_aliases(inverse=True)
   
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
def fn_bfile_info(bfile):
   # routine to get output fields from binary files:

   ###########################################################
   # get info like dimensions and variable names from .b file
   bid   = open(bfile,'r')
   lines = bid.readlines()
   bid.close()

   int_list = ['nx','ny','Nrecs','Norder'] # these are integers not floats
   binfo = {}  # dictionary with info about fields in corresponding .a file
   vlist = []  # list of variable names in corresponding .a file (in order)

   do_vlist = 0
   for lin in lines:
      ls = lin.split()
      if ls != []:
         # skip blank lines

         # if not at "Record number and name"
         if not(ls[0]=='Record' and ls[1]=='number'):
            key   = ls[1]
            if not do_vlist:
               if key in int_list:
                  val   = int(ls[0])
               else:
                  val   = float(ls[0])
               binfo.update({key:val})
            else:
               vlist.append(key)
         else:
            # have got down to list of variables
            do_vlist = 1

   if binfo['Nrecs']!=len(vlist):
      raise ValueError('Inconsistent number of records in file: '+bfile)

   return binfo,vlist
##############################################################

##############################################################
def fn_read_general_binary(afile):
   # routine to get output fields from binary files:

   ###########################################################
   # get dimensions and variable names from .b file
   bfile       = afile[:-2]+'.b'
   binfo,vlist = fn_bfile_info(bfile)

   nx    = binfo['nx']
   ny    = binfo['ny']
   if binfo['Norder']==1:
      order = 'fortran'
   else:
      order = 'C'
   ###########################################################

   ###########################################################
   # can now read data from .a file
   aid   = open(afile,'rb')

   out   = {}
   for key in vlist:
      out.update({key:get_array(aid,nx,ny,order=order)})

   aid.close()
   ###########################################################
   
   # outputs
   return out,binfo
##############################################################

##############################################################
def fn_check_out_bin(outdir):
   # routine to get output fields from binary files:
   afile       = outdir+'/wim_out.a'
   bfile       = outdir+'/wim_out.b'
   fields,info = fn_read_general_binary(afile)
   aliases     = key_aliases(inverse=True)
   
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
   out_fields  = {}

   keys  = ['dfloe','taux','tauy','Hs','Tp'] # can be got from s2.keys(),
   for n,key in enumerate(keys):
      out_fields.update({key:out_arrays[:,:,n]})
   
   # outputs
   return out_fields
##############################################################

##############################################################
def fn_check_prog(outdir,cts):
   # routine to get progress fields from binary files:
   # cts is a string eg '010' or '0010' corresponding to the time step
   if type(cts)==type(0):
      # convert from int to str of correct length
      import os
      fils  = os.listdir(outdir+'/binaries/prog/')
      n     = 0
      while '.swp' in fils[n] or '.DS_Store' in fils[n]:
          n = n+1
      cts0  = fils[n].strip('wim_prog')[:-2]
      fmt   = '%'+str(len(cts0))+'.'+str(len(cts0))+'d'
      cts   = fmt %(cts)
      # print(cts0,cts)

   afile       = outdir+'/binaries/prog/wim_prog'+cts+'.a'
   bfile       = outdir+'/binaries/prog/wim_prog'+cts+'.b'
   fields,info = fn_read_general_binary(afile)
   aliases     = key_aliases(inverse=True)
   
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

########################################################
def read_datfile(dname):
   print('\nOpening '+dname+'\n')
   fid = open(dname,'r')
   

   #####################################################
   # get header info
   in_hdr   = True
   Params   = None
   while in_hdr:
      lin0  = fid.readline()

      ##################################################
      if 'File info:' in lin0:
         # start to extract file info
         # (no of variables, length of records)
         info  = {}
         pos   = fid.tell()
         lin0  = fid.readline()
         while ':' in lin0:
            ss    = lin0.split()
            val   = int(ss[-1])
            vbl   = ss[1]
            info.update({vbl:val})
            #
            pos   = fid.tell()
            lin0  = fid.readline()

         # go back to start of line
         fid.seek(pos)
      ##################################################

      ##################################################
      elif 'Parameters:' in lin0:
         # start to extract parameters
         pos      = fid.tell()
         lin0     = fid.readline()
         Params   = {}
         while ':' in lin0:
            if '(' not in lin0:
               ss    = lin0.split()
               val   = float(ss[-1])
               obj   = var_info(val)
               vbl   = ss[1]
            else:
               i0    = lin0.index('(')
               i1    = lin0.index(')')
               unit  = lin0[i0+1:i1]
               val   = float(lin0[i1+1:].split()[-1])
               vbl   = lin0[:i0].split()[1]
               obj   = var_info(val,unit)

            # update dictionary
            Params.update({vbl:obj})

            # read next line
            pos   = fid.tell()
            lin0  = fid.readline()

         # go back to start of line
         fid.seek(pos)
      ##################################################

      ##################################################
      elif 'Variables:' in lin0:
         # get variable names and units
         Vlist    = []
         Vunits   = []
         for Nv in range(info['Nvars']):

            lin0  = fid.readline()
            if '(' not in lin0:
               vbl   = lin0.split()[1]
               unit  = ''
            else:
               ss    = lin0.split('(')
               vbl   = ss[0].split()[-1]
               unit  = ss[-1].split(')')[0]

            Vlist .append(vbl)
            Vunits.append(unit)
      ##################################################
         
      ##################################################
      elif lin0.split()!=[]:
         # check if 2nd-to-last line of header
         in_hdr   = not('#######' in lin0)
      ##################################################


   #####################################################
   # skip next line (blank) - then we are up to the data
   lin0     = fid.readline() 
   #####################################################


   #####################################################
   # get data
   arr      = []
   for i in range(info['Lvars']):
      line  = fid.readline().split()
      row   = []
      for thing in line:
         row.append(float(thing))
      arr.append(row)
   fid.close()
   #####################################################

   arr   = np.array(arr).transpose()
   out   = {}
   for n,vbl in enumerate(arr):
      obj   = var_info(vbl,Vunits[n])
      out  .update({Vlist[n]:obj})
   
   return out,Params
########################################################

########################################################
def write_datfile(dname,Vars,title=None,Params=None):
   print('\nWriting to '+dname+'\n')
   fid = open(dname,'w')

   #####################################################
   # write title
   if title is not None:
      if title[0]!='#':
         title = '# '+title
      fid.write(title+'\n\n')
   #####################################################

   #####################################################
   # Get variable info
   Varr     = []
   Vlines   = []
   for vbl in Vars.keys():
      line  = '# '+vbl+'('+Vars[vbl].unit+')\n'
      Vlines.append(line)
      Varr.append(Vars[vbl].data)
   Varr        = np.array(Varr).transpose()
   Lvars,Nvars = Varr.shape

   # Write file info:
   fid.write('# File info:\n')
   fid.write('# Nvars   : %d\n'     %(Nvars))
   fid.write('# Lvars   : %d\n\n'   %(Lvars))
   #####################################################

   #####################################################
   # write parameters
   if Params is not None:
      fid.write('# Parameters:\n')
      for vbl in Params.keys():
         if Params[vbl].unit is None:
            line  = '# '+vbl+' : '
         else:
            line  = '# '+vbl+'('+Params[vbl].unit+') : '
         line  = line+str(Params[vbl].value)+'\n'
         fid.write(line)
      fid.write('\n')
   #####################################################

   #####################################################
   # write variable info
   fid.write('# Variables:\n')
   for line in Vlines:
      fid.write(line)
   line  = 50*'#'
   fid.write(line+'\n\n')
   #####################################################

   #####################################################
   # write variables
   blk   = 4*' '
   Nv    = 0
   for row in Varr:
      line  = ''
      Nv    = Nv+1
      for val in row:
         line  = line+blk+str(val)

      line  = line[len(blk):]
      if Nv==Lvars:
         fid.write(line)
      else:
         fid.write(line+'\n')
   #####################################################

   # close file
   fid.close()
########################################################
