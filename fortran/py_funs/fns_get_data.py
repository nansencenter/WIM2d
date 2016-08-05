import numpy as np
import os,sys
import struct
import fns_plot_data as Fplt
import datetime as dtm

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


   data  = fid.read(rec_size) # no of bytes to read
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
            val   = ls[0]
            key   = ls[1]
            if not do_vlist:
               if key in int_list:
                  val   = int(val)
               elif ("T" in val) and ("Z" in val):
                  import datetime as dtm
                  val   = dtm.datetime.strptime(val,"%Y%m%dT%H%M%SZ")
               else:
                  val   = float(val)
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
   import os
   sz=os.path.getsize(afile)
   nv=len(vlist)
   fmt_size=int(sz/float(nv*nx*ny)) # 4 for single, 8 for double

   # print(sz,4*nx*ny*nv)
   # print(fmt_size,nx,ny,nv)

   aid   = open(afile,'rb')

   out   = {}
   for key in vlist:
      out.update({key:get_array(aid,nx,ny,order=order,fmt_size=fmt_size)})
      # print(key)
      # print(out[key].min(),out[key].max())

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

   keys  = ['dfloe','taux','tauy','Hs','Tp','mwd']
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
      # cts0  = fils[0].strip('wim_prog')[:-2]
      n     = 0
      while '.swp' in fils[n] or '.DS_Store' in fils[n]:
          n = n+1
      cts0  = fils[n].strip('wim_prog')[:-2]
      fmt   = '%'+str(len(cts0))+'.'+str(len(cts0))+'d'
      cts   = fmt %(cts)
      print(cts0,cts)

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


############################################################################
class wim_results:

   def __init__(self,outdir='.',grid_dir=None):

      self.bindir    = outdir+'/binaries'
      if not os.path.exists(self.bindir):
         raise ValueError(self.bindir+ ' does not exist')


      # =====================================================================
      # check for grid files
      self.grid_dir  = None
      if grid_dir is not None:
         lst   = os.listdir(grid_dir)
         for f in lst:
            if 'wim_grid' in f and '.a' in f:
               self.grid_dir  = grid_dir
               break
         
         if self.grid_dir is None:
            print('wim_grid.[a,b] not in '+grid_dir)
            raise ValueError('input "grid_dir" does not contain grid files')
         else:
            self.grid_dir  = grid_dir

      else:
         dirs  = [self.bindir,outdir+'/../grid']
         i     = 0

         while (self.grid_dir is None) and (i<len(dirs)):
            dir_i = dirs[i]
            lst   = os.listdir(dir_i)
            for f in lst:
               if 'wim_grid' in f and '.a' in f:
                  self.grid_dir  = dir_i
                  break
            i +=1

         if self.grid_dir is None:
            raise ValueError('wim_grid.[a,b] not found')
      # =====================================================================
      # print(self.grid_dir)
      # sys.exit()


      # =====================================================================
      # check for initial conditions
      binlist        = os.listdir(self.bindir)
      self.init_file = None

      for f in binlist:
         if ('wim_init' in f) and ('.a' in f):
            self.init_file = f
            break

      if self.init_file is not None:
         bfile = self.init_file.replace('.a','.b')
         info,self.variables_init   = fn_bfile_info(self.bindir+'/'+bfile)
         self.start_time  = info['t_out']
      # =====================================================================
      # print(self.init_file)
      # print(self.start_time)
      # print(self.variables_init)
      # sys.exit()


      # =====================================================================
      # check for final conditions
      self.out_file = None
      for f in binlist:
         if ('wim_out' in f) and ('.a' in f):
            self.out_file = f
            break

      if self.out_file is not None:
         bfile = self.out_file.replace('.a','.b')
         info,self.variables_out  = fn_bfile_info(self.bindir+'/'+bfile)
         self.finish_time  = info['t_out']
      # =====================================================================
      # print(self.out_file)
      # print(self.finish_time)
      # print(self.variables_out)
      # sys.exit()


      # =====================================================================
      # Check for progress files
      self.prog_dir  = self.bindir+'/prog'
      prog_files     = os.listdir(self.prog_dir)

      # find the .a files
      alist    = []
      steplist = []
      for pf in prog_files:
         if ('wim_prog' in pf) and ('.a'==pf[-2:]):
            alist.append(pf)
            stepno   = pf.strip('wim_prog').strip('.a') # step no eg 001 or date yyyymmddThhmmssZ
            if ('T' in stepno) and ('Z' in stepno):
               # date
               steplist.append(
                     dtm.datetime.strptime(stepno,'%Y%m%dT%H%M%SZ'))
            else:
               # integer step no
               steplist.append(int(stepno))

      # sort according to steplist:
      slist = sorted([(e,i) for i,e in enumerate(steplist)])
      self.prog_times   = [steplist[i] for e,i in slist]
      self.prog_files   = [alist[i] for e,i in slist]
      self.Nprog_files  = len(alist)
      # =====================================================================
      # print(self.prog_dir)
      # print(self.prog_files)
      # print(self.prog_times)
      # sys.exit()


      # =====================================================================
      # final checks
      if (self.init_file is None) and\
         (self.out_file is None)  and\
         (self.Nprog_files ==0 ):
         raise ValueError('No results in ' + outdir)

      # where to save figures
      self.figdir = outdir+'/figs'
      # =====================================================================

      return
   ##########################################################################


   ##########################################################################
   def get_grid(self):
      return fn_check_grid(self.grid_dir)
   ##########################################################################


   ##########################################################################
   def initial_fields(self):

      if self.init_file is None:
         raise ValueError('Initial conditions not outputted: wim_init*.[a,b]')
      else:
         afile = self.bindir+'/'+self.init_file
         return fn_read_general_binary(afile)
   ##########################################################################


   ##########################################################################
   def out_fields(self):

      if self.out_file is None:
         raise ValueError('Final conditions not outputted: wim_out*.[a,b]')
      else:
         afile = self.bindir+'/'+self.out_file
         return fn_read_general_binary(afile)
   ##########################################################################


   ##########################################################################
   def prog_fields(self,time_index=0):

      if self.Nprog_files==0:
         raise ValueError('Progress files not outputted: wim_prog*.[a,b]')
      else:
         afile = self.prog_dir+'/'+self.prog_files[time_index]
         return fn_read_general_binary(afile)
   ##########################################################################


   ##########################################################################
   def plot_initial(self):
      # Look at initial fields:

      if self.init_file is None:
         print('Initial conditions not outputted: wim_init*.[a,b]')
         print('- not plotting')
      else:
         print("Plotting initial conditions...")
         if not os.path.exists(self.figdir):
            os.mkdir(figdir)

         figdir1  = self.figdir+'/init/'
         if not os.path.exists(figdir1):
            os.mkdir(figdir1)

         grid_prams  = self.get_grid()
         afile       = self.bindir+'/'+self.init_file
         fields,info = fn_read_general_binary(afile)
         Fplt.fn_plot_gen(grid_prams,fields,figdir1)    # plot initial conditions

         print("Plots in "+figdir1)
         print(" ")
      return
   ##########################################################################


   ##########################################################################
   def plot_final(self):
      # Look at final fields:

      if self.out_file is None:
         print('Initial conditions not outputted: wim_out*.[a,b]')
         print('- not plotting')
      else:
         print("Plotting final conditions...")
         if not os.path.exists(self.figdir):
            os.mkdir(figdir)

         figdir1  = self.figdir+'/final/'
         if not os.path.exists(figdir1):
            os.mkdir(figdir1)

         grid_prams  = self.get_grid()
         afile       = self.bindir+'/'+self.out_file
         fields,info = fn_read_general_binary(afile)
         Fplt.fn_plot_gen(grid_prams,fields,figdir1)    # plot final conditions

         print("Plots in "+figdir1)
         print(" ")
      return
   ##########################################################################


   ################################################################
   def plot_prog(self):

      # =============================================================
      # Plot progress files (if they exist)
      if self.Nprog_files==0:
         print('No progress files wim*prog.[a,b] in '+
               self.bindir+'/prog')
         print('Not plotting')
         return
      else:
         print('\nPlotting progress files...\n')
      # =============================================================


      grid_prams  = self.get_grid()
      pdir        = self.prog_dir
      prog_files  = os.listdir(pdir)

      figdir3  = self.figdir
      if not os.path.exists(figdir3):
         os.mkdir(figdir3)
      figdir3  = self.figdir+'/prog'
      if not os.path.exists(figdir3):
         os.mkdir(figdir3)

      # =============================================================
      # determine the plotting limits
      if 0:
         # set colorbar axes manually
         zlims  = {'icec'  :[0,1],      \
                   'iceh'  :[0,5],      \
                   'Dmax'  :[0,300],    \
                   'tau_x' :[-.5,.5],  \
                   'tau_y' :[-.05,.05],\
                   'Hs'    :[0,4],        \
                   'Tp'    :[10,20],      \
                   'mwd'   :[-180,180]}

      elif 0:
         # let python choose
         # - different for each time step,
         #   so not good for a movie for example
         zlims  = {'icec'  :None,\
                   'iceh'  :None,\
                   'Dmax'  :None,\
                   'tau_x' :None,\
                   'tau_y' :None,\
                   'Hs'    :None,\
                   'Tp'    :None,\
                   'mwd'   :None}

      else:
         # set colorbar axes automatically
         zdef  = [1.e30,-1.e30]
         zlims  = {'icec' :1*zdef,\
                   'iceh' :1*zdef,\
                   'Dmax' :1*zdef,\
                   'tau_x':1*zdef,\
                   'tau_y':1*zdef,\
                   'Hs'   :1*zdef,\
                   'Tp'   :1*zdef,\
                   'mwd'  :1*zdef}

         # =============================================================
         # determine the plotting limits
         # by checking the arrays
         alist = self.prog_files
         tlist = []
         for pf in self.prog_files:
            afile = pdir+'/'+pf
            #
            fields,info = fn_read_general_binary(afile)
            for key in fields.keys():
               if key in zlims.keys():
                  zmin  = fields[key].min()
                  zmax  = fields[key].max()
                  if zmin<zlims[key][0]:
                     zlims[key][0]  = zmin
                  if zmax>zlims[key][1]:
                     zlims[key][1]  = zmax
            #
            afile = pdir+'/'+pf
            #
            fields,info = fn_read_general_binary(afile)
            for key in fields.keys():
               if key in zlims.keys():
                  zmin  = fields[key].min()
                  zmax  = fields[key].max()
                  if zmin<zlims[key][0]:
                     zlims[key][0]  = zmin
                  if zmax>zlims[key][1]:
                     zlims[key][1]  = zmax

            stepno   = pf.strip('wim_prog').strip('.a') # step no eg 001 or date yyyymmddThhmmssZ
            tlist.append(stepno)
            #
            afile = pdir+'/'+pf
            #
            fields,info = fn_read_general_binary(afile)
            for key in fields.keys():
               if key in zlims.keys():
                  zmin  = fields[key].min()
                  zmax  = fields[key].max()
                  if zmin<zlims[key][0]:
                     zlims[key][0]  = zmin
                  if zmax>zlims[key][1]:
                     zlims[key][1]  = zmax
         # =============================================================

      # =============================================================


      # =============================================================
      # do the plots
      for vbl in fields.keys():
         print('\n')
         print('Plotting results for '+vbl)
         print('\n')

         figdir3B    = figdir3+'/'+vbl
         if not os.path.exists(figdir3B):
            os.mkdir(figdir3B)

         for i,pf in enumerate(alist):
            afile       = pdir+'/'+pf
            print(afile)
            fields,info = fn_read_general_binary(afile)
            Fplt.fn_plot_gen(grid_prams,fields,figdir3B,\
                  zlims_in=zlims,text=tlist[i],vlist=[vbl])
      # =============================================================

      return
   # =============================================================
