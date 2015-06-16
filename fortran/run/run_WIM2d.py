################################################################
def do_run(RUN_OPT=0,in_fields=None,int_prams=None,real_prams=None):

   import numpy as np
   import os
   import sys

   dd   = os.path.abspath("..")
   sys.path.append(dd+"/bin")
   sys.path.append(dd+"/misc_py")

   import WIM2d_f2py as Mwim # fortran code compiled with f2py
   import fns_get_data as Fdat
   import fns_plot_data as Fplt

   run_dict = {0: 'old version (no in/out)',
               1: 'look at saved results of no in/out run',
               2: 'in/out',
               3: 'look at saved results of in/out run'}
   print("###################################################")
   print("Run option:")
   print(" "+str(RUN_OPT)+' : '+run_dict[RUN_OPT])
   print("###################################################")

   # check directories for outputs exist
   if (RUN_OPT < 2):
      outdir   = 'out'
   else:
      outdir   = 'out_io'

   figdir   = outdir+'/figs'
   dirs  = [outdir,outdir+'/log',
            outdir+'/binaries',outdir+'/binaries/prog',
            figdir,figdir+'/init',figdir+'/final']
   for j in range(0,len(dirs)):
      dirj  = dirs[j]
      if not os.path.exists(dirj):
         os.makedirs(dirj)

   if RUN_OPT%2 is 0:
      # if RUN_OPT is 0 or 2
      # clear out old progress files
      dd    = os.path.abspath(outdir+"/binaries/prog")
      files = os.listdir(dd)
      for f in files:
         os.remove(dd+"/"+f)
   else:
      # load out_fields
      print(" ")
      print("###################################################")
      print("Reading results from: "+outdir)
      print("###################################################")
      print(" ")
      out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   if RUN_OPT is 0:
      # run the "dumb" WIM
      print(" ")
      print("###################################################")
      print("Running WIM without input/output")
      print("###################################################")
      print(" ")

      Mwim.py_wim2d_run()
      
      print(" ")
      print("###################################################")
      print("Finished call to wim2d_run:")
      print("###################################################")
      print(" ")

      # load results from binary files
      out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   elif RUN_OPT is 2:
      # run wim2d with inputs and outputs

      ##########################################################
      SCATMOD     = 1
      ADV_DIM     = 2
      ADV_OPT     = 2
      CHECK_FINAL = 1
      CHECK_PROG  = 1
      CHECK_INIT  = 1
      DO_BREAKING = 1
      STEADY      = 1
      #
      int_prams_def  = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                                 CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                                 DO_BREAKING,STEADY])

      if int_prams is None:
         int_prams   = int_prams_def
      elif len(int_prams)!=len(int_prams_def):
         print('Length of int_prams = '+str(len(int_prams)))
         print('- should be: '+str(len(int_prams_def)))
         sys.exit('run_WIM2d.py, line 96')
      ##########################################################

      ##########################################################
      young          = 2.0e9     # Young's modulus [Pa]
      visc_rp        = 13.0      # Robinson-Palmer damping parameter [Pa/(m/s)] : ~13.0
      duration_hours = 17.77     # length of simulation [h]
      #
      duration       = duration_hours*60*60 # [s]
      real_prams_def = np.array([young,visc_rp,duration])

      if real_prams is None:
         real_prams  = real_prams_def
      elif len(real_prams)!=len(real_prams_def):
         print('Length of real_prams = '+str(len(real_prams)))
         print('- should be: '+str(len(real_prams_def)))
         sys.exit('run_WIM2d.py, line 116')
      ##########################################################

      ##########################################################
      if in_fields is not None:
         # 'in_fields' is given as input
         # - put data into 'in_arrays':
         keys  = ['icec','iceh','dfloe','Hs','Tp','mwd']

         nx,ny       = in_fields[keys[0]].shape
         in_arrays   = np.zeros((nx,ny,6))

         n  = 0
         for key in keys:
            in_arrays[:,:,n]  = in_fields[key]
            n                 = n+1

         del in_fields

      ##########################################################
      elif 1:
         # 'in_fields' not given as input
         # - read in inputs from saved files:
         in_dir                  = 'out/binaries'
         grid_prams              = Fdat.fn_check_grid(in_dir)
         ice_fields,wave_fields  = Fdat.fn_check_init(in_dir)

         # merge ice and wave fields:
         ice_fields.update(wave_fields)
         in_fields   = ice_fields
         del wave_fields

         # put data into 'in_arrays':
         keys  = ['icec','iceh','dfloe','Hs','Tp','mwd']

         nx,ny       = in_fields[keys[0]].shape
         in_arrays   = np.zeros((nx,ny,6))

         n  = 0
         for key in keys:
            in_arrays[:,:,n]  = in_fields[key]
            n                 = n+1

         del in_fields,ice_fields

      ##########################################################
      elif 1:
         # 'in_fields' not given as input
         # - specify 'in_arrays' manually

         in_dir      = 'out/binaries'
         grid_prams  = Fdat.fn_check_grid(in_dir)

         n        = grid_prams['nx']
         ny       = grid_prams['ny']
         X        = grid_prams['X']
         Y        = grid_prams['Y']
         LANDMASK = grid_prams['LANDMASK']
         xmin     = X.min()
         xmax     = X.max()

         # set ice conc/thickness
         ICE_MASK = np.zeros((nx,ny))
         ICE_MASK[np.logical_and(X>0.7*xmin,LANDMASK<1)]  = 1 # i>=24
         icec  = 0.75*ICE_MASK
         iceh  = 2.0*ICE_MASK
         dfloe = 250*ICE_MASK

         # set wave fields
         WAVE_MASK   = np.zeros((nx,ny))
         WAVE_MASK[X<xmin*0.8]   = 1   # i<=15
         Hs    = 2.0*WAVE_MASK
         Tp    = 12.0*WAVE_MASK
         mwd   = -90.0*WAVE_MASK

         in_arrays[:,:,0]  = icec
         in_arrays[:,:,1]  = iceh
         in_arrays[:,:,2]  = dfloe
         in_arrays[:,:,3]  = Hs
         in_arrays[:,:,4]  = Tp
         in_arrays[:,:,5]  = mwd
      ##########################################################

      # run the WIM
      print(" ")
      print("###################################################")
      print("Running WIM with input/output")
      print("###################################################")
      print(" ")
      
      out_arrays  = Mwim.py_wim2d_run_io(in_arrays,int_prams,real_prams)

      print(" ")
      print("###################################################")
      print("Finished call to wim2d_run_io:")
      print("###################################################")
      print(" ")

      # Dmax  = out_arrays[:,:,0]
      # Hs    = out_arrays[:,:,1]
      # Tp    = out_arrays[:,:,2]
      # tau_x = out_arrays[:,:,3]
      # tau_y = out_arrays[:,:,4]

      # convert out_arrays to Out_Fields object
      out_fields  = Fdat.fn_check_out_arr(out_arrays)
      del out_arrays
   #############################################################

   return out_fields,outdir
################################################################

################################################################
def do_run_vSdir(RUN_OPT=0,sdf_dir=None,wave_fields=None,ice_fields=None,int_prams=None,real_prams=None):

   import numpy as np
   import os
   import sys

   dd   = os.path.abspath("..")
   sys.path.append(dd+"/bin")
   sys.path.append(dd+"/py_funs")

   import WIM2d_f2py    as Mwim # fortran code compiled with f2py
   import fns_get_data  as Fdat
   import fns_plot_data as Fplt
   import fns_misc      as Fmisc

   run_dict = {0: 'run "vSdir", passing directional spectrum in/out',\
               1: 'look at saved results of "vSdir", passing directional spectrum in/out'}

   print("###################################################")
   print("Run option:")
   print(" "+str(RUN_OPT)+' : '+run_dict[RUN_OPT])
   print("###################################################")

   # check directories for outputs exist
   outdir   = 'out_2'
   figdir   = outdir+'/figs'
   dirs     = [outdir,
               outdir+'/log',
               outdir+'/binaries',
               outdir+'/binaries/prog',
               figdir,
               figdir+'/init',
               figdir+'/final']

   for j in range(0,len(dirs)):
      dirj  = dirs[j]
      if not os.path.exists(dirj):
         os.makedirs(dirj)

   if RUN_OPT%2 is 0:
      # if RUN_OPT is 0 or 2
      # clear out old progress files
      dd    = os.path.abspath(outdir+"/binaries/prog")
      files = os.listdir(dd)
      for f in files:
         os.remove(dd+"/"+f)
   else:
      # load out_fields
      print(" ")
      print("###################################################")
      print("Reading results from: "+outdir)
      print("###################################################")
      print(" ")
      out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   if RUN_OPT is 0:
      # run wim2d with inputs and outputs

      ##########################################################
      SCATMOD     = 1
      ADV_DIM     = 2
      ADV_OPT     = 2
      CHECK_FINAL = 1
      CHECK_PROG  = 1
      CHECK_INIT  = 1
      DO_BREAKING = 1
      STEADY      = 1
      #
      int_prams_def  = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                                 CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                                 DO_BREAKING,STEADY])

      if int_prams is None:
         int_prams   = int_prams_def
      elif len(int_prams)!=len(int_prams_def):
         print('Length of int_prams = '+str(len(int_prams)))
         print('- should be: '+str(len(int_prams_def)))
         sys.exit('run_WIM2d.py, line 305')
      ##########################################################

      ##########################################################
      young          = 2.0e9     # Young's modulus [Pa]
      visc_rp        = 13.0      # Robinson-Palmer damping parameter [Pa/(m/s)] : ~13.0
      duration_hours = 17.77     # length of simulation [h]
      #
      duration       = duration_hours*60*60 # [s]
      real_prams_def = np.array([young,visc_rp,duration])

      if real_prams is None:
         real_prams  = real_prams_def
      elif len(real_prams)!=len(real_prams_def):
         print('Length of real_prams = '+str(len(real_prams)))
         print('- should be: '+str(len(real_prams_def)))
         sys.exit('run_WIM2d.py, line 321')
      ##########################################################

      # dimensions of grid
      nx,ny,nfreq,ndir  = Mwim.wim2d_dimensions()
      freq_vec          = Mwim.wim2d_freq_vec(nfreq)

      ##########################################################
      if ice_fields is not None:
         # 'ice_fields' is given as input
         # - put data into 'ice_arrays':
         keys        = ['icec','iceh','dfloe']
         ice_arrays  = np.zeros((nx,ny,len(keys)))

         n  = 0
         for key in keys:
            ice_arrays[:,:,n] = ice_fields[key]
            n                 = n+1

         del ice_fields
      else:
         raise ValueError("No 'ice_fields' input")

      ##########################################################

      ##########################################################
      if sdf_dir is None:
         if wave_fields is not None:
            # 'wave_fields' is given as input
            # instead of sdf_dir
            Hs    = wave_fields['Hs']
            Tp    = wave_fields['Tp']
            mwd   = wave_fields['mwd']
            #
            Tmean    = np.mean(Tp [Hs>0])
            dir_mean = np.mean(mwd[Hs>0])

            ##########################################################
            if nfreq==1:
               # check only 1 freq
               stdev = np.std( Tp[Hs>0])
               if stdev/Tmean>1.e-3:

                  print('Expected one frequency, but')
                  print('std dev (Tp) = '+str(stdev)+'s')
                  print('mean    (Tp) = '+str(Tmean) +'s')
                  
                  raise ValueError('Tp array not consistent with 1 frequency')
            ##########################################################

            ##########################################################
            if ndir==1:
               # check only 1 dir
               stdev = np.std (mwd[Hs>0])
               if stdev/dir_mean>1.e-3:

                  print('Expected one direction, but')
                  print('std dev (mwd) = '+str(stdev)+'degrees')
                  print('mean    (mwd) = '+str(dir_mean) +'degrees')
                  
                  raise ValueError('mwd array not consistent with 1 direction')
            ##########################################################

            ################################################################
            sdf_dir  = np.zeros((nx,ny,ndir,nfreq))
            PI       = np.pi
            for i in range(nx):
               for j in range(ny):
                  if Hs[i,j]>0:

                     #########################################################
                     if nfreq==1:
                        # use single freq formula: S=A^2/2
                        sdf_dir[i,j,:,:]  = pow(Hs[i,j]/4.0,2)
                     else:
                        # use Bretschneider
                        for w in range(nfreq):
                           om                = 2*PI*freq_vec(w)
                           sdf_dir[i,j,:,w]  = Fmisc.SDF_Bretscneider(om)
                     #########################################################

                     #########################################################
                     if ndir>1:
                        theta_max   = 90.0
                        theta_min   = -270.0
                        dtheta      = (theta_min-theta_max)/float(ndir) #<0
                        
                        # modify spectrum depending on direction
                        for wth in range(ndir):
                           wavdir   = theta_max+wth*dtheta
                           chi      = PI/180.0*(wavdir-mwd[i,j])

                           if (np.cos(chi)>0.0):
                              theta_fac   = 2.0/PI*pow(np.cos(chi),2)
                           else:
                              theta_fac   = 0.0

                           sdf_dir[i,j,wth,:]   = sdf_dir[i,j,wth,:]*theta_fac
                     #########################################################

            del wave_fields # finished setting sdf_dir from wave_fields
         else:
            raise ValueError("No wave inputs (need 'sdf_dir' or 'wave_fields')")
      ############################################################################

      ##########################################################
      # now run the WIM
      print(" ")
      print("###################################################")
      print("Running WIM with sdf_dir inputted")
      print("###################################################")
      print(" ")
      
      print(Tmean,dir_mean)
      sdf_dir              = sdf_dir.reshape(sdf_dir.size,order='fortran')
      sdf_dir              = np.array(sdf_dir,dtype='float32')
      sdf_dir2,out_arrays  = Mwim.py_wim2d_run_vsdir(sdf_dir,ice_arrays,\
                                                     int_prams,real_prams,Tmean,dir_mean,\
                                                     ndir,nfreq)
      sdf_dir2             = sdf_dir2.reshape((nx,ny,ndir,nfreq),order='fortran')

      print(" ")
      print("###################################################")
      print("Finished call to wim2d_vSdir:")
      print("###################################################")
      print(" ")

      # Dmax  = out_arrays[:,:,0]
      # Hs    = out_arrays[:,:,1]
      # Tp    = out_arrays[:,:,2]
      # tau_x = out_arrays[:,:,3]
      # tau_y = out_arrays[:,:,4]
      # ?? mwd   = out_arrays[:,:,5] ??

      # convert out_arrays to Out_Fields object
      out_fields  = Fdat.fn_check_out_arr(out_arrays)
      del out_arrays
   #############################################################

   return out_fields,outdir
################################################################
