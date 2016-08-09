import numpy as np


################################################################
def do_run(in_fields=None,params_in={}):

   import numpy as np
   import os
   import sys

   w2d   = os.getenv('WIM2D_PATH')
   sys.path.append(w2d+"/py_funs")
   dd = os.path.abspath(".")
   sys.path.append(dd+"/bin")

   import WIM2d_f2py as Fwim # fortran code compiled with f2py
   import fns_get_data as Fdat
   import fns_plot_data as Fplt

   if in_fields is not None:
      RUN_OPT  = 1
   else:
      RUN_OPT  = 0
      if len(params_in)>0:
         import warnings
         warnings.warn('params_in input being ignored - no input arrays')

   run_dict = {0: 'old version (no in/out)',
               1: 'in/out'}

   print("###################################################")
   print("Run option:")
   print(" "+str(RUN_OPT)+' : '+run_dict[RUN_OPT])
   print("###################################################")

   # check directories for outputs exist
   if (RUN_OPT == 0):
      outdir   = 'out_py'
   else:
      outdir   = 'out_py_io_1'

   dirs     = [outdir,
               outdir+'/diagnostics',
               outdir+'/diagnostics/local',
               outdir+'/diagnostics/global',
               outdir+'/binaries',
               outdir+'/binaries/prog']

   # tell fortran where to save the outputs
   ifd_name = 'infile_dirs.txt'
   fid      = open(ifd_name,'w')
   fid.write('grid\n')
   fid.write(outdir+'\n')
   fid.close()

   for dirj in dirs:
      if not os.path.exists(dirj):
         os.makedirs(dirj)

   # clear out old progress files
   dd    = os.path.abspath(outdir+"/binaries/prog")
   files = os.listdir(dd)
   for f in files:
      os.remove(dd+"/"+f)

   #############################################################
   if RUN_OPT == 0:
      # run the non-IO WIM 
      # (can still specify parameters in infile_nonIO.txt)
      print(" ")
      print("###################################################")
      print("Running WIM without input/output")
      print("###################################################")
      print(" ")

      Fwim.py_wim2d_run()
      os.remove(ifd_name)

      print(" ")
      print("###################################################")
      print("Finished call to wim2d_run:")
      print("###################################################")
      print(" ")

      # load results from binary files
      return Fdat.wim_results(outdir)
      # out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   else:
      # run wim2d with inputs and outputs

      ##########################################################
      param_dict  = default_params()

      for key in params_in:
         if key not in param_dict:
            raise ValueError('Parameter '+key+' invalid')
         else:
            param_dict[key]   = params_in[key]

      param_vec   = param_dict2vec(param_dict)
      if 0:
         print(params_in)
         print(param_dict)
         print(param_vec)
         sys.exit()
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
      elif 0:
         # 'in_fields' not given as input
         # - read in inputs from saved files:
         in_dir                  = 'out_py/binaries'
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

      out_arrays  = Fwim.py_wim2d_run_io(in_arrays,param_vec)
      os.remove(ifd_name)

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
      # mwd   = out_arrays[:,:,4]

      # convert out_arrays to dictionary
      out_fields  = Fdat.fn_check_out_arr(out_arrays)
      del out_arrays
   #############################################################

   return out_fields,Fdat.wim_results(outdir)
################################################################

################################################################
def do_run_vSdir(sdf_dir=None,\
      wave_fields=None,ice_fields=None,\
      params_in={},mesh_e=None):

   import numpy as np
   import os
   import sys

   dd   = os.path.abspath("..")
   sys.path.append(dd+"/bin")
   sys.path.append(dd+"/py_funs")

   import WIM2d_f2py    as Fwim # fortran code compiled with f2py
   import fns_get_data  as Fdat
   import fns_plot_data as Fplt
   import fns_misc      as Fmisc

   RUN_OPT  = 0
   run_dict =\
      {0: 'run "vSdir", passing directional spectrum in/out'}

   print("###################################################")
   print("Run option:")
   print(" "+str(RUN_OPT)+' : '+run_dict[RUN_OPT])
   print("###################################################")

   TEST_IN_SPEC   = 0 # check Hs,Tp,mwd calc'd from input spec
   TEST_OUT_SPEC  = 0 # check Hs,Tp     calc'd from output spec

   # check directories for outputs exist
   if mesh_e is None:
      outdir = 'out_py_io_2'
   else:
      outdir = 'out_py_io_3'

   dirs  = [outdir,
            outdir+'/diagnostics',
            outdir+'/diagnostics/global',
            outdir+'/diagnostics/local',
            outdir+'/binaries',
            outdir+'/binaries/prog']

   # tell fortran where to save the outputs
   ifd_name = 'infile_dirs.txt'
   fid      = open(ifd_name,'w')
   fid.write('grid\n')
   fid.write(outdir+'\n')
   fid.close()

   for j in range(0,len(dirs)):
      dirj  = dirs[j]
      if not os.path.exists(dirj):
         os.makedirs(dirj)


   # clear out old progress files
   dd    = os.path.abspath(outdir+"/binaries/prog")
   files = os.listdir(dd)
   for f in files:
      os.remove(dd+"/"+f)

   #############################################################
   # run wim2d with inputs and outputs

   ##########################################################
   param_dict  = default_params()

   for key in params_in:
      if key not in param_dict:
         raise ValueError('Parameter '+key+'invalid')
      else:
         param_dict[key]   = params_in[key]

   param_vec   = param_dict2vec(param_dict)
   ##########################################################

   # dimensions of grid
   nx,ny,nfreq,ndir  = Fwim.wim2d_dimensions()
   freq_vec          = Fwim.wim2d_freq_vec(nfreq)

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
         Tmean           = np.mean(Tp [Hs>0])
         dir_mean        = np.mean(mwd[Hs>0])

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
                        om                = 2*PI*freq_vec[w]
                        sdf_dir[i,j,:,w]  = Fmisc.SDF_Bretschneider(om)
                  #########################################################

                  #########################################################
                  if ndir>1:
                     theta_max   = 90.0
                     theta_min   = -270.0
                     dtheta      = (theta_min-theta_max)/float(ndir) #<0
                     wavdir      = theta_max+np.arange(ndir)*dtheta

                     # modify spectrum depending on direction
                     D2R   = np.pi/180.
                     # test_sum = 0.
                     # test_sum_ = 0.
                     for wth in range(ndir):
                        theta_fac_   = Fmisc.theta_dirfrac(wavdir[wth]-.5*abs(dtheta),abs(dtheta),mwd[i,j])/abs(D2R*dtheta)
                        # chi      = PI/180.0*(wavdir[wth]-mwd[i,j])
                        # if (np.cos(chi)>0.0):
                        #    theta_fac   = 2.0/PI*pow(np.cos(chi),2)
                        # else:
                        #    theta_fac   = 0.0

                        sdf_dir[i,j,wth,:]   = sdf_dir[i,j,wth,:]*theta_fac_
                        # test_sum += theta_fac_
                        # test_sum_ += theta_fac_
                        # print(wavdir[wth],wavdir[wth]-.5*abs(dtheta),dtheta,mwd[i,j],theta_fac,theta_fac_)

                  # print(test_sum*abs(dtheta),test_sum_*abs(dtheta))
                  # sys.exit()
                  #########################################################

                  if 0:
                     # test spectrum is defined OK
                     sq  = np.squeeze(sdf_dir[i,j,:,0])
                     m0  = np.sum(sq)*abs((PI/180.)*dtheta)
                     print(i,j)
                     print(4*np.sqrt(m0),Hs[i,j])
                     #
                     import matplotlib.pyplot as plt
                     plt.plot(wavdir/180.,sq)
                     plt.show()
                     #
                     sys.exit()

         ###############################################################
         if TEST_IN_SPEC:
            # test integrals of spectrum are the same as intended:
            print('Testing input spectrum...')
            wave_fields2   = {'Hs':0,'Tp':0,'mwd':0}
            if nfreq==1:
               freq_vec_   = np.array([1./Tp_single_freq])
            else:
               freq_vec_   = freq_vec
            if ndir==1:
               wavdir_   = np.array([mwd_single_dir])
            else:
               wavdir_  = wavdir

            wave_fields2['Hs'],wave_fields2['Tp'],wave_fields2['mwd']   = \
                  Fmisc.spectrum_integrals(sdf_dir,freq_vec_,wavdir_)

            arrays   = [wave_fields2,wave_fields]
            keys     = arrays[0].keys()
            for key in keys:
               print('Comparing field: '+key)
               diff  = abs(arrays[0][key]-arrays[1][key])
               print('Maximum difference = '+str(np.max(diff))+'\n')
            sys.exit()
         ###############################################################

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

   sdf_dir  = sdf_dir.reshape(sdf_dir.size,order='fortran')
   sdf_dir  = np.array(sdf_dir,dtype='float32')
   
   print(param_vec)
   if mesh_e is None:
      sdf_dir2,out_arrays  = Fwim.py_wim2d_run_vsdir(
            sdf_dir,ice_arrays,param_vec,ndir,nfreq)
      os.remove(ifd_name)

      print(" ")
      print("###################################################")
      print("Finished call to wim2d_vSdir:")
      print("###################################################")
      print(" ")

   else:

      #################################################################
      # mesh_e=dictionary with keys ['x','y','conc','thick','Nfloes','broken'] 
      # > convert this to array
      nmesh_vars  = len(mesh_e)
      mesh_arr    = mesh_dict2arr(mesh_e,inverse=False)
      nmesh_e     = len(mesh_arr[:,0])
      # print(mesh_arr.transpose())
      # sys.exit()

      # reshape array to a vector
      mesh_arr    = mesh_arr.reshape((nmesh_e*nmesh_vars,),order='fortran')
      #################################################################


      #################################################################
      # run WIM, with interpolation of Nfloes onto the mesh
      out  = Fwim.py_wim2d_run_vsdir_mesh(
            sdf_dir,ice_arrays,mesh_arr,param_vec,ndir,nfreq,nmesh_e)
      sdf_dir2,out_arrays,mesh_arr  = out
      os.remove(ifd_name)

      print(" ")
      print("###################################################")
      print("Finished call to wim2d_vSdir_mesh:")
      print("###################################################")
      print(" ")
      #################################################################

      #################################################################
      mesh_arr = mesh_arr.reshape((nmesh_e,nmesh_vars),order='fortran')
      mesh_e   = mesh_dict2arr(mesh_arr,inverse=True)
      # print(mesh_arr.transpose())
      #################################################################

   sdf_dir2 = sdf_dir2.reshape((nx,ny,ndir,nfreq),order='fortran')


   # Dmax  = out_arrays[:,:,0]
   # Hs    = out_arrays[:,:,1]
   # Tp    = out_arrays[:,:,2]
   # tau_x = out_arrays[:,:,3]
   # tau_y = out_arrays[:,:,4]
   # mwd   = out_arrays[:,:,5]

   # convert out_arrays to dictionary
   out_fields  = Fdat.fn_check_out_arr(out_arrays)
   del out_arrays

   ###############################################################
   if TEST_OUT_SPEC:
      # test integrals of spectrum are the same as intended:
      print('Testing integrals of output spectrum...')
      wave_fields2   = {'Hs':0,'Tp':0,'mwd':0}
      if nfreq==1:
         Tp_single_freq = param_dict['T_init']
         freq_vec_   = np.array([1./Tp_single_freq])
      else:
         freq_vec_   = freq_vec

      if ndir==1:
         mwd_single_dir = param_dict['mwd_init']
         wavdir_   = np.array([mwd_single_dir])
      else:
         wavdir_  = wavdir

      wave_fields2['Hs'],\
            wave_fields2['Tp'],\
            wave_fields2['mwd']   = \
               Fmisc.spectrum_integrals(\
                  sdf_dir2,freq_vec_,wavdir_)

      arrays   = [wave_fields2,out_fields]
      keys     = ['Hs','Tp']
      # keys     = arrays[0].keys()
      for key in keys:
         print('Comparing field: '+key)
         diff  = abs(arrays[0][key]-arrays[1][key])
         print('Maximum difference = '+str(np.max(diff))+'\n')

      if 1:
         # plot Tp for visual comparison
         # - small differences?
         gf    = Fdat.fn_check_grid('inputs')
         labs  = ['$x$, km','$y$, km', '$T_p$, s']
         #
         f  = Fplt.cmap_3d(\
               gf['X']/1.e3,gf['Y']/1.e3,wave_fields2['Tp'],labs)
         f.show()
         #
         f  = Fplt.cmap_3d(\
               gf['X']/1.e3,gf['Y']/1.e3,out_fields['Tp'],labs)
         f.show()

      if 0:
         # plot Hs for visual comparison
         gf    = Fdat.fn_check_grid('inputs')
         labs  = ['$x$, km','$y$, km', '$H_s$, m']
         #
         f  = Fplt.cmap_3d(\
               gf['X']/1.e3,gf['Y']/1.e3,wave_fields2['Hs'],labs)
         f.show()
         #
         f  = Fplt.cmap_3d(\
               gf['X']/1.e3,gf['Y']/1.e3,out_fields['Hs'],labs)
         f.show()

      sys.exit()
   elif 0:
      # test out_fields are the same as in the binary files
      print('Testing output fields against binary files...')
      arrays      = 2*[1]
      arrays[0]   = out_fields
      arrays[1]   = Fdat.fn_check_out_bin('out_2/binaries')
      keys        = arrays[0].keys()
      for key in keys:
         print('Comparing field: '+key)
         diff  = abs(arrays[0][key]-arrays[1][key])
         print('Maximum difference = '+str(np.max(diff))+'\n')
      sys.exit()
   ###############################################################

   #############################################################

   return out_fields,Fdat.wim_results(outdir),mesh_e
################################################################


###################################################
def default_params(convert=False):
   
   param_dict  = {}
   ###################################################
   # default integer parameters
   param_dict.update({'SCATMOD'        : 0})
   param_dict.update({'ADV_DIM'        : 2})
   param_dict.update({'ADV_OPT'        : 2})
   param_dict.update({'DO_CHECK_FINAL' : 1})
   param_dict.update({'DO_CHECK_PROG'  : 1})
   param_dict.update({'DO_CHECK_INIT'  : 1})
   param_dict.update({'STEADY'         : 1})
   param_dict.update({'DO_BREAKING'    : 1})
   param_dict.update({'DO_ATTEN'       : 1})
   ###################################################


   ###################################################
   # default real parameters:
   param_dict.update({'young'    : 5.49e9})
   param_dict.update({'visc_rp'  : 13.0})
   duration_hours = 6.
   param_dict.update({'duration' : duration_hours*60*60})
   param_dict.update({'CFL'      : 0.7})
   ###################################################


   ###################################################
   # other integer params
   param_dict.update({"BRK_OPT"     :1})
   param_dict.update({"FSD_OPT"     :1})
   param_dict.update({"REF_HS_ICE"  :1})
   param_dict.update({"USE_ICE_VEL" :0})
   ###################################################


   ###################################################
   # init cons
   param_dict.update({"Hs_init"     :3.})
   param_dict.update({"T_init"      :12.})
   param_dict.update({"mwd_init"    :-90.})
   param_dict.update({"conc_init"   :.7})
   param_dict.update({"h_init"      :1.})
   param_dict.update({"Dmax_init"   :300.})
   ###################################################


   ###################################################
   # start time
   from datetime import datetime as dtm
   param_dict.update({'start_time':dtm(2015,1,1)})
   ###################################################


   ###################################################
   # diagnostics
   param_dict.update({"itest"    :-1})
   param_dict.update({"jtest"    :-1})
   param_dict.update({"dumpfreq" :10})
   ###################################################

   if convert:
      return param_dict2vec(param_dict)
   else:
      return param_dict
################################################################


# ================================================================
def param_dict2vec(param_dict):
   param_vec   = []

   # old int_prams
   param_vec.append(param_dict["SCATMOD"])
   param_vec.append(param_dict["ADV_DIM"])
   param_vec.append(param_dict["ADV_OPT"])
   param_vec.append(param_dict["BRK_OPT"])
   param_vec.append(param_dict["STEADY"])
   param_vec.append(param_dict["DO_ATTEN"])
   param_vec.append(param_dict["DO_CHECK_FINAL"])
   param_vec.append(param_dict["DO_CHECK_PROG"])
   param_vec.append(param_dict["DO_CHECK_INIT"])

   # old real_prams
   param_vec.append(param_dict["young"])
   param_vec.append(param_dict["visc_rp"])
   param_vec.append(param_dict["duration"])
   param_vec.append(param_dict["CFL"])

   # other integers
   param_vec.append(param_dict["FSD_OPT"])
   param_vec.append(param_dict["REF_HS_ICE"])
   param_vec.append(param_dict["USE_ICE_VEL"])

   # init cons
   param_vec.append(param_dict["Hs_init"])
   param_vec.append(param_dict["T_init"])
   param_vec.append(param_dict["mwd_init"])
   param_vec.append(param_dict["conc_init"])
   param_vec.append(param_dict["h_init"])
   param_vec.append(param_dict["Dmax_init"])

   # start time
   from datetime import datetime as dtm
   t0       = param_dict['start_time']
   reftime  = dtm(1900,1,1)
   tdiff    = t0-reftime
 
   model_day      = tdiff.days
   model_seconds  = tdiff.total_seconds()-24*3600*model_day
   param_vec.append(model_day)
   param_vec.append(model_seconds)

   # diagnostics
   param_vec.append(param_dict["itest"])
   param_vec.append(param_dict["jtest"])
   param_vec.append(param_dict["dumpfreq"])

   return np.array(param_vec)
# ================================================================


# ================================================================
def mesh_dict2arr(mesh_e,inverse=False):
   # mesh_e=dictionary with keys:
   # NB order important
   keys  = ['x','y','conc','thick','Nfloes','broken'] 
   if not inverse:
      out   = []
      for key in keys:
         out.append(1*mesh_e[key])
      out   = np.array(out).transpose()
   else:
      out   = {}
      for i,key in enumerate(keys):
         out.update({key:1*mesh_e[:,i]})

   return out
# ================================================================
