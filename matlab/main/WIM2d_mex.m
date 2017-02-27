function [out_fields,wave_stuff,mesh_e] =...
   WIM2d_mex(params_in,gridprams,ice_fields,wave_fields,wave_stuff,mesh_e)

%duration,ice_prams into params_in
%% params_in.MEX_OPT=1:
%% * inputs : params_in,gridprams,ice_fields,wave_fields
%% * outputs: out_fields
%%
%% params_in.MEX_OPT=2:
%% * inputs : params_in,gridprams,ice_fields,wave_fields,...
%%            wave_stuff (to give Sdir as input)
%% * outputs: out_fields,Sdir
%%
%% params_in.MEX_OPT=3:
%% * inputs : params_in,gridprams,ice_fields,wave_fields,...
%%            wave_stuff (to give Sdir as input),
%%            mesh_e (mesh quantities for interpolation)
%% * outputs: out_fields,Sdir,...
%%            mesh_e (output mesh quantities)
%%
%% ============================================================
%% Inputs:
%%
%% params_in = structure eg:
%%           int_prams: 9x1 vector
%%          real_prams: 4x1 vector
%%             MEX_OPT: 3
%%              DODISP: 1
%%
%% gridprams = structure eg:
%%      cice: [51x51 double]
%%      hice: [51x51 double]
%%      Dmax: [51x51 double]
%%
%%
%% ice_fields  = structure eg:
%%      cice: [51x51 double]
%%      hice: [51x51 double]
%%      Dmax: [51x51 double]
%%
%%
%% wave_fields = structure eg
%%           Hs: [150x10 double]
%%           Tp: [150x10 double]
%%          mwd: [150x10 double]
%%  STEADY_MASK: [150x10 logical]
%%
%%
%% wave_stuff = structure eg
%%        nfreq: 2
%%         ndir: 16
%%         freq: [23x1 double]
%%         dirs: [16x1 double]
%%     dir_spec: [150x20x16x23 double]
%%
%%
%% optional:
%% - give coordinates of FEM mesh centres
%% - do breaking on mesh as well to try to reduce numerical diffusion
%%   caused by interpolation between grid and mesh
%% mesh_e = structure eg
%%         xe: [760x1 double]
%%         ye: [760x1 double]
%%          c: [760x1 double]
%%      thick: [760x1 double]
%%       Dmax: [760x1 double]
%%     broken: [760x1 double]
%% ============================================================

% %%check params_in has the needed fields
% check_params_in_mex(params_in);
RMFORT6  = 1;

%% check if we want to do breaking on the mesh also
if ~exist('mesh_e','var')
   mesh_e  = [];
end
if ~exist('wave_stuff','var')
   wave_stuff  = [];
end

params_vec  = get_param_vec(params_in);
if params_in.MEX_OPT==1

   if params_in.DO_DISP; disp(' ');
      disp('*****************************************************************');
      disp('Running fortran code with mex function: run_WIM2d_io_mex_v2');
      disp('*****************************************************************');
      disp(' ');
   end

   %% input 2d arrays
   in_arrays         = zeros(gridprams.nx,gridprams.ny,6);
   in_arrays(:,:,1)  = ice_fields.cice;
   in_arrays(:,:,2)  = ice_fields.hice;
   in_arrays(:,:,3)  = ice_fields.Dmax;
   in_arrays(:,:,4)  = wave_fields.Hs;
   in_arrays(:,:,5)  = wave_fields.Tp;
   in_arrays(:,:,6)  = wave_fields.mwd;
   %for j=1:6
   %   disp([min(in_arrays(:)),max(in_arrays(:))]);
   %end
   %pause


   %% make the call!
   tic;
   out_arrays  = WIM2d_run_io_mex_v2(in_arrays(:),params_vec);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp','mwd'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end

   % delete annoying file
   if RMFORT6
      !rm -f fort.6
   end
   return%%MEX_OPT==1

elseif params_in.MEX_OPT==2

   if params_in.DO_DISP; disp(' ');
      disp('*****************************************************************');
      disp('Running fortran code with mex function: run_WIM2d_io_mex_vSdir');
      disp('*****************************************************************');
      disp(' ');
   end

   %% input 2d arrays
   in_arrays         = zeros(gridprams.nx,gridprams.ny,3);
   in_arrays(:,:,1)  = ice_fields.cice;
   in_arrays(:,:,2)  = ice_fields.hice;
   in_arrays(:,:,3)  = ice_fields.Dmax;

   % defined in case only one period or dirn
   T_init   = 1/max(wave_stuff.freq)
   dir_init = max(wave_stuff.dirs)

   %% make the call!
   tic;
   shp   = size(wave_stuff.dir_spec);
   [wave_stuff.dir_spec,out_arrays] =...
      WIM2d_run_io_mex_vSdir(...
         wave_stuff.dir_spec(:),in_arrays(:),params_vec);
   wave_stuff.dir_spec  = reshape(wave_stuff.dir_spec,shp);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp','mwd'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end
  
   % delete annoying file
   if RMFORT6
      !rm -f fort.6
   end
   return%%MEX_OPT==2

elseif params_in.MEX_OPT==3

   if params_in.DO_DISP; disp(' ');
      disp('*****************************************************************');
      disp('Running fortran code with mex function: run_WIM2d_io_mex_vSdir_mesh');
      disp('*****************************************************************');
      disp(' ');
   end

   %% input 2d arrays
   in_arrays         = zeros(gridprams.nx,gridprams.ny,3);
   in_arrays(:,:,1)  = ice_fields.cice;
   in_arrays(:,:,2)  = ice_fields.hice;
   in_arrays(:,:,3)  = ice_fields.Dmax;

   % defined in case only one period or dirn
   T_init   = 1/max(wave_stuff.freq);
   dir_init = max(wave_stuff.dirs);

   TEST_MESH_INTERP  = 0;
   FAKE_MESH         = 0;

   %% =============================================================================
   if TEST_MESH_INTERP & ~exist('mesh_e','var')
      %%create test mesh inputs:
      FAKE_MESH      = 1
      gridprams.x0   = min(gridprams.X(:));
      gridprams.y0   = min(gridprams.Y(:));
      %%
      xm0         = (gridprams.x0+gridprams.dx/2)+(0:gridprams.nx-2)*gridprams.dx;
      nmesh_e     = length(xm0);
      nmesh_vars  = 6;
      mesh_e      = zeros(nmesh_e,nmesh_vars);
      mesh_e(:,1) = xm0.';
      if gridprams.ny==1
         nmy   = 1;
      else
         nmy   = ceil(gridprams.ny/2);
      end
      mesh_e(:,2) = gridprams.Y(1,nmy);

      Nfloes      = 0*ice_fields.Dmax(:,nmy);
      jp          = find(ice_fields.Dmax(:,nmy)>0);
      Nfloes(jp)  = ice_fields.cice(jp,nmy)./ice_fields.Dmax(jp,nmy).^2;
      PP          = {ice_fields.cice(:,nmy),ice_fields.hice(:,nmy),...
                        Nfloes,0*ice_fields.Dmax(:,nmy)};
      for j=1:4
         mesh_e(:,j+2)  = avg(PP{j});
         if 0
            subplot(2,2,j);
            plot(xm0/1e3,mesh_e(:,j+2));
            hold on;
            plot(gridprams.X(:,1)/1e3,PP{j},'--g');
         end
      end
      clear PP jp;
      fnames   = {'xe','ye','c','thick','Nfloes','broken'};
      Mesh_e   = mesh_e;
      clear mesh_e;
      for j=1:nmesh_vars
         mesh_e.(fnames{j})   = Mesh_e(:,j);
      end
      clear fnames Mesh_e
   end
   %% =============================================================================

   %% get mesh variables
   nmesh_e     = length(mesh_e.xe);
   nmesh_vars  = length(fieldnames(mesh_e));
   mesh_arr    = [mesh_e.xe,mesh_e.ye,mesh_e.c,mesh_e.thick,mesh_e.Dmax,mesh_e.broken];
   %mesh0 = mesh_arr;

   %% make the call!
   tic;
   shp   = size(wave_stuff.dir_spec);
   [wave_stuff.dir_spec,out_arrays,mesh_arr] =...
      WIM2d_run_io_mex_vSdir_mesh(...
         wave_stuff.dir_spec(:),in_arrays(:),mesh_arr(:),...
         params_vec,nmesh_e);
   wave_stuff.dir_spec  = reshape(wave_stuff.dir_spec,shp);
   toc;


   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp','mwd'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end

   taux_max = max(out_fields.tau_x(:))

   %nmesh_e,nmesh_vars
   %[length(mesh_arr),nmesh_e*nmesh_vars]
   mesh_arr = reshape(mesh_arr,[nmesh_e,nmesh_vars]);
   %save mesh mesh_e mesh0 mesh_arr out_fields;

   %% recreate mesh_e
   mesh_e.Dmax    = mesh_arr(:,5);
   mesh_e.broken  = mesh_arr(:,6);

   if 0
      %look at Nfloes where ice is
      mesh_arr(mesh_out(:,5)>0,5)
   elseif 0
      %%look at where breaking occurred, next to thickness
      %% (proxy for original thickness)
      mesh_arr(:,[4,6])
   elseif 0
      %%look at difference between 1st 4 col's (should be ~0)
      mesh_arr(:,1:4)-mesh_e(:,1:4)
   elseif TEST_MESH_INTERP

      if FAKE_MESH==1
         figure(101);
         Nfloes_mesh    = mesh_arr(:,5);
         Dmax_mesh      = 0*mesh_arr(:,5);
         jp             = find(Nfloes_mesh>0);
         Dmax_mesh(jp)  = sqrt(mesh_arr(jp,3)./Nfloes_mesh(jp));
         plot(xm0/1e3,Dmax_mesh);
         hold on;
         plot(gridprams.X(:,nmy)/1e3,out_fields.Dmax(:,nmy),'--g');
         legend('mesh','grid')
         fn_fullscreen;
      end

      if gridprams.ny>1
         X_ = gridprams.X(:,1)/1e3;
         Y_ = gridprams.Y(1,:)/1e3;

         figure(102); fn_pcolor( X_,Y_,out_fields.Dmax ); colorbar; caxis([0 300]);
         figure(103); fn_pcolor( X_,Y_,out_fields.Hs ); colorbar;
         %figure(104); fn_pcolor( X_,Y_,out_fields.tau_x ); colorbar;
         %figure(105); fn_pcolor( X_,Y_,ice_fields.cice ); colorbar;
         %figure(106); fn_pcolor( X_,Y_,ice_fields.hice ); colorbar;
         %figure(107); fn_pcolor( X_,Y_,gridprams.LANDMASK ); colorbar;
      end

      disp(' ');
      Flds  = fieldnames(ice_fields);
      for j=1:length(Flds)
         v  = Flds{j};
         V  = ice_fields.(v);
         disp(['Range on grid of ',v,' = [',num2str(min(V(:))),',',num2str(max(V(:))),']']);
      end

      disp(' ');
      Flds  = fieldnames(out_fields);
      for j=1:length(Flds)
         v  = Flds{j};
         V  = out_fields.(v);
         disp(['Range on grid of ',v,' = [',num2str(min(V(:))),',',num2str(max(V(:))),']']);
      end

      disp(' ');
      Flds  = fieldnames(mesh_e);
      for j=1:length(Flds)
         v  = Flds{j};
         if strcmp(v,'xe')==0&strcmp(v,'ye')==0
            V  = mesh_e.(v);
            disp(['Range on mesh of ',v,' = [',num2str(min(V(:))),',',num2str(max(V(:))),']']);
         end
      end
      disp(' ');

      %broken   = mesh_arr(:,6)
      %error('Finished test of mesh interpolation');
      %pause
   end

   % delete annoying file
   if RMFORT6
      !rm -f fort.6
   end

   return;%%MEX_OPT==3
end%%choose mex function
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=avg(x)
y=.5*(x(1:end-1)+x(2:end));
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params_vec = get_param_vec(params_mex)
%% 
%% INPUT:
%% params_mex = 
%%            SCATMOD: 3
%%            ADV_DIM: 2
%%            ADV_OPT: 2
%%      DO_CHECK_INIT: 1
%%      DO_CHECK_PROG: 1
%%     DO_CHECK_FINAL: 1
%%             STEADY: 1
%%            BRK_OPT: 1
%%           DO_ATTEN: 1
%%              young: 5.490000000000000e+09
%%            drag_rp: 13
%%            visc_ws: 0
%%           duration: 21600
%%                CFL: 0.700000000000000
%%            FSD_OPT: 1
%%         REF_Hs_ICE: 0
%%        USE_ICE_VEL: 0
%%     TAKE_MAX_WAVES: 0
%%            Hs_init: 3
%%             T_init: 12
%%           dir_init: -90
%%          conc_init: 0.700000000000000
%%             h_init: 1
%%          Dmax_init: 300
%%               Dmin: 20
%%                 xi: 2
%%          fragility: .9
%%            Dthresh: 200
%%           cice_min: 0.05
%%          model_day: 42003
%%      model_seconds: 0
%%              itest: 25
%%              jtest: 5
%%           dumpfreq: 10
%%            MEX_OPT: 1
%%            DO_DISP: 1
%%
%% OUTPUT:
%% vector with some of these fields
%% - for ordering see:
%%   - fortran/infiles/infile_nonIO.txt
%%   - read_params_vec subroutine in fortran/src/main/mod_WIM2d_run.F

%% ================================================
%% old int_prams:
fields   = {};
n        = 9;
fields(end+1:end+n)  = {...
            'SCATMOD',...
            'ADV_DIM',...
            'ADV_OPT',...
            'BRK_OPT',...
            'STEADY',...
            'DO_ATTEN',...
            'DO_CHECK_INIT',...
            'DO_CHECK_PROG',...
            'DO_CHECK_FINAL'};
%% ================================================


%% ================================================
%% old real_prams:
n  = 6;
fields(end+1:end+n)  = {...
            'young',...
            'drag_rp',...
            'visc_ws',...
            'cohesion',...
            'friction',...
            'CFL'};
%% ================================================


%% ================================================
%% other integer parameters (in params):
n  = 4;
fields(end+1:end+n)  = {...
            'FSD_OPT',...
            'REF_Hs_ICE',...
            'USE_ICE_VEL',...
            'TAKE_MAX_WAVES'};
%% ================================================


%% ================================================
%% initial conditions (in params):
n  = 6;
fields(end+1:end+n)  = {...
            'Hs_init',...
            'T_init',...
            'dir_init',...
            'conc_init',...
            'h_init',...
            'Dmax_init'};
%% ================================================


%% ================================================
%% FSD info (in params):
n  = 5;
fields(end+1:end+n)  = {...
            'Dmin',...
            'xi',...
            'fragility',...
            'Dthresh',...
            'cice_min'};
%% ================================================


%% ================================================
%% duration & start time:
n  = 3;
fields(end+1:end+n)  = {...
            'duration',...
            'model_day',... % day relative to 1900-1-1
            'model_seconds'};
%% ================================================


%% ================================================
%% diagnostics:
n  = 3;
fields(end+1:end+n)  = {...
            'itest',...
            'jtest',...
            'dumpfreq'};
%% ================================================

Ni          = length(fields);
params_vec  = zeros(Ni,1);
for j=1:Ni
   cmd   = ['params_vec(',num2str(j),') = params_mex.',fields{j},';'];
   eval(cmd);
   if params_mex.DO_DISP==1
      disp([cmd(1:end-1),' = ',num2str(params_mex.(fields{j}))]);
   end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
