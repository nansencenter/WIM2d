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
%%          h: [760x1 double]
%%     Nfloes: [760x1 double]
%% ============================================================

% %%check params_in has the needed fields
% check_params_in_mex(params_in);

%% check if we want to do breaking on the mesh also
if ~exist('mesh_e','var')
   mesh_e  = [];
end
if ~exist('wave_stuff','var')
   wave_stuff  = [];
end

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
   out_arrays  = WIM2d_run_io_mex_v2(in_arrays(:),...
                  params_in.int_prams,params_in.real_prams);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end

   % delete annoying file
   !rm -f fort.6
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
         wave_stuff.dir_spec(:),in_arrays(:),...
         params_in.int_prams,params_in.real_prams,T_init,dir_init);
   wave_stuff.dir_spec  = reshape(wave_stuff.dir_spec,shp);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end
  
   % delete annoying file
   !rm -f fort.6
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
   if TEST_MESH_INTERP
      %%test mesh inputs:
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
         if 1
            subplot(2,2,j);
            plot(xm0/1e3,mesh_e(:,j+2));
            hold on;
            plot(gridprams.X(:,1)/1e3,PP{j},'--g');
         end
      end
      clear PP jp;
   end

   %% get mesh variables
   nmesh_e     = length(mesh_e.xe);
   nmesh_vars  = length(fieldnames(mesh_e));
   mesh_arr    = [mesh_e.xe,mesh_e.ye,mesh_e.c,mesh_e.h,mesh_e.Nfloes,mesh_e.broken];
   
   %% make the call!
   tic;
   shp   = size(wave_stuff.dir_spec);
   [wave_stuff.dir_spec,out_arrays,mesh_arr] =...
      WIM2d_run_io_mex_vSdir_mesh(...
         wave_stuff.dir_spec(:),in_arrays(:),mesh_arr(:),...
         params_in.int_prams,params_in.real_prams,T_init,dir_init,nmesh_e);
   wave_stuff.dir_spec  = reshape(wave_stuff.dir_spec,shp);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[gridprams.nx,gridprams.ny,Nout]);
   for j=1:Nout
      out_fields.(fldnames{j})   = out_arrays(:,:,j);
   end

   %nmesh_e,nmesh_vars
   %[length(mesh_arr),nmesh_e*nmesh_vars]
   mesh_arr = reshape(mesh_arr,[nmesh_e,nmesh_vars]);
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
      figure(101);
      Nfloes_mesh    = mesh_arr(:,5);
      Dmax_mesh      = 0*mesh_arr(:,5);
      jp             = find(Nfloes_mesh>0);
      Dmax_mesh(jp)  = sqrt(mesh_arr(jp,3)./Nfloes_mesh(jp));
      plot(xm0/1e3,Dmax_mesh);
      hold on;
      plot(gridprams.X(:,nmy)/1e3,out_fields.Dmax(:,nmy),'--g');
      fn_fullscreen;

      figure(102);
      pcolor(gridprams.X/1e3,gridprams.Y/1e3,out_fields.Dmax);
      colorbar;
      caxis([0 300]);

      error('Finished test of mesh interpolation');
   end

   %% recreate mesh_e
   fields   = {'Nfloes','broken'};
   for j=1:2
      fld            = fields{j};
      mesh_e.(fld)   = mesh_arr(:,4+j);
   end

   % delete annoying file
   !rm -f fort.6

   return;%%MEX_OPT==3
end%%choose mex function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=avg(x)
y=.5*(x(1:end-1)+x(2:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
