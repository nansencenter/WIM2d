%% /Users/twilliams/GITHUB-REPOSITORIES/WIM2d/fortran/run/run_WIM2d.m
%% Author: Timothy Williams
%% Date:   20141211, 17:49:41 CET
function out_fields  = run_WIM2d_io_mex(ice_fields,wave_fields,int_prams,real_prams);

%% extra parameters
DO_PLOT  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make directories in which to put outputs:
ifil  = 'infile.txt';
if exist(ifil)
   fid   = fopen(ifil);
   dd    = textscan(fid,'%s');
   fclose(fid);
   indir    = dd{1}{1};
   outdir   = dd{1}{2};
else
   indir    = 'out_io';
   outdir   = 'out_io';
end

dirs  = {outdir,
         [outdir,'/log'],
         [outdir,'/binaries'],
         [outdir,'/binaries/prog']};

for j=1:length(dirs)
   dd = dirs{j};
   if ~exist(dd)
      mkdir(dd);
   end
end

if length([dirs{4},'/*.a'])>0
   %%clean prog directory if necessary
   eval(['!rm ',dirs{4},'/*']);
   eval(['!rm -r ',outdir,'/figs/prog*']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('wave_fields','var')
   %% use default inputs
   %% - eg outputs from non-io version
   if 0
      %% run non-io version to compare to
      %% and get eg inputs
      %% (can also run python or shell versions
      %% to get the needed files)
      run_WIM2d_mex;
   end

   %%get init cons from non-io run;
   outdir0  = 'out';
   [grid_prams,ice_fields,wave_fields] = fn_check_init(outdir0);
end

if ~exist('int_prams','var')
   %% integer parameters
   SCATMOD     = 1;
   ADV_DIM     = 2;
   ADV_OPT     = 2;
   CHECK_FINAL = 1;
   CHECK_PROG  = 0;
   CHECK_INIT  = 1;
   DO_BREAKING = 1;
   STEADY      = 1;
   %%
   disp('run_WIM2d_io_mex: using default for int_prams:')
   int_prams   = [SCATMOD,ADV_DIM,ADV_OPT,...
                  CHECK_FINAL,CHECK_PROG,CHECK_INIT,...
                  DO_BREAKING,STEADY]
end

if ~exist('real_prams','var')
   %% real parameters
   young          = 5.49e9;
   visc_rp        = 13;
   duration_hours = 24;
   duration       = duration_hours*60*60;%s
   disp('run_WIM2d_io_mex: using default for real_prams:')
   real_prams     = [young,visc_rp,duration]
end

%% Call to mex-function:
%  [Dmax,Hs,Tp,taux,tauy]=...
%     WIM2d_run_io_mex(c,h,Dmax,Hs,Tp,mwd);
tic;
disp('Calling mex function WIM2d_run_io_mex...')

[out_fields.Dmax,out_fields.tau_x,out_fields.tau_y,...
   out_fields.Hs,out_fields.Tp]=...
      WIM2d_run_io_mex(ice_fields.cice,ice_fields.hice,ice_fields.Dmax,...
         wave_fields.Hs,wave_fields.Tp,wave_fields.mwd,int_prams,real_prams);

toc

%% delete annnoying file sometimes caused by "print*" commands in fortran code
if exist('fort.6','file')
   !rm fort.6
end

if DO_PLOT>0
   if ~exist('grid_prams','var');
      grid_prams  = fn_get_grid(indir);
   end

   disp('Plotting results...');
   figure(103);
   fn_plot_final(grid_prams,out_fields);
   fn_fullscreen;

   %% plot progress files and make movie if they exist
   DD = dir([outdir,'/binaries/prog/wim_prog*a']);
   if length(DD)>0 & DO_PLOT==2
      %% if some progress files exist,
      %% plot them and make movies
      MKMOV = 1;% 1, make movies; 0, don't
      cmd   = ['!./plot_prog.sh ',num2str(MKMOV)];
      eval(cmd);
   end
end

if 0
   %%test i-o:
   %%(compile with call to mex_io_gate_test)
   disp('Testing results...');
   d1 = norm(abs(out_fields.Dmax-(ice_fields.cice).^2))
   d2 = norm(abs(out_fields.Hs-(ice_fields.hice).^2))
   d3 = norm(abs(out_fields.Tp-(ice_fields.Dmax).^2))
   d4 = norm(abs(out_fields.taux-(wave_fields.Hs).^2-(wave_fields.Tp).^2-(wave_fields.mwd).^2))
   if 0
      {out_fields.tauy,int_prams}
      d5 = norm(abs(out_fields.tauy-(int_prams).^2))
   else
      d5 = norm(abs(out_fields.tauy-(wave_fields.Tp).^2))
   end
end
