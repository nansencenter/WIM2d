%% /Users/twilliams/GITHUB-REPOSITORIES/WIM2d/fortran/run/run_WIM2d.m
%% Author: Timothy Williams
%% Date:   20141211, 17:49:41 CET
function out_fields  = run_WIM2d_io_mex(ice_fields,wave_fields);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make directories in which to put outputs:
dirs  = {'out_io',
         'out_io/log',
         'out_io/binaries',
         'out_io/binaries/prog'};

for j=1:length(dirs)
   dd = dirs{j};
   if ~exist(dd)
      mkdir(dd);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
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
   %%
   SOLVER      = 1;
   ADV_DIM     = 2;
   int_prams   = [SOLVER,ADV_DIM];
end

%% Call to mex-function:
%  [Dmax,Hs,Tp,taux,tauy]=...
%     WIM2d_run_io_mex(c,h,Dmax,Hs,Tp,mwd);
disp('Calling mex function WIM2d_run_io_mex...')
[out_fields.Dmax,out_fields.Hs,out_fields.Tp,...
   out_fields.taux,out_fields.tauy]=...
      WIM2d_run_io_mex(ice_fields.cice,ice_fields.hice,ice_fields.Dmax,...
         wave_fields.Hs,wave_fields.Tp,wave_fields.mwd,int_prams);
if 0
   disp('Plotting results...');
   fn_plot_final(out_fields);
else
   %%test i-o:
   %%(compile with call to mex_io_gate_test)
   disp('Testing results...');
   d1 = norm(abs(out_fields.Dmax-(ice_fields.cice).^2))
   d2 = norm(abs(out_fields.Hs-(ice_fields.hice).^2))
   d3 = norm(abs(out_fields.Tp-(ice_fields.Dmax).^2))
   d4 = norm(abs(out_fields.taux-(wave_fields.Hs).^2-(wave_fields.Tp).^2-(wave_fields.mwd).^2))
   {out_fields.tauy,int_prams}
   d5 = norm(abs(out_fields.tauy-(int_prams).^2))
end
