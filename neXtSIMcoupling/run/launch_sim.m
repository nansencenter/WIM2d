DO_COMPILE     = 0;
USE_WIM        = 0;%%use waves
PLOT_STEPS     = 1;%%plot all steps after run
DIAGNOSTICS    = 1;%%diagnostics at end of run
test_and_exit  = 0;%%if 1, exit after 1 call to WIM (if USE_WIM==1)

if 0
   %% idealised domain (simplesquare)
   %% - doesn't need forcing files on johansen
   test_i   =  2;
elseif 0
   %% idealised domain
   %% - doesn't need forcing files on johansen
   test_i   =  16;
elseif 1
   %% idealised domain
   %% - doesn't need forcing files on johansen
   %% - use WIM
   test_i   = 16;
end


% --------------
% 1. compile mex files
if DO_COMPILE
   compile_mex_files;
end

% --------------
% 2. create the standard simul_in
[saved_simul_in,simul_in_name,domain,resol]=create_simul_in(test_i);

if USE_WIM
   create_simul_in_wim(saved_simul_in);
end

if 0
   % Rewrite the simul in file
   % (use restart)
   step_restart           = 12;
   simul_in.use_simul_out = 1;
   simul_in.step_nb       = step_restart;
   save(saved_simul_in,'simul_in')
end

if 1
   % Rewrite the simul in file
   simul_in = load(saved_simul_in);
   simul_in = simul_in.simul_in;

   if 1
      %change length of simulation
      days_in_sec                = 24*3600;
      simul_in.duration          = 2*days_in_sec;
   end

   if 0
      %change wind
      simul_in.constant_u  = -10;
      simul_in.constant_v  = 0;
   end

   if 1
      %change ice
      simul_in.init_concentration   = .7; %Initial ice 
   end

   if USE_WIM==1
      %add waves
      simul_in.wim.use_wim          = 1;
      simul_in.wim.MEX_OPT          = 1;
      simul_in.wim.DAMAGE_OPT       = 1;
      simul_in.wim.wim_break_damage = 0.999;
      simul_in.wim.coupling_option  = 2;
      simul_in.wim.test_and_exit    = test_and_exit;
      simul_in.wim.coupling_freq    = 20*simul_in.timestep;
      %simul_in.wim.coupling_freq    = 500*simul_in.timestep; %%long enough to get breaking in 1 call

      if strfind(simul_in.domain,'wim_grid')
         simul_in.wim.init_waves = 1;
      elseif strfind(simul_in.domain,'squaresmall')
         simul_in.wim.init_waves = 0;
      end

      simul_in.wim.init.Hs          = 4;
      simul_in.wim.init.Tp          = 12;
      simul_in.wim.init.mwd         = -90;
      simul_in.wim.int_prams.STEADY = 1;

      if simul_in.wim.MEX_OPT == 0
         simul_in.wim.single_freq   = 1;
         if simul_in.wim.single_freq==1;
            %% single frequency
            %% - not used if using mex function
            %% - need to pre-compile sizes into fortran code
            %%   (read wave_info.h)
            simul_in.wim.init.ndir  = 16;
            simul_in.wim.init.nfreq = 1;
            simul_in.wim.init.Tmin  = simul_in.wim.init.Tp;
            simul_in.wim.init.Tmax  = simul_in.wim.init.Tp;
         else
            %% wave spectrum
            %% - not used if using mex function
            %% - need to pre-compile sizes into fortran code
            %%   (read wave_info.h)
            simul_in.wim.init.ndir  = 16;
            simul_in.wim.init.nfreq = 25;
            simul_in.wim.init.Tmin  = 1/.042;
            simul_in.wim.init.Tmax  = 1/.4;
         end
      end
   end
   save(saved_simul_in,'simul_in');
   clear simul_in;
end


% --------------
% 3. run the simulation

if test_and_exit==1
   neXtSIM(saved_simul_in,1)
   return;
end

profile on
neXtSIM(saved_simul_in,1)
profile off
profsave(profile('info'),['test_',num2str(test_i),'_profile_results'])

% --------------
% 4. Plots of the scalar variables
load(saved_simul_in);
meshfile = getfield(simul_in,'meshfile');
domain   = getfield(simul_in,'domain');

% steps to be loaded
from_step = 0 ;
saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in.simul_in_name '_step' num2str(from_step) '.mat']
plot_param_v2('c',saved_simul_out,domain,'rev_gris',[0 1],[],'png')
plot_param_v2('h',saved_simul_out,domain,'jet',[0 5],[],'png')

to_step   = length(dir(['simul_out_' meshfile(1:end-4) '_' simul_in.simul_in_name '_step*.mat']))-1 ;
saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in.simul_in_name '_step' num2str(to_step) '.mat']
plot_param_v2('c',saved_simul_out,domain,'rev_gris',[0 1],[],'png')
plot_param_v2('h',saved_simul_out,domain,'jet',[0 5],[],'png')
plot_param_v2('log1md',saved_simul_out,domain,'jet',[-3.7 0],[],'png')

% --------------
% 5. Diagnostics
if DIAGNOSTICS
   filename=['diagnostics_' simul_in.meshfile(1:end-4) '_' simul_in.simul_in_name '.txt'];
   plot_diagnostics(filename);
end

% --------------
% 6. post_processing the output
%  [order]=main_script(simul_in_name,domain,resol)

% --------------
% 7. Plots of the dynamics
%  masked=0;
%  multifractal=0;

%  defo_name=['defo_mod_' simul_in_name '.mat'];
%  disp(defo_name)

%  plots_def_distrib({defo_name},{simul_in_name},order+2);
%  plots_def_multiscale({defo_name},{simul_in_name},[1,1],order+2,domain,multifractal);
%  plots_def_map({defo_name},{simul_in_name},domain,masked);

%  plots_vel_distrib({defo_name},{simul_in_name},[1:2:20],masked);
%  plots_vel_map({defo_name},{simul_in_name},domain,masked);

%%clean up directory
outdir   = setup_outdir(test_i,'mv');

if PLOT_STEPS & USE_WIM
   close all;
   figure(101);
   plot_steps_grid(outdir);
   plot_steps_mesh(outdir);

   %%gifs on grid
   cmd   = ['!../tools/make_gifs.sh ',outdir,' 0'];
   disp(cmd);
   eval(cmd);

   %%gifs on mesh
   cmd   = ['!../tools/make_gifs.sh ',outdir,' 1'];
   disp(cmd);
   eval(cmd);
elseif PLOT_STEPS
   close all;
   figure(101);
   plot_steps_mesh(outdir);

   %%gifs on mesh
   cmd   = ['!../tools/make_gifs.sh ',outdir,' 1'];
   disp(cmd);
   eval(cmd);
end
