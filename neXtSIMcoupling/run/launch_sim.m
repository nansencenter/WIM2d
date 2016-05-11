DO_COMPILE  = 0;
USE_WIM     = 0;

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
   USE_WIM  = 1;
end


% --------------
% 1. compile mex files
if DO_COMPILE
   compile_mex_files;
end

% --------------
% 2. create the standard simul_in
[saved_simul_in,simul_in_name,domain,resol]=create_simul_in(test_i);

WAVES = 1;
if WAVES
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
      simul_in.duration          = 1*days_in_sec;
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

   if WAVES==1
      %add waves
      simul_in.wim.use_wim          = 1;
      simul_in.wim.MEX_OPT          = 1;
      simul_in.wim.DAMAGE_OPT       = 1;
      simul_in.wim.wim_break_damage = .99;
      simul_in.wim.coupling_option  = 1;
      simul_in.wim.test_and_exit    = 0;

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
   test_and_exit  = simul_in.wim.test_and_exit;
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
filename=['diagnostics_' simul_in.meshfile(1:end-4) '_' simul_in.simul_in_name '.txt'];
import_data=importdata(filename)

list_absciss_to_plot={'p','p','p'};
list_ordinates_to_plot={{'NnR','NeR'},{'total_volume'},{'variation_total_volume','variation_total_volume_thermo','variation_total_volume_regrid','variation_total_volume_transport'}}

for k=1:length(list_absciss_to_plot)
  absciss_to_plot=list_absciss_to_plot{k};
  ordinates_to_plot=list_ordinates_to_plot{k}
  figure
  
  j=find(strcmp(import_data.textdata,absciss_to_plot));
  if(isempty(j))
      error('no absciss to plot')
  else
      absciss=import_data.data(:,j);
  end
  
  ColorOrder=get(gca,'ColorOrder');
  
  for i=1:length(ordinates_to_plot)
      j=find(strcmp(import_data.textdata,ordinates_to_plot{i}));
      if(~isempty(j))
          ordinate=import_data.data(:,j);
          plot(absciss,ordinate,'color',ColorOrder(i,:)); hold on;
      end
  end
  
  xbound=xlim;
  ybound=ylim;
  for i=1:length(ordinates_to_plot)
      
      x_fact=0.1;
      y_fact=0.1+0.05*(i-1);
      text(x_fact*xbound(2)+(1-x_fact)*xbound(1),y_fact*ybound(2)+(1-y_fact)*ybound(1),ordinates_to_plot{i},'color',ColorOrder(i,:),'fontsize',15)
  end
  saveas(gcf, [ filename(1:end-4) '_diag' num2str(k)], 'fig')
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
setup_outdir(test_i,'mv');
