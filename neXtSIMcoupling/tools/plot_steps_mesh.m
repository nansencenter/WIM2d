function plot_steps_mesh(rootdir);

if ~exist('rootdir','var');
   %% location of outputs
   run_no   = 2;
   if 0
      %%johansen
      rootdir  = '/Volumes/sim/tim';
   else
      %%external hard disk
      rootdir  = '/Volumes/Tim_Ext_HD2/WORK'
   end
   rootdir  = [rootdir,'/Model-Results/neXtSIM/Oban-test16/run',num2str(run_no)];
end
OVER_WRITE  = 0;
DO_RM       = 1;%don't leave files in working folder after plotting

outdir   = [rootdir,'/simul_out_steps_mat']
indir    = [rootdir,'/simul_in'];
figdir   = [rootdir,'/figs'];
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/simul_out_steps_mesh'];
eval(['!mkdir -p ',figdir]);

%%variables to plot
twlim = .5*[-1 1];%%tau_x range
vbls  = {'Dmax' ,'taux_waves','tauy_waves','c'       ,'h'  ,'log1md','Nfloes'};
cmaps = {'jet'  ,'jet'       ,'jet'       ,'rev_gris','jet','jet'   ,'jet'   };
lims  = {[0,300],twlim       ,.2*twlim    ,[0,1]     ,[0 2],[-3.7 0],[0,250] };

% steps to be loaded
dir0  = dir([outdir,'/simul_out_*_step*.mat']);
N0    = length(dir0)-1;
f0    = strsplit(dir0(1).name,'0.mat');
f0    = f0{1};%start of files
fmt   = ['%',sprintf( '%d.%dd', length(num2str(N0)), length(num2str(N0)) )];

%%get simul_in file
simul_name     = strsplit(dir0(1).name,'_step0.mat');
simul_name     = strsplit(simul_name{1},'simul_out_');
simul_name     = simul_name{2};
saved_simul_in = ['simul_in_',simul_name,'.mat'];
cmd            = ['!cp ',[indir,'/',saved_simul_in],' .'];
eval(cmd);

%% ==================================================
%% shorten list of var's
%% check simul_in to see if waves are present
shorten  = 1;
jkeep    = 1:7;
simul_in = load(saved_simul_in);
simul_in = simul_in.simul_in;
if isfield(simul_in,'wim')
   if simul_in.wim.use_wim==0
      jkeep    = [4,5,6];
      shorten  = 0;
   else
      jkeep    = [1,4,5,6,7];
      %jkeep    = [1,7];
      shorten  = 0;
   end
else
   jkeep    = [4,5,6];
   shorten  = 0;
end
clear simul_in;

if shorten==1
   %%shorten more manually
   jkeep = [1,4,5,6];
end
vbls  = vbls(jkeep);
cmaps = cmaps(jkeep);
lims  = lims(jkeep);
Nv = length(vbls);
%% ==================================================

domain   = 'wim_ideal';
for n=0:N0
   saved_simul_out   = [f0,num2str(n),'.mat'];
   saved_simul_out0  = [outdir,'/',saved_simul_out];

   for k=1:Nv
      vbl   = vbls{k};
      cmap  = cmaps{k};
      lim   = lims{k};
      %%
      ss    = strsplit(saved_simul_out,'step');
      ss2   = strsplit(ss{2},'.mat');
      m     = str2num(ss2{1});

      %% initial/final filenames
      figname0 = [f0,num2str(m),'_',vbl,'.png'];    %%initial
      figname  = [vbl,'/',f0,num2str(m,fmt),'.png'];%%final
      fig_full = [figdir,'/',figname];

      if ~(exist(fig_full)&~OVER_WRITE)
         eval(['!cp ',saved_simul_out0,' ',saved_simul_out]);
         plot_param_v2(vbl,saved_simul_out,domain,cmap,lim,[],'png')
         %%
         eval(['!mkdir -p ',figdir,'/',vbl]);
         eval(['!mv ',figname0,' ',fig_full]);
         disp(['saved to ',fig_full]);
         close;
         %%
         if DO_RM
            eval(['!rm ',saved_simul_out]);
         end%want to delete simul_out from current dir
      end%check if fig is present already
   end%loop over variables
end%loop over time steps

if DO_RM
   eval(['!rm ',saved_simul_in]);
end
