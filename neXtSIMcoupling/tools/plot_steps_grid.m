%% location of outputs
if 0
   %%johansen
   rootdir  = '/Volumes/sim/tim/Model-Results/nextsim/test15_waves';
else
   %%external hard disk
   rootdir  = '/Volumes/Tim_Ext_HD2/WORK/Model-Results/neXtSIM/Oban-test16/run1'
end

outdir   = [rootdir,'/simul_out_steps_mat'];
figdir   = [rootdir,'/figs'];
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/simul_out_steps_grid'];
eval(['!mkdir -p ',figdir]);

%%variables to plot
if 0
   twlim = .5*[-1 1];%%tau_x range
   vbls  = {'Dmax' ,'Hs' ,'Tp'  ,'taux_waves','tauy_waves'};
   cmaps = {'jet'  ,'jet','jet' ,'jet'       ,'jet'       };
   lims  = {[0,300],[0 6],[0 20],twlim       ,.1*twlim    };
else
   twlim = .5*[-1 1];%%tau_x range
   vbls  = {'Dmax' ,'taux_waves','tauy_waves'};
   cmaps = {'jet'  ,'jet'       ,'jet'       };
   lims  = {[0,300],twlim       ,.1*twlim    };
end
Nv    = length(vbls);

% steps to be loaded
dir0  = dir([outdir,'/simul_out_*_step*.mat']);
N0    = length(dir0)-1;
f0    = strsplit(dir0(1).name,'0.mat');
f0    = f0{1};%start of files
fmt   = ['%',sprintf( '%d.%dd', length(num2str(N0)), length(num2str(N0)) )];

for n=0:N0
   saved_simul_out0  = [outdir,'/',dir0(n+1).name];
   saved_simul_out   = dir0(n+1).name
   eval(['!cp ',saved_simul_out0,' ',saved_simul_out]);
   %%
   for k=1:Nv
      vbl   = vbls{k};
      cmap  = cmaps{k};
      lim   = lims{k};
      plot_param_grid(vbl,saved_simul_out,domain,cmap,lim,[],'png');
      %pause
      %%
      ss       = strsplit(saved_simul_out,'step');
      ss2      = strsplit(ss{2},'.mat');
      m        = str2num(ss2{1});
      figname0 = [f0,num2str(m),'_',vbl,'.png'];
      figname  = [vbl,'/',f0,num2str(m,fmt),'.png']
      %%
      eval(['!mkdir -p ',figdir,'/',vbl]);
      eval(['!mv ',figname0,' ',figdir,'/',figname]);
      close;
   end
   %pause

   %%
   eval(['!rm ',saved_simul_out]);
end
