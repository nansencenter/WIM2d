%% location of outputs
rootdir  = '/Volumes/sim/tim/Model-Results/nextsim/test15_waves';
outdir   = [rootdir,'/simul_out_steps_mat'];
figdir   = [rootdir,'/figs'];
eval(['!mkdir -p ',figdir]);
figdir   = [figdir,'/simul_out_steps_grid'];
eval(['!mkdir -p ',figdir]);

%%variables to plot
if 1
   twlim = .5*[-1 1];%%tau_x range
   vbls  = {'c','h','log1md','Dmax','taux_waves','tauy_waves'};
   cmaps = {'rev_gris','jet','jet','jet','jet','jet'};
   lims  = {[0,1],[0 2],[-3.7 0],[0,300],twlim,.2*twlim};
elseif 1
   twlim = .5*[-1 1];%%tau_x range
   vbls  = {'tauy_waves'};
   cmaps = {'jet'};
   lims  = {.2*twlim};
end
Nv    = length(vbls);

% steps to be loaded
dir0  = dir([outdir,'/simul_out_*_step*.mat']);
N0    = length(dir0)-1;
f0    = strsplit(dir0(1).name,'0.mat');
f0    = f0{1};%start of files
fmt   = ['%',sprintf( '%d.%dd', length(num2str(N0)), length(num2str(N0)) )];

for n=0:N0
   disp([n,N0]);
   %%
   saved_simul_out0  = [outdir,'/',dir0(n+1).name];
   saved_simul_out   = dir0(n+1).name
   eval(['!cp ',saved_simul_out0,' ',saved_simul_out]);
   %%
   for k=1:Nv
      vbl   = vbls{k};
      cmap  = cmaps{k};
      lim   = lims{k};
      plot_param_v2(vbl,saved_simul_out,domain,cmap,lim,[],'png')
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
   %%
   eval(['!rm ',saved_simul_out]);
end


%plot_param_v2('h',saved_simul_out,domain,'jet',[0 5],[],'pdf')
%
%to_step   = length(dir(['simul_out_' meshfile(1:end-4) '_' simul_in.simul_in_name '_step*.mat']))-1 ;
%saved_simul_out=['simul_out_' meshfile(1:end-4) '_' simul_in.simul_in_name '_step' num2str(to_step) '.mat']
%plot_param_v2('c',saved_simul_out,domain,'rev_gris',[0 1],[],'pdf')
%plot_param_v2('h',saved_simul_out,domain,'jet',[0 5],[],'pdf')
%plot_param_v2('log1md',saved_simul_out,domain,'jet',[-3.7 0],[],'pdf')

