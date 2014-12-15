%% test_WIM2d_F.m
%% Author: Timothy Williams
%% Date:   20141110, 08:37:21 CET

function test_WIM2d_F()
%clear;

SV_FIG   = 1;
IO_OPT   = 0;

if IO_OPT==0
   %% outputs from 'dumb' code with no in/out
   outdir   = '../out';
   matdir   = 'out';%%where to put outputs from this program
else
   %% outputs from code with in/out
   %% - called from python
   outdir   = '../out_io';
   matdir   = 'out_io';%%where to put outputs from this program
end
if ~exist(matdir)
   mkdir(matdir);
end

%%check initialisation
[grid_prams,ice_fields,wave_fields] = fn_check_init(outdir);
if 0
   nn = 1:26;
   tst_masks   = [nn',wave_fields.WAVE_MASK(nn,1)+ice_fields.ICE_MASK(nn,1)]
end

if 1%%plot and save initial conditions
   fig_dir  = [matdir,'/init_cons/'];
   if ~exist(fig_dir)
      mkdir(fig_dir);
   end
   figure(1),clf;
   fn_fullscreen;
   fn_plot_ice(grid_prams,ice_fields);
   saveas(gcf,[fig_dir,'ice.png']);
   %%
   figure(2),clf;
   fn_fullscreen;
   fn_plot_waves(grid_prams,wave_fields);
   saveas(gcf,[fig_dir,'waves.png']);
   %%
   cmd   = ['!cp ',outdir,'/binaries/wim_grid.* ',fig_dir]
   eval(cmd);
   cmd   = ['!cp ',outdir,'/binaries/wim_init.* ',fig_dir]
   eval(cmd);
end

progdir        = [outdir,'/binaries/prog/'];
D              = dir([progdir,'wim_prog*.a']);
nm             = D(1).name;
n0             = str2num(nm(9:11));
nm             = D(2).name;
nstep          = str2num(nm(9:11))-n0;
nm             = D(end).name;
nt             = str2num(nm(9:11))
%%
binary_final   = [outdir,'/binaries/wim_out'];

nvec  = (n0:nstep:nt);

%%further reduce freq of plotting
nstep_   = 50;
nf       = floor(nstep_/nstep);
nvec     = nvec(1:nf:end);

for r = 1:length(nvec)
   n  = nvec(r);
   disp([n,nt]);
   %%
   figure(3),clf;
   fn_fullscreen;
   out_fields  = fn_plot_prog(grid_prams,n,outdir);
   drawnow;
   %GEN_pause;
end

out_fields  = fn_plot_final(grid_prams,outdir);

if 1
   figure(4),clf;
   plot(grid_prams.X(:,1)/1e3,out_fields.Hs(:,1),'-k');
   set(gca,'yscale','log');
   GEN_proc_fig('{\itx}, km','{\itH}_s, m');
   %%
   xmin  = min(grid_prams.X(:))/1e3;%%km
   xmax  = max(grid_prams.X(:))/1e3;%%km
   ymax  = max(grid_prams.Y(:))/1e3;%%km
   axis([xmin,xmax,1e-3,1e1]);
end

%%these parameters determine where to save figure
SOLVER   = out_fields.SOLVER;
nw       = out_fields.n_wave_freq;
ndir     = out_fields.n_wavdir;
disp(out_fields);

if 1
   %% NB this definition won't
   %% work for all configurations
   dx       = grid_prams.dx;
   D_j      = out_fields.Dmax(:,1);
   MIZ_MASK = ((D_j>0)&(D_j<250));
   Wmiz     = sum(MIZ_MASK)*dx/1e3;
   %%
   disp(' ');
   disp(['MIZ width = ',num2str(Wmiz),' km']);
end

taux_max = max(out_fields.tau_x(:));
tauy_max = max(out_fields.tau_y(:));
disp(['max tau_x = ',num2str(taux_max),' Pa']);
disp(['max tau_y = ',num2str(tauy_max),' Pa']);
disp(' ');

if SV_FIG

   %%determine where to save files from parameters
   if nw==1
      if SOLVER==1
         fig_dir  = [matdir,'/isotropic_1freq'];  %%use this for monochromatic wave
      elseif SOLVER==0
         fig_dir  = [matdir,'out/simple_1freq'];  %%use this for monochromatic wave
      end
   else
      if SOLVER==1
         fig_dir  = [matdir,'out/isotropic_spec'];  %%use this for spectrum
      elseif SOLVER==0
         fig_dir  = [matdir,'out/simple_spec'];  %%use this for spectrum
      end
   end
   if ~exist(fig_dir)
      mkdir(fig_dir)
   end

   %%make subdirectories to separate file types
   Dirs  = {[fig_dir,'/binary/'],
            [fig_dir,'/log/'],
            [fig_dir,'/fig/'],
            [fig_dir,'/png/'],
            [fig_dir,'/att_fig/'],
            [fig_dir,'/att_png/']};
   for j=1:length(Dirs)
      if ~exist(Dirs{j})
         mkdir(Dirs{j});
      end
   end

   %%save binary file
   nd3   = num2str(ndir,'%3.3d');
   fn3   = [Dirs{1},'wim_out',nd3];
   cmd   = ['!cp ',binary_final,'.a ',fn3,'.a']
   eval(cmd);
   cmd   = ['!cp ',binary_final,'.b ',fn3,'.b']
   eval(cmd);

   %%save log file
   fil   = [Dirs{2},'wim2d_',nd3,'.log'];
   cmd   = ['!cp ',outdir,'/log/wim2d.log ',fil];
   eval(cmd);

   %%save main figures
   figure(3);
   saveas(gcf,[Dirs{3},'wim_final',nd3,'.fig']);
   saveas(gcf,[Dirs{4},'wim_final',nd3,'.png']);

   %%save plots of Hs
   figure(4);
   pos   = [0.13   0.121428571428571   0.775   0.803571428571429];
   set(gca,'position',pos);
   saveas(gcf,[Dirs{5},'wim_final',nd3,'_atten.fig']);
   saveas(gcf,[Dirs{6},'wim_final',nd3,'_atten.png']);
end
