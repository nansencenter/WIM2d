%% fig_Boltzmann_Steady.m
%% script to get, plot results from fn_Boltzmann_Steady.m (from Luke Bennetts)
%% Author: Timothy Williams
%% Date: 20150507
function fig_Boltzmann_steady(wth_km)

if ~exist('wth_km','var')
   wth  = 100e3;%m
   wstr = num2str(wth/1e3);% for fig names
elseif wth_km=='inf'
   wth  = 'inf';
   wstr = 'inf';% for fig names
else
   wth  = wth_km*1e3;%m
   wstr = num2str(wth_km);% for fig names
end

period   = 12;
fortyp   = 'freq';
freq     = 1/period;
conc     = 1;
alp_scat = 1.42798715E-04; %m^{-1}
absorb   = 0;
%%
lam0 = freq;
% ISO  = 0; %% use full elastic plate calculation
ISO  = alp_scat;%%isotropic scattering

RIGID = 5.49;%GPa - Young's modulus
Param = ParamDef_Default(RIGID);
%%
Ni         = 1;
Vert_Modes = 200;
Nangles    = 100;
Param = ModParam_def(Param,1,Vert_Modes,0,0,Nangles)
   %% NB Param.MIZ_length is only used for plotting
   %% - wth is the thing used for the calculations
%%
COMM  = 1;

%PLOT  = 0;%% No plotting inside fn_Boltzmann_Steady.m
PLOT  = 2;%% Plot inside - and compare to external energy calculation

outputs  = 'eigen-info';

figure(1);
clf;
out = fn_Boltzmann_Steady(fortyp, lam0, conc, wth, ISO, absorb,...
                          Param,outputs,COMM,PLOT);
legend_labs = {};
col_v2 = '-k';
if PLOT==2
   legend_labs{end+1} = 'SSB (v1)';
   col_v2 = '--r';
end

for j=1:length(out)
   disp(out(j).name)
   disp(out(j).value)
   if ~isempty(strfind(out(j).name,'eigen-info'))
      eigen_info = out(j).value;
   end
end

if 1
   Hs_inc  = 3;
   x_edge  = -220e3;
   if wth == 'inf'
      xvec_ss = x_edge+linspace(0,200e3,400).';
   else
      xvec_ss = x_edge+linspace(0,wth,400).';
   end
   Ndir    = length(eigen_info.angles);%%now different to Nangles - why?
   %%
   [Edir,j_inc] = fn_Boltzmann_calcEdir(eigen_info,xvec_ss-xvec_ss(1),Hs_inc);
   dtheta       = 2*pi/Ndir;
   m0           = dtheta*sum(Edir,2);
   Hs_ss        = real(4*sqrt(m0));
   m0_inc       = dtheta*sum(Edir(1,j_inc),2);%[Hs_inc,4*sqrt(m0_inc)]

   if 0
      %%test normalisation of Edir
      m0_inc   = dtheta*sum(Edir(:,incs),2);
      [m0_inc(1),m0(1)]
      [4*sqrt(m0_inc(1)),Hs_ss(1)]
      return
   end

   legend_labs{end+1} = 'SSB (v2)';
   xp_km = (xvec_ss-xvec_ss(1))/1e3;
   figure(1);%%ratio
   hold on;
   plot(xp_km,Hs_ss/Hs_inc,col_v2);
   xl = xlabel('x, km');
   yl = ylabel('H_s/H_{s,inc}');
   set(xl,'fontname','Times','fontsize',14);
   set(yl,'fontname','Times','fontsize',14);
   %%
   figure(2);%%absolute value of Hs (m)
   plot(xp_km,Hs_ss,'k');
   xl = xlabel('x, km');
   yl = ylabel('H_s, m');
   set(xl,'fontname','Times','fontsize',14);
   set(yl,'fontname','Times','fontsize',14);

   %% print to dat-file:
   outdir = 'out';
   if ~exist(outdir,'dir')
      mkdir(outdir)
   end
   dfile_mat   = [outdir,'/test_steady_mat.dat'];
   disp(' ');
   disp('Saving steady-state results from fn_Boltzmann_[Steady/calcEdir].m to:');
   disp(dfile_mat);
   disp(' ');
   fid   = fopen(dfile_mat,'w');
   fprintf(fid,'%s\n','# Steady-state results');
   fprintf(fid,'%s\n','# >> fig_Boltzmann_Steady.m -> fn_Boltzmann_Steady.m -> fn_Boltzmann_calcEdir.m');
   fprintf(fid,'%s\n','# x (m), Hs (m)');
   fprintf(fid,'%s\n','###############################################');
   fprintf(fid,'%s\n',' ');

   %% print data:
   for loop_x=1:length(xvec_ss)
      fprintf(fid,'%f%s%f\n',xvec_ss(loop_x),'    ',Hs_ss(loop_x));
   end
   fclose(fid);
end

%% add plots made by python/fortran code
w2d_path       = getenv('WIM2D_PATH');
figdir         = [w2d_path,'/fortran/run/fig_scripts/figs/TC2S'];
file_name{1}   = [figdir,'/test_steady2.dat'];%%steady state
file_name{2}   = [figdir,'/test_steady2_FT.dat'];%%steady state (FT)
file_name{3}   = [figdir,'/test_steady1.dat'];%%time-dependant code after a long time
cols           = {'--c','--m','-b'};
leg_text       = {'SSB (TW1)','SSB (TW2)','Time-dep'};

DO_INIT = 1;

%% add legends and save fig's
figure(1);
box on;
if length(legend_labs)>0
   legend(legend_labs);
end
figname  = [outdir,'/fig_eg_HsVsX_ratio_W',wstr,'.png'];
disp(['Saving to ',figname]);
saveas(gcf,figname);

figure(2);
box on;
legend(legend_labs{end});
figname  = [outdir,'/fig_eg_HsVsX_abs_W',wstr,'.png'];
disp(['Saving to ',figname]);
saveas(gcf,figname);
