%% fig_Boltzmann_Steady.m
%% script to get, plot results from fn_Boltzmann_Steady.m (from Luke Bennetts)
%% Author: Timothy Williams
%% Date: 20150507

period   = 12;
fortyp   = 'freq';
freq     = 1/period;
conc     = 1;
wth      = 100e3;%m
alp_scat = 1.42798715E-04; %m^{-1}
absorb   = 0;
%%
lam0  = freq;
ISO   = alp_scat;

RIGID = 5.49;%GPa - Young's modulus
Param = ParamDef_Default(RIGID);
%%
Ni          = 1;
Vert_Modes  = 200;
Nangles     = 100;
Param = ModParam_def(Param,1,Vert_Modes,0,0,Nangles);
%%
COMM  = 1;

%PLOT  = 0;%% No plotting inside fn_Boltzmann_Steady.m
PLOT  = 2;%% Plot inside

outputs  = 'eigen-info';

figure(1);
clf;
out = fn_Boltzmann_Steady(fortyp, lam0, conc, wth, ISO, absorb,...
                          Param,outputs,COMM,PLOT);
nplots   = 0;
if PLOT==2
   legend_labs = {'SSB (LB)'};
   nplots      = 1;
end

for j=1:length(out)
   disp(out(j).name)
   disp(out(j).value)
   if ~isempty(strfind(out(j).name,'eigen-info'))
      eigen_info  = out(j).value;
   end
end

if 1
   nplots   = nplots+1;
   Hs_inc   = 3;
   x_edge   = -220e3;
   xvec_ss  = x_edge+linspace(0,wth,400).';
   Ndir     = length(eigen_info.angles);%%now different to Nangles - why?
   %%
   [Edir,j_inc]   = fn_Boltzmann_calcEdir(eigen_info,xvec_ss-xvec_ss(1),Hs_inc);
   dtheta         = 2*pi/Ndir;
   m0             = dtheta*sum(Edir,2);
   Hs_ss          = real(4*sqrt(m0));

   if 0
      %%test normalisation of Edir
      m0_inc   = dtheta*sum(Edir(:,incs),2);
      [m0_inc(1),m0(1)]
      [4*sqrt(m0_inc(1)),Hs_ss(1)]
      return
   end

   legend_labs{end+1}   = 'SSB (LB2)';
   figure(1);%%ratio
   hold on;
   plot((xvec_ss-xvec_ss(1))/1e3,Hs_ss/Hs_ss(1),'--r');
   %%
   figure(2);%%absolute value of Hs (m)
   plot(xvec_ss/1e3,Hs_ss,'k');
   xl = xlabel('x, km');
   yl = ylabel('H_s, m');
   set(xl,'fontname','Times','fontsize',14);
   set(yl,'fontname','Times','fontsize',14);
end

%% add plots made by python/fortran code
w2d_path       = getenv('WIM2D_PATH');
figdir         = [w2d_path,'/fortran/run/fig_scripts/figs/TC2S'];
file_name{1}   = [figdir,'/test_steady2.dat'];%%steady state
file_name{2}   = [figdir,'/test_steady2_FT.dat'];%%steady state (FT)
file_name{3}   = [figdir,'/test_steady1.dat'];%%time-dependant code after a long time
cols           = {'--c','--m','-b'};
leg_text       = {'SSB (TW1)','SSB (TW2)','Time-dep'};

DO_INIT  = 1;
for j=1:3
   fname = file_name{j}
   if exist(fname)
      legend_labs{end+1}   = leg_text{j};
      nplots   = nplots+1;
      fid      = fopen(fname);

      %% read past header
      found_hashes   = 0;
      while ~found_hashes
         C  = fgets(fid);
         if length(C)>=5
            if strcmp('#####',C(1:5))
               found_hashes   = 1;

               %% next line will be blank, so skip it
               %% - rest of file is 2 columns of data
               C  = fgets(fid);
            end
         end
      end

      %%read data
      C  = textscan(fid,'%f%f');
      fclose(fid);

      x  = C{1};
      Hs = C{2};
      if DO_INIT
         x0       = x(1);
         DO_INIT  = 0;
      end
      jz = find(x>=x0,1,'first');
      y0 = Hs(jz);
      %%
      figure(1);
      hold on;
      plot((x-x0)/1e3,Hs/y0,cols{j});
      hold off;
      xlim([0,wth/1e3]);
      %%
      figure(2);
      hold on;
      plot(x/1e3,Hs,cols{j});
      xl = xlabel('x, km');
      yl = ylabel('H_s, m');
      set(xl,'fontname','Times','fontsize',14);
      set(yl,'fontname','Times','fontsize',14);
      hold off;
   end
end

%% add legends and save fig's
if ~exist('out','dir')
   !mkdir out
end

figure(1);
box on;
legend(legend_labs{1:nplots});
figname  = 'out/fig_eg_HsVsX_ratio.png';
disp(['Saving to ',figname]);
saveas(gcf,figname);

if nplots>1
   figure(2);
   box on;
   if PLOT==2
      legend(legend_labs{2:nplots});
   else
      legend(legend_labs{1:nplots});
   end
   figname  = 'out/fig_eg_HsVsX_abs.png';
   disp(['Saving to ',figname]);
   saveas(gcf,figname);
end
