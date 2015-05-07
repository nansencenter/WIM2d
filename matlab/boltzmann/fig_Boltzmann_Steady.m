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

out = fn_Boltzmann_Steady(fortyp, lam0, conc, wth, ISO, absorb);
for j=1:length(out)
   disp(out(j).name)
   disp(out(j).value)
end

%% add plots made by python/fortran code
figdir         = '/Users/timill/GITHUB-REPOSITORIES/WIM2d/fortran/run/fig_scripts/figs/TC2S';
file_name{1}   = [figdir,'/test_steady1.dat'];%%time-dependant code after a long time
file_name{2}   = [figdir,'/test_steady2.dat'];%%steady state
file_name{3}   = [figdir,'/test_steady2_FT.dat'];%%steady state (FT)
cols           = {'-b','--r','--c'};
legend_labs    = {'SSB (LB)','SSB (TW1)','SSB (TW2)','Time-dep'};

nplots   = 1;
for j=3:-1:1
   fname = file_name{j};
   if exist(fname)
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
      if nplots==2
         x0 = x(1);
      end
      jz = find(x>=x0,1,'first');
      y0 = Hs(jz);
      %%
      figure(1);
      hold on;
      plot((x-x0),Hs/y0,cols{j});
      hold off;
      xlim([0,wth]);
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
   legend(legend_labs{2:nplots});
   figname  = 'out/fig_eg_HsVsX_abs.png';
   disp(['Saving to ',figname]);
   saveas(gcf,figname);
end
