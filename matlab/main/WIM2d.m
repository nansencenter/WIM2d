function [ice_fields,wave_fields,ice_prams,grid_prams,Dmax_all,brkcrt] =...
   WIM_2D(grid_prams,wave_prams,ice_prams)
% clear;

DO_SAVE     = 0;
DO_PLOT     = 1;  %% change this to 0
                  %% if graphics aren't supported;
USE_ICE_VEL = 0   %% if 0, approx ice group vel by water group vel;  
DO_ATTEN    = 1   %% if 0, just advect waves
                  %%  without attenuation;
DO_BREAKING = 1   %% if 0, turn off breaking for testing
STEADY      = 1   %% Steady-state solution: top-up waves inside wave mask
SOLVER      = 1   %% 0: old way; 1: scatter E isotropically

OPT      = 1;%%ice-water-land configuration;
PLOT_OPT = 2;%%plot option

CHK_ATTEN   = 0;%%check by running with old attenuation

%if DO_SAVE
%   filename=['wim2d_out',num2str(Tm),...
%            's_',num2str(mwd_dim),...
%            'deg_spread','.mat'];
%else
%   disp('not saving');
%   %disp('not saving - push any key');
%   %pause;
%end

%% set attenuation model;
%% also give progress report every 'reps' time
%%  steps;
reps  = 1;
format long

disp('Initialization')

%% Environment
% homedir = '/home/nersc/dany';
% wavedir = [homedir,'/waves'];
% mizdir  = [homedir,'/validation/miz'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%check what inputs we are given;
HAVE_GRID   = exist('grid_prams','var');
HAVE_ICE    = exist('ice_fields','var');
if HAVE_ICE
   if ~isfield(ice_fields,'cice')
      ice_prams   = ice_fields;
      HAVE_ICE    = 0;
   else
      ice_prams   = struct('c'        ,'given',...
                           'h'        ,'given',...
                           'bc_opt'   ,0,...
                           'young_opt',1);
      ice_prams   = fn_fill_iceprams(ice_prams);
   end
end
HAVE_WAVES  = exist('wave_prams','var');
if HAVE_WAVES
   %%check if full info is given or just partial info;
   if isfield(wave_prams,'dir_spec')
      wave_fields = wave_prams;
   else
      HAVE_WAVES  = 0;
   end
end
HAVE3       = HAVE_GRID+HAVE_ICE+HAVE_WAVES;
if (HAVE3>0)&(HAVE3<3)
   disp(' ');
   disp('*********************************************************************');
   disp('Please specify all 3 of ''grid_prams'',''ice_fields'',''wave_fields''');
   disp('(or specify none, and let configuration be chosen through local variable ''OPT'')');
   disp(['HAVE_GRID   = ',num2str(HAVE_GRID )]);
   disp(['HAVE_ICE    = ',num2str(HAVE_ICE  )]);
   disp(['HAVE_WAVES  = ',num2str(HAVE_WAVES)]);
   disp('*********************************************************************');
   disp(' ');
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Grid
if HAVE_GRID==0

   %%can pass in all/none of grid_prams,ice_fields,wave_fields
   %%but not some
   if exist('ice_fields','var')|exist('wave_fields','var')
      disp('WIM2d.m');
      disp('Please specify all 3 of ''grid_prams'',''ice_fields'',''wave_fields''');
      disp('(at least ''grid_prams'' not given)');
      return;
   end

   if 1
      nx = 150;
      ny = 50;
   else
      nx = 149;
      ny = 151;
   end
   dx = 4000; % m
   dy = 4000; % m
   %%
   grid_prams  = struct('nx',nx,'ny',ny,...
                        'dx',dx,'dy',dy);
   grid_prams  = get_grid(grid_prams,OPT);
   %% grid_prams  = structure, eg:
   %%        nx: 51
   %%        ny: 51
   %%        dx: 4000
   %%        dy: 4000
   %%         X: [51x51 double]
   %%         Y: [51x51 double]
   %%  LANDMASK: [51x51 double]
   %%      scuy: [51x51 double]
   %%      scvx: [51x51 double]
   %%      scp2: [51x51 double]
   %%     scp2i: [51x51 double]
end
nx       = grid_prams.nx;
ny       = grid_prams.ny;
dy       = grid_prams.dy;
dy       = grid_prams.dy;
X        = grid_prams.X;
Y        = grid_prams.Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice/water;
if HAVE_ICE==0
   h           = 2;
   c           = 0.75;
   bc_opt      = 0;%%breaking condition (0=beam;1=Marchenko)
   young_opt   = 0;%%young's modulus option
   ice_prams   = struct('c'         ,c,...
                        'h'         ,h,...
                        'bc_opt'    ,bc_opt,...
                        'young_opt' ,young_opt);
   %%
   [ice_fields,ice_prams]  = iceinit(ice_prams,grid_prams,OPT);
   %% ice_fields  = structure:
   %%      cice: [51x51 double]
   %%      hice: [51x51 double]
   %%      Dmax: [51x51 double]
   %%  WTR_MASK: [51x51 logical]
   %%  ICE_MASK: [51x51 double]

   %% ice_prams = structure eg:
   %%               c: 0.750000000000000
   %%               h: 2
   %%           young: 2.000000000000000e+09
   %%          bc_opt: 0
   %%          rhowtr: 1.025000000000000e+03
   %%          rhoice: 9.225000000000000e+02
   %%               g: 9.810000000000000
   %%         poisson: 0.300000000000000
   %%             vbf: 0.100000000000000
   %%              vb: 100
   %%         sigma_c: 2.741429878818372e+05
   %%        strain_c: 1.370714939409186e-04
   %%  flex_rig_coeff: 1.831501831501831e+08
   %%            Dmin: 20
   %%              xi: 2
   %%       fragility: 0.900000000000000
end

cice     = ice_fields.cice;
hice     = ice_fields.hice;
Dmax     = ice_fields.Dmax;
WTR_MASK = ice_fields.WTR_MASK;
ICE_MASK = ice_fields.ICE_MASK;

%% add wave stress computation
ice_fields.tau_x  = 0*cice;
ice_fields.tau_y  = 0*cice;

%% Model Parameters
rho_wtr  = ice_prams.rhowtr;   % Water density   [m/s^2]
g        = ice_prams.g;        % Gravity         [m/s^2]
strain_c = ice_prams.strain_c; % breaking strain [-] 
young    = ice_prams.young;    % Young's modulus [Pa]
visc_rp  = ice_prams.visc_rp;  % Robinson-Palmer coeff [Pa/(m/s)]

% Parameters for the floe size distribution
Dmin        = ice_prams.Dmin;      % [m]
xi          = ice_prams.xi;        % [-]
fragility   = ice_prams.fragility; % [-]

if 0
   fn_plot_ice(grid_prams,ice_fields);
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%waves
if HAVE_WAVES==0
   if ~exist('wave_prams','var');
      Hs0         = 2;
      Tp0         = 12;
      %Tp0         = 6;
      wave_prams  = struct('Hs',Hs0,...
                           'Tp',Tp0);
   end
   wave_fields = waves_init(grid_prams,wave_prams,ice_fields,OPT);
   WAVE_MASK   = wave_fields.WAVE_MASK;
   wave_stuff  = set_incident_waves(grid_prams,wave_fields);
end

nw       = wave_stuff.nfreq;     %% number of frequencies
om       = 2*pi*wave_stuff.freq; %% radial freq
ndir     = wave_stuff.ndir;      %% number of directions
wavdir   = wave_stuff.dirs;      %% wave from, degrees, clockwise
Sdir     = wave_stuff.dir_spec;  %% initial directional spectrum
if STEADY==1
   S_inc    = Sdir;
   theta    = -pi/180*(90+wavdir);
   [i1,j1]  = find(WAVE_MASK==1,1,'first');
   J_STEADY = find(S_inc(i1,j1,:,ceil(nw/2))>0);
   % squeeze(S_inc(i1,j1,:))
   % GEN_pause;
end
%%
T  = 2*pi./om;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Water wavelength and wave speed
%% is a function only of wave period
wlng     = g.*T.^2./(2.*pi);
ap       = sqrt(g.*wlng./(2.*pi)); % Phase speed
ag       = ap./2;                  % Group speed

ag_eff      = zeros(nx,ny,nw);
ap_eff      = zeros(nx,ny,nw);
wlng_ice    = zeros(nx,ny,nw);
disp_ratio  = ones (nx,ny,nw);
atten_nond  = zeros(nx,ny,nw);
damping     = zeros(nx,ny,nw);

for i = 1:nx
for j = 1:ny
   if ICE_MASK(i,j)==1
      %% if ice is present:
      %% get ice wavelengths, group velocities,
      %% displacement ratios (convert wtr disp to ice),
      %% attenuation;
      if DO_ATTEN==1
         [damping_rp,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(hice(i,j),om,young,visc_rp);
         %%
         if CHK_ATTEN==1
            %%check with old version
         end
         atten_nond(i,j,:) = alp_scat;
         damping(i,j,:)    = damping_rp;
      else
         [damping_rp,kice,kwtr,int_adm,NDprams] =...
            RT_param_outer(hice(i,j),om,young,visc_rp);
         modT  = 1;
      end

      if USE_ICE_VEL==0
         %%use wtr group vel;
         ag_eff(i,j,:)  = ag;
         ap_eff(i,j,:)  = ap;
      else
         %%TODO check if this is correct
         %%weighted avg of ice and wtr group vel;
         ag_ice         = GEN_get_ice_groupvel(hice(i,j),T,Inf,young);
         ag_eff(i,j,:)  = cice(i,j)*ag_ice+...
                           +(1-cice(i,j))*ag;

         %%weighted avg of ice and wtr phase vel;
         ap_ice         = om./k_ice;
         ap_eff(i,j,:)  = cice(i,j)*ap_ice+...
                           +(1-cice(i,j))*ap;
      end

      wlng_ice(i,j,:)   = 2*pi./kice;
      disp_ratio(i,j,:) = wlng.*(kice/2/pi).*modT;
   else
      ag_eff(i,j,:)     = ag;
      ap_eff(i,j,:)     = ap;
      wlng_ice(i,j,:)   = wlng;
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CFL   = .7;
amax  = max(ag_eff(:))
dt    = CFL*dx/max(ag_eff(:)); 

if 0
   nt    = 50;
   L     = nt*dt*amax
else
   L     = max(X(:))-min(X(:));
   amin  = min(ag_eff(:));
   uc    = amin+.7*(amax-amin);
   nt    = round( L/uc/dt );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display parameters
Info  = { '------------------------------------';
         ['sigma_c    = ' num2str(ice_prams.sigma_c) ' Pa'];
         ['fragility  = ' num2str(fragility)];
         ['strain_c   = ' num2str(strain_c)];
         ['h          = ' num2str(ice_prams.h) ' m const'];
         ['Tp         = ' num2str(wave_prams.Tp) ' s'];
         ['Hs         = ' num2str(wave_prams.Hs) ' m'];
         ['CFL        = ' num2str(CFL)];
         ['dt         = ' num2str(dt)];
         ['nt         = ' num2str(nt)];
         ['nfreq      = ' num2str(nw)];
         ['ndir       = ' num2str(ndir)];
         ['SOLVER     = ' num2str(SOLVER)];
         '------------------------------------';
         ' '};

disp(strvcat(Info));

%% Integration
t0       = now;   %%days
t0_fac   = 24*60; %%days to minutes
brkcrt  = zeros(nx,ny,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define weights for numerical quadrature;
if nw>1%% weights for integral over frequency
   %% (Simpson's rule);
   wt_simp            = 2+0*om;
   wt_simp([1 end])   = 1;
   wt_simp(2:2:end-1) = 4;

   %%NB om needs to be equally spaced;
   dom    = abs(om(2)-om(1));
   wt_om  = dom/3*wt_simp;

   %if 0
   %   Hs2      = 4*sqrt(wt_om*S00');
   %   testHs   = [Hs2,Hs0]
   %   pause;
   %end
else
   wt_om = 1;
end

% %% weights for integral over directions
% %% NB using radians for mwd;
% if ndir>1
%    wt_theta = ones(ndir,1)*(2*pi/ndir);
% else
%    wt_theta = 1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GET_OUT  = 1;
if GET_OUT
   Dmax_all         = zeros(nx,ny,floor(nt/reps));
   Dmax_all(:,:,1)  = Dmax;
end

if DO_PLOT
   %%
   figure(1),clf;
   fn_full_screen;
   fn_plot_ice(grid_prams,ice_fields);
   pause(0.1);
   %%
   figure(2),clf;
   fn_full_screen;
   Tc    = 12;%check this period
   jchq  = find(abs(T-Tc)==min(abs(T-Tc)));
   jdir  = round(ndir/2);
   {jdir,wave_stuff.dirs,Sdir}
   s1 = struct('dir',wave_stuff.dirs(jdir),...
               'period',Tc,...
               'Sdir',Sdir(:,:,jdir,jchq));
   %%
   if PLOT_OPT==1
      fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
   else
      fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
         Dmax,ice_fields.tau_y);
   end
   if OPT==1
      subplot(2,2,1);
      hold on;
      x0 = min(X(:));
      x1 = X(find(WAVE_MASK(:,1)==0,1,'first'),1);
      % {x0/1e3,x1/1e3}
      yc = .3*max(Y(:))/1e3;
      x_ = [min(X(:)),x0,x0,x1,x1,max(X(:))]/1e3;
      y_ = [0,0,yc*[1,1],0,0];
      plot(x_,y_,'k');
      plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc/Hs0,'--k');
      hold off;
   end
   %%
   clear s1;
   %X,Y
   %GEN_pause
   pause(0.1);
end

disp('beginning main integration...');

for n = 2:nt
   disp([n nt])

   %% spectral moments;
   mom0  = zeros(nx,ny);
   mom2  = zeros(nx,ny);

   %% wave stresses;
   tau_x = zeros(nx,ny);
   tau_y = zeros(nx,ny);

   %% variances of stress and strain;
   var_stress     = zeros(nx,ny);
   var_strain     = zeros(nx,ny);

   % %% test integrals;
   % var_boundary   = cell(1,length(Jy_boundary));
   % for r=1:length(Jy_boundary)
   %    var_boundary{r}   = 0;
   % end
   % var_boundary0  = var_boundary;  
   %%
   for jw   = 1:nw

      %% CALC DIMENSIONAL ATTEN COEFF;
      atten_dim   = 0*X;
      damp_dim    = 0*X;
      for i = 1:nx
      for j = 1:ny

         %%top-up waves in wave mask if STEADY==1
         %%(steady-state solution);
         if WAVE_MASK(i,j)>0 & STEADY==1
            Sdir(i,j,J_STEADY,:)  = S_inc(i,j,J_STEADY,:);
         end
         
         if ICE_MASK(i,j)>0 & DO_ATTEN==1
            Dave  = floe_scaling(fragility,xi,...
                     Dmin,Dmax(i,j));

            %% get expected no of floes met per unit
            %%  distance if travelling in a line;
            if Dmax(i,j) < 200
               c1d = cice(i,j)/Dave;
               %% floes per unit length;
            else
               c1d = cice(i)/Dmax(i,j);
               %% uniform lengths
            end

            %% ENERGY attenuation coeff;
            atten_dim(i,j) = atten_nond(i,j,jw)*c1d;%%scattering
            damp_dim(i,j)  = 2*damping(i,j,jw)*cice(i,j);%%damping
         end
      end% j
      end% i
      % max(ag_eff(:))
      % max(atten_dim(:))
      % max(damp_dim(:))
      %1e3*max(ag_eff(:))*max(atten_dim(:))
      % GEN_pause

      s1.ndir        = ndir;
      s1.wavdir      = wavdir;
      %s1.Sdir        = reshape( Sdir(:,:,jw,:), nx,ny,ndir);
      s1.Sdir        = Sdir(:,:,:,jw);
      s1.ag_eff      = ag_eff(:,:,jw);
      s1.atten_dim   = atten_dim;
      s1.damp_dim    = damp_dim;
      s1.ICE_MASK    = ice_fields.ICE_MASK;

      if ndir==1
         if SOLVER~=0
            disp('warning: changing SOLVER option as not enough directions');
            disp(['(ndir = ',num2str(ndir)]);
         end
         SOLVER   = 0;
      end

      if SOLVER==0
         %% Simple attenuation scheme - doesn't conserve scattered energy
         [Sdir(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
            adv_atten_simple(grid_prams,ice_prams,s1,dt);
         clear s1 S_out;
      elseif SOLVER==1
         %% same as SOLVER==0, but scattered energy
         %% is distributed isotropically
         [Sdir(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
            adv_atten_isotropic(grid_prams,ice_prams,s1,dt);
         clear s1 S_out;
      end

      %% integrate stress densities over frequency
      %% TODO: check if this is correct for ice-covered water
      tmp   = rho_wtr*g*tau_x_om./ap_eff(:,:,jw);  %%[Pa*s]
      tau_x = tau_x+wt_om(jw)*tmp;                 %%[Pa]
      tmp   = rho_wtr*g*tau_y_om./ap_eff(:,:,jw);  %%[Pa*s]
      tau_y = tau_y+wt_om(jw)*tmp;                 %%[Pa]
      clear tmp;
      %GEN_pause

      %% INTEGRALS FOR BREAKING PROB:
      for i = 1:nx
      for j = 1:ny

         %% INTEGRATE SPECTRUM OVER DIRECTION;
         %S_freq   = wt_theta'*squeeze(Sdir(i,j,jw,:));

         %% convert from water amp's to ice amp's;
         F     = disp_ratio(i,j,jw);%%|T| also
         k_ice = 2*pi/wlng_ice(i,j,jw);

         %% SPECTRAL MOMENTS;
         %%take abs as small errors can make S_freq negative
         mom0(i,j)   = mom0(i,j)+abs( wt_om(jw)*S_freq(i,j)*F^2 );
         mom2(i,j)   = mom2(i,j)+abs( wt_om(jw)*om(jw)^2*S_freq(i,j)*F^2 );

         if ICE_MASK(i,j)==1
            %% VARIANCE OF STRAIN;
            strain_density    = abs( S_freq(i,j)*F^2*...
                                  (k_ice^2*hice(i,j)/2)^2 );
            var_strain(i,j)   = var_strain(i,j)+...
                                 + wt_om(jw)*strain_density;
         end
      end%% end spatial loop x;
      end%% end spatial loop y;

   end%% end spectral loop;
   %mom0,mom2,wlng_ice,return

   % if 1%%plot sig wave height in each cell;
   %    for r = []%1:length(Jy_boundary)
   %       Hs
   %       [4*sqrt(var_boundary0{r});4*sqrt(var_boundary{r})]
   %    end
   %    %%
   %    Hs_all   = sqrt(mom0)*4;
   %    sig_eps  = sqrt(var_strain)*2;
   %    Tc       = 2*pi*sqrt(mom0./mom2);
   %    %%
   %    chq_Hs   = Hs_all(1:10,1:8)
   %    chq_Tc   = Tc(1:10,1:8)
   %    chq_eps  = sig_eps(1:10,1:8)
   %    chq_Dmax = Dmax(1:10,1:8)
   %    %%
   %    pcolor(X,Y,Dmax);
   %    pause;
   % end

   %%calc Hs, Tw (into wave_fields.Tp)
   wave_fields.Hs       = 4*sqrt(mom0);
   wave_fields.Tp       = 0*X;
   jnz                  = find(mom2>0);
   wave_fields.Tp(jnz)  = 2*pi*sqrt(mom0(jnz)./mom2(jnz));

   %% wave stresses
   ice_fields.tau_x  = tau_x;
   ice_fields.tau_y  = tau_y;

   %% FINALLY DO FLOE BREAKING;
   for i=1:nx
   for j=1:ny
      if ICE_MASK(i,j)==1 & mom0(i,j)>0
         %% only try breaking if ice is present
         %%  & some waves have arrived;

         %% significant strain amp
         sig_strain  = 2*sqrt(var_strain(i,j));

         %%  probability of critical strain
         %%  being exceeded from Rayleigh distribution;
         Pstrain  = exp( -strain_c^2/(2*var_strain(i,j)) );
         P_crit   = (1-DO_BREAKING)+exp(-1);%%this is critical prob if monochromatic wave

         %% FLOE BREAKING:
         BREAK_CRIT     = ( Pstrain>=P_crit );%%breaks if larger than this
         brkcrt(i,j,n)  = BREAK_CRIT;

         if BREAK_CRIT
            %% use crest period to work out wavelength
            %% - half this is max poss floe length;
            T_crit      = wave_fields.Tp(i,j);
            wlng_crest  = ...
               GEN_get_ice_wavelength(hice(i,j),T_crit,Inf,young);

            Dc = max(Dmin,wlng_crest/2);
            Dmax(i,j)   = min(Dc,Dmax(i,j));
         end%% end breaking action;

         if 0%i==11 & j==1
            BREAK_CRIT
            Dmax(i,j)
         end

      elseif WTR_MASK(i,j)==1%% only water present
         Dmax(i,j)   = 0;
      end
      
   end%% end spatial loop j in y;
   end%% end spatial loop i in x;

   ice_fields.Dmax   = Dmax;
   jmiz              = find((Dmax>0)&(Dmax<250));
   % {jmiz}

   %% progress report;
   if round(n/reps)==(n/reps)
      if GET_OUT
         Dmax_all(:,:,n/reps) = Dmax;
      end
      disp([num2str(n),...
            ' time steps done, out of ',...
               num2str(nt)]);
      %%
      t1    = now;
      disp(['time taken (mins): ',...
            num2str(t0_fac*(t1-t0))]);

      if DO_PLOT
         figure(2),clf;
         %%
         if PLOT_OPT==1
            s1 = struct('dir',wave_stuff.dirs(jdir),...
                        'period',Tc,...
                        'Sdir',Sdir(:,:,jdir,jchq));
            fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
         else
            fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
               Dmax,ice_fields.tau_y);
         end

         if OPT==1
            subplot(2,2,1);
            hold on;
            x0 = min(X(:))+uc*n*dt;
            x1 = X(find(WAVE_MASK(:,1)==0,1,'first'),1)+uc*n*dt;
            % {uc,x0/1e3,x1/1e3}
            yc = .3*max(Y(:))/1e3;
            x_ = [min(X(:)),x0,x0,x1,x1,max(X(:))]/1e3;
            y_ = [0,0,yc*[1,1],0,0];
            plot(x_,y_,'k');
            plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc/Hs0,'--k');
            hold off;
            %%
            subplot(2,2,2);
            hold on;
            plot(x_,y_,'k');
            plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc/Hs0,'--k');
            hold off;
         end

         clear s1;
         %GEN_pause;
         pause(0.1);
      end
      %%
      if 0
         Jlook = nx-16+(1:16);
         Dmax(edge:edge+9,Jlook)
      end
   end


   for jw=[]%11
   [n,T(jw)],%[ag_ice(jw,1:11,1)]
   testSatt=S(1:11,1,jw,jmwd)
   pause
   end

end%% end time loop

if DO_PLOT%%check exponential attenuation
   figure(2),clf;
   if PLOT_OPT==1
      s1 = struct('dir',wave_stuff.dirs(jdir),...
                  'period',Tc,...
                  'Sdir',Sdir(:,:,jdir,jchq));
      fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
   else
      fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
         Dmax,ice_fields.tau_y);
   end
   %%
   figure(3);
   plot(X(:,1)/1e3,wave_fields.Hs(:,1),'-k');
   set(gca,'yscale','log');
   GEN_proc_fig('{\itx}, km','{\itH}_s, m');

   if 1%%save figures
      fig_dir  = 'test_B/B';  %%use this for monochromatic wave
      %fig_dir  = 'test_B2/B'; %%use this for full freq spec

      figure(3);
      saveas(gcf,[fig_dir,num2str(ndir,'%3.3d'),'_atten.fig']);
      saveas(gcf,[fig_dir,num2str(ndir,'%3.3d'),'_atten.png']);
      %%
      figure(2);
      saveas(gcf,[fig_dir,num2str(ndir,'%3.3d'),'.fig']);
      saveas(gcf,[fig_dir,num2str(ndir,'%3.3d'),'.png']);
   end
end

%%display info again
disp(strvcat(Info));

%% save time-stepped Dmax;
if DO_SAVE
   save(filename,'out');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec(X,Y,Hs,Tw,Dmax,s1)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

subplot(2,2,1);
H  = pcolor(X/1e3,Y/1e3,Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

subplot(2,2,2);
H  = pcolor(X/1e3,Y/1e3,Dmax);
caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

subplot(2,2,3);
H  = pcolor(X/1e3,Y/1e3,Tw);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itT}_{\rm w}, m');
GEN_font(ttl);

subplot(2,2,4);
H  = pcolor(X/1e3,Y/1e3,s1.Sdir);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = [num2str(s1.period),'s, ',num2str(s1.dir),'^o'];
ttl   = title({'\itS, \rmm^2s',ttl});
GEN_font(ttl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec_2(X,Y,Hs,tau_x,Dmax,tau_y)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

subplot(2,2,1);
H  = pcolor(X/1e3,Y/1e3,Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

subplot(2,2,2);
H  = pcolor(X/1e3,Y/1e3,Dmax);
caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

subplot(2,2,3);
H  = pcolor(X/1e3,Y/1e3,tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

subplot(2,2,4);
H  = pcolor(X/1e3,Y/1e3,tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
