function [Dmiz,Lmiz,out] =...
   WIM_2D_ideal(Tm,mwd_dim,SHARP_DIST,AC_option)


%% DEPENDENCIES:
%% 1. 'ALP_FXN': provides attenuation coefficients
%% 2. 'SDF_FXN' SDF_PiersonMoscowitz.m,
%%   SDF_Bretschneider.m, SDF_jonswap.m
%% 3. 
%
DO_SAVE  = 0;
DO_PLOT  = 1;%% change this to 0
             %% if graphics aren't supported;
DD_elseoption  = 0
USE_ICE_WLNG   = 0%% if 0, approx ice wlng by water wlng;  
DO_ATTEN       = 1%% if 0, just advect waves
                  %%  without attenuation;

OPT   = 1;%%ice-water-land configuration;
if ~exist('SHARP_DIST')
   SHARP_DIST     = 1
end
%%
SINGLE_FREQ    = 1;
%%
if ~exist('Tm')
   Tm       = 10;%% peak period;
end
Hs = [];%% sig wave height 
        %%  ([]->Pierson-Moskowitz spectrum
if ~exist('mwd_dim')
   mwd_dim  = 180%% waves towards south;
end
mwd      = mwd_dim*pi/180;%% convert dirn to radians;
%%

if ~exist('AC_option')
   AC_option   = 6;
end

if DO_SAVE
   filename=['wim2d_out',num2str(Tm),...
            's_',num2str(mwd_dim),...
            'deg_spread',num2str(SHARP_DIST),'.mat'];
else
   disp('not saving');
   %disp('not saving - push any key');
   %pause;
end

%% set attenuation model;
%% also give progress report every 'reps' time
%%  steps;
reps  = 10;
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
HAVE_WAVES  = exist('wave_fields','var');
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

   nx = 51;
   ny = 51;
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
LANDMASK = grid_prams.LANDMASK;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice/water;
if HAVE_ICE==0
   h           = 2;
   c           = 0.75;
   ice_prams   = struct('c',c,'h',h);
   ice_fields  = iceinit(ice_prams,grid_prams,OPT);
   %% ice_fields  = structure:
   %%      cice: [51x51 double]
   %%      hice: [51x51 double]
   %%      Dmax: [51x51 double]
   %%  WTR_MASK: [51x51 logical]
   %%  ICE_MASK: [51x51 double]
   cice     = ice_fields.cice;
   hice     = ice_fields.hice;
   Dmax     = ice_fields.Dmax;
   WTR_MASK = ice_fields.WTR_MASK;
   ICE_MASK = ice_fields.ICE_MASK;
end
if 0
   fn_plot_ice(grid_prams,ice_fields);
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%waves
if HAVE_WAVES==0
   wave_prams  = struct('Hs',Hs,'Tp',Tp);
   wave_fields = waves_init(grid_prams,ice_fields,wave_prams,OPT);
end

wave_stuff  = get_wave_spec(wave_fields);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%waves

CFL   = .7;
dt    = CFL*dx/max(ag_eff); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waves

if 1%%default
   f      = 1/16;%0.042;% min freq/ resolution
   f1     = 1/2.5;%0.4;% max freq
   nw     = 21;% NB needs to be odd for Simpson's rule;
   freq   = linspace(f,f1,nw);
   T      = 1./freq;
end
om = 2*pi./T;
%%
ndir        = 8;
theta       = linspace(0,2*pi,ndir+1)';
theta(end)  = [];
dtheta      = theta(2);%% directional resolution of
                       %% spectral density function;

%%find closest direction to mwd;
err   = abs(theta-mwd);
jmwd  = find(err==min(err));
mwd   = theta(jmwd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% allocate memory for ice
%%  properties in each grid cell
%%  (constant in time);
cice        = zeros(ny,nx);
hice        = cice;
ag_ice      = zeros(ny,nx,nw);
wavlen_ice  = ag_ice;
%%
Dmax        = zeros(ny,nx);
S           = zeros(ny,nx,nw,ndir);

%% Model Parameters
rhoice  = 922.5;        % Ice density           [kg/m^3]
rho     = 1025.0;       % Water density         [kg/m^3]
rhoave  = 0.5.*(rho+rhoice);
g       = 9.81;         % Gravity               [m/s^2]
%%
young    = 5e9;              % Young's modulus       [Pa]
poisson  = .3;             % Poisson's ratio
lame_lam = young*poisson/(1+poisson)/(1-2*poisson);
lame_mu  = young/2/(1+poisson);
                        % Lame constants \lambda & \mu  [Pa]
flex_rig_coeff = young/12/(1-poisson^2);

salt    = 5;            % Ice salinity          [psu]
temp    = -10;          % Ice temperature       [oC]
% Brine volume (Frankenstein and Gardner 1967)
vb      = salt.*(49.185./abs(temp) + 0.532);  % [ppt]
% Flexural strength (Timco and O'Brien 1994)
sigma_c     = 1.76.*exp(-5.88.*sqrt(10.^-3.*vb)); % [MPa]
sigma_c     = 1.e6.*sigma_c;                      % [Pa]
strain_c2   = (1-poisson^2)/young*sigma_c;


% Fatigue (Langhorne et al. 1998)
mu      = 0.6;          %                       [-]
%  Endurance limit (Langhorne et al. 1998)
strain_c1   = 3.e-5;       %                    [-] 
strain_c    = min(strain_c1,strain_c2);

% Parameters for the floe size distribution
Dmin        = 20;           %                       [m]
xi          = 2;            %                       [-]
fragility   = 0.9;        %                       [-]


% disp('min wave speed = ')
% disp(min(ag))
% disp('max wave speed = ')
% disp(max(ag))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: move this to iceinit.m
%% ICE CONDITIONS:
%% set ice conc's, thicknesses & pre-compute
%%  ice wavelength, group vel, etc;
%% NB may need to move this inside the time loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Water wavelength and wave speed
%% is a function only of wave period
wlng     = g.*T.^2./(2.*pi);
ap       = sqrt(g.*wlng./(2.*pi));    % Phase speed
ag       = ap./2;                  % Group speed

ag_ice      = zeros(nw,ny,nx);
wlng_ice    = zeros(nw,ny,nx);
disp_ratio  = ones(nw,ny,nx);
atten_nond  = zeros(nw,ny,nx);
damping     = zeros(nw,ny,nx);

%% EXTERNAL CALL (GEN_get_ice_groupvel.m):
for i = 1:ny
   for j = 1:nx
      %% get ice wavelengths, group velocities,
      %%  and displacement ratios;
      if cice(i,j)>0 & USE_ICE_WLNG
         [ag_ice(:,i,j),wlng_ice(:,i,j)]  =...
            GEN_get_ice_groupvel(hice(i,j),T,Inf);
         disp_ratio(:,i,j) = wlng./wlng_ice(:,i,j);
      else
         ag_ice(:,i,j)     = ag;
         wlng_ice(:,i,j)   = wlng;
      end

      %% get damping and attenuation from scattering;
      if cice(i,j)>0 & DO_ATTEN
         if AC_option==0
            atten_nond(:,i,j) = alpha_km08(T,hice(i,j));
         elseif AC_option==6%%LB & damping;
            [atten_nond(:,i,j),damping(:,i,j)]  =...
               feval(ALPfxn,(2*pi)./T,hice(i));
         else
            atten_nond(:,i,j) =...
               feval(ALPfxn,(2*pi)./T,hice(i));
         end
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test waves don't cross more than 1 grid cell per time step;
if USE_ICE_WLNG
   ag_fac   = 1.2;
else
   ag_fac   = 1;
end
if ( max(ag_fac*ag)*dt >= dx )%% group vel increases with presence
                              %%  of ice so test 1.5*ag instead of
                              %%  ag itself;
    disp('***  Violation CFL   ***')
    disp(['*** Reduce time step by factor of ',...
           num2str(max(ag_fac*ag)*dt/dx) ,' ***'])
           dt,max(ag),max(ag_fac*ag)
           [dx,max(ag_fac*ag)*dt]
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incident wave & boundary conditions;

%% EXTERNAL CALL (SDF_*.m):
if SINGLE_FREQ
   jtp      = find(abs(T-Tm)==min(abs(T-Tm)));
   S00      = 0*T;
   S00(jtp) = 1/4;%% Hs=1;
   wt_om    = S00;
   if ~isempty(Hs)
      S00   = Hs^2/16*S00;
   end
   disp('warning using only a single frequency');
   disp('push any key');
   pause;
else
   if isempty(Hs)% Pierson-Moskowitz spectrum
      [S00,Hs] = SDF_PiersonMoscowitz(om,{Tm});
   else% Pierson-Moskowitz spectrum
      S00 = SDF_Bretschneider(om,{Tm,Hs});
   end
end

%% add angular factor
%% NB mwd is in radians;
dth   = theta-mwd;
if SHARP_DIST%% have a sharp directional spectrum
             %%  for testing;
   theta_Sfac        = zeros(ndir,1);
   theta_Sfac(jmwd)  = 1;
else
   theta_Sfac  = (1+cos(dth))/2/(pi/2);
      %%\int_0^pi[cos^2\theta]d\theta = pi/2
   for j=1:length(theta)
      th    = theta(j);
      dth0  = dth(j);
      %theta(j)\in[0,2*pi]
      %mwd\in[0,2*pi]
      if th<=mwd
         if ~((abs(dth0)<pi)|(abs(dth0+2*pi)<pi))
            theta_Sfac(j)  = 0;
         end
      else
         if ~((abs(dth0)<pi)|(abs(dth0-2*pi)<pi))
            theta_Sfac(j)  = 0;
         end
      end
   end
end

%% boundary conditions:

%% BUILT-IN CALL (squeeze.m):
for jw   = 1:nw
   for jth  = 1:ndir
      if ~TEST1D
         %%  set spectrum in water
         %%  to the incident wave spectrum;
         for r = 1:length(Jy0)
            S(Jy0{r},Jx0{r},jw,jth)  =...
               S00(jw)*theta_Sfac(jth);
         end
      else
         %%  set spectrum at upper boundary
         %%  to the incident wave spectrum;
         for r = 1:length(Jy_boundary)
            ji = Jy_boundary{r};
            jj = Jx_boundary{r};
            %%
            S(ji,jj,jw,jth)  =...
               S00(jw)*theta_Sfac(jth);
         end
      end
   end
   if 0%% test theta integration
      wt_theta = ones(ndir,1)*(2*pi/ndir);
      [S00(jw), wt_theta'*squeeze(S(ji,jj,jw,:))],pause
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Display parameters
disp('------------------------------------')
disp(['sigma_c    = ' num2str(sigma_c) ' Pa']);
disp(['fragility  = ' num2str(fragility)]);
disp(['fatigue    = ' num2str(mu)]);
disp(['strain_c   = ' num2str(strain_c)]);
disp(['h          = ' num2str(h) ' m const']);
disp(['Tm         = ' num2str(Tm) ' s']);
disp(['Hs         = ' num2str(Hs) ' m']);
disp(['AC_option  = ' num2str(AC_option)]);
disp('');

%% Integration
disp('Integrating ...')
t0       = 24*60*rem(now,1);
brkcrt  = zeros(ny,nx,nt);

%% define weights for numerical quadrature;
if ~SINGLE_FREQ%% weights for integral over frequency
   %% (Simpson's rule);
   wt_simp            = 2+0*T;
   wt_simp([1 end])   = 1;
   wt_simp(2:2:end-1) = 4;

   %%NB om needs to be equally spaced;
   dom    = abs(om(2)-om(1));
   wt_om  = dom/3*wt_simp;
   if 0
      Hs2      = 4*sqrt(wt_om*S00');
      testHs   = [Hs2,Hs]
      pause;
   end
end


%% weights for integral over directions
%% NB using radians for mwd;
if ~SHARP_DIST
   wt_theta = ones(ndir,1)*(2*pi/ndir);
else
   wt_theta       = zeros(ndir,1);
   wt_theta(jmwd) = 1;
end

GET_OUT  = 1;
if GET_OUT
   out         = zeros(ny,nx,floor(nt/reps));
   out(:,:,1)  = Dmax;
end

if DO_PLOT
   %%
   jpp   =...
      find(abs(T-Tm)==min(abs(T-Tm)));
   %%
   subplot(1,3,1);
   %GEN_plot_matrix(xgrid/1e3,ygrid/1e3,cice,[-1 1]);
   H  = pcolor(xgrid/1e3,ygrid/1e3,cice);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   GEN_proc_fig('\itx, \rmkm','\ity, \itkm');
   ttl   = title('Concentration');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,2);
   %GEN_plot_matrix(xgrid/1e3,ygrid/1e3,Dmax,[-1 500]);
   H  = pcolor(xgrid/1e3,ygrid/1e3,Dmax);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \itkm');
   ttl   = title('\itD_{\rmmax}');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,3);
   %GEN_plot_matrix(xgrid/1e3,ygrid/1e3,S(:,:,jpp,jmwd));
   H  = pcolor(xgrid/1e3,ygrid/1e3,S(:,:,jpp,jmwd));
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \itkm');
   colorbar;
   ttl   = title('\itS(T_\rmp,<\theta>), m^2s');
   GEN_font(ttl);
   %%
   disp('beginning main integration (push any key):');
   pause
end

for n = 2:nt
   disp([n nt])

   %% spectral moments;
   mom0  = zeros(ny,nx);
   mom2  = zeros(ny,nx);

   %% variances of stress and strain;
   var_stress     = zeros(ny,nx);
   var_strain     = zeros(ny,nx);

   %% test integrals;
   var_boundary   = cell(1,length(Jy_boundary));
   for r=1:length(Jy_boundary)
      var_boundary{r}   = 0;
   end
   var_boundary0  = var_boundary;  
   %%
   for jw   = 1:nw
      for jth  = 1:ndir
         %% boundary conditions:
         %%  set spectrum at edge of water boundary
         %%  to the incident wave spectrum;
         for r = 1:length(Jy_boundary)
            S(Jy_boundary{r},Jx_boundary{r},jw,jth)   =...
               S00(jw)*theta_Sfac(jth);
            var_boundary0{r}  = var_boundary0{r}+...
               wt_om(jw)*wt_theta(jth)*...
                  S(Jy_boundary{r},Jx_boundary{r},jw,jth);
         end

         %% advect waves;
         wave_props  = {{ag(jw),squeeze(ag_ice(jw,:,:))},...
                           180/pi*theta(jth)};
         grid_props  = {dt,dx,dy};
         %%
         [S(:,:,jw,jth),ag_eff]  =...
            FBwave_advect2D_v1(S(:,:,jw,jth),...
               wave_props,cice,grid_props);

         for r = 1:length(Jy_boundary)
            var_boundary{r}   = var_boundary{r} +...
               wt_om(jw)*wt_theta(jth)*...
                  S(Jy_boundary{r},Jx_boundary{r},jw,jth);
         end
      end%% end directional loop;


      for i = 1:ny
         for j = 1:nx

            %% ATTEN WAVES;
            if cice(i,j)>0
               Dave  = floe_scaling(fragility,xi,...
                        Dmin,Dmax(i,j));

               %% get expected no of floes met per unit
               %%  distance if travelling in a line;
               if Dmax(i,j) < 200
                  c1d = cice(i,j)/Dave;
                  %% floes per unit length;
               else
                  c1d = cice(i)/Dmax(i,j);
                  %% ??
               end

               %% ENERGY attenuation coeff;
               atten_dim   = atten_nond(jw,i,j)*c1d+...
                                 2*damping(jw,i,j)*cice(i,j);
               S(i,j,jw,:) = S(i,j,jw,:)*...
                  exp(-atten_dim*ag_eff(i,j)*dt);
            end

            %% INTEGRATE SPECTRUM OVER DIRECTION;
            Sint              = wt_theta'*squeeze(S(i,j,jw,:));

            %% convert from water amp's to ice amp's;
            F     = disp_ratio(jw,i,j);%%|T| also
            k_ice = 2*pi/wlng_ice(jw,i,j);

            %% SPECTRAL MOMENTS;
            mom0(i,j)   = mom0(i,j)+wt_om(jw)*Sint*F^2;
            mom2(i,j)   = mom2(i,j)+wt_om(jw)*om(jw)^2*Sint*F^2;

            %% VARIANCES OF STRESS AND STRAIN;
            strain_density    = Sint*F^2*...
                                  (k_ice^2*hice(i,j)/2)^2;
            var_strain(i,j)   = var_strain(i,j)+...
                                 + wt_om(jw)*strain_density;
         end%% end spatial loop x;
      end%% end spatial loop y;

   end%% end spectral loop;

   if 1%%plot sig wave height in each cell;
      for r = []%1:length(Jy_boundary)
         Hs
         [4*sqrt(var_boundary0{r});4*sqrt(var_boundary{r})]
      end
      %%
      Hs_all   = sqrt(mom0)*4;
      sig_eps  = sqrt(var_strain)*2;
      Tc       = 2*pi*sqrt(mom0./mom2);
      %%
      chq_Hs   = Hs_all(1:10,1:8)
      chq_Tc   = Tc(1:10,1:8)
      chq_eps  = sig_eps(1:10,1:8)
      chq_Dmax = Dmax(1:10,1:8)
      %%
      [Y,X]    = meshgrid(ygrid,xgrid);
%      subplot(1,2,1)
%      surf(X,Y,Hs_all');
%      set(gca,'CameraPosition',[mean(xgrid),mean(ygrid),...
%         Hs+1]);

%      subplot(1,2,2),
      surf(X,Y,Dmax');
      set(gca,'CameraPosition',[mean(xgrid),mean(ygrid),...
         500+1]);
%         colorbar;
      pause;
   end
   %% FINALLY DO FLOE BREAKING;
   for i=1:ny
      for j=1:nx
         if cice(i,j)>0 & mom0(i,j)>0
            %% only try breaking if ice is present
            %%  & some waves have arrived;
            Nwaves   = dt/2/pi*sqrt(mom2(i,j)/mom0(i,j));
            P_crit   = 1/Nwaves;
            T_crit   = 2*pi*sqrt(mom0(i,j)/mom2(i,j));

            %% significant strain
            sig_strain  = 4*sqrt(var_strain(i,j));

            %% from sig stress
            %%  probability of critical stress
            %%  being exceeded from Rayleigh distribution;
            Pstrain  = exp( -strain_c/(2*var_strain(i,j)) );
               %FB_wave_height_Raleigh_truncdist(strain_c,sig_strain);

            %% FLOE BREAKING:
            BREAK_CRIT     = ( Pstrain>=P_crit );%%breaks if larger than this
            brkcrt(i,j,n)  = BREAK_CRIT;

            if BREAK_CRIT
               %% use crest period to work out wavelength
               %% - half this is max poss floe length;
               if ~USE_ICE_WLNG
                  wlng_crest = g.*T_crit.^2./(2.*pi);
               else
                  wlng_crest =...
                     GEN_get_ice_wavelength(hice(i,j),T_crit); 
               end

               Dc = wlng_crest/2;
               if Dc >= Dmin & Dmax(i,j)>Dc
                  Dmax(i,j)   = Dc;
               end
            end%% end breaking action;

            if 0%i==11 & j==1
               BREAK_CRIT
               Dmax(i,j)
            end

         elseif cice(i,j)<=0%% only water present
            Dmax(i,j)   = 0;
         end
         
      end%% end spatial loop in x;
   end%% end spatial loop in y;

   %% progress report;
   if round(n/reps)==(n/reps)
      if GET_OUT
         out(:,:,n)  = Dmax;
      end
      disp([num2str(n),...
            ' time steps done, out of ',...
               num2str(nt)]);
      %%
      t1    = 24*60*rem(now,1);
      disp(['time taken (mins): ',...
            num2str(t1-t0)]);

      jmiz  = find(Dmax<500&Dmax>0);
      Dmiz  = max(Dmax(jmiz))
      Lmiz  = dx/1e3*length(jmiz);
      if DO_PLOT
         %%plot Dmax
         subplot(1,3,2);
         pcolor(xgrid/1e3,ygrid/1e3,Dmax);
         daspect([1 1 1]);
         GEN_proc_fig('\itx, \rmkm','\ity, \itkm');
         ttl   = title('\itD_{\rm max}, m');
         GEN_font(ttl);
         colorbar;

         %%plot peak period
         subplot(1,3,3);
         pcolor(xgrid/1e3,ygrid/1e3,S(:,:,jpp,jmwd));
         daspect([1 1 1]);
         GEN_proc_fig('\itx, \rmkm','\ity, \itkm');
         colorbar;
         ttl   = title('\itS(T_\rmp,<\theta>), m^2s');
         GEN_font(ttl);
         pause(.1);
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

%% save time-stepped Dmax;
if DO_SAVE
   save(filename,'out');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_ice(grid_prams,ice_fields);

X  = grid_prams.X;
Y  = grid_prams.Y;
%%
cice  = ice_fields.cice;
hice  = ice_fields.hice;
Dmax  = ice_fields.Dmax;

subplot(1,3,1);
%GEN_plot_matrix(X/1e3,Y/1e3,cice,[-1 1]);
H  = pcolor(X/1e3,Y/1e3,cice);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
colorbar;
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('Concentration');
GEN_font(ttl);
colorbar;
%%
subplot(1,3,3);
%GEN_plot_matrix(X/1e3,Y/1e3,Dmax,[-1 500]);
H  = pcolor(X/1e3,Y/1e3,Dmax);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('\itD_{\rmmax}');
GEN_font(ttl);
colorbar;
%%
subplot(1,3,2);
%GEN_plot_matrix(X/1e3,Y/1e3,S(:,:,jpp,jmwd));
H  = pcolor(X/1e3,Y/1e3,hice);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('\ith, \rmm');
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

